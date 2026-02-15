"""Validate AlphaGenome predictions against patient RNA-seq data.

Reads the scored predictions TSV (output of ``s4_score_variants.py``) and the
patient RNA-seq CSV, joins on Ensembl gene ID, and produces two CSV files:

1. **validation_table.csv** — per-variant data enriched with raw counts and
   FPKM, ready for plotting.
2. **validation_correlations.csv** — Pearson & Spearman correlation statistics
   for several predicted-vs-observed comparisons.

This script does **not** call the AlphaGenome API.

Usage
-----
    python s5_validate.py                                     # defaults
    python s5_validate.py --scored my_scored.tsv --rna my_rna.csv
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

from src.constants import (
    EXAMPLE_RNA_PATH,
    SCORED_VARIANTS,
    VALIDATION_TABLE,
    VALIDATION_CORRELATIONS,
    LOG_DIR,
)
from src import utils
from src.exceptions import PipelineInputError

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = str(LOG_DIR / "validation.log")
TPM_EXPRESSED_THRESHOLD = 1.0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate AlphaGenome predictions against patient RNA-seq.",
    )
    parser.add_argument(
        "--scored", default=SCORED_VARIANTS,
        help=f"Input scored-predictions TSV (default: {SCORED_VARIANTS}).",
    )
    parser.add_argument(
        "--rna", default=EXAMPLE_RNA_PATH,
        help=f"Patient RNA-seq CSV (default: {EXAMPLE_RNA_PATH}).",
    )
    parser.add_argument(
        "--table", default=VALIDATION_TABLE,
        help=f"Output per-variant CSV (default: {VALIDATION_TABLE}).",
    )
    parser.add_argument(
        "--correlations", default=VALIDATION_CORRELATIONS,
        help=f"Output correlations CSV (default: {VALIDATION_CORRELATIONS}).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_scored(path):
    """Load the scored variants TSV, coercing '.' to NaN."""
    df = pd.read_csv(path, sep="\t", na_values=".")
    logging.info("Loaded %d scored variant(s) from %s", len(df), path)
    return df


def load_rna(path):
    """Load the RNA-seq CSV and return a DataFrame indexed by stripped gene ID."""
    rna = pd.read_csv(path, comment="#")
    rna = rna[rna["gene_id"].str.startswith("ENSG", na=False)].copy()
    rna["gene_id_stripped"] = rna["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    rna = rna.set_index("gene_id_stripped")
    logging.info("Loaded RNA-seq data for %d genes from %s", len(rna), path)
    return rna


# ---------------------------------------------------------------------------
# Correlation analysis
# ---------------------------------------------------------------------------

def _corr_row(df, x_col, y_col, label, transform=None):
    """Compute one correlation and return a dict (one row of the output CSV).

    Parameters
    ----------
    transform : str or None
        If ``"log10"``, apply log₁₀(x + 1) to both columns before correlating.
    """
    mask = df[[x_col, y_col]].notna().all(axis=1)
    x = df.loc[mask, x_col].astype(float).values
    y = df.loc[mask, y_col].astype(float).values

    if transform == "log10":
        x = np.log10(x + 1)
        y = np.log10(y + 1)

    n = len(x)
    row = {
        "comparison": label,
        "x": x_col,
        "y": y_col,
        "transform": transform or "none",
        "n": n,
    }

    if n < 3:
        row.update(pearson_r=np.nan, pearson_p=np.nan,
                   spearman_rho=np.nan, spearman_p=np.nan)
        logging.warning("Only %d points for '%s' — skipping.", n, label)
    else:
        pr, pp = stats.pearsonr(x, y)
        sr, sp = stats.spearmanr(x, y)
        row.update(pearson_r=pr, pearson_p=pp,
                   spearman_rho=sr, spearman_p=sp)
    return row


def compute_all_correlations(df):
    """Return a list of correlation dicts for all comparisons of interest."""
    rows = []

    # Predicted fold-change vs observed TPM
    rows.append(_corr_row(df, "LOG2_FC", "OBSERVED_TPM",
                          "LOG2_FC vs OBSERVED_TPM"))
    rows.append(_corr_row(df, "LOG2_FC", "OBSERVED_TPM",
                          "LOG2_FC vs OBSERVED_TPM (log₁₀)",
                          transform="log10"))

    # Predicted fold-change vs raw counts (if available)
    if "unstranded" in df.columns:
        rows.append(_corr_row(df, "LOG2_FC", "unstranded",
                              "LOG2_FC vs raw_counts"))

    # Stratified: expressed genes only (TPM ≥ 1)
    expressed = df[df["OBSERVED_TPM"] >= TPM_EXPRESSED_THRESHOLD]
    if len(expressed) >= 3:
        rows.append(_corr_row(expressed, "LOG2_FC", "OBSERVED_TPM",
                              "LOG2_FC vs TPM (expressed only, TPM≥1)"))

    return rows


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def validate(
    scored_path=SCORED_VARIANTS,
    rna_path=EXAMPLE_RNA_PATH,
    table_path=VALIDATION_TABLE,
    correlations_path=VALIDATION_CORRELATIONS,
):
    """Join predictions with RNA-seq and write CSV outputs."""

    utils.validate_file(scored_path, "Scored variants file")

    # --- Load & join -------------------------------------------------------
    scored = load_scored(scored_path)
    rna = load_rna(rna_path)

    extra_cols = ["unstranded", "fpkm_unstranded", "tpm_unstranded"]
    extra = rna[[c for c in extra_cols if c in rna.columns]]

    df = scored.merge(extra, left_on="GENE_ID", right_index=True, how="left")

    # Cross-check TPM
    if "tpm_unstranded" in df.columns:
        tpm_match = np.allclose(
            df["OBSERVED_TPM"].fillna(0),
            df["tpm_unstranded"].fillna(0),
            atol=0.01,
        )
        logging.info("TPM cross-check: %s",
                      "PASS ✅" if tpm_match else "MISMATCH ⚠️")

    # --- Validation table CSV ----------------------------------------------
    out_cols = [
        "CHROM", "POS", "REF", "ALT", "GENE", "GENE_ID",
        "LOG2_FC", "STATUS",
        "VAF", "OBSERVED_TPM", "EXPRESSED", "NMD_FLAG", "VACCINE_PRIORITY",
        "unstranded", "fpkm_unstranded",
    ]
    available = [c for c in out_cols if c in df.columns]
    df[available].to_csv(table_path, index=False)
    logging.info("Validation table → %s (%d rows)", table_path, len(df))

    # --- Correlations CSV --------------------------------------------------
    corr_rows = compute_all_correlations(df)
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(correlations_path, index=False)
    logging.info("Correlations → %s (%d comparisons)", correlations_path, len(corr_df))

    # --- Summary ------------------------------------------------------------
    logging.info("="*60)
    logging.info("  Validation outputs")
    logging.info("="*60)
    logging.info("  Table:        %s  (%d variants)", table_path, len(df))
    logging.info("  Correlations: %s  (%d comparisons)", correlations_path, len(corr_df))
    logging.info("="*60)
    logging.info("\n%s", corr_df.to_string(index=False))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    try:
        validate(
            scored_path=args.scored,
            rna_path=args.rna,
            table_path=args.table,
            correlations_path=args.correlations,
        )
    except PipelineInputError as exc:
        logging.critical("%s", exc)
        sys.exit(1)
