"""Compare variant gene expression against GTEx lung tissue baselines.

Queries the GTEx Portal API (v2) for median gene-level TPM in normal lung
tissue, joins with the validation table, and classifies each gene's
silencing status:

- **Tissue-normal silence** — gene is lowly expressed in both the patient's
  tumour RNA-seq *and* in GTEx normal lung.  Likely reflects normal lung
  biology rather than a tumour-specific event.
- **Tumour-specific silencing** — gene is expressed in GTEx normal lung but
  silenced in the patient's tumour RNA-seq.  Suggests the tumour has lost
  expression of this gene.
- **Tumour over-expression** — gene is much more highly expressed in the
  tumour than in normal lung tissue.

This script does **not** call the AlphaGenome API.

Usage
-----
    python s6_gtex_baseline.py                            # defaults
    python s6_gtex_baseline.py --table output/my_table.csv
    python s6_gtex_baseline.py --tissue Lung              # default
"""

import argparse
import logging
import os
import sys
import time

import pandas as pd
import requests

from constants import VALIDATION_TABLE, GTEX_COMPARISON
import utils

# ---------------------------------------------------------------------------
# GTEx API configuration
# ---------------------------------------------------------------------------
GTEX_API_BASE = "https://gtexportal.org/api/v2"
GTEX_DATASET = "gtex_v8"
DEFAULT_TISSUE = "Lung"

# Thresholds
TPM_EXPRESSED_THRESHOLD = 1.0   # TPM ≥ 1 → "expressed"
OVEREXPR_FOLD = 4.0             # tumour/GTEx ratio above this → "over-expressed"

# Rate-limit: be polite to the public API
API_DELAY = 0.3  # seconds between requests

LOG_FILENAME = "log/gtex_baseline.log"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare variant gene expression against GTEx lung baselines.",
    )
    parser.add_argument(
        "--table", default=VALIDATION_TABLE,
        help=f"Input validation table CSV (default: {VALIDATION_TABLE}).",
    )
    parser.add_argument(
        "--output", "-o", default=GTEX_COMPARISON,
        help=f"Output GTEx comparison CSV (default: {GTEX_COMPARISON}).",
    )
    parser.add_argument(
        "--tissue", default=DEFAULT_TISSUE,
        help=f"GTEx tissue site detail ID (default: {DEFAULT_TISSUE}).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# GTEx API queries
# ---------------------------------------------------------------------------

def _get_json(url, params=None, retries=3):
    """GET request with retry logic."""
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, ValueError) as exc:
            logging.warning("GTEx API attempt %d/%d failed: %s",
                            attempt + 1, retries, exc)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    return None


def resolve_gencode_id(gene_id_stripped):
    """Look up the versioned GENCODE ID for a stripped Ensembl gene ID.

    The GTEx median-expression endpoint requires the *versioned* GENCODE ID
    (e.g. ``ENSG00000141736.13``).  This function queries the GTEx reference
    gene endpoint to resolve the version.
    """
    data = _get_json(
        f"{GTEX_API_BASE}/reference/gene",
        params={"geneId": gene_id_stripped, "datasetId": GTEX_DATASET},
    )
    if data and data.get("data"):
        return data["data"][0].get("gencodeId")
    return None


def query_median_expression(gencode_id, tissue):
    """Query GTEx for median TPM of a gene in the specified tissue.

    Returns the median TPM as a float, or ``None`` if not found.
    """
    data = _get_json(
        f"{GTEX_API_BASE}/expression/medianGeneExpression",
        params={
            "gencodeId": gencode_id,
            "datasetId": GTEX_DATASET,
            "tissueSiteDetailId": tissue,
        },
    )
    if data and data.get("data"):
        return data["data"][0].get("median")
    return None


def fetch_gtex_baselines(gene_ids, tissue):
    """Fetch GTEx median lung TPM for a list of stripped Ensembl gene IDs.

    Returns a dict ``{gene_id_stripped: median_tpm}``.
    """
    results = {}
    for gene_id in gene_ids:
        # Step 1: resolve versioned GENCODE ID
        gencode_id = resolve_gencode_id(gene_id)
        if not gencode_id:
            logging.warning("Gene %s not found in GTEx reference.", gene_id)
            results[gene_id] = None
            time.sleep(API_DELAY)
            continue

        # Step 2: query median expression
        median_tpm = query_median_expression(gencode_id, tissue)
        if median_tpm is not None:
            logging.info("  %s (%s) — GTEx %s median TPM = %.2f",
                         gene_id, gencode_id, tissue, median_tpm)
        else:
            logging.warning("  %s (%s) — no expression data for tissue %s",
                            gene_id, gencode_id, tissue)
        results[gene_id] = median_tpm
        time.sleep(API_DELAY)

    return results


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------

def classify_silencing(tumour_tpm, gtex_tpm):
    """Classify a gene's expression status relative to GTEx normal lung.

    Parameters
    ----------
    tumour_tpm : float
        Observed TPM in the patient's tumour RNA-seq.
    gtex_tpm : float or None
        Median TPM in GTEx normal lung tissue.

    Returns
    -------
    str
        One of:
        - ``"tissue-normal silence"`` — low in both tumour and GTEx
        - ``"tumour-specific silencing"`` — expressed in GTEx, silenced in tumour
        - ``"tumour over-expression"`` — tumour TPM >> GTEx TPM
        - ``"comparable"`` — similar expression in both
        - ``"no GTEx data"`` — GTEx lookup failed
    """
    if gtex_tpm is None:
        return "no GTEx data"

    tumour_expressed = tumour_tpm >= TPM_EXPRESSED_THRESHOLD
    gtex_expressed = gtex_tpm >= TPM_EXPRESSED_THRESHOLD

    if not tumour_expressed and not gtex_expressed:
        return "tissue-normal silence"
    if not tumour_expressed and gtex_expressed:
        return "tumour-specific silencing"

    # Both expressed — check for over-expression
    if gtex_tpm > 0 and tumour_tpm / gtex_tpm >= OVEREXPR_FOLD:
        return "tumour over-expression"

    return "comparable"


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def gtex_baseline(
    table_path=VALIDATION_TABLE,
    output_path=GTEX_COMPARISON,
    tissue=DEFAULT_TISSUE,
):
    """Load the validation table, query GTEx, and write the comparison CSV."""

    if not os.path.isfile(table_path):
        logging.critical("Validation table not found: %s", table_path)
        sys.exit(1)

    # --- Load validation table ---------------------------------------------
    df = pd.read_csv(table_path)
    logging.info("Loaded %d variants from %s", len(df), table_path)

    gene_ids = df["GENE_ID"].dropna().unique().tolist()
    logging.info("Querying GTEx %s for %d unique gene(s)...", tissue, len(gene_ids))

    # --- Fetch GTEx baselines ----------------------------------------------
    gtex_lookup = fetch_gtex_baselines(gene_ids, tissue)

    # --- Join & classify ---------------------------------------------------
    df["GTEX_LUNG_TPM"] = df["GENE_ID"].map(gtex_lookup)
    df["SILENCING_CLASS"] = df.apply(
        lambda row: classify_silencing(row["OBSERVED_TPM"], row["GTEX_LUNG_TPM"]),
        axis=1,
    )

    # Fold-change: tumour / GTEx
    df["TUMOUR_VS_GTEX_RATIO"] = df.apply(
        lambda row: (
            round(row["OBSERVED_TPM"] / row["GTEX_LUNG_TPM"], 2)
            if pd.notna(row["GTEX_LUNG_TPM"]) and row["GTEX_LUNG_TPM"] > 0
            else None
        ),
        axis=1,
    )

    # --- Write output ------------------------------------------------------
    out_cols = [
        "GENE", "GENE_ID", "OBSERVED_TPM", "GTEX_LUNG_TPM",
        "TUMOUR_VS_GTEX_RATIO", "SILENCING_CLASS",
        "LOG2_FC", "VAF", "NMD_FLAG", "VACCINE_PRIORITY",
    ]
    available = [c for c in out_cols if c in df.columns]
    df[available].to_csv(output_path, index=False)
    logging.info("GTEx comparison → %s (%d rows)", output_path, len(df))

    # --- Console summary ---------------------------------------------------
    print(f"\n{'='*70}")
    print("  GTEx Baseline Comparison")
    print(f"{'='*70}")
    print(f"  Tissue:  {tissue}")
    print(f"  Genes:   {len(gene_ids)}")
    print(f"  Output:  {output_path}")
    print(f"{'='*70}\n")

    summary = df[["GENE", "OBSERVED_TPM", "GTEX_LUNG_TPM",
                   "TUMOUR_VS_GTEX_RATIO", "SILENCING_CLASS"]].copy()
    summary.columns = ["Gene", "Tumour TPM", "GTEx Lung TPM",
                        "Tumour/GTEx", "Classification"]
    print(summary.to_string(index=False))
    print()

    # Classification counts
    print("Classification summary:")
    for cls, count in df["SILENCING_CLASS"].value_counts().items():
        print(f"  {cls}: {count}")
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    gtex_baseline(
        table_path=args.table,
        output_path=args.output,
        tissue=args.tissue,
    )
