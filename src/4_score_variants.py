"""Score raw AlphaGenome predictions with biological context.

Reads the raw predictions TSV produced by ``3_gene_expression_prediction.py``
and enriches each variant with:

* **LOG2_FC** — log₂(ALT_EXPR / REF_EXPR)
* **STATUS** — Gain_of_Expression / Loss_of_Expression / Neutral
* **VAF** — variant allele frequency from the VCF tumour sample
* **OBSERVED_TPM** — RNA-seq expression (TPM) for the gene
* **EXPRESSED** — whether the gene is expressed above a threshold
* **NMD_FLAG** — whether the variant triggers nonsense-mediated decay
* **VACCINE_PRIORITY** — composite HIGH / MEDIUM / LOW suitability score

This script is *cheap* — it never touches the AlphaGenome API, so you can
re-run it freely while tuning thresholds or adding new metrics.

Usage
-----
    python 4_score_variants.py                              # defaults
    python 4_score_variants.py --predictions raw.tsv --vcf input.vcf
    python 4_score_variants.py -o my_scored.tsv
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
from cyvcf2 import VCF

from constants import (
    HIGH_IMPACT_VCF,
    EXAMPLE_RNA_PATH,
    RAW_PREDICTIONS,
    SCORED_VARIANTS,
)
import utils

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = "log/score_variants.log"

# Expression fold-change thresholds
FC_GAIN_THRESHOLD = 1.0
FC_LOSS_THRESHOLD = -1.0

# Vaccine target scoring thresholds
TPM_EXPRESSED_THRESHOLD = 1.0   # TPM >= 1 → gene is expressed
VAF_CLONAL_THRESHOLD = 0.2     # VAF >= 0.2 → likely clonal variant

SCORED_HEADER = (
    "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID"
    "\tREF_EXPR\tALT_EXPR\tLOG2_FC\tSTATUS"
    "\tVAF\tOBSERVED_TPM\tEXPRESSED\tNMD_FLAG\tVACCINE_PRIORITY\n"
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Score raw AlphaGenome predictions with biological context.",
    )
    parser.add_argument(
        "--predictions", default=RAW_PREDICTIONS,
        help=f"Input raw-predictions TSV (default: {RAW_PREDICTIONS}).",
    )
    parser.add_argument(
        "--vcf", default=HIGH_IMPACT_VCF,
        help=f"Input VCF for VAF / NMD extraction (default: {HIGH_IMPACT_VCF}).",
    )
    parser.add_argument(
        "--rna", default=EXAMPLE_RNA_PATH,
        help=f"RNA-seq CSV for TPM lookup (default: {EXAMPLE_RNA_PATH}).",
    )
    parser.add_argument(
        "--output", "-o", default=SCORED_VARIANTS,
        help=f"Output scored TSV path (default: {SCORED_VARIANTS}).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

def compute_log2_fc(ref_sum, alt_sum):
    """Compute log₂ fold-change, guarding against division by zero."""
    if ref_sum == 0 and alt_sum == 0:
        return 0.0
    return float(np.log2((alt_sum + 1e-9) / (ref_sum + 1e-9)))


def classify(log2_fc):
    """Return a human-readable status label for a log₂ fold-change."""
    if log2_fc > FC_GAIN_THRESHOLD:
        return "Gain_of_Expression"
    if log2_fc < FC_LOSS_THRESHOLD:
        return "Loss_of_Expression"
    return "Neutral"


def vaccine_priority(log2_fc, vaf, tpm, nmd_flag):
    """Assign a composite vaccine-suitability priority.

    A *good* vaccine target is a clonal somatic variant in an expressed
    gene whose mutant allele is predicted to maintain (or increase)
    expression and is not subject to nonsense-mediated decay.

    Returns one of: HIGH, MEDIUM, LOW.
    """
    expressed = (not np.isnan(tpm)) and tpm >= TPM_EXPRESSED_THRESHOLD
    clonal = (not np.isnan(vaf)) and vaf >= VAF_CLONAL_THRESHOLD
    loss = log2_fc < FC_LOSS_THRESHOLD

    # Automatic LOW: NMD kills the transcript, or predicted expression loss
    if nmd_flag or loss:
        return "LOW"
    # HIGH: expressed gene, clonal variant, no predicted loss
    if expressed and clonal:
        return "HIGH"
    # Everything else
    return "MEDIUM"


# ---------------------------------------------------------------------------
# VCF index for VAF / NMD lookup
# ---------------------------------------------------------------------------

def _build_vcf_index(vcf_file):
    """Build a dict mapping (chrom, pos, ref, alt) → {vaf, nmd} from *vcf_file*.

    This lets us look up per-variant VCF fields efficiently without
    re-iterating the VCF for every prediction row.
    """
    index = {}
    if not os.path.isfile(vcf_file):
        logging.warning("VCF not found (%s) — VAF / NMD will be unavailable.", vcf_file)
        return index

    vcf = VCF(vcf_file)
    for variant in vcf:
        if variant.FILTER is not None:
            continue
        if not variant.ALT:
            continue
        key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
        index[key] = {
            "vaf": utils.get_vaf(variant, tumor_index=1),
            "nmd": utils.has_nmd(variant),
        }
    vcf.close()
    logging.info("VCF index: %d variant(s) loaded.", len(index))
    return index


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def score_variants(
    predictions_file=RAW_PREDICTIONS,
    vcf_file=HIGH_IMPACT_VCF,
    rna_file=EXAMPLE_RNA_PATH,
    output_file=SCORED_VARIANTS,
):
    """Read raw predictions, enrich with biological metrics, write scored TSV."""

    # --- Validate inputs ---------------------------------------------------
    if not os.path.isfile(predictions_file):
        logging.critical("Raw predictions file not found: %s", predictions_file)
        sys.exit(1)

    # --- Load supporting data ----------------------------------------------
    tpm_lookup = utils.load_tpm_lookup(rna_file)
    vcf_index = _build_vcf_index(vcf_file)

    # --- Read raw predictions ----------------------------------------------
    raw = pd.read_csv(predictions_file, sep="\t")
    logging.info("Loaded %d raw prediction(s) from %s", len(raw), predictions_file)

    if raw.empty:
        logging.warning("No predictions to score.")
        with open(output_file, "w") as fh:
            fh.write(SCORED_HEADER)
        return

    # --- Score each variant ------------------------------------------------
    rows = []
    for _, r in raw.iterrows():
        chrom = str(r["CHROM"])
        pos = int(r["POS"])
        ref = str(r["REF"])
        alt = str(r["ALT"])
        gene = str(r["GENE"])
        gene_id = str(r["GENE_ID"])
        ref_expr = float(r["REF_EXPR"])
        alt_expr = float(r["ALT_EXPR"])

        log2_fc = compute_log2_fc(ref_expr, alt_expr)
        status = classify(log2_fc)

        # VCF-derived metrics
        vcf_data = vcf_index.get((chrom, pos, ref, alt), {})
        vaf = vcf_data.get("vaf", float("nan"))
        nmd_flag = vcf_data.get("nmd", False)

        # RNA-seq TPM
        observed_tpm = tpm_lookup.get(gene_id, float("nan"))
        expressed = (not np.isnan(observed_tpm)) and observed_tpm >= TPM_EXPRESSED_THRESHOLD

        priority = vaccine_priority(log2_fc, vaf, observed_tpm, nmd_flag)

        tpm_str = f"{observed_tpm:.2f}" if not np.isnan(observed_tpm) else "."
        vaf_str = f"{vaf:.3f}" if not np.isnan(vaf) else "."

        rows.append(
            f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{gene_id}"
            f"\t{ref_expr:.6f}\t{alt_expr:.6f}\t{log2_fc:.4f}\t{status}"
            f"\t{vaf_str}\t{tpm_str}\t{expressed}\t{nmd_flag}\t{priority}"
        )

    # --- Write output ------------------------------------------------------
    with open(output_file, "w") as fh:
        fh.write(SCORED_HEADER)
        for row in rows:
            fh.write(row + "\n")

    logging.info("Wrote %d scored variant(s) to %s", len(rows), output_file)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    score_variants(
        predictions_file=args.predictions,
        vcf_file=args.vcf,
        rna_file=args.rna,
        output_file=args.output,
    )
