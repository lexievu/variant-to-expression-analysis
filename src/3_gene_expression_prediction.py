"""AlphaGenome expression-impact prediction for somatic variants.

Reads a filtered VCF (by default the output of ``2_vcf_filter.py``), queries
the AlphaGenome API for each variant, and writes a TSV of predicted expression
changes (log₂ fold-change).

Robustness features
-------------------
* Exponential-backoff retry on transient API / network errors.
* Configurable rate-limiting delay between API calls.
* Full stack-trace logging on every failure.
* Checkpoint / resume — skips variants already present in the output TSV.
* Graceful file handling via context managers.
* Input validation (API key, VCF path).
* Progress counter so you can monitor long runs.

Usage
-----
    python 3_gene_expression_prediction.py                     # defaults
    python 3_gene_expression_prediction.py --vcf my.vcf -o out.tsv
    python 3_gene_expression_prediction.py --resume            # pick up where you left off
"""

import argparse
import logging
import os
import sys
import time
import traceback

import numpy as np
from cyvcf2 import VCF
from alphagenome.models import dna_client
from alphagenome.data import genome
from dotenv import load_dotenv

from constants import HIGH_IMPACT_VCF, EXAMPLE_RNA_PATH
import utils

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = "log/gene_expression_prediction.log"
DEFAULT_VCF = HIGH_IMPACT_VCF
DEFAULT_OUTPUT = "output/alphagenome_hits.tsv"
DEFAULT_TISSUE = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

# Retry / rate-limit settings
MAX_RETRIES = 3
RETRY_BASE_DELAY = 2.0          # seconds; doubles each retry
RATE_LIMIT_DELAY = 0.5          # seconds between successive API calls

# Expression fold-change thresholds
FC_GAIN_THRESHOLD = 1.0
FC_LOSS_THRESHOLD = -1.0

# Vaccine target scoring thresholds
TPM_EXPRESSED_THRESHOLD = 1.0   # TPM >= 1 → gene is expressed
VAF_CLONAL_THRESHOLD = 0.2     # VAF >= 0.2 → likely clonal variant

DEFAULT_RNA = EXAMPLE_RNA_PATH

HEADER = (
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
        description="Predict expression impact of somatic variants via AlphaGenome.",
    )
    parser.add_argument(
        "--vcf", default=DEFAULT_VCF,
        help=f"Input VCF of filtered variants (default: {DEFAULT_VCF}).",
    )
    parser.add_argument(
        "--output", "-o", default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--tissue", default=DEFAULT_TISSUE,
        help=f"UBERON ontology ID for tissue type (default: {DEFAULT_TISSUE}).",
    )
    parser.add_argument(
        "--rna", default=DEFAULT_RNA,
        help=f"Path to RNA-seq CSV for TPM lookup (default: {DEFAULT_RNA}).",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume a previous run — skip variants already in the output file.",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_checkpoint(output_file):
    """Return a set of (chrom, pos, ref, alt) tuples already in *output_file*.

    If the file doesn't exist or is empty, returns an empty set.
    """
    done = set()
    if not os.path.isfile(output_file):
        return done
    try:
        with open(output_file) as fh:
            for line in fh:
                if line.startswith("CHROM"):
                    continue  # skip header
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    done.add((parts[0], int(parts[1]), parts[2], parts[3]))
    except Exception:
        logging.warning(
            "Could not read checkpoint file %s — starting fresh.\n%s",
            output_file, traceback.format_exc(),
        )
    logging.info("Checkpoint: %d variant(s) already processed.", len(done))
    return done


def _predict_with_retry(model, interval, ag_variant, tissue_id,
                        max_retries=MAX_RETRIES,
                        base_delay=RETRY_BASE_DELAY):
    """Call *model.predict_variant* with exponential-backoff retry.

    Returns the prediction outputs on success, or raises the last
    exception after all retries are exhausted.
    """
    last_exc = None
    for attempt in range(1, max_retries + 1):
        try:
            return model.predict_variant(
                interval=interval,
                variant=ag_variant,
                ontology_terms=[tissue_id],
                requested_outputs=[dna_client.OutputType.RNA_SEQ],
            )
        except Exception as exc:
            last_exc = exc
            if attempt < max_retries:
                delay = base_delay * (2 ** (attempt - 1))
                logging.warning(
                    "Attempt %d/%d failed (%s). Retrying in %.1fs…",
                    attempt, max_retries, exc, delay,
                )
                time.sleep(delay)
            else:
                logging.error(
                    "All %d attempts failed. Last error: %s",
                    max_retries, exc,
                )
    raise last_exc  # type: ignore[misc]


def _classify(log2_fc):
    """Return a human-readable status label for a log₂ fold-change."""
    if log2_fc > FC_GAIN_THRESHOLD:
        return "Gain_of_Expression"
    if log2_fc < FC_LOSS_THRESHOLD:
        return "Loss_of_Expression"
    return "Neutral"


def _compute_log2_fc(ref_sum, alt_sum):
    """Compute log₂ fold-change, guarding against division by zero."""
    if ref_sum == 0 and alt_sum == 0:
        return 0.0
    return float(np.log2((alt_sum + 1e-9) / (ref_sum + 1e-9)))


def _vaccine_priority(log2_fc, vaf, tpm, nmd_flag):
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
# Main pipeline
# ---------------------------------------------------------------------------

def run_predictions(
    vcf_file=DEFAULT_VCF,
    output_file=DEFAULT_OUTPUT,
    tissue_id=DEFAULT_TISSUE,
    rna_file=DEFAULT_RNA,
    resume=False,
    dna_length=DNA_SEQUENCE_LENGTH,
):
    """Iterate over variants in *vcf_file* and write AlphaGenome predictions."""

    # --- Validate inputs ---------------------------------------------------
    load_dotenv()
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        logging.critical("ALPHAGENOME_API_KEY not set in environment or .env")
        sys.exit(1)

    if not os.path.isfile(vcf_file):
        logging.critical("VCF file not found: %s", vcf_file)
        sys.exit(1)

    # --- Checkpoint --------------------------------------------------------
    already_done = set()
    file_mode = "w"
    if resume:
        already_done = _load_checkpoint(output_file)
        if already_done:
            file_mode = "a"  # append to existing results

    # --- TPM lookup --------------------------------------------------------
    tpm_lookup = utils.load_tpm_lookup(rna_file)

    # --- Initialise model & VCF --------------------------------------------
    logging.info("Input VCF:  %s", vcf_file)
    logging.info("Output TSV: %s", output_file)
    logging.info("Tissue:     %s", tissue_id)
    logging.info("Resume:     %s", resume)
    logging.info("Connecting to AlphaGenome…")
    model = dna_client.create(api_key)
    vcf = VCF(vcf_file)

    count_total = 0
    count_skipped = 0
    count_saved = 0
    count_errors = 0
    start_time = time.time()

    with open(output_file, file_mode) as outfile:
        if file_mode == "w":
            outfile.write(HEADER)

        for variant in vcf:
            if variant.FILTER is not None:
                continue

            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF

            if not variant.ALT:
                logging.warning("Skipping %s:%d — no ALT allele", chrom, pos)
                continue
            alt = variant.ALT[0]
            count_total += 1

            # --- Checkpoint: skip if already done --------------------------
            variant_key = (chrom, pos, ref, alt)
            if variant_key in already_done:
                count_skipped += 1
                logging.debug("Skipping %s:%d (already processed)", chrom, pos)
                continue

            gene_name = utils.get_gene_name(variant)
            gene_id = utils.get_gene_id(variant)  # version-stripped

            # --- Extract additional metrics --------------------------------
            vaf = utils.get_vaf(variant, tumor_index=1)
            nmd_flag = utils.has_nmd(variant)
            observed_tpm = tpm_lookup.get(gene_id, float('nan'))
            expressed = (not np.isnan(observed_tpm)) and observed_tpm >= TPM_EXPRESSED_THRESHOLD

            ag_variant = genome.Variant(chrom, pos, ref, alt)
            start = max(1, pos - dna_length // 2)
            end = pos + dna_length // 2
            interval = genome.Interval(chrom, start, end)

            try:
                # --- Predict (with retry) ----------------------------------
                outputs = _predict_with_retry(
                    model, interval, ag_variant, tissue_id,
                )

                # --- Score -------------------------------------------------
                ref_sum = float(np.sum(outputs.reference.rna_seq.values))
                alt_sum = float(np.sum(outputs.alternate.rna_seq.values))
                log2_fc = _compute_log2_fc(ref_sum, alt_sum)
                status = _classify(log2_fc)
                priority = _vaccine_priority(log2_fc, vaf, observed_tpm, nmd_flag)

                tpm_str = f"{observed_tpm:.2f}" if not np.isnan(observed_tpm) else "."
                vaf_str = f"{vaf:.3f}" if not np.isnan(vaf) else "."

                line = (
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_name}\t{gene_id}"
                    f"\t{ref_sum:.2f}\t{alt_sum:.2f}\t{log2_fc:.4f}\t{status}"
                    f"\t{vaf_str}\t{tpm_str}\t{expressed}\t{nmd_flag}\t{priority}\n"
                )
                outfile.write(line)
                outfile.flush()
                count_saved += 1

                logging.info(
                    "[%d/%d] %s:%d %s — %s (log2FC=%.4f, VAF=%s, TPM=%s, NMD=%s, priority=%s)",
                    count_saved, count_total, chrom, pos, gene_name,
                    status, log2_fc, vaf_str, tpm_str, nmd_flag, priority,
                )

            except Exception:
                count_errors += 1
                logging.error(
                    "Failed permanently at %s:%d (%s/%s). Traceback:\n%s",
                    chrom, pos, ref, alt, traceback.format_exc(),
                )

            # --- Rate-limit between API calls ------------------------------
            time.sleep(RATE_LIMIT_DELAY)

    vcf.close()

    # --- Summary -----------------------------------------------------------
    elapsed = time.time() - start_time
    logging.info(
        "Done in %.1fs — %d total variants, %d saved, %d skipped (checkpoint), %d errors.",
        elapsed, count_total, count_saved, count_skipped, count_errors,
    )
    if count_errors:
        logging.warning(
            "%d variant(s) failed after all retries. See log for stack traces.",
            count_errors,
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    run_predictions(
        vcf_file=args.vcf,
        output_file=args.output,
        tissue_id=args.tissue,
        rna_file=args.rna,
        resume=args.resume,
    )