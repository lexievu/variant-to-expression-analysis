"""AlphaGenome expression prediction — raw API output.

Reads a filtered VCF (by default the output of ``s2_vcf_filter.py``), queries
the AlphaGenome API for each variant, and writes a **raw predictions TSV**
containing only the expensive API outputs.

Downstream scoring (VAF, TPM, NMD, vaccine priority) is handled by
``s4_score_variants.py``, which can be re-run cheaply without touching the API.

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
    python s3_gene_expression_prediction.py                     # defaults
    python s3_gene_expression_prediction.py --vcf my.vcf -o out.tsv
    python s3_gene_expression_prediction.py --resume            # pick up where you left off
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

from constants import HIGH_IMPACT_VCF, RAW_PREDICTIONS
import utils

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = "log/gene_expression_prediction.log"
DEFAULT_VCF = HIGH_IMPACT_VCF
DEFAULT_OUTPUT = RAW_PREDICTIONS
DEFAULT_TISSUE = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

# Retry / rate-limit settings
MAX_RETRIES = 3
RETRY_BASE_DELAY = 2.0          # seconds; doubles each retry
RATE_LIMIT_DELAY = 0.5          # seconds between successive API calls

RAW_HEADER = "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tREF_EXPR\tALT_EXPR\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Query AlphaGenome for expression predictions (raw output).",
    )
    parser.add_argument(
        "--vcf", default=DEFAULT_VCF,
        help=f"Input VCF of filtered variants (default: {DEFAULT_VCF}).",
    )
    parser.add_argument(
        "--output", "-o", default=DEFAULT_OUTPUT,
        help=f"Output raw-predictions TSV path (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--tissue", default=DEFAULT_TISSUE,
        help=f"UBERON ontology ID for tissue type (default: {DEFAULT_TISSUE}).",
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


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_predictions(
    vcf_file=DEFAULT_VCF,
    output_file=DEFAULT_OUTPUT,
    tissue_id=DEFAULT_TISSUE,
    resume=False,
    dna_length=DNA_SEQUENCE_LENGTH,
):
    """Query AlphaGenome for each variant and write raw predictions."""

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
            outfile.write(RAW_HEADER)

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

            ag_variant = genome.Variant(chrom, pos, ref, alt)
            start = max(1, pos - dna_length // 2)
            end = pos + dna_length // 2
            interval = genome.Interval(chrom, start, end)

            try:
                # --- Predict (with retry) ----------------------------------
                outputs = _predict_with_retry(
                    model, interval, ag_variant, tissue_id,
                )

                # --- Write raw API output only -----------------------------
                ref_sum = float(np.sum(outputs.reference.rna_seq.values))
                alt_sum = float(np.sum(outputs.alternate.rna_seq.values))

                line = (
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_name}\t{gene_id}"
                    f"\t{ref_sum:.6f}\t{alt_sum:.6f}\n"
                )
                outfile.write(line)
                outfile.flush()
                count_saved += 1

                logging.info(
                    "[%d/%d] %s:%d %s — ref=%.6f alt=%.6f",
                    count_saved, count_total, chrom, pos, gene_name,
                    ref_sum, alt_sum,
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
        resume=args.resume,
    )