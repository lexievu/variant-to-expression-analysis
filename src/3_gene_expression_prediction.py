"""AlphaGenome expression-impact prediction for somatic variants.

Reads a filtered VCF, queries the AlphaGenome API for each variant, and
writes a TSV of predicted expression changes (log₂ fold-change).

Robustness features
-------------------
* Exponential-backoff retry on transient API / network errors.
* Configurable rate-limiting delay between API calls.
* Full stack-trace logging on every failure.
* Graceful file handling via context managers.
* Input validation (API key, VCF path).
* Progress counter so you can monitor long runs.
"""

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

from constants import QUICK_10
import utils

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
LOG_FILENAME = "log/gene_expression_prediction.log"
VCF_FILE = QUICK_10
OUTPUT_FILENAME = "output/alphagenome_hits.tsv"
TISSUE_ID = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

# Retry / rate-limit settings
MAX_RETRIES = 3                 # attempts per variant (1 initial + retries)
RETRY_BASE_DELAY = 2.0          # seconds; doubles each retry
RATE_LIMIT_DELAY = 0.5          # seconds between successive API calls

# Expression fold-change thresholds
FC_GAIN_THRESHOLD = 1.0         # log2 FC > 1  → Gain_of_Expression
FC_LOSS_THRESHOLD = -1.0        # log2 FC < -1 → Loss_of_Expression

HEADER = (
    "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID"
    "\tREF_EXPR\tALT_EXPR\tLOG2_FC\tSTATUS\n"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
                    "All %d attempts failed for variant. Last error: %s",
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


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_predictions(
    vcf_file=VCF_FILE,
    output_file=OUTPUT_FILENAME,
    tissue_id=TISSUE_ID,
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

    # --- Initialise model & VCF --------------------------------------------
    logging.info("Connecting to AlphaGenome…")
    model = dna_client.create(api_key)
    vcf = VCF(vcf_file)

    count_processed = 0
    count_saved = 0
    count_errors = 0
    start_time = time.time()

    with open(output_file, "w") as outfile:
        outfile.write(HEADER)

        for variant in vcf:
            if variant.FILTER is not None:
                continue

            # --- Prepare inputs --------------------------------------------
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF

            if not variant.ALT:
                logging.warning("Skipping %s:%d — no ALT allele", chrom, pos)
                continue
            alt = variant.ALT[0]

            gene_name = utils.get_gene_name(variant)
            gene_id = utils.get_gene_id(variant)  # version-stripped

            ag_variant = genome.Variant(chrom, pos, ref, alt)
            start = max(1, pos - dna_length // 2)
            end = pos + dna_length // 2
            interval = genome.Interval(chrom, start, end)

            count_processed += 1

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

                line = (
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_name}\t{gene_id}"
                    f"\t{ref_sum:.2f}\t{alt_sum:.2f}\t{log2_fc:.4f}\t{status}\n"
                )
                outfile.write(line)
                outfile.flush()  # persist each result immediately
                count_saved += 1

                logging.info(
                    "[%d/%d] %s:%d %s — %s (log2FC=%.4f)",
                    count_saved, count_processed, chrom, pos, gene_name,
                    status, log2_fc,
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
        "Done in %.1fs — processed %d variants, saved %d results, %d errors.",
        elapsed, count_processed, count_saved, count_errors,
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
    run_predictions()