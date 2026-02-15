"""AlphaGenome expression prediction — per-gene log₂ fold-change via GeneMaskLFCScorer.

Reads a filtered VCF (by default the output of ``s2_vcf_filter.py``), queries
the AlphaGenome API using ``score_variant`` with ``GeneMaskLFCScorer`` to
compute **per-gene log₂ fold-change** (masking to exon bins only), and writes
a raw predictions TSV filtered to the target gene from the VCF CSQ annotation.

This replaces the previous whole-window summing strategy which diluted the
target gene's signal across 21–99 neighbouring genes (see docs/GENE_DILUTION.md).

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
import pandas as pd
from cyvcf2 import VCF
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
from dotenv import load_dotenv

from src.constants import HIGH_IMPACT_VCF, RAW_PREDICTIONS, LOG_DIR, DOTENV_PATH
from src import utils
from src.exceptions import PipelineInputError

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = str(LOG_DIR / "gene_expression_prediction.log")
DEFAULT_VCF = HIGH_IMPACT_VCF
DEFAULT_OUTPUT = RAW_PREDICTIONS
DEFAULT_TISSUE = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

# Retry / rate-limit settings
MAX_RETRIES = 3
RETRY_BASE_DELAY = 2.0          # seconds; doubles each retry
RATE_LIMIT_DELAY = 0.5          # seconds between successive API calls

RAW_HEADER = "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tLOG2_FC\n"


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


def _score_with_retry(model, interval, ag_variant,
                      max_retries=MAX_RETRIES,
                      base_delay=RETRY_BASE_DELAY):
    """Call *model.score_variant* with ``GeneMaskLFCScorer`` and retry.

    Returns a tidy ``pd.DataFrame`` (via ``tidy_scores``) on success,
    or raises the last exception after all retries are exhausted.
    """
    scorer = variant_scorers.GeneMaskLFCScorer(
        requested_output=dna_client.OutputType.RNA_SEQ,
    )
    last_exc = None
    for attempt in range(1, max_retries + 1):
        try:
            scores = model.score_variant(
                interval=interval,
                variant=ag_variant,
                variant_scorers=[scorer],
            )
            return variant_scorers.tidy_scores(scores)
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
    load_dotenv(DOTENV_PATH)
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        raise PipelineInputError(
            "ALPHAGENOME_API_KEY not set in environment or .env"
        )

    utils.validate_file(vcf_file, "Input VCF")

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
                # --- Score with GeneMaskLFCScorer (with retry) -------------
                scores_df = _score_with_retry(
                    model, interval, ag_variant,
                )

                # --- Extract target gene LFC -------------------------------
                log2_fc = float("nan")
                if scores_df is not None and not scores_df.empty:
                    target = scores_df[
                        scores_df["gene_name"] == gene_name
                    ]
                    if target.empty and gene_id:
                        # Fallback: match on ENSEMBL ID prefix
                        target = scores_df[
                            scores_df["gene_id"].str.startswith(
                                gene_id.split(".")[0]
                            )
                        ]
                    if not target.empty:
                        log2_fc = float(target.iloc[0]["raw_score"])
                    else:
                        logging.warning(
                            "Gene %s (%s) not found in %d GeneMask results "
                            "at %s:%d",
                            gene_name, gene_id,
                            len(scores_df), chrom, pos,
                        )

                line = (
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_name}\t{gene_id}"
                    f"\t{log2_fc:.6f}\n"
                )
                outfile.write(line)
                outfile.flush()
                count_saved += 1

                logging.info(
                    "[%d/%d] %s:%d %s — log2_fc=%.6f",
                    count_saved, count_total, chrom, pos, gene_name,
                    log2_fc,
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
    try:
        run_predictions(
            vcf_file=args.vcf,
            output_file=args.output,
            tissue_id=args.tissue,
            resume=args.resume,
        )
    except PipelineInputError as exc:
        logging.critical("%s", exc)
        sys.exit(1)