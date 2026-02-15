"""Extract somatic variants from a TCGA paired tumour/normal VCF.

Reads a GATK MuTect2 VCF containing both NORMAL and TUMOR samples and
writes a tab-separated file listing every **somatic** variant — i.e. one
where the normal sample is homozygous reference (0/0) and the tumour
sample carries at least one alternate allele.

No quality or impact filtering is applied at this stage; the purpose is
purely to separate tumour-specific mutations from germline variants.
Down-stream filtering is handled by ``s2_vcf_filter.py``.

Output columns
--------------
CHROM, POS, REF, ALT, Normal_GT, Tumor_GT, Type, Change

Usage
-----
    python -m src.s1_initial_vcf_processing                    # defaults
    python -m src.s1_initial_vcf_processing --vcf path/to.vcf  # custom input
    python -m src.s1_initial_vcf_processing -o my_output.txt   # custom output
"""

from __future__ import annotations

import argparse
import logging

from cyvcf2 import VCF

from src import utils
from src.constants import DATA_PATH, LOG_DIR, OUTPUT_DIR

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_OUTPUT = str(OUTPUT_DIR / "somatic_variants.txt")
LOG_FILENAME = str(LOG_DIR / "initial_vcf_processing.log")

HEADER = "CHROM\tPOS\tREF\tALT\tNormal_GT\tTumor_GT\tType\tChange\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (defaults to ``sys.argv[1:]`` when *None*).

    Returns:
        Namespace with *vcf* and *output* attributes.
    """
    parser = argparse.ArgumentParser(
        description="Extract somatic variants from a paired tumour/normal VCF.",
    )
    parser.add_argument(
        "--vcf",
        default=DATA_PATH,
        help=f"Path to the input VCF file (default: {DATA_PATH}).",
    )
    parser.add_argument(
        "--output", "-o",
        default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT}).",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_genotype(alleles: list[int]) -> str:
    """Convert a cyvcf2 genotype list ``[a1, a2, phased]`` to ``'a1/a2'``."""
    return f"{alleles[0]}/{alleles[1]}"


def is_somatic(normal_gt: list[int], tumor_gt: list[int]) -> bool:
    """Return True if the variant is somatic.

    A variant is somatic when the normal sample is homozygous reference
    (both alleles == 0) and the tumour sample carries at least one
    alternate allele (> 0).

    Args:
        normal_gt: Genotype list ``[allele1, allele2, phased]`` for the
            NORMAL sample.
        tumor_gt: Genotype list ``[allele1, allele2, phased]`` for the
            TUMOR sample.
    """
    normal_is_ref = normal_gt[0] == 0 and normal_gt[1] == 0
    tumor_has_alt = tumor_gt[0] > 0 or tumor_gt[1] > 0
    return normal_is_ref and tumor_has_alt


def _find_sample_indices(vcf: VCF) -> tuple[int, int]:
    """Return ``(normal_idx, tumor_idx)`` from the VCF sample list.

    Raises:
        ValueError: If NORMAL or TUMOR sample names are not found.
    """
    samples = list(vcf.samples)
    try:
        normal_idx = samples.index("NORMAL")
    except ValueError as exc:
        raise ValueError(
            f"'NORMAL' sample not found in VCF samples: {samples}"
        ) from exc
    try:
        tumor_idx = samples.index("TUMOR")
    except ValueError as exc:
        raise ValueError(
            f"'TUMOR' sample not found in VCF samples: {samples}"
        ) from exc
    return normal_idx, tumor_idx


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def extract_somatic_variants(vcf_path: str, output_path: str) -> dict:
    """Read *vcf_path* and write somatic variants to *output_path*.

    Returns:
        A summary dict with keys ``total``, ``somatic``, and ``skipped``.
    """
    utils.validate_file(vcf_path, label="Input VCF")
    logging.info("Opening VCF: %s", vcf_path)

    vcf = VCF(vcf_path)
    normal_idx, tumor_idx = _find_sample_indices(vcf)
    logging.info(
        "Sample indices — NORMAL: %d, TUMOR: %d", normal_idx, tumor_idx,
    )

    total = 0
    somatic = 0
    skipped = 0

    with open(output_path, "w") as fh:
        fh.write(HEADER)

        for variant in vcf:
            total += 1

            normal_gt = variant.genotypes[normal_idx]
            tumor_gt = variant.genotypes[tumor_idx]

            if is_somatic(normal_gt, tumor_gt):
                somatic += 1
                line = (
                    f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t"
                    f"{','.join(variant.ALT)}\t"
                    f"{_format_genotype(normal_gt)}\t"
                    f"{_format_genotype(tumor_gt)}\n"
                )
                fh.write(line)
            else:
                skipped += 1

            if total % 5000 == 0:
                logging.info(
                    "Processed %d variants (%d somatic so far)…", total, somatic,
                )

    vcf.close()

    stats = {"total": total, "somatic": somatic, "skipped": skipped}
    _log_summary(stats, output_path)
    return stats


def _log_summary(stats: dict, output_path: str) -> None:
    """Write a human-readable summary to the log."""
    logging.info("=" * 60)
    logging.info("INITIAL VCF PROCESSING SUMMARY")
    logging.info("=" * 60)
    logging.info("Total variants scanned:   %d", stats["total"])
    logging.info("Somatic variants written:  %d", stats["somatic"])
    logging.info("Skipped (germline/other):  %d", stats["skipped"])
    logging.info("Output file: %s", output_path)
    logging.info("=" * 60)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    utils.setup_logging(LOG_FILENAME)
    logging.info("Config: %s", vars(args))
    extract_somatic_variants(args.vcf, args.output)
