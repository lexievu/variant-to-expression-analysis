"""Validate AlphaGenome gene-mask coordinates against UCSC Genome Browser.

Queries AlphaGenome's ``score_variant`` with ``GeneMaskLFCScorer`` and the
UCSC REST API (``wgEncodeGencodeBasicV47`` track on hg38) for the same
genomic regions, then writes a comparison TSV so you can verify that
AlphaGenome's internal GENCODE gene models match the public reference.

The output contains one row per gene × variant, with columns from both
sources side-by-side (gene boundaries, strand, gene type).

Usage
-----
    python -m src.a_gene_masking_validation                # all variants in high-impact VCF
    python -m src.a_gene_masking_validation --vcf my.vcf   # custom VCF
    python -m src.a_gene_masking_validation --skip-api      # UCSC-only (no AlphaGenome cost)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
import traceback
from urllib.error import URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from dotenv import load_dotenv

from src.constants import (
    DOTENV_PATH,
    HIGH_IMPACT_VCF,
    LOG_DIR,
    OUTPUT_DIR,
)
from src import utils
from src.exceptions import PipelineInputError

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LOG_FILENAME = str(LOG_DIR / "gene_masking_validation.log")
DEFAULT_VCF = HIGH_IMPACT_VCF
DEFAULT_OUTPUT = str(OUTPUT_DIR / "gene_masking_validation.tsv")
DEFAULT_TISSUE = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576
UCSC_API_BASE = "https://api.genome.ucsc.edu"
UCSC_TRACK = "wgEncodeGencodeBasicV47"  # GENCODE Basic v47 on hg38

# Rate-limit
UCSC_DELAY = 1.0          # seconds between UCSC API calls (per their policy)
AG_DELAY = 0.5            # seconds between AlphaGenome API calls
AG_MAX_RETRIES = 3
AG_RETRY_BASE_DELAY = 2.0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Validate AlphaGenome gene masks against UCSC GENCODE.",
    )
    parser.add_argument(
        "--vcf", default=DEFAULT_VCF,
        help=f"Input VCF of filtered variants (default: {DEFAULT_VCF}).",
    )
    parser.add_argument(
        "--output", "-o", default=DEFAULT_OUTPUT,
        help=f"Output comparison TSV (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--tissue", default=DEFAULT_TISSUE,
        help=f"UBERON ontology ID for tissue type (default: {DEFAULT_TISSUE}).",
    )
    parser.add_argument(
        "--skip-api", action="store_true",
        help="Skip AlphaGenome API calls (produce UCSC-only output).",
    )
    return parser.parse_args(argv)


# ===================================================================
# UCSC Genome Browser helpers
# ===================================================================

def fetch_ucsc_genes(chrom: str, start: int, end: int,
                     track: str = UCSC_TRACK) -> list[dict]:
    """Fetch gene annotations from the UCSC REST API for a genomic region.

    Uses GENCODE Basic (``wgEncodeGencodeBasicV47``) which provides
    Ensembl transcript IDs, gene symbols, exon coordinates, and gene type.

    Args:
        chrom: Chromosome name (e.g. ``chr17``).
        start: 0-based start coordinate.
        end: 1-based end coordinate.
        track: UCSC track name.

    Returns:
        List of dicts, one per transcript, with keys from the UCSC track
        schema (``name``, ``chrom``, ``strand``, ``txStart``, ``txEnd``,
        ``cdsStart``, ``cdsEnd``, ``exonStarts``, ``exonEnds``, ``name2``,
        ``geneType``, etc.).
    """
    url = (
        f"{UCSC_API_BASE}/getData/track"
        f"?genome=hg38;track={track};chrom={chrom}"
        f";start={start};end={end}"
    )
    logging.info("UCSC query: %s", url)

    last_exc = None
    for attempt in range(1, 4):
        try:
            with urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            if "statusCode" in data and data["statusCode"] >= 400:
                logging.warning("UCSC API error: %s", data.get("error", ""))
                return []
            return data.get(track, [])
        except (URLError, OSError, json.JSONDecodeError) as exc:
            last_exc = exc
            if attempt < 3:
                delay = 2 ** attempt
                logging.warning(
                    "UCSC attempt %d/3 failed (%s). Retrying in %ds…",
                    attempt, exc, delay,
                )
                time.sleep(delay)
    logging.error("All UCSC retries exhausted: %s", last_exc)
    return []


def parse_exon_coords(transcript: dict) -> list[tuple[int, int]]:
    """Extract (start, end) tuples for every exon in a UCSC transcript.

    ``exonStarts`` / ``exonEnds`` are comma-separated strings of 0-based
    coordinates from the GENCODE genePred format.
    """
    starts = [int(s) for s in transcript["exonStarts"].rstrip(",").split(",") if s]
    ends = [int(s) for s in transcript["exonEnds"].rstrip(",").split(",") if s]
    return list(zip(starts, ends))


def ucsc_genes_to_records(transcripts: list[dict],
                          target_gene_name: str | None = None,
                          ) -> list[dict]:
    """Collapse UCSC transcripts into per-gene summaries.

    Groups by ``name2`` (gene symbol) and computes the gene-level
    envelope (min txStart, max txEnd) and the union of all CDS exon
    coordinates.

    If *target_gene_name* is given, only that gene is returned.
    """
    from collections import defaultdict
    by_gene: dict[str, dict] = defaultdict(lambda: {
        "transcripts": [],
        "tx_starts": [],
        "tx_ends": [],
        "cds_starts": [],
        "cds_ends": [],
        "strands": set(),
        "gene_types": set(),
        "all_exons": [],
    })

    for tx in transcripts:
        gene_name = tx.get("name2", "")
        if target_gene_name and gene_name != target_gene_name:
            continue
        g = by_gene[gene_name]
        g["transcripts"].append(tx["name"])  # ENST ID
        g["tx_starts"].append(tx["txStart"])
        g["tx_ends"].append(tx["txEnd"])
        g["cds_starts"].append(tx["cdsStart"])
        g["cds_ends"].append(tx["cdsEnd"])
        g["strands"].add(tx["strand"])
        g["gene_types"].add(tx.get("geneType", ""))
        g["all_exons"].extend(parse_exon_coords(tx))

    records = []
    for gene_name, g in by_gene.items():
        # Union of exon intervals (sorted, not merged — for transparency)
        exons_sorted = sorted(set(g["all_exons"]))
        records.append({
            "ucsc_gene_name": gene_name,
            "ucsc_strand": ",".join(sorted(g["strands"])),
            "ucsc_gene_type": ",".join(sorted(g["gene_types"])),
            "ucsc_tx_start": min(g["tx_starts"]),
            "ucsc_tx_end": max(g["tx_ends"]),
            "ucsc_cds_start": min(g["cds_starts"]),
            "ucsc_cds_end": max(g["cds_ends"]),
            "ucsc_n_transcripts": len(g["transcripts"]),
            "ucsc_n_exons_union": len(exons_sorted),
            "ucsc_transcript_ids": ";".join(g["transcripts"]),
            "ucsc_exon_starts": ";".join(str(s) for s, _ in exons_sorted),
            "ucsc_exon_ends": ";".join(str(e) for _, e in exons_sorted),
        })
    return records


# ===================================================================
# AlphaGenome helpers
# ===================================================================

def query_alphagenome(model, chrom, pos, ref, alt, tissue_id,
                      dna_length=DNA_SEQUENCE_LENGTH):
    """Run ``score_variant`` with ``GeneMaskLFCScorer`` and return per-gene metadata.

    Returns a list of dicts with the gene-level information AlphaGenome
    uses internally for gene-mask scoring, plus the raw LFC score.
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    scorer = variant_scorers.GeneMaskLFCScorer(
        requested_output=dna_client.OutputType.RNA_SEQ,
    )

    ag_variant = genome.Variant(chrom, pos, ref, alt)
    start = max(1, pos - dna_length // 2)
    end = pos + dna_length // 2
    interval = genome.Interval(chrom, start, end)

    last_exc = None
    for attempt in range(1, AG_MAX_RETRIES + 1):
        try:
            scores = model.score_variant(
                interval=interval,
                variant=ag_variant,
                variant_scorers=[scorer],
            )
            break
        except Exception as exc:
            last_exc = exc
            if attempt < AG_MAX_RETRIES:
                delay = AG_RETRY_BASE_DELAY * (2 ** (attempt - 1))
                logging.warning(
                    "AlphaGenome attempt %d/%d failed (%s). Retry in %.1fs…",
                    attempt, AG_MAX_RETRIES, exc, delay,
                )
                time.sleep(delay)
    else:
        raise last_exc  # type: ignore[misc]

    # Convert AnnData to a tidy DataFrame
    df = variant_scorers.scores_to_dataframe(scores, match_gene_strand=False)

    records = []
    if df is not None and not df.empty:
        for _, row in df.iterrows():
            records.append({
                "ag_gene_id": row.get("gene_id", ""),
                "ag_gene_name": row.get("gene_name", ""),
                "ag_gene_type": row.get("gene_type", ""),
                "ag_gene_strand": row.get("gene_strand", ""),
                "ag_track_strand": row.get("track_strand", ""),
                "ag_lfc_score": row.get("raw_score", np.nan),
                "ag_scored_interval": row.get("scored_interval", ""),
            })
    return records


# ===================================================================
# Main comparison pipeline
# ===================================================================

def run_validation(
    vcf_file=DEFAULT_VCF,
    output_file=DEFAULT_OUTPUT,
    tissue_id=DEFAULT_TISSUE,
    skip_api=False,
    dna_length=DNA_SEQUENCE_LENGTH,
):
    """Compare AlphaGenome gene masks to UCSC GENCODE for every variant."""

    load_dotenv(DOTENV_PATH)
    utils.validate_file(vcf_file, "Input VCF")

    model = None
    if not skip_api:
        api_key = os.getenv("ALPHAGENOME_API_KEY")
        if not api_key:
            raise PipelineInputError(
                "ALPHAGENOME_API_KEY not set. Use --skip-api for UCSC-only mode."
            )
        from alphagenome.models import dna_client
        logging.info("Connecting to AlphaGenome…")
        model = dna_client.create(api_key)

    vcf = VCF(vcf_file)
    all_rows: list[dict] = []
    variant_count = 0

    for variant in vcf:
        if variant.FILTER is not None:
            continue

        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        if not variant.ALT:
            continue
        alt = variant.ALT[0]
        gene_name = utils.get_gene_name(variant)
        gene_id = utils.get_gene_id(variant)
        variant_count += 1

        logging.info(
            "[%d] Processing %s:%d %s>%s  gene=%s (%s)",
            variant_count, chrom, pos, ref, alt, gene_name, gene_id,
        )

        # ---- UCSC lookup --------------------------------------------------
        window_start = max(0, pos - dna_length // 2)
        window_end = pos + dna_length // 2

        ucsc_transcripts = fetch_ucsc_genes(chrom, window_start, window_end)
        time.sleep(UCSC_DELAY)

        # Get summaries for the target gene and all genes in the window
        target_ucsc = ucsc_genes_to_records(ucsc_transcripts, target_gene_name=gene_name)
        all_ucsc = ucsc_genes_to_records(ucsc_transcripts)

        # ---- AlphaGenome lookup -------------------------------------------
        ag_records: list[dict] = []
        if model is not None:
            try:
                ag_records = query_alphagenome(
                    model, chrom, pos, ref, alt, tissue_id, dna_length,
                )
                time.sleep(AG_DELAY)
            except Exception:
                logging.error(
                    "AlphaGenome failed for %s:%d. Traceback:\n%s",
                    chrom, pos, traceback.format_exc(),
                )

        # ---- Build comparison rows ----------------------------------------
        variant_info = {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "vcf_gene_name": gene_name,
            "vcf_gene_id": gene_id,
            "ucsc_genes_in_window": len(all_ucsc),
        }

        if ag_records:
            # One row per AlphaGenome gene, left-joined with UCSC target gene
            for ag in ag_records:
                row = {**variant_info, **ag}
                # Try to match this AG gene to the UCSC data by gene name
                matching_ucsc = [
                    u for u in all_ucsc
                    if u["ucsc_gene_name"] == ag.get("ag_gene_name")
                ]
                if matching_ucsc:
                    row.update(matching_ucsc[0])
                else:
                    # Fill UCSC columns with empty
                    for key in _ucsc_empty_record():
                        row.setdefault(key, "")
                all_rows.append(row)
        elif target_ucsc:
            # No AlphaGenome data — output UCSC only
            for u in target_ucsc:
                row = {**variant_info}
                row.update(u)
                for key in _ag_empty_record():
                    row.setdefault(key, "")
                all_rows.append(row)
        else:
            # Neither source had data
            row = {**variant_info}
            row.update(_ucsc_empty_record())
            row.update(_ag_empty_record())
            all_rows.append(row)

    vcf.close()

    # ---- Write output -----------------------------------------------------
    if not all_rows:
        logging.warning("No variants processed — nothing to write.")
        return

    df = pd.DataFrame(all_rows)
    # Reorder columns for readability
    col_order = _column_order()
    present = [c for c in col_order if c in df.columns]
    extra = [c for c in df.columns if c not in col_order]
    df = df[present + extra]
    df.to_csv(output_file, sep="\t", index=False)
    logging.info("Wrote %d rows to %s", len(df), output_file)

    # ---- Print summary to console -----------------------------------------
    _print_summary(df)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ucsc_empty_record() -> dict:
    return {
        "ucsc_gene_name": "", "ucsc_strand": "", "ucsc_gene_type": "",
        "ucsc_tx_start": "", "ucsc_tx_end": "",
        "ucsc_cds_start": "", "ucsc_cds_end": "",
        "ucsc_n_transcripts": "", "ucsc_n_exons_union": "",
        "ucsc_transcript_ids": "", "ucsc_exon_starts": "", "ucsc_exon_ends": "",
    }


def _ag_empty_record() -> dict:
    return {
        "ag_gene_id": "", "ag_gene_name": "", "ag_gene_type": "",
        "ag_gene_strand": "", "ag_track_strand": "",
        "ag_lfc_score": "", "ag_scored_interval": "",
    }


def _column_order() -> list[str]:
    """Preferred column ordering for the output TSV."""
    return [
        # Variant identity
        "chrom", "pos", "ref", "alt",
        "vcf_gene_name", "vcf_gene_id",
        # Summary
        "ucsc_genes_in_window",
        # AlphaGenome gene-mask output
        "ag_gene_id", "ag_gene_name", "ag_gene_type",
        "ag_gene_strand", "ag_track_strand", "ag_lfc_score",
        "ag_scored_interval",
        # UCSC GENCODE reference
        "ucsc_gene_name", "ucsc_strand", "ucsc_gene_type",
        "ucsc_tx_start", "ucsc_tx_end",
        "ucsc_cds_start", "ucsc_cds_end",
        "ucsc_n_transcripts", "ucsc_n_exons_union",
        "ucsc_transcript_ids",
        "ucsc_exon_starts", "ucsc_exon_ends",
    ]


def _print_summary(df: pd.DataFrame) -> None:
    """Log a human-readable comparison summary."""
    logging.info("=" * 72)
    logging.info("GENE MASKING VALIDATION SUMMARY")
    logging.info("=" * 72)

    n_variants = df[["chrom", "pos"]].drop_duplicates().shape[0]
    logging.info("Variants processed:  %d", n_variants)

    if "ag_gene_name" in df.columns:
        ag_genes = df[df["ag_gene_name"].astype(str).ne("")]["ag_gene_name"].nunique()
        logging.info("AlphaGenome genes:   %d", ag_genes)

    if "ucsc_gene_name" in df.columns:
        ucsc_genes = df[df["ucsc_gene_name"].astype(str).ne("")]["ucsc_gene_name"].nunique()
        logging.info("UCSC GENCODE genes:  %d", ucsc_genes)

    # Check gene-name agreement
    if "ag_gene_name" in df.columns and "ucsc_gene_name" in df.columns:
        matched = df[
            (df["ag_gene_name"].astype(str).ne(""))
            & (df["ucsc_gene_name"].astype(str).ne(""))
        ]
        if not matched.empty:
            agree = (matched["ag_gene_name"] == matched["ucsc_gene_name"]).sum()
            logging.info("Name matches:        %d/%d", agree, len(matched))

    # Check strand agreement
    if "ag_gene_strand" in df.columns and "ucsc_strand" in df.columns:
        both = df[
            (df["ag_gene_strand"].astype(str).ne(""))
            & (df["ucsc_strand"].astype(str).ne(""))
        ]
        if not both.empty:
            strand_agree = (both["ag_gene_strand"] == both["ucsc_strand"]).sum()
            logging.info("Strand matches:      %d/%d", strand_agree, len(both))

    # Neighbor gene dilution — how many genes per window?
    if "ucsc_genes_in_window" in df.columns:
        per_variant = df.groupby(["chrom", "pos"])["ucsc_genes_in_window"].first()
        logging.info("Genes per 1 MB window (UCSC):")
        logging.info(
            "  min=%d  median=%.0f  max=%d",
            per_variant.min(), per_variant.median(), per_variant.max(),
        )
        logging.info("  → This is why whole-window summing dilutes the target gene signal.")

    logging.info("=" * 72)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    try:
        run_validation(
            vcf_file=args.vcf,
            output_file=args.output,
            tissue_id=args.tissue,
            skip_api=args.skip_api,
        )
    except PipelineInputError as exc:
        logging.critical("%s", exc)
        sys.exit(1)
