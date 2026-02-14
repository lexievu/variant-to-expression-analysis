"""Indirect comparison of AlphaGenome gene masks vs UCSC exon coordinates.

Performs three independent data-collection steps (each cached to disk so
that API calls are never repeated) and one offline comparison step:

1. **AlphaGenome raw tracks** — call ``predict_variant`` and save per-variant
   REF/ALT RNA-seq arrays to ``output/ag_raw_tracks/<key>.npz``.
2. **AlphaGenome gene-mask scores** — call ``score_variant`` with
   ``GeneMaskLFCScorer`` and save to ``output/ag_genemask_scores.tsv``.
3. **Offline comparison** — apply the UCSC exon mask (from step 1 of
   ``a_gene_masking_validation.py``) to the saved raw tracks, compute a
   manual LFC, and compare to the ``GeneMaskLFCScorer`` score.

Requires ``output/gene_masking_validation.tsv`` from
``a_gene_masking_validation.py --skip-api`` (UCSC exon data).

Usage
-----
    python -m src.b_gene_masking_comparison                  # full run (API)
    python -m src.b_gene_masking_comparison --skip-api       # offline only
    python -m src.b_gene_masking_comparison --compare-only   # step 3 only
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
import traceback
from pathlib import Path

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
LOG_FILENAME = str(LOG_DIR / "gene_masking_comparison.log")
DEFAULT_VCF = HIGH_IMPACT_VCF
DEFAULT_UCSC_TSV = str(OUTPUT_DIR / "gene_masking_validation.tsv")
DEFAULT_GENEMASK_TSV = str(OUTPUT_DIR / "ag_genemask_scores.tsv")
DEFAULT_COMPARISON_TSV = str(OUTPUT_DIR / "gene_masking_comparison.tsv")
RAW_TRACKS_DIR = OUTPUT_DIR / "ag_raw_tracks"
DEFAULT_TISSUE = "UBERON:0002048"  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

# Rate-limit
AG_DELAY = 0.5            # seconds between AlphaGenome API calls
AG_MAX_RETRIES = 3
AG_RETRY_BASE_DELAY = 2.0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Compare AlphaGenome GeneMaskLFC to manual UCSC-exon "
                    "masking of raw prediction tracks.",
    )
    parser.add_argument(
        "--vcf", default=DEFAULT_VCF,
        help=f"Input VCF of filtered variants (default: {DEFAULT_VCF}).",
    )
    parser.add_argument(
        "--ucsc-tsv", default=DEFAULT_UCSC_TSV,
        help=f"UCSC validation TSV from a_gene_masking_validation.py "
             f"(default: {DEFAULT_UCSC_TSV}).",
    )
    parser.add_argument(
        "--genemask-output", default=DEFAULT_GENEMASK_TSV,
        help=f"GeneMaskLFCScorer cached TSV (default: {DEFAULT_GENEMASK_TSV}).",
    )
    parser.add_argument(
        "--comparison-output", default=DEFAULT_COMPARISON_TSV,
        help=f"Final comparison TSV (default: {DEFAULT_COMPARISON_TSV}).",
    )
    parser.add_argument(
        "--tissue", default=DEFAULT_TISSUE,
        help=f"UBERON ontology ID for tissue type (default: {DEFAULT_TISSUE}).",
    )
    parser.add_argument(
        "--skip-api", action="store_true",
        help="Skip AlphaGenome API calls (steps 1 & 2). "
             "Run comparison from whatever is already cached.",
    )
    parser.add_argument(
        "--compare-only", action="store_true",
        help="Skip all API calls; run only offline comparison (step 3) "
             "from previously cached files.",
    )
    return parser.parse_args(argv)


# ===================================================================
# Variant key helpers
# ===================================================================

def _variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Return a filesystem-safe unique key for a variant."""
    return f"{chrom}_{pos}_{ref}_{alt}"


# ===================================================================
# AlphaGenome retry helper
# ===================================================================

def _ag_retry(fn, *args, **kwargs):
    """Call *fn* with exponential-backoff retry."""
    last_exc = None
    for attempt in range(1, AG_MAX_RETRIES + 1):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            last_exc = exc
            if attempt < AG_MAX_RETRIES:
                delay = AG_RETRY_BASE_DELAY * (2 ** (attempt - 1))
                logging.warning(
                    "AlphaGenome attempt %d/%d failed (%s). Retry in %.1fs…",
                    attempt, AG_MAX_RETRIES, exc, delay,
                )
                time.sleep(delay)
    raise last_exc  # type: ignore[misc]


# ===================================================================
# Step 1 — predict_variant → raw tracks to disk
# ===================================================================

def save_raw_tracks(model, chrom, pos, ref, alt, tissue_id,
                    dna_length=DNA_SEQUENCE_LENGTH,
                    out_dir=RAW_TRACKS_DIR) -> Path:
    """Call ``predict_variant`` and cache REF/ALT RNA-seq arrays to disk.

    Saves a ``.npz`` file containing:
        * ``ref_values`` — numpy array of shape ``(bins, tracks)``
        * ``alt_values`` — numpy array of shape ``(bins, tracks)``
        * ``resolution`` — int, bin size in base pairs
        * ``interval_start``, ``interval_end`` — genomic coords of the window

    Returns the path to the saved ``.npz`` file.
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    ag_variant = genome.Variant(chrom, pos, ref, alt)
    start = max(1, pos - dna_length // 2)
    end = pos + dna_length // 2
    interval = genome.Interval(chrom, start, end)

    outputs = _ag_retry(
        model.predict_variant,
        interval=interval,
        variant=ag_variant,
        ontology_terms=[tissue_id],
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
    )

    ref_td = outputs.reference.rna_seq
    alt_td = outputs.alternate.rna_seq

    os.makedirs(out_dir, exist_ok=True)
    key = _variant_key(chrom, pos, ref, alt)
    npz_path = Path(out_dir) / f"{key}.npz"
    np.savez_compressed(
        npz_path,
        ref_values=ref_td.values,
        alt_values=alt_td.values,
        resolution=np.array([ref_td.resolution]),
        interval_start=np.array([start]),
        interval_end=np.array([end]),
    )
    logging.info("Saved raw tracks to %s  (shape=%s, resolution=%d)",
                 npz_path, ref_td.values.shape, ref_td.resolution)
    return npz_path


def load_raw_tracks(chrom, pos, ref, alt, tracks_dir=RAW_TRACKS_DIR):
    """Load cached raw tracks from disk.  Returns dict or None."""
    key = _variant_key(chrom, pos, ref, alt)
    npz_path = Path(tracks_dir) / f"{key}.npz"
    if not npz_path.is_file():
        return None
    data = np.load(npz_path)
    return {
        "ref_values": data["ref_values"],
        "alt_values": data["alt_values"],
        "resolution": int(data["resolution"][0]),
        "interval_start": int(data["interval_start"][0]),
        "interval_end": int(data["interval_end"][0]),
    }


# ===================================================================
# Step 2 — score_variant (GeneMaskLFCScorer) → TSV
# ===================================================================

def query_genemask_scores(model, chrom, pos, ref, alt, tissue_id,
                          dna_length=DNA_SEQUENCE_LENGTH) -> list[dict]:
    """Call ``score_variant`` with ``GeneMaskLFCScorer`` and return records."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    scorer = variant_scorers.GeneMaskLFCScorer(
        requested_output=dna_client.OutputType.RNA_SEQ,
    )

    ag_variant = genome.Variant(chrom, pos, ref, alt)
    start = max(1, pos - dna_length // 2)
    end = pos + dna_length // 2
    interval = genome.Interval(chrom, start, end)

    scores = _ag_retry(
        model.score_variant,
        interval=interval,
        variant=ag_variant,
        variant_scorers=[scorer],
    )

    df = variant_scorers.tidy_scores(scores)
    records = []
    if df is not None and not df.empty:
        for _, row in df.iterrows():
            records.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
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
# Steps 1 & 2 — AlphaGenome API calls (with caching)
# ===================================================================

def run_alphagenome_calls(vcf_file, genemask_output, tissue_id,
                          dna_length=DNA_SEQUENCE_LENGTH):
    """Call predict_variant + score_variant for each variant, saving to disk.

    * Raw tracks are saved per-variant as ``.npz`` files in
      ``output/ag_raw_tracks/``.
    * GeneMaskLFCScorer results are appended to *genemask_output* TSV.

    Already-cached variants (existing ``.npz`` / TSV rows) are skipped.
    """
    load_dotenv(DOTENV_PATH)
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        raise PipelineInputError(
            "ALPHAGENOME_API_KEY not set. Use --skip-api for offline mode."
        )
    utils.validate_file(vcf_file, "Input VCF")

    from alphagenome.models import dna_client
    logging.info("Connecting to AlphaGenome…")
    model = dna_client.create(api_key)

    # Load already-saved genemask rows so we can skip them
    existing_gm_keys: set[str] = set()
    if os.path.isfile(genemask_output):
        try:
            prev = pd.read_csv(genemask_output, sep="\t")
            for _, r in prev.iterrows():
                existing_gm_keys.add(
                    _variant_key(r["chrom"], int(r["pos"]),
                                 r["ref"], r["alt"]))
        except Exception:
            pass

    vcf = VCF(vcf_file)
    genemask_rows: list[dict] = []
    variant_count = 0

    for variant in vcf:
        if variant.FILTER is not None:
            continue
        chrom, pos, ref = variant.CHROM, variant.POS, variant.REF
        if not variant.ALT:
            continue
        alt = variant.ALT[0]
        variant_count += 1
        key = _variant_key(chrom, pos, ref, alt)

        # --- Step 1: raw tracks -------------------------------------------
        cached = load_raw_tracks(chrom, pos, ref, alt)
        if cached is not None:
            logging.info("[%d] Raw tracks already cached: %s", variant_count, key)
        else:
            logging.info("[%d] predict_variant  %s:%d %s>%s",
                         variant_count, chrom, pos, ref, alt)
            try:
                save_raw_tracks(model, chrom, pos, ref, alt, tissue_id,
                                dna_length)
                time.sleep(AG_DELAY)
            except Exception:
                logging.error("predict_variant failed for %s:\n%s",
                              key, traceback.format_exc())

        # --- Step 2: GeneMask scores --------------------------------------
        if key in existing_gm_keys:
            logging.info("[%d] GeneMask scores already cached: %s",
                         variant_count, key)
        else:
            logging.info("[%d] score_variant (GeneMask)  %s:%d %s>%s",
                         variant_count, chrom, pos, ref, alt)
            try:
                records = query_genemask_scores(model, chrom, pos, ref, alt,
                                                tissue_id, dna_length)
                genemask_rows.extend(records)
                time.sleep(AG_DELAY)
            except Exception:
                logging.error("score_variant failed for %s:\n%s",
                              key, traceback.format_exc())

    vcf.close()

    # Append new genemask rows
    if genemask_rows:
        df_new = pd.DataFrame(genemask_rows)
        write_header = not os.path.isfile(genemask_output)
        df_new.to_csv(genemask_output, sep="\t", index=False,
                      mode="a", header=write_header)
        logging.info("Steps 1–2 done — appended %d GeneMask rows to %s",
                     len(df_new), genemask_output)
    else:
        logging.info("Steps 1–2 done — no new variants to process.")


# ===================================================================
# Step 3 — Offline comparison (no API calls)
# ===================================================================

def build_exon_mask(exon_starts_str: str, exon_ends_str: str,
                    interval_start: int, interval_end: int,
                    resolution: int) -> np.ndarray:
    """Build a boolean bin mask from UCSC exon coordinates.

    Args:
        exon_starts_str: Semicolon-separated 0-based exon start positions.
        exon_ends_str:   Semicolon-separated 0-based exon end positions.
        interval_start:  Genomic start of the prediction window.
        interval_end:    Genomic end of the prediction window.
        resolution:      Bin size in base pairs.

    Returns:
        1-D boolean array of length ``(interval_end - interval_start) //
        resolution`` where True = bin overlaps at least one exon.
    """
    n_bins = (interval_end - interval_start) // resolution
    mask = np.zeros(n_bins, dtype=bool)

    if not exon_starts_str or not exon_ends_str:
        return mask

    starts = [int(s) for s in exon_starts_str.split(";") if s]
    ends = [int(e) for e in exon_ends_str.split(";") if e]

    for exon_start, exon_end in zip(starts, ends):
        # Convert genomic coords → bin indices (clip to window)
        rel_start = max(0, exon_start - interval_start)
        rel_end = min(interval_end - interval_start, exon_end - interval_start)
        if rel_end <= rel_start:
            continue
        bin_start = rel_start // resolution
        bin_end = (rel_end + resolution - 1) // resolution  # round up
        bin_end = min(bin_end, n_bins)
        mask[bin_start:bin_end] = True

    return mask


def compute_manual_lfc(ref_values: np.ndarray, alt_values: np.ndarray,
                       mask: np.ndarray) -> float:
    """Compute log₂(sum(ALT_masked)) − log₂(sum(REF_masked)).

    Matches the aggregation strategy of ``GeneMaskLFCScorer``
    (``DIFF_LOG2_SUM``).  If either masked sum is ≤ 0 or mask is empty,
    returns NaN.
    """
    if not mask.any():
        return float("nan")

    # Sum across masked bins and all tracks (to match GeneMask behaviour)
    ref_masked = ref_values[mask].sum()
    alt_masked = alt_values[mask].sum()

    if ref_masked <= 0 or alt_masked <= 0:
        return float("nan")

    return float(np.log2(alt_masked) - np.log2(ref_masked))


def run_comparison(ucsc_tsv, genemask_tsv, comparison_output,
                   tracks_dir=RAW_TRACKS_DIR):
    """Merge UCSC exon masks with raw tracks; compare to GeneMask scores.

    Reads cached files only — no API calls.
    """
    if not os.path.isfile(ucsc_tsv):
        logging.error("UCSC TSV not found: %s", ucsc_tsv)
        logging.error("Run  python -m src.a_gene_masking_validation --skip-api  first.")
        return

    ucsc_df = pd.read_csv(ucsc_tsv, sep="\t")
    logging.info("Loaded %d UCSC rows from %s", len(ucsc_df), ucsc_tsv)

    # Load GeneMask scores if available
    gm_df = None
    if os.path.isfile(genemask_tsv):
        gm_df = pd.read_csv(genemask_tsv, sep="\t")
        logging.info("Loaded %d GeneMask rows from %s", len(gm_df), genemask_tsv)
    else:
        logging.warning("GeneMask TSV not found: %s  (AG columns will be empty)",
                        genemask_tsv)

    rows: list[dict] = []

    for _, ucsc_row in ucsc_df.iterrows():
        chrom = ucsc_row["chrom"]
        pos = int(ucsc_row["pos"])
        ref = ucsc_row["ref"]
        alt = ucsc_row["alt"]
        gene_name = ucsc_row["vcf_gene_name"]
        gene_id = ucsc_row["vcf_gene_id"]

        base: dict = {
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "vcf_gene_name": gene_name, "vcf_gene_id": gene_id,
            "ucsc_genes_in_window": ucsc_row.get("ucsc_genes_in_window", ""),
            "ucsc_gene_name": ucsc_row.get("ucsc_gene_name", ""),
            "ucsc_strand": ucsc_row.get("ucsc_strand", ""),
            "ucsc_gene_type": ucsc_row.get("ucsc_gene_type", ""),
            "ucsc_tx_start": ucsc_row.get("ucsc_tx_start", ""),
            "ucsc_tx_end": ucsc_row.get("ucsc_tx_end", ""),
            "ucsc_n_exons_union": ucsc_row.get("ucsc_n_exons_union", ""),
        }

        # --- Manual mask LFC from raw tracks --------------------------------
        tracks = load_raw_tracks(chrom, pos, ref, alt, tracks_dir)
        if tracks is not None:
            mask = build_exon_mask(
                str(ucsc_row.get("ucsc_exon_starts", "")),
                str(ucsc_row.get("ucsc_exon_ends", "")),
                tracks["interval_start"],
                tracks["interval_end"],
                tracks["resolution"],
            )
            manual_lfc = compute_manual_lfc(
                tracks["ref_values"], tracks["alt_values"], mask,
            )

            # Naive whole-window LFC for comparison
            ref_total = tracks["ref_values"].sum()
            alt_total = tracks["alt_values"].sum()
            if ref_total > 0 and alt_total > 0:
                whole_window_lfc = float(np.log2(alt_total) - np.log2(ref_total))
            else:
                whole_window_lfc = float("nan")

            base["manual_ucsc_lfc"] = manual_lfc
            base["whole_window_lfc"] = whole_window_lfc
            base["n_masked_bins"] = int(mask.sum())
            base["n_total_bins"] = len(mask)
            base["pct_bins_masked"] = f"{100 * mask.sum() / len(mask):.2f}"
            base["resolution"] = tracks["resolution"]
        else:
            base["manual_ucsc_lfc"] = ""
            base["whole_window_lfc"] = ""
            base["n_masked_bins"] = ""
            base["n_total_bins"] = ""
            base["pct_bins_masked"] = ""
            base["resolution"] = ""

        # --- GeneMask LFC from score_variant --------------------------------
        if gm_df is not None and not gm_df.empty:
            match = gm_df[
                (gm_df["chrom"] == chrom)
                & (gm_df["pos"] == pos)
                & (gm_df["ref"] == ref)
                & (gm_df["alt"] == alt)
            ]
            target_match = match[match["ag_gene_name"] == gene_name]
            if not target_match.empty:
                row_gm = target_match.iloc[0]
                base["ag_gene_id"] = row_gm.get("ag_gene_id", "")
                base["ag_gene_name"] = row_gm.get("ag_gene_name", "")
                base["ag_gene_strand"] = row_gm.get("ag_gene_strand", "")
                base["ag_gene_type"] = row_gm.get("ag_gene_type", "")
                base["ag_lfc_score"] = row_gm.get("ag_lfc_score", "")
            else:
                for k in ("ag_gene_id", "ag_gene_name", "ag_gene_strand",
                          "ag_gene_type", "ag_lfc_score"):
                    base[k] = ""
        else:
            for k in ("ag_gene_id", "ag_gene_name", "ag_gene_strand",
                      "ag_gene_type", "ag_lfc_score"):
                base[k] = ""

        # --- LFC difference (if both available) -----------------------------
        ag_val = pd.to_numeric(base.get("ag_lfc_score"), errors="coerce")
        manual_val = pd.to_numeric(base.get("manual_ucsc_lfc"), errors="coerce")
        if pd.notna(ag_val) and pd.notna(manual_val):
            base["lfc_diff"] = float(ag_val - manual_val)
        else:
            base["lfc_diff"] = ""

        rows.append(base)

    if not rows:
        logging.warning("No comparison rows produced.")
        return

    df = pd.DataFrame(rows)

    # Reorder columns for readability
    col_order = [
        "chrom", "pos", "ref", "alt",
        "vcf_gene_name", "vcf_gene_id",
        "ucsc_genes_in_window",
        # Gene metadata comparison
        "ucsc_gene_name", "ucsc_strand", "ucsc_gene_type",
        "ag_gene_id", "ag_gene_name", "ag_gene_strand", "ag_gene_type",
        # Score comparison — the key columns
        "ag_lfc_score", "manual_ucsc_lfc", "lfc_diff",
        "whole_window_lfc",
        # Mask details
        "ucsc_tx_start", "ucsc_tx_end", "ucsc_n_exons_union",
        "n_masked_bins", "n_total_bins", "pct_bins_masked", "resolution",
    ]
    present = [c for c in col_order if c in df.columns]
    extra = [c for c in df.columns if c not in col_order]
    df = df[present + extra]

    df.to_csv(comparison_output, sep="\t", index=False)
    logging.info("Step 3 done — wrote %d rows to %s", len(df), comparison_output)
    _print_summary(df)


# ---------------------------------------------------------------------------
# Summary printing
# ---------------------------------------------------------------------------

def _print_summary(df: pd.DataFrame) -> None:
    """Log a human-readable comparison summary."""
    logging.info("=" * 72)
    logging.info("GENE MASKING COMPARISON SUMMARY")
    logging.info("=" * 72)

    n_variants = df[["chrom", "pos"]].drop_duplicates().shape[0]
    logging.info("Variants processed:  %d", n_variants)

    # --- Gene metadata agreement ---
    if "ag_gene_name" in df.columns and "ucsc_gene_name" in df.columns:
        both = df[
            df["ag_gene_name"].astype(str).ne("")
            & df["ucsc_gene_name"].astype(str).ne("")
        ]
        if not both.empty:
            name_agree = (both["ag_gene_name"] == both["ucsc_gene_name"]).sum()
            logging.info("Gene name matches:   %d/%d", name_agree, len(both))

    if "ag_gene_strand" in df.columns and "ucsc_strand" in df.columns:
        both = df[
            df["ag_gene_strand"].astype(str).ne("")
            & df["ucsc_strand"].astype(str).ne("")
        ]
        if not both.empty:
            strand_agree = (both["ag_gene_strand"] == both["ucsc_strand"]).sum()
            logging.info("Strand matches:      %d/%d", strand_agree, len(both))

    # --- Score comparison (AG GeneMask vs manual UCSC mask) ---
    has_ag = "ag_lfc_score" in df.columns
    has_manual = "manual_ucsc_lfc" in df.columns

    if has_ag and has_manual:
        score_df = df[
            pd.to_numeric(df["ag_lfc_score"], errors="coerce").notna()
            & pd.to_numeric(df["manual_ucsc_lfc"], errors="coerce").notna()
        ].copy()
        if not score_df.empty:
            ag = pd.to_numeric(score_df["ag_lfc_score"])
            manual = pd.to_numeric(score_df["manual_ucsc_lfc"])
            diff = (ag - manual).abs()
            corr = ag.corr(manual) if len(ag) > 1 else float("nan")
            logging.info("GeneMaskLFC vs Manual-UCSC-masked LFC (n=%d):",
                         len(score_df))
            logging.info("  Pearson r:         %.4f", corr)
            logging.info("  Mean |diff|:       %.6f", diff.mean())
            logging.info("  Max  |diff|:       %.6f", diff.max())
            if pd.notna(corr) and corr > 0.95:
                logging.info("  Strong agreement — gene masks are "
                             "functionally equivalent.")
            elif pd.notna(corr) and corr > 0.8:
                logging.warning("  Moderate agreement — masks may differ "
                                "at edges.")
            elif pd.notna(corr):
                logging.warning("  Weak agreement — investigate mask "
                                "differences.")

    # --- Signal amplification: masked vs whole-window ---
    if has_manual:
        lfc_df = df[
            pd.to_numeric(df.get("manual_ucsc_lfc", pd.Series()),
                          errors="coerce").notna()
            & pd.to_numeric(df.get("whole_window_lfc", pd.Series()),
                            errors="coerce").notna()
        ].copy()
        if not lfc_df.empty:
            logging.info("Signal amplification (|masked LFC| vs "
                         "|whole-window LFC|):")
            for _, r in lfc_df.iterrows():
                m = abs(float(r["manual_ucsc_lfc"]))
                w = abs(float(r["whole_window_lfc"]))
                gene = r.get("vcf_gene_name", "?")
                ratio = m / w if w > 1e-12 else float("inf")
                logging.info("  %s  masked=%.6f  whole=%.6f  ratio=%.1f×",
                             gene, m, w, ratio)

    # --- Dilution stats ---
    if "ucsc_genes_in_window" in df.columns:
        per_v = df.groupby(["chrom", "pos"])["ucsc_genes_in_window"].first()
        logging.info("Genes per 1 MB window (UCSC):")
        logging.info("  min=%s  median=%.0f  max=%s",
                     per_v.min(), per_v.median(), per_v.max())
        logging.info("  This is why whole-window summing dilutes the "
                     "target gene signal.")

    # --- Mask coverage ---
    if "pct_bins_masked" in df.columns:
        pcts = pd.to_numeric(df["pct_bins_masked"], errors="coerce").dropna()
        if not pcts.empty:
            logging.info("Target gene exon coverage (%% of 1 MB window "
                         "bins masked):")
            logging.info("  min=%.2f%%  median=%.2f%%  max=%.2f%%",
                         pcts.min(), pcts.median(), pcts.max())

    logging.info("=" * 72)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run(
    vcf_file=DEFAULT_VCF,
    ucsc_tsv=DEFAULT_UCSC_TSV,
    genemask_output=DEFAULT_GENEMASK_TSV,
    comparison_output=DEFAULT_COMPARISON_TSV,
    tissue_id=DEFAULT_TISSUE,
    skip_api=False,
    compare_only=False,
):
    """Run all steps of the gene-masking comparison."""
    load_dotenv(DOTENV_PATH)

    # Steps 1 & 2 — AlphaGenome API (skippable)
    if not compare_only and not skip_api:
        logging.info("=" * 30 + " STEPS 1–2: AlphaGenome API " + "=" * 30)
        run_alphagenome_calls(vcf_file, genemask_output, tissue_id)

    # Step 3 — Offline comparison
    logging.info("=" * 30 + " STEP 3: Offline comparison " + "=" * 30)
    run_comparison(ucsc_tsv, genemask_output, comparison_output)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    utils.setup_logging(LOG_FILENAME)
    args = parse_args()
    logging.info("Config: %s", vars(args))
    try:
        run(
            vcf_file=args.vcf,
            ucsc_tsv=args.ucsc_tsv,
            genemask_output=args.genemask_output,
            comparison_output=args.comparison_output,
            tissue_id=args.tissue,
            skip_api=args.skip_api,
            compare_only=args.compare_only,
        )
    except PipelineInputError as exc:
        logging.critical("%s", exc)
        sys.exit(1)
