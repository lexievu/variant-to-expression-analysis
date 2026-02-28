"""Tests for src/s4_score_variants.py — scoring and vaccine-priority logic."""

import math
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from src import s4_score_variants as score_mod
from src.exceptions import PipelineInputError


# ===================================================================
# compute_log2_fc
# ===================================================================

class TestComputeLog2Fc:
    """Test the log₂ fold-change calculation."""

    def test_equal_values_returns_zero(self):
        assert score_mod.compute_log2_fc(100.0, 100.0) == pytest.approx(0.0, abs=1e-6)

    def test_doubling_returns_one(self):
        # alt = 2 * ref → log₂(2) = 1.0
        assert score_mod.compute_log2_fc(100.0, 200.0) == pytest.approx(1.0, abs=1e-4)

    def test_halving_returns_minus_one(self):
        assert score_mod.compute_log2_fc(200.0, 100.0) == pytest.approx(-1.0, abs=1e-4)

    def test_both_zero_returns_zero(self):
        assert score_mod.compute_log2_fc(0.0, 0.0) == 0.0

    def test_ref_zero_alt_positive(self):
        # With the epsilon guard: log₂((alt + 1e-9) / 1e-9) ≈ log₂(alt/1e-9)
        result = score_mod.compute_log2_fc(0.0, 100.0)
        assert result > 30  # very large positive FC

    def test_alt_zero_ref_positive(self):
        result = score_mod.compute_log2_fc(100.0, 0.0)
        assert result < -30  # very large negative FC

    def test_large_values(self):
        result = score_mod.compute_log2_fc(1e6, 1e6)
        assert result == pytest.approx(0.0, abs=1e-6)

    def test_small_difference(self):
        """Typical case: near-neutral prediction."""
        result = score_mod.compute_log2_fc(106168.0, 106162.0)
        assert abs(result) < 0.001


# ===================================================================
# classify
# ===================================================================

class TestClassify:
    """Test fold-change classification into Gain / Loss / Neutral."""

    def test_neutral_zero(self):
        assert score_mod.classify(0.0) == "Neutral"

    def test_neutral_small_positive(self):
        assert score_mod.classify(0.5) == "Neutral"

    def test_neutral_small_negative(self):
        assert score_mod.classify(-0.5) == "Neutral"

    def test_neutral_at_threshold(self):
        """Exactly at the threshold (1.0) — should still be Neutral (not >)."""
        assert score_mod.classify(1.0) == "Neutral"

    def test_gain_above_threshold(self):
        assert score_mod.classify(1.1) == "Gain_of_Expression"

    def test_gain_large(self):
        assert score_mod.classify(5.0) == "Gain_of_Expression"

    def test_loss_below_threshold(self):
        assert score_mod.classify(-1.1) == "Loss_of_Expression"

    def test_loss_at_threshold(self):
        """Exactly -1.0 — should be Neutral (not <)."""
        assert score_mod.classify(-1.0) == "Neutral"

    def test_loss_large(self):
        assert score_mod.classify(-5.0) == "Loss_of_Expression"


# ===================================================================
# vaccine_priority
# ===================================================================

class TestVaccinePriority:
    """Test the composite vaccine-target priority scoring."""

    # --- HIGH priority: expressed, clonal, no NMD, no loss ---
    def test_high_expressed_clonal_no_nmd(self):
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.3, tpm=10.0, nmd_flag=False
        ) == "HIGH"

    def test_high_with_gain(self):
        """Gain of expression should still be HIGH if expressed + clonal."""
        assert score_mod.vaccine_priority(
            log2_fc=2.0, vaf=0.5, tpm=50.0, nmd_flag=False
        ) == "HIGH"

    # --- LOW priority: NMD flag ---
    def test_low_nmd_overrides_everything(self):
        """NMD always forces LOW, even if expressed and clonal."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.5, tpm=100.0, nmd_flag=True
        ) == "LOW"

    # --- LOW priority: expression loss ---
    def test_low_expression_loss(self):
        """Predicted loss of expression → LOW."""
        assert score_mod.vaccine_priority(
            log2_fc=-1.5, vaf=0.5, tpm=10.0, nmd_flag=False
        ) == "LOW"

    def test_low_nmd_and_loss(self):
        """Both NMD and loss → still LOW (not double-counted)."""
        assert score_mod.vaccine_priority(
            log2_fc=-2.0, vaf=0.5, tpm=10.0, nmd_flag=True
        ) == "LOW"

    # --- MEDIUM priority: various partial conditions ---
    def test_medium_expressed_but_subclonal(self):
        """Expressed but subclonal (VAF < 0.2) → MEDIUM."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.1, tpm=10.0, nmd_flag=False
        ) == "MEDIUM"

    def test_medium_clonal_but_not_expressed(self):
        """Clonal but gene not expressed (TPM < 1) → MEDIUM."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.5, tpm=0.5, nmd_flag=False
        ) == "MEDIUM"

    def test_medium_nan_vaf(self):
        """Missing VAF → not clonal → MEDIUM."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=float('nan'), tpm=10.0, nmd_flag=False
        ) == "MEDIUM"

    def test_medium_nan_tpm(self):
        """Missing TPM → not expressed → MEDIUM."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.5, tpm=float('nan'), nmd_flag=False
        ) == "MEDIUM"

    # --- Edge: threshold boundaries ---
    def test_vaf_exactly_at_threshold(self):
        """VAF = 0.2 exactly → clonal (≥), so HIGH if expressed."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.2, tpm=5.0, nmd_flag=False
        ) == "HIGH"

    def test_tpm_exactly_at_threshold(self):
        """TPM = 1.0 exactly → expressed (≥), so HIGH if clonal."""
        assert score_mod.vaccine_priority(
            log2_fc=0.0, vaf=0.3, tpm=1.0, nmd_flag=False
        ) == "HIGH"

    def test_loss_at_boundary(self):
        """log2_fc = -1.0 exactly → NOT loss (must be <), so not forced LOW."""
        assert score_mod.vaccine_priority(
            log2_fc=-1.0, vaf=0.5, tpm=10.0, nmd_flag=False
        ) == "HIGH"


# ===================================================================
# parse_args (CLI from s4_score_variants.py)
# ===================================================================

class TestScoreParseArgs:
    """Test CLI argument parsing for the scoring script."""

    def test_defaults(self):
        args = score_mod.parse_args([])
        from src.constants import RAW_PREDICTIONS, HIGH_IMPACT_VCF, EXAMPLE_RNA_PATH, SCORED_VARIANTS
        assert args.predictions == RAW_PREDICTIONS
        assert args.vcf == HIGH_IMPACT_VCF
        assert args.rna == EXAMPLE_RNA_PATH
        assert args.output == SCORED_VARIANTS

    def test_custom_predictions(self):
        args = score_mod.parse_args(['--predictions', 'my_raw.tsv'])
        assert args.predictions == 'my_raw.tsv'

    def test_custom_output_short(self):
        args = score_mod.parse_args(['-o', 'out.tsv'])
        assert args.output == 'out.tsv'

    def test_custom_vcf_and_rna(self):
        args = score_mod.parse_args(['--vcf', 'v.vcf', '--rna', 'r.csv'])
        assert args.vcf == 'v.vcf'
        assert args.rna == 'r.csv'


# ===================================================================
# score_variants — end-to-end pipeline
# ===================================================================

class TestScoreVariantsPipeline:
    """Test the full scoring pipeline with temporary files."""

    def _write_raw(self, tmp_dir):
        """Create a minimal raw predictions TSV."""
        path = os.path.join(tmp_dir, "raw.tsv")
        with open(path, "w") as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tLOG2_FC\tACTIVE_EXPR\n")
            f.write("chr1\t100\tA\tT\tTP53\tENSG00000141510\t-1.000000\t55.000000\n")
            f.write("chr7\t200\tG\tC\tEGFR\tENSG00000146648\t1.000000\t120.000000\n")
        return path

    def _write_rna(self, tmp_dir):
        """Create a minimal RNA-seq CSV."""
        path = os.path.join(tmp_dir, "rna.csv")
        with open(path, "w") as f:
            f.write("gene_id,tpm_unstranded\n")
            f.write("ENSG00000141510.18,42.5\n")
            f.write("ENSG00000146648.12,0.3\n")
        return path

    def test_pipeline_produces_output(self):
        with tempfile.TemporaryDirectory() as td:
            raw_path = self._write_raw(td)
            rna_path = self._write_rna(td)
            out_path = os.path.join(td, "scored.tsv")

            score_mod.score_variants(
                predictions_file=raw_path,
                vcf_file="/nonexistent/skip.vcf",  # VCF index gracefully missing
                rna_file=rna_path,
                output_file=out_path,
            )

            assert os.path.isfile(out_path)
            df = pd.read_csv(out_path, sep="\t", na_values=".")
            assert len(df) == 2
            # TP53: LOG2_FC = -1.0 (from raw predictions) → Neutral at boundary
            assert df.iloc[0]["GENE"] == "TP53"
            assert df.iloc[0]["LOG2_FC"] == pytest.approx(-1.0, abs=0.01)
            # EGFR: LOG2_FC = 1.0 (from raw predictions) → Neutral at boundary
            assert df.iloc[1]["GENE"] == "EGFR"
            assert df.iloc[1]["LOG2_FC"] == pytest.approx(1.0, abs=0.01)

    def test_pipeline_with_tpm_lookup(self):
        """TPM values should be looked up and reflected in output."""
        with tempfile.TemporaryDirectory() as td:
            raw_path = self._write_raw(td)
            rna_path = self._write_rna(td)
            out_path = os.path.join(td, "scored.tsv")

            score_mod.score_variants(
                predictions_file=raw_path,
                vcf_file="/nonexistent/skip.vcf",
                rna_file=rna_path,
                output_file=out_path,
            )

            df = pd.read_csv(out_path, sep="\t", na_values=".")
            # TP53 is expressed (42.5 TPM ≥ 1.0)
            assert df.iloc[0]["OBSERVED_TPM"] == pytest.approx(42.5)
            assert df.iloc[0]["EXPRESSED"] == True
            # EGFR is not expressed (0.3 TPM < 1.0)
            assert df.iloc[1]["OBSERVED_TPM"] == pytest.approx(0.3)
            assert df.iloc[1]["EXPRESSED"] == False

    def test_empty_predictions_writes_header_only(self):
        with tempfile.TemporaryDirectory() as td:
            raw_path = os.path.join(td, "raw.tsv")
            with open(raw_path, "w") as f:
                f.write("CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tLOG2_FC\tACTIVE_EXPR\n")
            rna_path = self._write_rna(td)
            out_path = os.path.join(td, "scored.tsv")

            score_mod.score_variants(
                predictions_file=raw_path,
                vcf_file="/nonexistent/skip.vcf",
                rna_file=rna_path,
                output_file=out_path,
            )

            assert os.path.isfile(out_path)
            with open(out_path) as f:
                lines = f.readlines()
            # Header only, no data rows
            assert len(lines) == 1
            assert lines[0].startswith("CHROM")

    def test_missing_predictions_raises(self):
        with pytest.raises(PipelineInputError, match="Raw predictions file"):
            score_mod.score_variants(
                predictions_file="/nonexistent/raw.tsv",
            )
