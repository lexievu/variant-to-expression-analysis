"""Tests for src/s3_gene_expression_prediction.py — checkpoint, retry, CLI, validation."""

import os
import tempfile
import textwrap
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from src import s3_gene_expression_prediction as pred_mod
from src.exceptions import PipelineInputError


# ===================================================================
# parse_args
# ===================================================================

class TestPredParseArgs:
    """Test CLI argument parsing for the prediction script."""

    def test_defaults(self):
        args = pred_mod.parse_args([])
        assert args.vcf == pred_mod.DEFAULT_VCF
        assert args.output == pred_mod.DEFAULT_OUTPUT
        assert args.tissue == pred_mod.DEFAULT_TISSUE
        assert args.resume is False

    def test_custom_vcf(self):
        args = pred_mod.parse_args(["--vcf", "my.vcf"])
        assert args.vcf == "my.vcf"

    def test_custom_output_short(self):
        args = pred_mod.parse_args(["-o", "out.tsv"])
        assert args.output == "out.tsv"

    def test_custom_tissue(self):
        args = pred_mod.parse_args(["--tissue", "UBERON:0000955"])
        assert args.tissue == "UBERON:0000955"

    def test_resume_flag(self):
        args = pred_mod.parse_args(["--resume"])
        assert args.resume is True


# ===================================================================
# _load_checkpoint
# ===================================================================

class TestLoadCheckpoint:
    """Test the checkpoint/resume reader."""

    def test_missing_file_returns_empty(self):
        result = pred_mod._load_checkpoint("/nonexistent/path.tsv")
        assert result == set()

    def test_reads_existing_variants(self):
        content = textwrap.dedent("""\
            CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tLOG2_FC\tACTIVE_EXPR
            chr1\t12345\tA\tT\tTP53\tENSG00000141510\t-1.000000\t50.000000
            chr7\t55249063\tG\tC\tEGFR\tENSG00000146648\t0.500000\t120.000000
        """)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write(content)
            path = f.name
        try:
            result = pred_mod._load_checkpoint(path)
            assert ("chr1", 12345, "A", "T") in result
            assert ("chr7", 55249063, "G", "C") in result
            assert len(result) == 2
        finally:
            os.unlink(path)

    def test_skips_header_line(self):
        content = "CHROM\tPOS\tREF\tALT\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write(content)
            path = f.name
        try:
            result = pred_mod._load_checkpoint(path)
            assert len(result) == 0
        finally:
            os.unlink(path)

    def test_empty_file_returns_empty(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            path = f.name
        try:
            result = pred_mod._load_checkpoint(path)
            assert result == set()
        finally:
            os.unlink(path)


# ===================================================================
# _score_with_retry
# ===================================================================

class TestScoreWithRetry:
    """Test the exponential-backoff retry wrapper."""

    def test_success_on_first_attempt(self):
        model = MagicMock()
        expected = MagicMock()
        model.score_variant.return_value = expected

        with patch("src.s3_gene_expression_prediction.variant_scorers") as mock_vs:
            mock_vs.GeneMaskLFCScorer.return_value = MagicMock()
            mock_vs.GeneMaskActiveScorer.return_value = MagicMock()
            mock_vs.tidy_scores.return_value = expected
            result = pred_mod._score_with_retry(
                model, "interval", "variant",
                max_retries=3, base_delay=0.0,
            )
        assert result is expected
        assert model.score_variant.call_count == 1

    def test_success_after_retry(self):
        model = MagicMock()
        expected = MagicMock()
        model.score_variant.side_effect = [
            RuntimeError("transient"),
            expected,
        ]

        with patch("src.s3_gene_expression_prediction.variant_scorers") as mock_vs:
            mock_vs.GeneMaskLFCScorer.return_value = MagicMock()
            mock_vs.GeneMaskActiveScorer.return_value = MagicMock()
            mock_vs.tidy_scores.return_value = expected
            result = pred_mod._score_with_retry(
                model, "interval", "variant",
                max_retries=3, base_delay=0.0,
            )
        assert result is expected
        assert model.score_variant.call_count == 2

    def test_raises_after_all_retries_exhausted(self):
        model = MagicMock()
        model.score_variant.side_effect = RuntimeError("permanent")

        with patch("src.s3_gene_expression_prediction.variant_scorers") as mock_vs:
            mock_vs.GeneMaskLFCScorer.return_value = MagicMock()
            mock_vs.GeneMaskActiveScorer.return_value = MagicMock()
            with pytest.raises(RuntimeError, match="permanent"):
                pred_mod._score_with_retry(
                    model, "interval", "variant",
                    max_retries=2, base_delay=0.0,
                )
        assert model.score_variant.call_count == 2


# ===================================================================
# run_predictions — input validation
# ===================================================================

class TestRunPredictionsValidation:
    """Test input validation in run_predictions."""

    def test_missing_api_key_raises(self):
        with patch.dict(os.environ, {}, clear=True), \
             patch("src.s3_gene_expression_prediction.load_dotenv"):
            with pytest.raises(PipelineInputError, match="ALPHAGENOME_API_KEY"):
                pred_mod.run_predictions(vcf_file="dummy.vcf")

    def test_missing_vcf_raises(self):
        with patch.dict(os.environ, {"ALPHAGENOME_API_KEY": "fake"}, clear=False), \
             patch("src.s3_gene_expression_prediction.load_dotenv"):
            with pytest.raises(PipelineInputError, match="Input VCF"):
                pred_mod.run_predictions(vcf_file="/nonexistent/file.vcf")
