"""Tests for src/s6_gtex_baseline.py — classification, API helpers, CLI, pipeline."""

import os
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from src import s6_gtex_baseline as gtex_mod
from src.exceptions import PipelineInputError


# ===================================================================
# parse_args
# ===================================================================

class TestGtexParseArgs:
    """Test CLI argument parsing."""

    def test_defaults(self):
        args = gtex_mod.parse_args([])
        from src.constants import VALIDATION_TABLE, GTEX_COMPARISON
        assert args.table == VALIDATION_TABLE
        assert args.output == GTEX_COMPARISON
        assert args.tissue == "Lung"

    def test_custom_table(self):
        args = gtex_mod.parse_args(["--table", "my_table.csv"])
        assert args.table == "my_table.csv"

    def test_custom_output_short(self):
        args = gtex_mod.parse_args(["-o", "out.csv"])
        assert args.output == "out.csv"

    def test_custom_tissue(self):
        args = gtex_mod.parse_args(["--tissue", "Brain_Cortex"])
        assert args.tissue == "Brain_Cortex"

# ===================================================================
# classify_silencing
# ===================================================================

class TestClassifySilencing:
    """Test expression classification against GTEx baselines."""

    # --- No GTEx data ---
    def test_none_gtex_returns_no_data(self):
        assert gtex_mod.classify_silencing(10.0, None) == "no GTEx data"

    # --- Tissue-normal silence: both < 1 TPM ---
    def test_both_low(self):
        assert gtex_mod.classify_silencing(0.1, 0.2) == "tissue-normal silence"

    def test_both_zero(self):
        assert gtex_mod.classify_silencing(0.0, 0.0) == "tissue-normal silence"

    # --- Tumour-specific silencing: tumour low, GTEx high ---
    def test_tumour_silenced_gtex_expressed(self):
        assert gtex_mod.classify_silencing(0.3, 5.0) == "tumour-specific silencing"

    def test_tumour_zero_gtex_expressed(self):
        assert gtex_mod.classify_silencing(0.0, 10.0) == "tumour-specific silencing"

    # --- Tumour over-expression: tumour ≥ 4× GTEx ---
    def test_overexpression(self):
        # 20 / 4 = 5.0 → ratio = 5.0 ≥ 4.0
        assert gtex_mod.classify_silencing(20.0, 4.0) == "tumour over-expression"

    def test_overexpression_exact_threshold(self):
        # 4.0 / 1.0 = 4.0 → exactly at threshold → over-expression
        assert gtex_mod.classify_silencing(4.0, 1.0) == "tumour over-expression"

    def test_overexpression_large_ratio(self):
        assert gtex_mod.classify_silencing(100.0, 2.0) == "tumour over-expression"

    # --- Comparable: both expressed, tumour < 4× GTEx ---
    def test_comparable(self):
        assert gtex_mod.classify_silencing(5.0, 4.0) == "comparable"

    def test_comparable_equal(self):
        assert gtex_mod.classify_silencing(10.0, 10.0) == "comparable"

    def test_comparable_tumour_slightly_lower(self):
        assert gtex_mod.classify_silencing(3.0, 10.0) == "comparable"

    # --- Edge: threshold boundaries ---
    def test_tumour_at_threshold_gtex_below(self):
        """Tumour = 1.0 (expressed), GTEx = 0.5 (not expressed)."""
        # tumour_expressed=True, gtex_expressed=False
        # Falls to the over-expression / comparable branch
        # tumour / gtex = 1.0 / 0.5 = 2.0 < 4.0 → comparable
        assert gtex_mod.classify_silencing(1.0, 0.5) == "comparable"

    def test_tumour_below_threshold_gtex_at_threshold(self):
        """Tumour = 0.9 (not expressed), GTEx = 1.0 (expressed)."""
        assert gtex_mod.classify_silencing(0.9, 1.0) == "tumour-specific silencing"

    def test_gtex_zero_tumour_expressed(self):
        """GTEx = 0, tumour expressed — gtex_expressed is False."""
        # Both expressed check: tumour yes, gtex no
        # Falls to: tumour_expressed=True, gtex_expressed=False
        # Then: gtex_tpm > 0 is False → skip over-expression → comparable
        assert gtex_mod.classify_silencing(10.0, 0.0) == "comparable"


# ===================================================================
# _get_json
# ===================================================================

class TestGetJson:
    """Test the HTTP GET helper with retries."""

    @patch("src.s6_gtex_baseline.requests.get")
    def test_successful_request(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"data": [{"id": 1}]}
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = gtex_mod._get_json("http://example.com/api")
        assert result == {"data": [{"id": 1}]}
        mock_get.assert_called_once()

    @patch("src.s6_gtex_baseline.requests.get")
    def test_returns_none_after_all_retries(self, mock_get):
        import requests
        mock_get.side_effect = requests.RequestException("timeout")

        result = gtex_mod._get_json("http://example.com/api", retries=2)
        assert result is None
        assert mock_get.call_count == 2

    @patch("src.s6_gtex_baseline.requests.get")
    def test_retries_on_failure_then_succeeds(self, mock_get):
        import requests
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"data": []}
        mock_resp.raise_for_status.return_value = None

        mock_get.side_effect = [
            requests.RequestException("transient"),
            mock_resp,
        ]

        result = gtex_mod._get_json("http://example.com/api", retries=3)
        assert result == {"data": []}
        assert mock_get.call_count == 2


# ===================================================================
# resolve_gencode_id
# ===================================================================

class TestResolveGencodeId:
    """Test GENCODE ID resolution from GTEx reference API."""

    @patch("src.s6_gtex_baseline._get_json")
    def test_found(self, mock_json):
        mock_json.return_value = {
            "data": [{"gencodeId": "ENSG00000141510.17"}]
        }
        result = gtex_mod.resolve_gencode_id("ENSG00000141510")
        assert result == "ENSG00000141510.17"

    @patch("src.s6_gtex_baseline._get_json")
    def test_not_found_returns_none(self, mock_json):
        mock_json.return_value = {"data": []}
        assert gtex_mod.resolve_gencode_id("ENSG_FAKE") is None

    @patch("src.s6_gtex_baseline._get_json")
    def test_api_failure_returns_none(self, mock_json):
        mock_json.return_value = None
        assert gtex_mod.resolve_gencode_id("ENSG00000141510") is None


# ===================================================================
# query_median_expression
# ===================================================================

class TestQueryMedianExpression:
    """Test median expression lookup from GTEx."""

    @patch("src.s6_gtex_baseline._get_json")
    def test_found(self, mock_json):
        mock_json.return_value = {
            "data": [{"median": 12.5}]
        }
        result = gtex_mod.query_median_expression("ENSG00000141510.17", "Lung")
        assert result == 12.5

    @patch("src.s6_gtex_baseline._get_json")
    def test_no_data_returns_none(self, mock_json):
        mock_json.return_value = {"data": []}
        assert gtex_mod.query_median_expression("ENSG.17", "Lung") is None

    @patch("src.s6_gtex_baseline._get_json")
    def test_api_failure_returns_none(self, mock_json):
        mock_json.return_value = None
        assert gtex_mod.query_median_expression("ENSG.17", "Lung") is None


# ===================================================================
# fetch_gtex_baselines
# ===================================================================

class TestFetchGtexBaselines:
    """Test batch GTEx lookup orchestration."""

    @patch("src.s6_gtex_baseline.time.sleep")  # skip delays
    @patch("src.s6_gtex_baseline.query_median_expression")
    @patch("src.s6_gtex_baseline.resolve_gencode_id")
    def test_successful_lookup(self, mock_resolve, mock_expr, mock_sleep):
        mock_resolve.return_value = "ENSG00000141510.17"
        mock_expr.return_value = 42.0

        result = gtex_mod.fetch_gtex_baselines(["ENSG00000141510"], "Lung")
        assert result == {"ENSG00000141510": 42.0}

    @patch("src.s6_gtex_baseline.time.sleep")
    @patch("src.s6_gtex_baseline.resolve_gencode_id")
    def test_unresolvable_gene(self, mock_resolve, mock_sleep):
        mock_resolve.return_value = None

        result = gtex_mod.fetch_gtex_baselines(["ENSG_FAKE"], "Lung")
        assert result == {"ENSG_FAKE": None}

    @patch("src.s6_gtex_baseline.time.sleep")
    @patch("src.s6_gtex_baseline.query_median_expression")
    @patch("src.s6_gtex_baseline.resolve_gencode_id")
    def test_no_expression_data(self, mock_resolve, mock_expr, mock_sleep):
        mock_resolve.return_value = "ENSG00000141510.17"
        mock_expr.return_value = None

        result = gtex_mod.fetch_gtex_baselines(["ENSG00000141510"], "Lung")
        assert result == {"ENSG00000141510": None}


# ===================================================================
# gtex_baseline — end-to-end pipeline
# ===================================================================

class TestGtexBaselinePipeline:
    """Test the full pipeline with mocked API calls."""

    def _write_table(self, tmp_dir):
        path = os.path.join(tmp_dir, "validation_table.csv")
        df = pd.DataFrame({
            "GENE": ["TP53", "EGFR"],
            "GENE_ID": ["ENSG00000141510", "ENSG00000146648"],
            "OBSERVED_TPM": [42.5, 0.3],
            "LOG2_FC": [-1.0, 1.0],
            "VAF": [0.3, 0.5],
            "NMD_FLAG": [False, False],
            "VACCINE_PRIORITY": ["HIGH", "MEDIUM"],
        })
        df.to_csv(path, index=False)
        return path

    @patch("src.s6_gtex_baseline.fetch_gtex_baselines")
    def test_pipeline_creates_output(self, mock_fetch):
        mock_fetch.return_value = {
            "ENSG00000141510": 15.0,
            "ENSG00000146648": 8.0,
        }
        with tempfile.TemporaryDirectory() as td:
            table_path = self._write_table(td)
            out_path = os.path.join(td, "gtex.csv")

            gtex_mod.gtex_baseline(
                table_path=table_path,
                output_path=out_path,
                tissue="Lung",
            )

            assert os.path.isfile(out_path)
            df = pd.read_csv(out_path)
            assert len(df) == 2
            assert "SILENCING_CLASS" in df.columns
            assert "GTEX_LUNG_TPM" in df.columns
            assert "TUMOUR_VS_GTEX_RATIO" in df.columns

    @patch("src.s6_gtex_baseline.fetch_gtex_baselines")
    def test_pipeline_classification(self, mock_fetch):
        mock_fetch.return_value = {
            "ENSG00000141510": 15.0,   # TP53: tumour 42.5, GTEx 15 → ratio 2.83 → comparable
            "ENSG00000146648": 8.0,    # EGFR: tumour 0.3, GTEx 8 → tumour-specific silencing
        }
        with tempfile.TemporaryDirectory() as td:
            table_path = self._write_table(td)
            out_path = os.path.join(td, "gtex.csv")

            gtex_mod.gtex_baseline(
                table_path=table_path,
                output_path=out_path,
            )

            df = pd.read_csv(out_path)
            classes = dict(zip(df["GENE"], df["SILENCING_CLASS"]))
            assert classes["TP53"] == "comparable"
            assert classes["EGFR"] == "tumour-specific silencing"

    def test_missing_table_raises(self):
        with pytest.raises(PipelineInputError, match="Validation table"):
            gtex_mod.gtex_baseline(table_path="/nonexistent/table.csv")
