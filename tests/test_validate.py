"""Tests for src/s5_validate.py — correlation and data-loading logic."""

import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from importlib import import_module
validate_mod = import_module('s5_validate')


# ===================================================================
# _corr_row
# ===================================================================

class TestCorrRow:
    """Test the single-correlation helper."""

    def _make_df(self, x, y):
        return pd.DataFrame({"x": x, "y": y})

    def test_perfect_positive_correlation(self):
        df = self._make_df([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])
        row = validate_mod._corr_row(df, "x", "y", "test")
        assert row["n"] == 5
        assert row["pearson_r"] == pytest.approx(1.0, abs=1e-6)
        assert row["spearman_rho"] == pytest.approx(1.0, abs=1e-6)
        assert row["comparison"] == "test"

    def test_perfect_negative_correlation(self):
        df = self._make_df([1, 2, 3, 4, 5], [10, 8, 6, 4, 2])
        row = validate_mod._corr_row(df, "x", "y", "neg")
        assert row["pearson_r"] == pytest.approx(-1.0, abs=1e-6)

    def test_fewer_than_three_returns_nan(self):
        df = self._make_df([1, 2], [3, 4])
        row = validate_mod._corr_row(df, "x", "y", "small")
        assert row["n"] == 2
        assert np.isnan(row["pearson_r"])
        assert np.isnan(row["spearman_rho"])

    def test_nan_values_excluded(self):
        """NaN rows should be dropped before correlating."""
        df = self._make_df([1, np.nan, 3, 4, 5], [2, 4, np.nan, 8, 10])
        row = validate_mod._corr_row(df, "x", "y", "nan_test")
        assert row["n"] == 3  # only rows 0, 3, 4 survive

    def test_log10_transform(self):
        df = self._make_df([10, 100, 1000], [20, 200, 2000])
        row = validate_mod._corr_row(df, "x", "y", "log_test", transform="log10")
        assert row["transform"] == "log10"
        assert row["pearson_r"] == pytest.approx(1.0, abs=1e-4)

    def test_transform_none(self):
        df = self._make_df([1, 2, 3], [4, 5, 6])
        row = validate_mod._corr_row(df, "x", "y", "no_tf")
        assert row["transform"] == "none"

    def test_all_nan_returns_nan(self):
        df = self._make_df([np.nan, np.nan], [np.nan, np.nan])
        row = validate_mod._corr_row(df, "x", "y", "all_nan")
        assert row["n"] == 0
        assert np.isnan(row["pearson_r"])


# ===================================================================
# load_scored
# ===================================================================

class TestLoadScored:
    """Test the TSV loader for scored predictions."""

    def test_loads_tsv(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("GENE\tLOG2_FC\tVAF\n")
            f.write("TP53\t-1.2\t0.3\n")
            f.write("ERBB2\t0.5\t.\n")
            path = f.name
        try:
            df = validate_mod.load_scored(path)
            assert len(df) == 2
            assert df.iloc[0]["GENE"] == "TP53"
            # '.' should be converted to NaN
            assert np.isnan(df.iloc[1]["VAF"])
        finally:
            os.unlink(path)


# ===================================================================
# load_rna
# ===================================================================

class TestLoadRna:
    """Test the RNA-seq CSV loader."""

    def test_loads_and_indexes_by_stripped_id(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("gene_id,gene_name,tpm_unstranded\n")
            f.write("ENSG00000000003.15,TSPAN6,42.5\n")
            f.write("ENSG00000000005.6,TNMD,0.0\n")
            path = f.name
        try:
            df = validate_mod.load_rna(path)
            assert "ENSG00000000003" in df.index
            assert "ENSG00000000005" in df.index
        finally:
            os.unlink(path)

    def test_skips_non_ensembl_rows(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("gene_id,gene_name,tpm_unstranded\n")
            f.write("ENSG00000000003.15,TSPAN6,10.0\n")
            f.write("N_unmapped,,0\n")
            path = f.name
        try:
            df = validate_mod.load_rna(path)
            assert len(df) == 1
        finally:
            os.unlink(path)


# ===================================================================
# parse_args
# ===================================================================

class TestValidateParseArgs:
    """Test CLI argument parsing for the validate script."""

    def test_defaults(self):
        args = validate_mod.parse_args([])
        from constants import SCORED_VARIANTS, EXAMPLE_RNA_PATH
        assert args.scored == SCORED_VARIANTS
        assert args.rna == EXAMPLE_RNA_PATH

    def test_custom_scored(self):
        args = validate_mod.parse_args(['--scored', 'my_scored.tsv'])
        assert args.scored == 'my_scored.tsv'

    def test_custom_rna(self):
        args = validate_mod.parse_args(['--rna', 'my_rna.csv'])
        assert args.rna == 'my_rna.csv'
