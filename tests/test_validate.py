"""Tests for src/s5_validate.py — correlation and data-loading logic."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from src import s5_validate as validate_mod
from src.exceptions import PipelineInputError


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
        from src.constants import SCORED_VARIANTS, EXAMPLE_RNA_PATH
        assert args.scored == SCORED_VARIANTS
        assert args.rna == EXAMPLE_RNA_PATH

    def test_custom_scored(self):
        args = validate_mod.parse_args(['--scored', 'my_scored.tsv'])
        assert args.scored == 'my_scored.tsv'

    def test_custom_rna(self):
        args = validate_mod.parse_args(['--rna', 'my_rna.csv'])
        assert args.rna == 'my_rna.csv'


# ===================================================================
# compute_all_correlations
# ===================================================================

class TestComputeAllCorrelations:
    """Test the batch correlation builder."""

    def _make_df(self, n=10):
        """Create a DataFrame with the columns expected by compute_all_correlations."""
        rng = np.random.default_rng(42)
        return pd.DataFrame({
            "LOG2_FC": rng.uniform(-2, 2, n),
            "OBSERVED_TPM": rng.uniform(0, 50, n),
            "ACTIVE_EXPR": rng.uniform(10, 200, n),
        })

    def test_returns_at_least_two_comparisons(self):
        df = self._make_df()
        rows = validate_mod.compute_all_correlations(df)
        # LOG2_FC vs TPM + LOG2_FC vs TPM (log₁₀) = 2 minimum
        assert len(rows) >= 2

    def test_includes_log10_transforms(self):
        df = self._make_df()
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert any("log" in l.lower() for l in labels)

    def test_with_unstranded_column(self):
        df = self._make_df()
        df["unstranded"] = np.random.default_rng(0).integers(100, 10000, len(df))
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert any("raw_counts" in l for l in labels)

    def test_expressed_stratification(self):
        """When enough expressed genes exist, an extra stratified row appears."""
        df = self._make_df(20)
        df["OBSERVED_TPM"] = np.random.default_rng(7).uniform(2, 100, 20)
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert any("expressed only" in l for l in labels)

    def test_too_few_expressed_skips_stratification(self):
        df = self._make_df(5)
        df["OBSERVED_TPM"] = 0.1  # all below threshold
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert not any("expressed only" in l for l in labels)

    def test_active_expr_correlations_included(self):
        """When ACTIVE_EXPR column is present, extra correlations are computed."""
        df = self._make_df(10)
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert any("ACTIVE_EXPR" in l for l in labels)
        assert any("ACTIVE_EXPR" in l and "log" in l.lower() for l in labels)

    def test_active_expr_absent_no_crash(self):
        """When ACTIVE_EXPR column is missing, correlations still work."""
        rng = np.random.default_rng(42)
        df = pd.DataFrame({
            "LOG2_FC": rng.uniform(-2, 2, 10),
            "OBSERVED_TPM": rng.uniform(0, 50, 10),
        })
        rows = validate_mod.compute_all_correlations(df)
        labels = [r["comparison"] for r in rows]
        assert not any("ACTIVE_EXPR" in l for l in labels)
        assert len(rows) >= 2  # LOG2_FC correlations still present


# ===================================================================
# validate — end-to-end pipeline
# ===================================================================

class TestValidatePipeline:
    """Test the full validation pipeline with temporary files."""

    def _write_scored(self, tmp_dir):
        path = os.path.join(tmp_dir, "scored.tsv")
        with open(path, "w") as f:
            f.write(
                "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID"
                "\tLOG2_FC\tACTIVE_EXPR\tSTATUS"
                "\tVAF\tOBSERVED_TPM\tEXPRESSED\tNMD_FLAG\tVACCINE_PRIORITY\n"
            )
            f.write(
                "chr1\t100\tA\tT\tTP53\tENSG00000141510"
                "\t-1.0\t55.0\tNeutral"
                "\t0.3\t42.5\tTrue\tFalse\tHIGH\n"
            )
        return path

    def _write_rna(self, tmp_dir):
        path = os.path.join(tmp_dir, "rna.csv")
        with open(path, "w") as f:
            f.write("gene_id,gene_name,tpm_unstranded,unstranded,fpkm_unstranded\n")
            f.write("ENSG00000141510.18,TP53,42.5,5000,12.3\n")
        return path

    def test_pipeline_creates_both_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            scored_path = self._write_scored(td)
            rna_path = self._write_rna(td)
            table_path = os.path.join(td, "table.csv")
            corr_path = os.path.join(td, "corr.csv")

            validate_mod.validate(
                scored_path=scored_path,
                rna_path=rna_path,
                table_path=table_path,
                correlations_path=corr_path,
            )

            assert os.path.isfile(table_path)
            assert os.path.isfile(corr_path)

            table_df = pd.read_csv(table_path)
            assert len(table_df) == 1
            assert "GENE" in table_df.columns

            corr_df = pd.read_csv(corr_path)
            assert len(corr_df) >= 1
            assert "pearson_r" in corr_df.columns

    def test_missing_scored_raises(self):
        with pytest.raises(PipelineInputError, match="Scored variants file"):
            validate_mod.validate(scored_path="/nonexistent/scored.tsv")
