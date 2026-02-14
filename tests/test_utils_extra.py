"""Tests for untested functions in src/utils.py: get_vaf, has_nmd, load_tpm_lookup."""

import os
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import utils


# ---------------------------------------------------------------------------
# Minimal mock that mimics the cyvcf2 Variant interface for format() / INFO
# ---------------------------------------------------------------------------

class MockVariantVAF:
    """Mock variant with FORMAT/AF and INFO/CSQ fields."""

    def __init__(self, af_values=None, csq=None):
        """
        af_values: 2-d list  [[normal_af], [tumour_af]] mimicking format('AF').
        csq: CSQ string for INFO.get('CSQ').
        """
        self._af = af_values
        self._csq = csq
        self.INFO = self

    def format(self, field):
        if field == 'AF' and self._af is not None:
            return self._af
        return None

    def get(self, field):
        if field == 'CSQ':
            return self._csq
        return None


# ===================================================================
# get_vaf
# ===================================================================

class TestGetVaf:
    """Test variant allele frequency extraction."""

    def test_typical_tumour_vaf(self):
        v = MockVariantVAF(af_values=[[0.0], [0.35]])
        assert utils.get_vaf(v) == pytest.approx(0.35)

    def test_explicit_tumour_index(self):
        v = MockVariantVAF(af_values=[[0.0], [0.42]])
        assert utils.get_vaf(v, tumor_index=1) == pytest.approx(0.42)

    def test_normal_index(self):
        v = MockVariantVAF(af_values=[[0.01], [0.42]])
        assert utils.get_vaf(v, tumor_index=0) == pytest.approx(0.01)

    def test_no_af_returns_nan(self):
        v = MockVariantVAF(af_values=None)
        assert np.isnan(utils.get_vaf(v))

    def test_index_error_returns_nan(self):
        """Only one sample but requesting index 1 → NaN."""
        v = MockVariantVAF(af_values=[[0.3]])
        assert np.isnan(utils.get_vaf(v, tumor_index=1))

    def test_zero_vaf(self):
        v = MockVariantVAF(af_values=[[0.0], [0.0]])
        assert utils.get_vaf(v) == pytest.approx(0.0)

    def test_high_vaf(self):
        v = MockVariantVAF(af_values=[[0.0], [1.0]])
        assert utils.get_vaf(v) == pytest.approx(1.0)


# ===================================================================
# has_nmd
# ===================================================================

class TestHasNmd:
    """Test NMD transcript detection."""

    def test_nmd_present(self):
        v = MockVariantVAF(csq="A|NMD_transcript_variant|MODIFIER|GENE1|ENSG001")
        assert utils.has_nmd(v) is True

    def test_no_nmd(self):
        v = MockVariantVAF(csq="A|missense_variant|HIGH|GENE1|ENSG001")
        assert utils.has_nmd(v) is False

    def test_no_csq(self):
        v = MockVariantVAF(csq=None)
        assert utils.has_nmd(v) is False

    def test_nmd_in_second_transcript(self):
        csq = "A|missense_variant|HIGH|G|E,A|NMD_transcript_variant|MODIFIER|G|E"
        v = MockVariantVAF(csq=csq)
        assert utils.has_nmd(v) is True

    def test_nmd_substring_matches(self):
        """'NMD' is checked with an `in` test, so partial matches count."""
        v = MockVariantVAF(csq="A|NMD_something|MODIFIER|G|E")
        assert utils.has_nmd(v) is True

    def test_empty_consequence_field(self):
        v = MockVariantVAF(csq="A||HIGH|GENE1|ENSG001")
        assert utils.has_nmd(v) is False

    def test_short_fields_no_crash(self):
        """Only one field — should not crash, just return False."""
        v = MockVariantVAF(csq="A")
        assert utils.has_nmd(v) is False


# ===================================================================
# load_tpm_lookup
# ===================================================================

class TestLoadTpmLookup:
    """Test the RNA-seq CSV → gene-ID-to-TPM lookup builder."""

    def _make_csv(self, tmp_dir, rows, header="gene_id,tpm_unstranded"):
        """Write a CSV and return its path."""
        path = os.path.join(tmp_dir, "rna.csv")
        with open(path, "w") as f:
            f.write(header + "\n")
            for row in rows:
                f.write(row + "\n")
        return path

    def test_basic_lookup(self):
        with tempfile.TemporaryDirectory() as td:
            csv = self._make_csv(td, [
                "ENSG00000000003.15,42.5",
                "ENSG00000000005.6,0.0",
            ])
            lookup = utils.load_tpm_lookup(csv)
        assert lookup["ENSG00000000003"] == pytest.approx(42.5)
        assert lookup["ENSG00000000005"] == pytest.approx(0.0)
        assert len(lookup) == 2

    def test_strips_version(self):
        with tempfile.TemporaryDirectory() as td:
            csv = self._make_csv(td, ["ENSG00000000003.15,10.0"])
            lookup = utils.load_tpm_lookup(csv)
        assert "ENSG00000000003" in lookup
        assert "ENSG00000000003.15" not in lookup

    def test_skips_non_ensembl(self):
        with tempfile.TemporaryDirectory() as td:
            csv = self._make_csv(td, [
                "ENSG00000000003.15,10.0",
                "N_unmapped,0",
                "N_ambiguous,0",
            ])
            lookup = utils.load_tpm_lookup(csv)
        assert len(lookup) == 1
        assert "N_unmapped" not in lookup

    def test_missing_file_returns_empty(self):
        lookup = utils.load_tpm_lookup("/nonexistent/rna.csv")
        assert lookup == {}

    def test_comment_header_skipped(self):
        """RNA CSVs start with a # comment line — verify it's ignored."""
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "rna.csv")
            with open(path, "w") as f:
                f.write("# some comment metadata\n")
                f.write("gene_id,tpm_unstranded\n")
                f.write("ENSG00000000003.1,5.5\n")
            lookup = utils.load_tpm_lookup(path)
        assert lookup["ENSG00000000003"] == pytest.approx(5.5)

    def test_empty_csv_returns_empty(self):
        with tempfile.TemporaryDirectory() as td:
            csv = self._make_csv(td, [])
            lookup = utils.load_tpm_lookup(csv)
        assert lookup == {}
