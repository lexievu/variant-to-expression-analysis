"""Tests for src/utils.py — CSQ annotation parsing utilities."""

import sys
import os
import pytest

# Add src/ to path so we can import utils directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import utils


# ---------------------------------------------------------------------------
# Helpers — lightweight mock that mimics cyvcf2 Variant.INFO.get('CSQ')
# ---------------------------------------------------------------------------

class MockVariant:
    """Minimal stand-in for a cyvcf2 Variant with a CSQ INFO field."""

    def __init__(self, csq_string=None):
        self._csq = csq_string
        self.INFO = self  # allows variant.INFO.get('CSQ')

    def get(self, key, default=None):
        if key == 'CSQ':
            return self._csq
        return default


# ---------------------------------------------------------------------------
# Fixtures — reusable CSQ strings
# ---------------------------------------------------------------------------

SINGLE_TRANSCRIPT = "T|missense_variant|MODERATE|KRAS|ENSG00000133703.13|extra"

MULTI_TRANSCRIPT = (
    "T|stop_gained|HIGH|TP53|ENSG00000141510.18|extra,"
    "T|missense_variant|MODERATE|TP53|ENSG00000141510.18|extra,"
    "A|synonymous_variant|LOW|BRCA1|ENSG00000012048.23|extra"
)

SHORT_FIELDS = "T|missense_variant|MODERATE"  # only 3 fields, < 5

NO_VERSION = "T|missense_variant|MODERATE|EGFR|ENSG00000146648|extra"


# ===================================================================
# get_gene_name
# ===================================================================

class TestGetGeneName:
    def test_returns_symbol(self):
        v = MockVariant(SINGLE_TRANSCRIPT)
        assert utils.get_gene_name(v) == "KRAS"

    def test_no_csq_returns_dot(self):
        v = MockVariant(None)
        assert utils.get_gene_name(v) == "."

    def test_short_fields_returns_symbol_if_index3_exists(self):
        # 3 pipe-delimited fields → index 3 does NOT exist
        v = MockVariant(SHORT_FIELDS)
        assert utils.get_gene_name(v) == "."

    def test_multi_transcript_returns_first_symbol(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        # split('|')[3] operates on the full string before any comma split,
        # so it returns the SYMBOL from the first transcript
        assert utils.get_gene_name(v) == "TP53"

    def test_empty_csq_string(self):
        v = MockVariant("")
        # Empty string is truthy-ish for `if csq:` — but split gives [""]
        # index 3 will raise IndexError → returns "."
        assert utils.get_gene_name(v) == "."


# ===================================================================
# get_gene_id
# ===================================================================

class TestGetGeneId:
    def test_returns_stripped_id(self):
        v = MockVariant(SINGLE_TRANSCRIPT)
        assert utils.get_gene_id(v) == "ENSG00000133703"

    def test_returns_full_id_when_strip_disabled(self):
        v = MockVariant(SINGLE_TRANSCRIPT)
        assert utils.get_gene_id(v, strip_version=False) == "ENSG00000133703.13"

    def test_no_csq_returns_dot(self):
        v = MockVariant(None)
        assert utils.get_gene_id(v) == "."

    def test_no_version_suffix(self):
        v = MockVariant(NO_VERSION)
        # No dot in gene_id → returned as-is
        assert utils.get_gene_id(v) == "ENSG00000146648"

    def test_short_fields_returns_dot(self):
        v = MockVariant(SHORT_FIELDS)
        assert utils.get_gene_id(v) == "."

    def test_empty_csq(self):
        v = MockVariant("")
        assert utils.get_gene_id(v) == "."


# ===================================================================
# parse_csq
# ===================================================================

class TestParseCsq:
    def test_single_transcript(self):
        v = MockVariant(SINGLE_TRANSCRIPT)
        result = utils.parse_csq(v)
        assert len(result) == 1
        entry = result[0]
        assert entry['allele'] == 'T'
        assert entry['consequence'] == 'missense_variant'
        assert entry['impact'] == 'MODERATE'
        assert entry['gene_name'] == 'KRAS'
        assert entry['gene_id'] == 'ENSG00000133703.13'
        assert entry['gene_id_stripped'] == 'ENSG00000133703'

    def test_multi_transcript(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        result = utils.parse_csq(v)
        assert len(result) == 3
        assert result[0]['impact'] == 'HIGH'
        assert result[1]['impact'] == 'MODERATE'
        assert result[2]['gene_name'] == 'BRCA1'

    def test_no_csq_returns_empty(self):
        v = MockVariant(None)
        assert utils.parse_csq(v) == []

    def test_short_fields_skipped(self):
        v = MockVariant(SHORT_FIELDS)
        assert utils.parse_csq(v) == []

    def test_mixed_valid_and_short(self):
        mixed = SINGLE_TRANSCRIPT + ",T|too_short"
        v = MockVariant(mixed)
        result = utils.parse_csq(v)
        # Only the first transcript has ≥ 5 fields
        assert len(result) == 1
        assert result[0]['gene_name'] == 'KRAS'

    def test_gene_id_without_version(self):
        v = MockVariant(NO_VERSION)
        result = utils.parse_csq(v)
        assert result[0]['gene_id'] == 'ENSG00000146648'
        assert result[0]['gene_id_stripped'] == 'ENSG00000146648'

    def test_empty_csq_string(self):
        v = MockVariant("")
        assert utils.parse_csq(v) == []


# ===================================================================
# parse_csq_field (from 2_vcf_filter.py) — integration-style tests
# ===================================================================

# 2_vcf_filter no longer has top-level side effects (logging is deferred
# to setup_logging()), so we can import directly.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from importlib import import_module
vcf_filter = import_module('2_vcf_filter')


class TestParseCsqField:
    """Test the filter-level function that wraps utils.parse_csq
    with impact and gene-set checks."""

    def test_keeps_high_impact_in_gene_set(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        genes = {'ENSG00000141510'}  # TP53
        result = vcf_filter.parse_csq_field(v, genes, {'HIGH'})
        assert len(result) == 1
        assert result[0]['impact'] == 'HIGH'
        assert result[0]['gene_id'] == 'ENSG00000141510'

    def test_skips_when_gene_not_in_set(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        genes = {'ENSG_NOT_PRESENT'}
        assert vcf_filter.parse_csq_field(v, genes, {'HIGH'}) == []

    def test_skips_moderate_when_only_high(self):
        v = MockVariant(SINGLE_TRANSCRIPT)  # MODERATE impact
        genes = {'ENSG00000133703'}
        assert vcf_filter.parse_csq_field(v, genes, {'HIGH'}) == []

    def test_keeps_moderate_when_requested(self):
        v = MockVariant(SINGLE_TRANSCRIPT)
        genes = {'ENSG00000133703'}
        result = vcf_filter.parse_csq_field(v, genes, {'HIGH', 'MODERATE'})
        assert len(result) == 1
        assert result[0]['consequence'] == 'missense_variant'

    def test_high_and_moderate_together(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        genes = {'ENSG00000141510'}  # TP53 — has both HIGH and MODERATE transcripts
        result = vcf_filter.parse_csq_field(v, genes, {'HIGH', 'MODERATE'})
        assert len(result) == 2
        impacts = {r['impact'] for r in result}
        assert impacts == {'HIGH', 'MODERATE'}

    def test_empty_gene_set(self):
        v = MockVariant(MULTI_TRANSCRIPT)
        assert vcf_filter.parse_csq_field(v, set(), {'HIGH'}) == []

    def test_no_csq(self):
        v = MockVariant(None)
        assert vcf_filter.parse_csq_field(v, {'ENSG00000141510'}, {'HIGH'}) == []


# ===================================================================
# parse_args (CLI from 2_vcf_filter.py)
# ===================================================================

class TestParseArgs:
    """Test CLI argument parsing."""

    def test_defaults(self):
        args = vcf_filter.parse_args([])
        from constants import DATA_PATH, EXAMPLE_RNA_PATH, DEFAULT_IMPACT_LEVELS
        assert args.vcf == DATA_PATH
        assert args.rna == EXAMPLE_RNA_PATH
        assert args.output == vcf_filter.DEFAULT_OUTPUT
        assert args.impact == DEFAULT_IMPACT_LEVELS

    def test_custom_impact_high_only(self):
        args = vcf_filter.parse_args(['--impact', 'HIGH'])
        assert args.impact == {'HIGH'}

    def test_custom_impact_comma_separated(self):
        args = vcf_filter.parse_args(['--impact', 'HIGH,MODERATE'])
        assert args.impact == {'HIGH', 'MODERATE'}

    def test_impact_case_insensitive(self):
        args = vcf_filter.parse_args(['--impact', 'high,moderate'])
        assert args.impact == {'HIGH', 'MODERATE'}

    def test_custom_output(self):
        args = vcf_filter.parse_args(['-o', 'my_output.vcf'])
        assert args.output == 'my_output.vcf'

    def test_custom_vcf_and_rna(self):
        args = vcf_filter.parse_args(['--vcf', 'a.vcf', '--rna', 'b.csv'])
        assert args.vcf == 'a.vcf'
        assert args.rna == 'b.csv'
