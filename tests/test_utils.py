"""Tests for src/utils.py — CSQ annotation parsing utilities."""

import sys
import os
import logging
import tempfile
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

    def test_multi_transcript_quirk_no_comma_split(self):
        # get_gene_name does a raw split('|') without first splitting by comma.
        # This works because index 3 always falls within the first transcript's
        # fields. Document this behaviour so a future refactor doesn't break it.
        csq = "A|cons1|HIGH|GENE_A|ENSG_A|x,B|cons2|LOW|GENE_B|ENSG_B|y"
        v = MockVariant(csq)
        assert utils.get_gene_name(v) == "GENE_A"

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

    def test_multi_transcript_returns_first_id(self):
        # Like get_gene_name, get_gene_id does a raw split('|') without
        # first splitting by comma, so it returns the Gene (index 4) from
        # the first transcript.
        csq = "A|cons1|HIGH|GENE_A|ENSG_FIRST.5|x,B|cons2|LOW|GENE_B|ENSG_SECOND.3|y"
        v = MockVariant(csq)
        assert utils.get_gene_id(v) == "ENSG_FIRST"
        assert utils.get_gene_id(v, strip_version=False) == "ENSG_FIRST.5"


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

    def test_skips_low_impact(self):
        """LOW impact should be excluded even when gene is in the valid set."""
        csq = "T|synonymous_variant|LOW|BRCA1|ENSG00000012048.23|extra"
        v = MockVariant(csq)
        genes = {'ENSG00000012048'}
        assert vcf_filter.parse_csq_field(v, genes, {'HIGH'}) == []
        assert vcf_filter.parse_csq_field(v, genes, {'HIGH', 'MODERATE'}) == []


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

    def test_whitespace_in_impact(self):
        """Spaces around commas should be stripped."""
        args = vcf_filter.parse_args(['--impact', 'HIGH, MODERATE'])
        assert args.impact == {'HIGH', 'MODERATE'}

    def test_trailing_comma_produces_empty_element(self):
        """A trailing comma will produce an empty string element.
        Document this edge-case (not necessarily desirable, but current behaviour)."""
        args = vcf_filter.parse_args(['--impact', 'HIGH,'])
        assert '' in args.impact  # empty string from trailing comma

    def test_unknown_impact_accepted(self):
        """parse_args does not validate impact values — any string is accepted."""
        args = vcf_filter.parse_args(['--impact', 'BANANA'])
        assert args.impact == {'BANANA'}


# ===================================================================
# setup_logging (from utils.py)
# ===================================================================

class TestSetupLogging:
    """Test the shared logging configuration helper."""

    def _reset_root_logger(self):
        """Remove all handlers from the root logger between tests."""
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)
            h.close()
        root.setLevel(logging.WARNING)  # Reset to default level

    def setup_method(self):
        self._saved_handlers = logging.getLogger().handlers[:]
        self._reset_root_logger()

    def teardown_method(self):
        self._reset_root_logger()
        # Restore any handlers that existed before (e.g. pytest's LogCaptureHandler)
        root = logging.getLogger()
        for h in self._saved_handlers:
            root.addHandler(h)

    @pytest.mark.no_capture
    def test_creates_handlers(self):
        """setup_logging should attach at least a StreamHandler."""
        # Ensure no handlers exist
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)

        with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as f:
            log_path = f.name
        try:
            utils.setup_logging(log_path)
            handler_types = {type(h) for h in root.handlers}
            assert logging.FileHandler in handler_types
            assert logging.StreamHandler in handler_types
        finally:
            os.unlink(log_path)

    @pytest.mark.no_capture
    def test_idempotent(self):
        """Calling setup_logging twice should not duplicate handlers."""
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)

        with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as f:
            log_path = f.name
        try:
            utils.setup_logging(log_path)
            count_after_first = len(root.handlers)
            utils.setup_logging(log_path)
            count_after_second = len(root.handlers)
            # Guard clause `if not logger.handlers` prevents duplication
            assert count_after_second == count_after_first
        finally:
            os.unlink(log_path)

    @pytest.mark.no_capture
    def test_custom_level(self):
        """setup_logging should respect the level parameter."""
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)

        with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as f:
            log_path = f.name
        try:
            utils.setup_logging(log_path, level=logging.DEBUG)
            assert root.level == logging.DEBUG
        finally:
            os.unlink(log_path)


# ===================================================================
# load_rna_gene_ids (from 2_vcf_filter.py)
# ===================================================================

class TestLoadRnaGeneIds:
    """Test the RNA gene-id loading function from the filter script."""

    def _write_csv(self, lines):
        """Helper to write lines to a temp CSV file."""
        f = tempfile.NamedTemporaryFile(
            mode='w', suffix='.csv', delete=False
        )
        f.write('\n'.join(lines) + '\n')
        f.close()
        return f.name

    def test_strips_version_suffixes(self):
        path = self._write_csv([
            'gene_id,gene_name',
            'ENSG00000141510.18,TP53',
            'ENSG00000133703.13,KRAS',
        ])
        try:
            gene_ids = vcf_filter.load_rna_gene_ids(path)
            assert 'ENSG00000141510' in gene_ids
            assert 'ENSG00000133703' in gene_ids
            # Version-suffixed form should NOT be present
            assert 'ENSG00000141510.18' not in gene_ids
        finally:
            os.unlink(path)

    def test_header_row_excluded(self):
        """The first row is used as column headers by pandas, so it should
        not appear as a gene_id value."""
        path = self._write_csv([
            'gene_id,gene_name',
            'ENSG00000141510.18,TP53',
        ])
        try:
            gene_ids = vcf_filter.load_rna_gene_ids(path)
            assert 'gene_id' not in gene_ids
        finally:
            os.unlink(path)

    def test_includes_non_ensembl_rows(self):
        """Non-Ensembl IDs (no dot) should pass through as-is."""
        path = self._write_csv([
            'gene_id,gene_name',
            'SOME_CUSTOM_ID,MYSTERY',
        ])
        try:
            gene_ids = vcf_filter.load_rna_gene_ids(path)
            assert 'SOME_CUSTOM_ID' in gene_ids
        finally:
            os.unlink(path)

    def test_missing_file_returns_empty(self):
        """A non-existent file should return an empty set (error is logged)."""
        gene_ids = vcf_filter.load_rna_gene_ids('/tmp/this_file_does_not_exist_12345.csv')
        assert gene_ids == set()
