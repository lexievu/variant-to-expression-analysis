"""Tests for src/s1_initial_vcf_processing.py — initial somatic extraction."""

import logging
import os
import tempfile
import textwrap

import pytest

from src.s1_initial_vcf_processing import (
    HEADER,
    _format_genotype,
    _log_summary,
    is_somatic,
    parse_args,
)


# ===================================================================
# parse_args
# ===================================================================

class TestParseArgs:
    def test_defaults(self):
        args = parse_args([])
        assert args.vcf  # non-empty default
        assert args.output.endswith("somatic_variants.txt")

    def test_custom_vcf(self):
        args = parse_args(["--vcf", "/tmp/my.vcf"])
        assert args.vcf == "/tmp/my.vcf"

    def test_custom_output_long(self):
        args = parse_args(["--output", "/tmp/out.txt"])
        assert args.output == "/tmp/out.txt"

    def test_custom_output_short(self):
        args = parse_args(["-o", "/tmp/out.txt"])
        assert args.output == "/tmp/out.txt"


# ===================================================================
# _format_genotype
# ===================================================================

class TestFormatGenotype:
    def test_homozygous_ref(self):
        assert _format_genotype([0, 0, False]) == "0/0"

    def test_heterozygous(self):
        assert _format_genotype([0, 1, False]) == "0/1"

    def test_homozygous_alt(self):
        assert _format_genotype([1, 1, False]) == "1/1"

    def test_multi_allelic(self):
        assert _format_genotype([1, 2, False]) == "1/2"

    def test_phased_ignored(self):
        # The phased flag (index 2) should not affect the output
        assert _format_genotype([0, 1, True]) == "0/1"


# ===================================================================
# is_somatic
# ===================================================================

class TestIsSomatic:
    """Somatic = normal is 0/0 AND tumour has at least one ALT allele."""

    def test_classic_somatic_het(self):
        assert is_somatic([0, 0, False], [0, 1, False]) is True

    def test_classic_somatic_hom(self):
        assert is_somatic([0, 0, False], [1, 1, False]) is True

    def test_germline_both_het(self):
        """Both samples carry the variant → germline, not somatic."""
        assert is_somatic([0, 1, False], [0, 1, False]) is False

    def test_germline_both_hom(self):
        assert is_somatic([1, 1, False], [1, 1, False]) is False

    def test_normal_het_tumor_ref(self):
        """Normal has variant, tumour doesn't — not somatic."""
        assert is_somatic([0, 1, False], [0, 0, False]) is False

    def test_both_ref(self):
        """Neither sample has the variant — not somatic."""
        assert is_somatic([0, 0, False], [0, 0, False]) is False

    def test_multi_allelic_tumor(self):
        """Tumour carries allele 2 (multi-allelic ALT)."""
        assert is_somatic([0, 0, False], [0, 2, False]) is True

    def test_tumor_hom_alt2(self):
        assert is_somatic([0, 0, False], [2, 2, False]) is True

    def test_normal_negative_alleles(self):
        """cyvcf2 uses -1 for missing genotypes — should not count as 0."""
        assert is_somatic([-1, -1, False], [0, 1, False]) is False


# ===================================================================
# _log_summary
# ===================================================================

class TestLogSummary:
    def test_logs_all_stats(self, caplog):
        stats = {"total": 100, "somatic": 42, "skipped": 58}
        with caplog.at_level(logging.INFO):
            _log_summary(stats, "/tmp/out.txt")

        combined = caplog.text
        assert "100" in combined
        assert "42" in combined
        assert "58" in combined
        assert "/tmp/out.txt" in combined

    def test_logs_no_print(self, caplog):
        """Verify summary goes through logging, not print."""
        stats = {"total": 10, "somatic": 3, "skipped": 7}
        with caplog.at_level(logging.INFO):
            _log_summary(stats, "/dev/null")
        # At least the header + 3 stat lines + footer
        assert len(caplog.records) >= 5


# ===================================================================
# HEADER constant
# ===================================================================

class TestHeader:
    def test_header_columns(self):
        cols = HEADER.strip().split("\t")
        assert cols == [
            "CHROM", "POS", "REF", "ALT",
            "Normal_GT", "Tumor_GT", "Type", "Change",
        ]

    def test_header_ends_with_newline(self):
        assert HEADER.endswith("\n")


# ===================================================================
# extract_somatic_variants — integration with a minimal VCF
# ===================================================================

class TestExtractSomaticVariants:
    """Integration tests using a real (tiny) VCF written to a temp file.

    We create a minimal VCF with two samples (NORMAL, TUMOR), one somatic
    and one germline variant, then verify the output file.
    """

    MINIMAL_VCF = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
        chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
        chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1
        chr2\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t1/1
    """)

    @pytest.fixture()
    def vcf_and_output(self, tmp_path):
        """Write a minimal VCF and return (vcf_path, output_path)."""
        vcf_path = str(tmp_path / "test.vcf")
        with open(vcf_path, "w") as fh:
            fh.write(self.MINIMAL_VCF)
        output_path = str(tmp_path / "somatic.txt")
        return vcf_path, output_path

    def test_counts(self, vcf_and_output):
        from src.s1_initial_vcf_processing import extract_somatic_variants

        vcf_path, output_path = vcf_and_output
        stats = extract_somatic_variants(vcf_path, output_path)

        assert stats["total"] == 3
        assert stats["somatic"] == 2  # chr1:100 and chr2:300
        assert stats["skipped"] == 1  # chr1:200 (germline)

    def test_output_file_contents(self, vcf_and_output):
        from src.s1_initial_vcf_processing import extract_somatic_variants

        vcf_path, output_path = vcf_and_output
        extract_somatic_variants(vcf_path, output_path)

        with open(output_path) as fh:
            lines = fh.readlines()

        # Header + 2 somatic rows
        assert len(lines) == 3
        assert lines[0] == HEADER

        # First somatic variant: chr1:100 A>G, normal 0/0, tumor 0/1
        fields = lines[1].strip().split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "100"
        assert fields[2] == "A"
        assert fields[3] == "G"
        assert fields[4] == "0/0"
        assert fields[5] == "0/1"

    def test_missing_vcf_raises(self, tmp_path):
        from src.s1_initial_vcf_processing import extract_somatic_variants
        from src.exceptions import PipelineInputError

        with pytest.raises(PipelineInputError, match="not found"):
            extract_somatic_variants(
                str(tmp_path / "nonexistent.vcf"),
                str(tmp_path / "out.txt"),
            )

    def test_no_print_calls(self, vcf_and_output, capsys):
        """Ensure no output is sent to stdout/stderr via print()."""
        from src.s1_initial_vcf_processing import extract_somatic_variants

        vcf_path, output_path = vcf_and_output
        extract_somatic_variants(vcf_path, output_path)

        captured = capsys.readouterr()
        assert captured.out == ""


# ===================================================================
# _find_sample_indices
# ===================================================================

class TestFindSampleIndices:
    """Test the sample-index lookup helper."""

    def test_standard_order(self):
        from src.s1_initial_vcf_processing import _find_sample_indices

        class FakeVCF:
            samples = ["NORMAL", "TUMOR"]

        normal_idx, tumor_idx = _find_sample_indices(FakeVCF())
        assert normal_idx == 0
        assert tumor_idx == 1

    def test_reversed_order(self):
        from src.s1_initial_vcf_processing import _find_sample_indices

        class FakeVCF:
            samples = ["TUMOR", "NORMAL"]

        normal_idx, tumor_idx = _find_sample_indices(FakeVCF())
        assert normal_idx == 1
        assert tumor_idx == 0

    def test_missing_normal_raises(self):
        from src.s1_initial_vcf_processing import _find_sample_indices

        class FakeVCF:
            samples = ["SAMPLE_A", "TUMOR"]

        with pytest.raises(ValueError, match="NORMAL"):
            _find_sample_indices(FakeVCF())

    def test_missing_tumor_raises(self):
        from src.s1_initial_vcf_processing import _find_sample_indices

        class FakeVCF:
            samples = ["NORMAL", "SAMPLE_B"]

        with pytest.raises(ValueError, match="TUMOR"):
            _find_sample_indices(FakeVCF())
