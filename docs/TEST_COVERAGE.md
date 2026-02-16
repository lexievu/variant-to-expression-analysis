# Test Coverage

Summary of the project's automated test suite.  
Run the suite with `python -m pytest tests/ -v` and measure coverage with:

```bash
python -m coverage run --source=src -m pytest tests/ -q
python -m coverage report --show-missing
```

---

## Coverage by Module

| Module | Stmts | Miss | Cover | Notes |
|--------|------:|-----:|------:|-------|
| `__init__.py` | 0 | 0 | 100% | Package marker |
| `constants.py` | 16 | 0 | 100% | All path constants exercised |
| `exceptions.py` | 1 | 0 | 100% | `PipelineInputError` |
| `utils.py` | 82 | 0 | 100% | CSQ parsing, gene ID lookup, validation helpers |
| `s1_initial_vcf_processing.py` | 71 | 5 | 93% | Somatic extraction, sample-index detection, summary logging |
| `s4_score_variants.py` | 102 | 19 | 81% | Pipeline, TPM lookup, NMD, priority labels |
| `s5_validate.py` | 90 | 8 | 91% | Correlation computation, pipeline output |
| `s6_gtex_baseline.py` | 108 | 8 | 93% | API helpers, classification, pipeline |
| `s2_vcf_filter.py` | 69 | 32 | 54% | CLI & argument parsing covered; `filter_vcf()` requires real VCF fixtures |
| `s3_gene_expression_prediction.py` | 144 | 77 | 47% | Parse args, checkpoint, retry covered; inner AlphaGenome loop excluded |
| **TOTAL** | **683** | **149** | **78%** | |

---

## Test Files

| File | Tests | Covers |
|------|------:|--------|
| `test_initial_vcf_processing.py` | ~30 | s1 `is_somatic`, `extract_somatic_variants`, `_find_sample_indices`, `_log_summary`, `parse_args` |
| `test_utils.py` | ~60 | `get_gene_name`, `get_gene_id`, `parse_csq`, `parse_csq_field`, `parse_args`, `setup_logging`, `load_rna_gene_ids` |
| `test_utils_extra.py` | ~25 | `get_vaf`, `has_nmd`, `load_tpm_lookup`, `validate_gene_id`, `validate_file` |
| `test_prediction.py` | ~14 | s3 `parse_args`, `load_checkpoint`, `_score_with_retry` (mocked API), `run_predictions` input validation |
| `test_score_variants.py` | ~20 | s4 `classify`, `vaccine_priority`, pipeline integration (temp files), missing-file error |
| `test_validate.py` | ~25 | s5 `compute_all_correlations` (Pearson/Spearman, stratified, edge cases), pipeline integration, missing-file error |
| `test_gtex_baseline.py` | ~30 | s6 `classify_silencing`, `parse_args`, `_get_json` retries, `resolve_gencode_id`, `query_median_expression`, `fetch_gtex_baselines`, pipeline integration |
| `conftest.py` | — | Shared fixture: log-capture |

**Total: 210 tests — all passing.**

---

## What Is *Not* Covered (and Why)

### `s2_vcf_filter.py` — `filter_vcf()` (lines 84–122)

The core VCF filtering loop uses `cyvcf2.VCF` and `cyvcf2.Writer`, which
operate on real VCF file handles.  Unit-testing this without shipping a test
VCF fixture is impractical and brittle.  The underlying parsing logic
(`parse_csq_field`, `parse_csq`) **is** fully tested.

### `s3_gene_expression_prediction.py` — inner prediction loop (lines 176–274)

This loop calls `model.score_variant` with `GeneMaskLFCScorer` to make API
predictions.  Fully mocking the AlphaGenome model yields low-value, high-
maintenance tests.  Argument parsing, checkpoint resume, retry logic, and
input validation **are** covered.

### `__main__` entry points (all pipeline scripts)

Each script's `if __name__ == "__main__":` block is a thin wrapper that calls
the already-tested pipeline function.  These are excluded by convention.

---

## Running Coverage Locally

```bash
conda activate biotech_challenge

# Quick test run
python -m pytest tests/ -q

# Full coverage report
python -m coverage run --source=src -m pytest tests/ -q
python -m coverage report --show-missing

# HTML report (opens in browser)
python -m coverage html
open htmlcov/index.html
```
