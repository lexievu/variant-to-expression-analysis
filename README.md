# Variant-to-Expression Analysis

Evaluating [AlphaGenome](https://deepmind.google/technologies/alphagenome/)'s ability to predict gene expression changes from cancer-associated somatic DNA variants, using TCGA lung adenocarcinoma (LUAD) whole-genome sequencing data.

---

## Purpose

Somatic mutations in cancer can alter gene expression in ways that drive tumor progression. This pipeline asks: **can a deep learning DNA model (AlphaGenome) accurately predict the expression impact of real cancer variants?**

Starting from a TCGA LUAD paired tumor/normal VCF, the pipeline:

1. **Extracts** somatic variants from paired tumour/normal VCF (normal 0/0, tumour carries ALT).
2. **Filters** for high-confidence somatic variants (PASS, tumour ALT allele, configurable VEP impact level, expressed in patient RNA-seq).
3. Queries the **AlphaGenome API** using `GeneMaskLFCScorer` to predict per-gene exon-masked log₂ fold-change for each variant.
4. **Scores** raw predictions with biological context — VAF, RNA-seq TPM, NMD, and a composite vaccine-priority label (LOG2_FC passed through from step 3).
5. **Validates** predictions against real patient RNA-seq data via Pearson/Spearman correlations (GENCODE v36, version-stripped gene IDs).
6. **Compares** tumour expression to GTEx normal-tissue baselines to classify silencing status.

---

## Pipeline Overview

```
TCGA LUAD VCF (paired tumor/normal)
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s1_initial_vcf_processing.py                   │
│  Extract somatic variants (normal=0/0)          │
│  → output/somatic_variants.txt                  │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s2_vcf_filter.py                               │
│  PASS + ALT allele + VEP impact + RNA match     │
│  → output/high_impact_variants.vcf              │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s3_gene_expression_prediction.py               │
│  AlphaGenome API → GeneMaskLFCScorer per-gene   │
│  exon-masked LOG2_FC. Retry, checkpoint/resume  │
│  → output/raw_predictions.tsv                   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s4_score_variants.py                           │
│  LOG2_FC (passthrough), VAF, TPM, NMD, priority │
│  Re-runnable without API calls                  │
│  → output/scored_variants.tsv                   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s5_validate.py                                 │
│  Pearson & Spearman correlation vs patient RNA  │
│  → output/validation_table.csv                  │
│  → output/validation_correlations.csv           │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s6_gtex_baseline.py                            │
│  GTEx API → normal-tissue expression baselines  │
│  Silencing classification per gene              │
│  → output/gtex_comparison.csv                   │
└─────────────────────────────────────────────────┘
```

---

## Repository Structure

```
├── data/
│   ├── Example_RNA.csv                  # Patient RNA-seq data (GENCODE v36 gene IDs + TPM)
│   └── VCF_File/                        # TCGA LUAD somatic VCF (MuTect2)
├── docs/
│   ├── GENE_DILUTION.md                 # Gene dilution problem & GeneMaskLFCScorer validation
│   ├── METRICS.md                       # Scoring metric definitions
│   ├── PIPELINE_REPORT.md               # Comprehensive reader-friendly pipeline report
│   ├── RESULTS.md                       # Initial prediction results
│   ├── TEST_COVERAGE.md                 # Test coverage summary
│   ├── UCSC_CONTEXT.md                  # UCSC Genome Browser context note
│   └── VALIDATION.md                    # RNA-seq correlation analysis
├── notebooks/
│   └── prediction_vs_rnaseq.ipynb       # Prediction vs RNA-seq exploration
├── output/                              # Pipeline outputs (git-ignored data files)
├── src/
│   ├── __init__.py                      # Package marker
│   ├── constants.py                     # Shared paths & configuration constants
│   ├── exceptions.py                    # PipelineInputError exception
│   ├── utils.py                         # CSQ parsing, gene ID lookup, validation helpers
│   ├── s1_initial_vcf_processing.py     # Somatic variant extraction (normal=0/0)
│   ├── s2_vcf_filter.py                 # Variant filtering (CLI: --impact, --output)
│   ├── s3_gene_expression_prediction.py # AlphaGenome predictions (GeneMaskLFCScorer)
│   ├── s4_score_variants.py             # Biological scoring (VAF, TPM, NMD, priority)
│   ├── s5_validate.py                   # Correlation analysis against patient RNA-seq
│   └── s6_gtex_baseline.py              # GTEx normal-tissue expression comparison
├── tests/                               # 210 unit & integration tests
│   ├── conftest.py                      # Shared fixtures
│   ├── test_initial_vcf_processing.py   # s1 somatic extraction tests
│   ├── test_utils.py                    # Core utility tests
│   ├── test_utils_extra.py              # Validation helper & VAF/NMD tests
│   ├── test_prediction.py               # s3 prediction module tests
│   ├── test_score_variants.py           # s4 scoring tests
│   ├── test_validate.py                 # s5 validation tests
│   └── test_gtex_baseline.py            # s6 GTEx baseline tests
├── .github/workflows/
│   └── tests.yml                        # CI: pytest on push/PR to main (micromamba)
├── environment.yml                      # Conda environment spec (pinned deps)
├── pyproject.toml                       # PEP 621 package metadata & build config
├── .env.example                         # Template for API key configuration
├── log/                                 # Run logs
└── PLANNING.md                          # Execution roadmap & progress tracker
```

---

## Getting Started

### Prerequisites

- Python 3.11+
- A valid **AlphaGenome API key** (stored in a `.env` file at the project root)

### Installation

```bash
# Clone the repository
git clone https://github.com/<your-username>/variant-to-expression-analysis.git
cd variant-to-expression-analysis

# Create the conda environment (installs all pinned deps + editable package)
conda env create -f environment.yml
conda activate biotech_challenge
```

Alternatively, install with pip alone:

```bash
pip install -e ".[dev]"
```

### Configuration

Copy the example env file and add your API key:

```bash
cp .env.example .env
# Edit .env and set ALPHAGENOME_API_KEY=your_real_key
```

### Running the Pipeline

Activate the environment and run scripts sequentially:

```bash
conda activate biotech_challenge

# Step 1: Extract somatic variants from paired VCF
python src/s1_initial_vcf_processing.py                   # → output/somatic_variants.txt

# Step 2: Filter variants (defaults to HIGH impact; use --impact to change)
python src/s2_vcf_filter.py                              # → output/high_impact_variants.vcf
python src/s2_vcf_filter.py --impact HIGH,MODERATE        # Include missense variants
python src/s2_vcf_filter.py --impact HIGH -o custom.vcf   # Custom output path

# Step 3: Predict expression via AlphaGenome (expensive — uses API)
python src/s3_gene_expression_prediction.py               # → output/raw_predictions.tsv
python src/s3_gene_expression_prediction.py --resume      # Resume interrupted run

# Step 4: Score variants with biological context (cheap — no API)
python src/s4_score_variants.py                           # → output/scored_variants.tsv

# Step 5: Validate predictions against patient RNA-seq
python src/s5_validate.py                                 # → output/validation_table.csv
                                                          # → output/validation_correlations.csv

# Step 6: Compare with GTEx normal-tissue baselines
python src/s6_gtex_baseline.py                            # → output/gtex_comparison.csv
python src/s6_gtex_baseline.py --tissue Lung              # default tissue
```

### Running Tests

```bash
# Run all 210 tests
python -m pytest tests/ -v

# With coverage report
python -m coverage run --source=src -m pytest tests/ -q
python -m coverage report --show-missing
```

Tests also run automatically on every push/PR to `main` via GitHub Actions.  
See [docs/TEST_COVERAGE.md](docs/TEST_COVERAGE.md) for a detailed coverage breakdown.

> **Tip:** The prediction script uses a small variant subset by default to keep API costs low during development. Adjust `constants.py` to process more variants.

---

## Key Details

| Parameter | Value |
|---|---|
| Cancer type | Lung Adenocarcinoma (LUAD) |
| Variant caller | GATK MuTect2 |
| Genome build | GRCh38 / hg38 |
| AlphaGenome tissue | Lung (`UBERON:0002048`) |
| Prediction window | 1,048,576 bp (1 MB) centered on variant |
| Expression metric | Per-gene exon-masked log₂ fold-change (GeneMaskLFCScorer) |
| Gain threshold | log₂ FC > 1.0 |
| Loss threshold | log₂ FC < −1.0 |

For a detailed explanation of all scoring metrics (VAF, TPM, NMD, vaccine priority), see [docs/METRICS.md](docs/METRICS.md).
For initial prediction results and interpretation (**model predictions only — not yet validated**), see [docs/RESULTS.md](docs/RESULTS.md).
For the RNA-seq correlation analysis, see [docs/VALIDATION.md](docs/VALIDATION.md).
For why UCSC Genome Browser context retrieval is not needed, see [docs/UCSC_CONTEXT.md](docs/UCSC_CONTEXT.md).
For the test suite coverage summary, see [docs/TEST_COVERAGE.md](docs/TEST_COVERAGE.md).

---

## License

See [LICENSE](LICENSE) for details.
