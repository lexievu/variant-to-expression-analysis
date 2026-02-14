# Variant-to-Expression Analysis

Evaluating [AlphaGenome](https://deepmind.google/technologies/alphagenome/)'s ability to predict gene expression changes from cancer-associated somatic DNA variants, using TCGA lung adenocarcinoma (LUAD) whole-genome sequencing data.

---

## Purpose

Somatic mutations in cancer can alter gene expression in ways that drive tumor progression. This pipeline asks: **can a deep learning DNA model (AlphaGenome) accurately predict the expression impact of real cancer variants?**

Starting from a TCGA LUAD paired tumor/normal VCF, the pipeline:

1. **Filters** for high-confidence somatic variants (PASS, configurable VEP impact level, expressed in patient RNA-seq).
2. Queries the **AlphaGenome API** to predict reference vs. alternate RNA-seq expression across a 1 MB window around each variant.
3. **Scores** raw predictions with biological context — log₂ fold-change, VAF, RNA-seq TPM, NMD, and a composite vaccine-priority label.
4. Links predictions to patient RNA-seq data via **harmonised Ensembl gene IDs** (GENCODE v36, version-stripped).

The results can be benchmarked against real patient RNA-seq data and GTEx expression baselines.

---

## Pipeline Overview

```
TCGA LUAD VCF (paired tumor/normal)
        │
        ▼
┌─────────────────────────────────────────────────┐
│  s2_vcf_filter.py                                │
│  PASS + VEP impact (configurable) + RNA match   │
│  → output/high_impact_variants.vcf              │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s3_gene_expression_prediction.py                │
│  AlphaGenome API → raw ref/alt expression sums  │
│  Retry, rate-limit, checkpoint/resume           │
│  → output/raw_predictions.tsv                   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  s4_score_variants.py                            │
│  log₂ FC, VAF, TPM, NMD, vaccine priority      │
│  Re-runnable without API calls                  │
│  → output/scored_variants.tsv                   │
└─────────────────────────────────────────────────┘
```

---

## Repository Structure

```
├── data/
│   ├── Example_RNA.csv                  # Patient RNA-seq data (GENCODE v36 gene IDs + TPM)
│   └── VCF_File/                        # TCGA LUAD somatic VCF (MuTect2)
├── notebooks/
│   ├── analysis.ipynb                   # Exploratory analysis & API prototyping
│   └── alphagenome_input.csv            # Exported variant list for API queries
├── output/
│   ├── high_impact_variants.vcf         # Filtered variants (PASS + impact + RNA match)
│   ├── raw_predictions.tsv              # AlphaGenome raw expression sums
│   ├── scored_variants.tsv              # Scored variants (FC, VAF, TPM, NMD, priority)
│   └── ...                              # Intermediate outputs
├── src/
│   ├── s2_vcf_filter.py                  # Variant filtering (CLI: --impact, --output)
│   ├── s3_gene_expression_prediction.py  # AlphaGenome raw predictions (API, retry, resume)
│   ├── s4_score_variants.py              # Biological scoring (VAF, TPM, NMD, priority)
│   ├── constants.py                     # Shared configuration (UPPERCASE constants)
│   └── utils.py                         # CSQ parsing, gene ID extraction, logging setup
├── tests/
│   └── test_utils.py                    # Unit & integration tests (44 tests)
├── .github/workflows/
│   └── tests.yml                        # CI: pytest on push/PR to main (micromamba)
├── environment.yml                      # Conda environment spec (bioconda + conda-forge)
├── log/                                 # Filter run logs
└── PLANNING.md                          # Execution roadmap & progress tracker
```

---

## Getting Started

### Prerequisites

- Python 3.9+
- A valid **AlphaGenome API key** (stored in a `.env` file at the project root)

### Installation

```bash
# Clone the repository
git clone https://github.com/<your-username>/variant-to-expression-analysis.git
cd variant-to-expression-analysis

# Create the conda environment
conda create -n biotech_challenge python=3.10
conda activate biotech_challenge
pip install cyvcf2 numpy pandas python-dotenv alphagenome
```

### Configuration

Create a `.env` file in the `src/` directory:

```
ALPHAGENOME_API_KEY=your_api_key_here
```

### Running the Pipeline

Activate the environment and run scripts sequentially:

```bash
conda activate biotech_challenge

# Step 1: Filter variants (defaults to HIGH impact; use --impact to change)
python src/s2_vcf_filter.py                              # → output/high_impact_variants.vcf
python src/s2_vcf_filter.py --impact HIGH,MODERATE        # Include missense variants
python src/s2_vcf_filter.py --impact HIGH -o custom.vcf   # Custom output path

# Step 2: Predict expression via AlphaGenome (expensive — uses API)
python src/s3_gene_expression_prediction.py               # → output/raw_predictions.tsv
python src/s3_gene_expression_prediction.py --resume      # Resume interrupted run

# Step 3: Score variants with biological context (cheap — no API)
python src/s4_score_variants.py                           # → output/scored_variants.tsv
```

### Running Tests

```bash
python -m pytest tests/ -v
```

Tests also run automatically on every push/PR to `main` via GitHub Actions.

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
| Expression metric | log₂(alt_sum / ref_sum) fold-change |
| Gain threshold | log₂ FC > 1.0 |
| Loss threshold | log₂ FC < −1.0 |

For a detailed explanation of all scoring metrics (VAF, TPM, NMD, vaccine priority), see [docs/METRICS.md](docs/METRICS.md).
For initial prediction results and interpretation (**model predictions only — not yet validated**), see [docs/RESULTS.md](docs/RESULTS.md).
For the RNA-seq correlation analysis, see [docs/VALIDATION.md](docs/VALIDATION.md).
For why UCSC Genome Browser context retrieval is not needed, see [docs/UCSC_CONTEXT.md](docs/UCSC_CONTEXT.md).

---

## License

See [LICENSE](LICENSE) for details.
