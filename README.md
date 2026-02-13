# Variant-to-Expression Analysis

Evaluating [AlphaGenome](https://deepmind.google/technologies/alphagenome/)'s ability to predict gene expression changes from cancer-associated somatic DNA variants, using TCGA lung adenocarcinoma (LUAD) whole-genome sequencing data.

---

## Purpose

Somatic mutations in cancer can alter gene expression in ways that drive tumor progression. This pipeline asks: **can a deep learning DNA model (AlphaGenome) accurately predict the expression impact of real cancer variants?**

Starting from a TCGA LUAD paired tumor/normal VCF, the pipeline:

1. Identifies **somatic variants** (mutations present only in the tumor).
2. Applies **quality filters** (GATK MuTect2 PASS, TLOD ≥ 6.3, tumor depth ≥ 20) to retain high-confidence calls.
3. Queries the **AlphaGenome API** to predict reference vs. alternate RNA-seq expression across a 1 MB window around each variant.
4. Computes **log₂ fold-change** and classifies each variant as *Gain of Expression*, *Loss of Expression*, or *Neutral*.

The results can be benchmarked against real patient RNA-seq data and GTEx expression baselines.

---

## Pipeline Overview

```
TCGA LUAD VCF (paired tumor/normal)
        │
        ▼
┌──────────────────────────┐
│  1. vcf_processing.py    │  Extract somatic variants (normal 0/0, tumor variant)
│     → somatic_variants.txt│
│     → gene_names.txt     │
└──────────┬───────────────┘
           │
           ▼
┌──────────────────────────┐
│  2. vcf_high_confidence.py│  Filter: PASS + TLOD ≥ 6.3 + tumor depth ≥ 20
│     → high_confidence.vcf │
└──────────┬───────────────┘
           │
           ▼
┌───────────────────────────────────┐
│  3. gene_expression_prediction.py │  AlphaGenome API → ref/alt RNA-seq predictions
│     → alphagenome_hits.tsv        │  Log₂ FC classification
└───────────────────────────────────┘
```

An **alternative filtering path** is available via `2_vcf_filter.py`, which cross-references VEP consequence annotations (CSQ field) against an RNA gene list to retain only HIGH-impact variants with measurable expression.

---

## Repository Structure

```
├── data/
│   ├── Example_RNA.csv                  # Reference RNA gene list for filtering
│   └── VCF_File/                        # TCGA LUAD somatic VCF (MuTect2)
├── notebooks/
│   ├── analysis.ipynb                   # Exploratory analysis & API prototyping
│   └── alphagenome_input.csv            # Exported variant list for API queries
├── output/
│   ├── somatic_variants.txt             # Stage 1 output
│   ├── high_confidence.vcf              # Stage 2 output (full)
│   ├── high_confidence_10.vcf           # Stage 2 output (10-variant subset)
│   ├── high_impact_variants.vcf         # Alternative filter output
│   ├── alphagenome_hits.tsv             # Stage 3 predictions
│   └── gene_names.txt                   # Unique gene symbols from somatic variants
├── src/
│   ├── vcf_processing.py                # Stage 1: somatic variant extraction
│   ├── vcf_high_confidence.py           # Stage 2: quality-based filtering
│   ├── 2_vcf_filter.py                  # Stage 2 (alt): impact + RNA cross-ref filter
│   ├── 3_gene_expression_prediction.py  # Stage 3: AlphaGenome expression prediction
│   ├── global_variables.py              # Shared file paths
│   └── utils.py                         # Helper functions (gene name parsing)
└── log/                                 # Filter run logs
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

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install cyvcf2 numpy pandas python-dotenv alphagenome
```

### Configuration

Create a `.env` file in the project root:

```
ALPHAGENOME_API_KEY=your_api_key_here
```

### Running the Pipeline

Execute the scripts sequentially from the project root:

```bash
python src/vcf_processing.py              # → output/somatic_variants.txt
python src/vcf_high_confidence.py          # → output/high_confidence_10.vcf
python src/3_gene_expression_prediction.py # → output/alphagenome_hits.tsv
```

> **Tip:** The pipeline uses a small variant subset (`high_confidence_10.vcf`) by default to keep API costs low during development. Adjust `global_variables.py` to process more variants.

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

---

## License

See [LICENSE](LICENSE) for details.
