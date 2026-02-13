# Variant-to-Expression Analysis

Evaluating [AlphaGenome](https://deepmind.google/technologies/alphagenome/)'s ability to predict gene expression changes from cancer-associated somatic DNA variants, using TCGA lung adenocarcinoma (LUAD) whole-genome sequencing data.

---

## Purpose

Somatic mutations in cancer can alter gene expression in ways that drive tumor progression. This pipeline asks: **can a deep learning DNA model (AlphaGenome) accurately predict the expression impact of real cancer variants?**

Starting from a TCGA LUAD paired tumor/normal VCF, the pipeline:

1. **Filters** for high-confidence somatic variants (PASS, HIGH-impact VEP consequence, expressed in patient RNA-seq).
2. Queries the **AlphaGenome API** to predict reference vs. alternate RNA-seq expression across a 1 MB window around each variant.
3. Computes **log₂ fold-change** and classifies each variant as *Gain of Expression*, *Loss of Expression*, or *Neutral*.
4. Links predictions to patient RNA-seq data via **harmonised Ensembl gene IDs** (GENCODE v36, version-stripped).

The results can be benchmarked against real patient RNA-seq data and GTEx expression baselines.

---

## Pipeline Overview

```
TCGA LUAD VCF (paired tumor/normal)
        │
        ▼
┌───────────────────────────────────────────────┐
│  2_vcf_filter.py                              │
│  PASS + VEP HIGH-impact + RNA-seq gene match  │
│  → output/high_impact_variants.vcf             │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│  3_gene_expression_prediction.py               │
│  AlphaGenome API → ref/alt RNA-seq predictions  │
│  Log₂ FC classification + harmonised gene IDs   │
│  → output/alphagenome_hits.tsv                  │
└───────────────────────────────────────────────┘
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
│   ├── high_impact_variants.vcf         # Filtered variants (PASS + HIGH impact + RNA match)
│   ├── alphagenome_hits.tsv             # AlphaGenome expression predictions
│   └── ...                              # Intermediate outputs
├── src/
│   ├── 2_vcf_filter.py                  # Variant filtering (PASS + VEP impact + RNA cross-ref)
│   ├── 3_gene_expression_prediction.py   # AlphaGenome expression prediction
│   ├── constants.py                     # Shared file paths
│   └── utils.py                         # CSQ annotation parsing (gene name, gene ID, full parse)
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

python src/2_vcf_filter.py                 # → output/high_impact_variants.vcf
python src/3_gene_expression_prediction.py  # → output/alphagenome_hits.tsv
```

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

---

## License

See [LICENSE](LICENSE) for details.
