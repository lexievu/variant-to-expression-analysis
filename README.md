# Variant-to-Expression Analysis

Evaluating [AlphaGenome](https://deepmind.google/technologies/alphagenome/)'s ability to predict gene expression changes from cancer-associated DNA mutations, using TCGA lung adenocarcinoma (LUAD) whole-genome sequencing data.

---

## Background — Key Concepts

If you're coming from a biology background, here's a quick refresher on the specialist terms used throughout this project:

| Term | What it means |
|---|---|
| **Somatic variant / mutation** | A DNA change that has occurred in a body (somatic) cell — not inherited. Cancer is driven by somatic mutations that accumulate in tumour cells. |
| **VCF file** | Variant Call Format — a standard text file that lists every DNA position where a patient's genome differs from the reference human genome. Think of it as a spreadsheet of mutations. |
| **Gene expression** | How actively a gene is being "read" (transcribed) by the cell. Measured as the amount of mRNA produced. Higher expression → more mRNA → typically more protein. |
| **RNA-seq / TPM** | RNA sequencing measures the mRNA in a sample. TPM (Transcripts Per Million) is a normalised unit that lets us compare expression levels across genes and samples. |
| **Log₂ fold-change (LOG2_FC)** | A way of expressing how much gene expression has changed. A LOG2_FC of +1 means expression has *doubled*; −1 means it has *halved*; 0 means no change. |
| **AlphaGenome** | A deep-learning AI model (by Google DeepMind) that reads a stretch of DNA sequence and predicts how much each nearby gene will be expressed. We use it to ask: "if we introduce this cancer mutation into the DNA, how does the model think expression will change?" |
| **GeneMaskLFCScorer** | A scoring method within AlphaGenome that masks (covers) a gene's exons to isolate the effect of a single mutation on that gene's predicted expression, returning a log₂ fold-change. |
| **GeneMaskActiveScorer** | A companion scorer that returns the absolute expression level (the higher of REF and ALT signals across exons), used as a predicted-expression proxy to compare against observed TPM. |
| **VEP impact** | The Variant Effect Predictor classifies each mutation's likely effect on protein function as LOW, MODERATE, HIGH, or MODIFIER. We focus on HIGH-impact variants (e.g. those that introduce a premature stop codon). |
| **VAF** | Variant Allele Frequency — the fraction of DNA molecules in the tumour sample that carry the mutation (0–1). A VAF of 0.5 means roughly half the tumour cells have it. |
| **NMD** | Nonsense-Mediated Decay — a cellular quality-control mechanism that destroys mRNAs containing premature stop codons, effectively silencing the gene. |
| **GTEx** | Genotype-Tissue Expression project — a public database of gene expression measured in healthy human tissues. We use it as a "normal" baseline to see if a tumour gene is abnormally high or low. |
| **TCGA** | The Cancer Genome Atlas — a large public dataset of cancer genomics data. Our input VCF comes from the LUAD (lung adenocarcinoma) cohort. |

---

## Purpose

Somatic mutations in cancer can change how much a gene is expressed (turned on or off), and those expression changes can drive tumour growth. This pipeline asks a simple but important question:

> **Can a deep-learning AI model (AlphaGenome) accurately predict the expression impact of real cancer mutations?**

To answer this, we take a real lung cancer patient's mutation data and:

1. **Extract** somatic (tumour-only) mutations by comparing the patient's tumour DNA against their matched normal (healthy) DNA.
2. **Filter** down to a small set of high-confidence, high-impact mutations — those most likely to affect gene expression.
3. **Predict** expression changes using AlphaGenome: for each mutation, the model estimates how much the nearby gene's expression would increase or decrease (reported as a log₂ fold-change).
4. **Score** each variant with additional biological context — e.g. what fraction of tumour cells carry it (VAF), how much mRNA the gene actually produces (TPM), and whether the mutation triggers mRNA destruction (NMD).
5. **Validate** the AI's predictions by comparing them to the patient's *actual* RNA-seq expression data, using statistical correlation.
6. **Compare** the tumour's expression to healthy lung tissue (from the GTEx database) to see which genes are abnormally silenced or overexpressed.

---

## Pipeline Overview

The diagram below shows how data flows through the six scripts. Each box is one step; the arrow shows what feeds into the next.

```
TCGA LUAD VCF (paired tumor/normal)
A file listing every mutation found in a lung
cancer patient's tumour vs their healthy tissue.
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 1 — s1_initial_vcf_processing.py          │
│  Pull out somatic mutations: keep only those    │
│  present in the tumour but absent from normal   │
│  tissue (normal genotype = 0/0).                │
│  → output/somatic_variants.txt                  │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 2 — s2_vcf_filter.py                      │
│  Apply quality filters:                         │
│   • Passed variant-caller QC (PASS flag)        │
│   • Tumour carries the alternative allele       │
│   • Predicted HIGH impact on protein (VEP)      │
│   • Gene is actually expressed in RNA-seq       │
│  → output/high_impact_variants.vcf              │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 3 — s3_gene_expression_prediction.py      │
│  Send each mutation to the AlphaGenome AI and   │
│  ask: "how would this mutation change the       │
│  nearby gene's expression?"                     │
│  Returns a LOG2_FC (log₂ fold-change) and an    │
│  ACTIVE_EXPR (absolute expression level) per    │
│  gene, using both scorers in a single API call. │
│  → output/raw_predictions.tsv                   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 4 — s4_score_variants.py                  │
│  Annotate each variant with extra biology:      │
│   • VAF (how common is it in the tumour?)       │
│   • TPM (how active is the gene in RNA-seq?)    │
│   • NMD (does the mutation trigger mRNA decay?) │
│   • Vaccine-priority flag                       │
│  No API calls — can be re-run cheaply.          │
│  → output/scored_variants.tsv                   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 5 — s5_validate.py                        │
│  Compare the AI's predictions to the patient's  │
│  real RNA-seq data using correlation statistics  │
│  (Pearson & Spearman).                          │
│  → output/validation_table.csv                  │
│  → output/validation_correlations.csv           │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│  Step 6 — s6_gtex_baseline.py                   │
│  Compare tumour gene expression to healthy lung │
│  tissue (GTEx database) and flag genes that are │
│  abnormally silenced or overexpressed.          │
│  → output/gtex_comparison.csv                   │
└─────────────────────────────────────────────────┘
```

---

## Repository Structure

```
├── data/
│   ├── Example_RNA.csv                  # Patient RNA-seq data (gene IDs + expression in TPM)
│   └── VCF_File/                        # TCGA lung cancer mutation file (from MuTect2 variant caller)
├── docs/
│   ├── GENE_DILUTION.md                 # Why GeneMaskLFCScorer is needed (gene dilution problem)
│   ├── METRICS.md                       # Definitions of every scoring metric
│   ├── PIPELINE_REPORT.md               # Reader-friendly summary of the full pipeline
│   ├── RESULTS.md                       # Initial prediction results
│   ├── TEST_COVERAGE.md                 # Test coverage summary
│   ├── UCSC_CONTEXT.md                  # Note on UCSC Genome Browser context
│   └── VALIDATION.md                    # RNA-seq correlation analysis
├── notebooks/
│   └── prediction_vs_rnaseq.ipynb       # Interactive exploration: predictions vs real expression
├── output/                              # Pipeline outputs (auto-generated, not stored in git)
├── src/                                 # Source code for each pipeline step
│   ├── constants.py                     # Shared settings (file paths, thresholds, tissue IDs)
│   ├── exceptions.py                    # Custom error types
│   ├── utils.py                         # Helper functions (parsing, gene ID lookup, etc.)
│   ├── s1_initial_vcf_processing.py     # Step 1 — somatic mutation extraction
│   ├── s2_vcf_filter.py                 # Step 2 — variant filtering (quality + impact)
│   ├── s3_gene_expression_prediction.py # Step 3 — AlphaGenome expression prediction
│   ├── s4_score_variants.py             # Step 4 — biological scoring (VAF, TPM, NMD)
│   ├── s5_validate.py                   # Step 5 — correlation with real RNA-seq data
│   └── s6_gtex_baseline.py              # Step 6 — comparison to healthy tissue (GTEx)
├── tests/                               # 210 automated tests to catch bugs
│   └── ...
├── environment.yml                      # Conda environment specification (pinned dependencies)
├── pyproject.toml                       # Python package metadata
└── PLANNING.md                          # Project roadmap & progress tracker
```

---

## Getting Started

### Prerequisites

- **Python 3.11 or newer** — the programming language the pipeline is written in.
- A valid **AlphaGenome API key** — required to query the AI model. Store it in a `.env` file at the project root (see Configuration below).

### Installation

```bash
# 1. Download ("clone") the repository to your computer
git clone https://github.com/<your-username>/variant-to-expression-analysis.git
cd variant-to-expression-analysis

# 2. Create a conda environment with all the required software libraries
conda env create -f environment.yml
conda activate biotech_challenge
```

If you prefer pip over conda:

```bash
pip install -e ".[dev]"
```

### Configuration

```bash
# Copy the template environment file and fill in your API key
cp .env.example .env
# Then open .env in a text editor and set:
#   ALPHAGENOME_API_KEY=your_real_key
```

### Running the Pipeline

Run each step in order — the output of one step feeds into the next:

```bash
conda activate biotech_challenge

# Step 1: Extract somatic (tumour-only) mutations
python src/s1_initial_vcf_processing.py                   # → output/somatic_variants.txt

# Step 2: Filter to high-confidence, high-impact mutations
python src/s2_vcf_filter.py                              # → output/high_impact_variants.vcf
python src/s2_vcf_filter.py --impact HIGH,MODERATE        # Optionally include moderate-impact mutations too
python src/s2_vcf_filter.py --impact HIGH -o custom.vcf   # Or write to a custom output file

# Step 3: Get AlphaGenome's expression predictions (uses the API — can be slow/expensive)
python src/s3_gene_expression_prediction.py               # → output/raw_predictions.tsv
python src/s3_gene_expression_prediction.py --resume      # Resume if the run was interrupted

# Step 4: Add biological context scores (fast — no API calls needed)
python src/s4_score_variants.py                           # → output/scored_variants.tsv

# Step 5: Validate predictions against the patient's real RNA-seq data
python src/s5_validate.py                                 # → output/validation_table.csv
                                                          # → output/validation_correlations.csv

# Step 6: Compare tumour expression to healthy tissue (GTEx)
python src/s6_gtex_baseline.py                            # → output/gtex_comparison.csv
python src/s6_gtex_baseline.py --tissue Lung              # (Lung is the default tissue)
```

### Running Tests

```bash
# Run all 210 automated tests
python -m pytest tests/ -v

# With a coverage report (shows which lines of code are tested)
python -m coverage run --source=src -m pytest tests/ -q
python -m coverage report --show-missing
```

Tests also run automatically on every push/PR to `main` via GitHub Actions.  
See [docs/TEST_COVERAGE.md](docs/TEST_COVERAGE.md) for a detailed coverage breakdown.

> **Tip:** The prediction script uses a small variant subset by default to keep API costs low during development. Adjust `constants.py` to process more variants.

---

## Key Details

| Parameter | Value | What it means |
|---|---|---|
| Cancer type | Lung Adenocarcinoma (LUAD) | A common subtype of non-small-cell lung cancer |
| Variant caller | GATK MuTect2 | The software that identified mutations by comparing tumour vs normal DNA |
| Genome build | GRCh38 / hg38 | The version of the human reference genome used as a baseline |
| AlphaGenome tissue | Lung (`UBERON:0002048`) | The AI model is told to make lung-specific predictions |
| Prediction window | 1,048,576 bp (~1 million bases) | How much DNA context the model reads around each mutation |
| Expression metric | Per-gene exon-masked log₂ fold-change + absolute expression level | The predicted change in gene expression (LOG2_FC) and the predicted absolute level (ACTIVE_EXPR = max of REF/ALT mean across exons) |
| Gain threshold | LOG2_FC > 1.0 | Expression at least *doubled* — flagged as a gain |
| Loss threshold | LOG2_FC < −1.0 | Expression at least *halved* — flagged as a loss |

### Key biological finding

All 8 HIGH-impact variants were predicted **Neutral** (LOG2_FC ≈ 0). This is the expected result: these are protein-disrupting mutations (frameshift, stop_gained, splice-site) that alter the protein, not the DNA regulatory landscape. AlphaGenome models transcription from sequence features (promoters, enhancers, splice signals), so coding-region disruptions don't change predicted transcription. Any expression reduction from these variants would occur post-transcriptionally via NMD — outside the model's scope. See [docs/RESULTS.md](docs/RESULTS.md) for the full interpretation and implications for vaccine target selection.

### Further reading

- [docs/METRICS.md](docs/METRICS.md) — detailed explanation of every scoring metric (VAF, TPM, NMD, vaccine priority).
- [docs/RESULTS.md](docs/RESULTS.md) — initial prediction results and interpretation (**model predictions only — not yet validated**).
- [docs/VALIDATION.md](docs/VALIDATION.md) — RNA-seq correlation analysis (how well do the AI's predictions match reality?).
- [docs/GENE_DILUTION.md](docs/GENE_DILUTION.md) — the gene dilution problem and why GeneMaskLFCScorer is needed.
- [docs/UCSC_CONTEXT.md](docs/UCSC_CONTEXT.md) — why UCSC Genome Browser context retrieval is not needed.
- [docs/TEST_COVERAGE.md](docs/TEST_COVERAGE.md) — test suite coverage summary.

---

## License

See [LICENSE](LICENSE) for details.
