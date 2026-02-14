# Scoring Metrics for Cancer Vaccine Target Prioritisation

This document explains the metrics used in `s4_score_variants.py` to evaluate whether a somatic variant is likely to produce a well-expressed neoantigen suitable for cancer vaccine targeting.

> **Architecture note:** Raw expression predictions are produced by `s3_gene_expression_prediction.py` (AlphaGenome API). All scoring and metric calculation is performed by `s4_score_variants.py`, which can be re-run without API calls.

> The core question: **"Will this mutation produce enough mutant protein on the tumour cell surface for the immune system to recognise?"**

---

## Output Columns

The prediction script outputs a TSV with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | str | Chromosome |
| `POS` | int | Genomic position (1-based) |
| `REF` | str | Reference allele |
| `ALT` | str | Alternate allele |
| `GENE` | str | Gene symbol (from VEP CSQ) |
| `GENE_ID` | str | Ensembl gene ID (version-stripped) |
| `REF_EXPR` | float | AlphaGenome predicted expression — reference allele |
| `ALT_EXPR` | float | AlphaGenome predicted expression — alternate allele |
| `LOG2_FC` | float | log₂(ALT_EXPR / REF_EXPR) fold-change |
| `STATUS` | str | Gain_of_Expression / Loss_of_Expression / Neutral |
| `VAF` | float | Tumour variant allele frequency |
| `OBSERVED_TPM` | float | Patient RNA-seq expression (TPM) |
| `EXPRESSED` | bool | Whether OBSERVED_TPM ≥ 1.0 |
| `NMD_FLAG` | bool | Whether any transcript triggers nonsense-mediated decay |
| `VACCINE_PRIORITY` | str | Composite score: HIGH / MEDIUM / LOW |

---

## Metric Rationale

### 1. AlphaGenome Predicted Fold-Change (`LOG2_FC`, `STATUS`)

**What it measures:** Whether the somatic mutation is predicted to alter gene expression relative to the reference genome.

**Why it matters:** A mutation that causes a predicted loss of expression (log₂FC < −1.0) suggests the mutant allele will not be transcribed at meaningful levels — the neoantigen won't be produced. Conversely, neutral or gain-of-expression predictions suggest the mutant protein will be made.

**Thresholds:**
- `LOG2_FC > 1.0` → **Gain_of_Expression** (predicted ≥2× increase)
- `LOG2_FC < −1.0` → **Loss_of_Expression** (predicted ≥2× decrease)
- Otherwise → **Neutral**

**Limitations:** AlphaGenome predicts expression from DNA sequence features alone. It does not account for epigenetic silencing, post-transcriptional regulation, or tumour microenvironment effects.

---

### 2. Variant Allele Frequency (`VAF`)

**What it measures:** The fraction of tumour DNA reads carrying the mutant allele, extracted from the VCF `FORMAT/AF` field (tumour sample).

**Why it matters:** VAF is a proxy for **clonality**. A high-VAF variant (≥ 0.2) is present in most tumour cells, meaning the neoantigen will be broadly presented across the tumour. Low-VAF variants (< 0.1) may exist only in a small subclone — even if well-expressed, the immune response would only target a minority of tumour cells.

**Threshold:** `VAF ≥ 0.2` is considered clonal for prioritisation purposes.

**Source:** Extracted from `FORMAT/AF[TUMOR]` in MuTect2 VCFs (sample index 1).

---

### 3. Observed Tumour TPM (`OBSERVED_TPM`, `EXPRESSED`)

**What it measures:** The patient's actual RNA-seq expression level for the mutated gene, measured in Transcripts Per Million (TPM).

**Why it matters:** This is the ground-truth expression measurement from the patient's tumour. A gene with TPM < 1 is effectively silent — regardless of what AlphaGenome predicts, the protein is not being made and no neoantigen will be presented on the cell surface. This is the single most important metric for filtering out non-viable targets.

**Threshold:** `TPM ≥ 1.0` is considered expressed. This is a common RNA-seq threshold used in the TCGA consortium and GTEx to distinguish expressed from silenced genes.

**Source:** Joined from `Example_RNA.csv` (`tpm_unstranded` column) on version-stripped Ensembl gene ID.

---

### 4. Nonsense-Mediated Decay Flag (`NMD_FLAG`)

**What it measures:** Whether VEP annotates any transcript of the variant with `NMD_transcript_variant`.

**Why it matters:** Nonsense-mediated decay is a cellular quality-control mechanism that degrades mRNAs containing premature stop codons. Stop-gained and frameshift variants in non-final exons typically trigger NMD, meaning the mutant mRNA is destroyed before translation. This is one of the most common reasons a HIGH-impact variant fails to produce a neoantigen.

Critically, NMD is **invisible to AlphaGenome** — the model predicts expression from DNA sequence features and does not model mRNA surveillance. A variant might have a neutral AlphaGenome prediction but still be silenced by NMD.

**Source:** Parsed from the VEP CSQ `Consequence` field (index 1). Any transcript containing the string `NMD` triggers the flag.

---

### 5. Composite Vaccine Priority (`VACCINE_PRIORITY`)

**What it measures:** A combined assessment of all four metrics above, distilled into a single actionable label.

**Scoring logic:**

```
IF NMD_FLAG = True OR STATUS = Loss_of_Expression:
    → LOW   (transcript likely degraded or silenced)

ELSE IF EXPRESSED = True AND VAF ≥ 0.2:
    → HIGH  (clonal variant in an expressed gene, no predicted loss)

ELSE:
    → MEDIUM (partial evidence — investigate further)
```

**Interpretation:**

| Priority | Meaning | Action |
|----------|---------|--------|
| **HIGH** | Clonal variant in an expressed gene with no evidence of silencing. Best candidates for vaccine targeting. | Prioritise for further validation (immunogenicity assays, MHC binding prediction). |
| **MEDIUM** | Some positive signals but incomplete evidence. Gene may be lowly expressed, variant may be subclonal, or TPM data is unavailable. | Review manually; consider additional evidence (e.g., protein-level data, MHC binding). |
| **LOW** | Variant is likely non-functional as a vaccine target due to NMD degradation or predicted expression loss. | Deprioritise unless strong orthogonal evidence exists. |

---

## Example Interpretation

| Variant | VAF | TPM | NMD | log₂FC | Priority | Reasoning |
|---------|-----|-----|-----|--------|----------|-----------|
| TP53 stop_gained | 0.41 | 35.2 | No | −0.2 | **HIGH** | Clonal, well-expressed, no NMD or loss |
| KRAS missense | 0.08 | 12.1 | No | 0.1 | **MEDIUM** | Expressed but subclonal (VAF < 0.2) |
| BRCA1 frameshift | 0.35 | 8.4 | Yes | 0.0 | **LOW** | NMD will degrade the transcript |
| ELFN1 stop_gained | 0.42 | 0.0 | No | −1.5 | **LOW** | Gene not expressed + predicted loss |

---

## Limitations & Future Directions

1. **No MHC binding prediction.** Even a well-expressed neoantigen is useless if it cannot be presented by the patient's HLA alleles. Tools like NetMHCpan could be integrated as a downstream filter.

2. **No GTEx baseline.** We do not yet compare tumour TPM against normal tissue expression (GTEx lung). A gene silenced in both tumour and normal lung is a tissue-level absence, not a tumour-specific event.

3. **AlphaGenome tissue specificity.** Predictions use UBERON:0002048 (Lung) but the model's tissue-specific accuracy for cancer samples has not been independently validated.

4. **Single-sample evidence.** All metrics derive from a single TCGA patient. Vaccine prioritisation in practice requires population-level recurrence data and multi-omic confirmation.
