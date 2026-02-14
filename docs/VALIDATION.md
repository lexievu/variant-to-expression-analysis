# Validation: AlphaGenome Predictions vs. Patient RNA-seq

> **⚠️ These are AlphaGenome model predictions, not validated experimental results.** This report compares the model's predicted expression values against observed RNA-seq data from the same patient. It is a consistency check, not an independent validation of model accuracy.

---

## Data Summary

| Metric | Value |
|--------|-------|
| Patient | TCGA-05-4384 (LUAD) |
| Total variants scored | 8 |
| Expressed genes (TPM ≥ 1) | 7 |
| Silenced genes (TPM < 1) | 1 (ELFN1-AS1) |
| Prediction status | All 8 Neutral |
| Vaccine priority | 3 HIGH, 3 MEDIUM, 2 LOW |

---

## Per-Variant Table

| Gene | REF_EXPR | ALT_EXPR | log₂FC | TPM | Raw Counts | VAF | NMD | Priority |
|------|----------|----------|--------|-----|------------|-----|-----|----------|
| FAM107A | 106,169 | 106,163 | −0.0001 | 10.14 | 608 | 0.236 | No | HIGH |
| LAMC3 | 141,537 | 141,972 | +0.0044 | 25.02 | 1,411 | 0.456 | No | HIGH |
| DOT1L | 257,204 | 257,515 | +0.0017 | 6.31 | 853 | 0.351 | No | HIGH |
| TTC7A | 134,680 | 134,496 | −0.0020 | 26.33 | 2,186 | 0.250 | Yes | LOW |
| TMTC1 | 34,362 | 34,393 | +0.0013 | 4.90 | 367 | 0.171 | No | MEDIUM |
| MMP25 | 218,183 | 218,130 | −0.0004 | 4.47 | 244 | 0.135 | No | MEDIUM |
| ERBB2 | 272,124 | 271,699 | −0.0023 | 222.22 | 11,791 | 0.360 | Yes | LOW |
| ELFN1-AS1 | 70,284 | 70,282 | −0.0000 | 0.44 | 38 | 0.415 | No | MEDIUM |

Full data: [`output/validation_table.csv`](../output/validation_table.csv)

---

## Correlation Analysis

We computed Pearson and Spearman correlations between AlphaGenome's predicted expression values and the patient's actual RNA-seq measurements. The question: **does the model's predicted expression magnitude track real-world gene expression levels?**

### Results

| Comparison | Transform | n | Pearson r | p-value | Spearman ρ | p-value |
|------------|-----------|---|-----------|---------|------------|---------|
| ALT_EXPR vs TPM | none | 8 | +0.548 | 0.160 | +0.476 | 0.233 |
| ALT_EXPR vs TPM | log₁₀ | 8 | +0.511 | 0.195 | +0.476 | 0.233 |
| REF_EXPR vs TPM | none | 8 | +0.550 | 0.158 | +0.476 | 0.233 |
| REF_EXPR vs TPM | log₁₀ | 8 | +0.512 | 0.195 | +0.476 | 0.233 |
| LOG2_FC vs TPM | none | 8 | −0.472 | 0.237 | −0.357 | 0.385 |
| ALT_EXPR vs raw counts | none | 8 | +0.546 | 0.161 | +0.548 | 0.160 |
| ALT_EXPR vs raw counts | log₁₀ | 8 | +0.430 | 0.288 | +0.548 | 0.160 |
| ALT_EXPR vs TPM (expressed only) | none | 7 | +0.521 | 0.230 | +0.286 | 0.535 |

Full data: [`output/validation_correlations.csv`](../output/validation_correlations.csv)

### Interpretation

**No correlations reach statistical significance** (all p > 0.05). However, the sample size (n = 8) gives very low statistical power — with only 8 data points, a correlation would need r ≈ 0.71 to reach significance at α = 0.05.

That said, several patterns are worth noting:

1. **Moderate positive trend (r ≈ 0.55).** Both ALT_EXPR and REF_EXPR show a consistent positive Pearson correlation with observed TPM. Genes with higher predicted expression sums tend to have higher real-world TPM. This is directionally encouraging but far from conclusive.

2. **REF ≈ ALT.** The REF and ALT correlations are nearly identical (r = 0.550 vs 0.548), which is consistent with the model predicting essentially no expression change for any of these variants (all log₂FC < ±0.005).

3. **LOG2_FC negatively correlated with TPM (r = −0.47).** This is weak and non-significant, but hints that higher-expressed genes tend to have slightly more negative predicted fold-changes. This could be noise, or it could reflect a subtle model behaviour where abundant transcripts are predicted to be more sensitive to disruption.

4. **ERBB2 dominates.** With TPM = 222 and raw counts = 11,791, ERBB2 is a high-leverage point. Removing it would likely change the correlations substantially. Any interpretation must acknowledge that a single outlier drives much of the signal.

---

## GTEx Normal Lung Baseline Comparison

To distinguish tumour-specific expression changes from normal lung biology, we compared observed tumour TPM against **GTEx v8 median TPM** for normal lung tissue (queried via the [GTEx Portal API](https://gtexportal.org/api/v2)).

| Gene | Tumour TPM | GTEx Lung TPM | Tumour/GTEx | Classification |
|------|-----------|---------------|-------------|----------------|
| TTC7A | 26.33 | 30.76 | 0.86 | comparable |
| FAM107A | 10.14 | 100.48 | 0.10 | comparable |
| ELFN1-AS1 | 0.44 | 0.08 | 5.34 | tissue-normal silence |
| LAMC3 | 25.02 | 42.06 | 0.59 | comparable |
| TMTC1 | 4.90 | 13.61 | 0.36 | comparable |
| MMP25 | 4.47 | 24.08 | 0.19 | comparable |
| ERBB2 | 222.22 | 47.78 | 4.65 | tumour over-expression |
| DOT1L | 6.31 | 28.86 | 0.22 | comparable |

**Classification criteria:**
- **Tissue-normal silence** — TPM < 1 in both tumour and GTEx (gene is normally low in lung)
- **Tumour-specific silencing** — TPM < 1 in tumour but ≥ 1 in GTEx (potential loss-of-expression hit)
- **Tumour over-expression** — tumour TPM ≥ 4× GTEx median
- **Comparable** — similar expression levels in both

### Key findings

1. **ELFN1-AS1** (the only silenced gene) is classified as **tissue-normal silence** — it has a median of just 0.08 TPM in GTEx normal lung. Its low expression in the tumour (0.44 TPM) is consistent with normal lung biology, not a tumour-specific silencing event.

2. **ERBB2** shows **tumour over-expression** (4.65× normal). This is consistent with known ERBB2/HER2 amplification in a subset of lung adenocarcinomas. Despite being flagged as LOW vaccine priority (due to NMD), its over-expression makes it a well-established therapeutic target.

3. **FAM107A** has a tumour/GTEx ratio of just 0.10 — it is expressed at 10.14 TPM in the tumour but 100.48 TPM in normal lung. While both are above the expressed threshold, this 10-fold *under*-expression in the tumour may reflect tumour-specific down-regulation (FAM107A is a known tumour suppressor).

4. **No tumour-specific silencing** was detected — no gene crossed from expressed in GTEx to silenced in the tumour. This limits the set of candidate loss-of-expression hits.

Full data: [`output/gtex_comparison.csv`](../output/gtex_comparison.csv)
Script: [`src/s6_gtex_baseline.py`](../src/s6_gtex_baseline.py)

---

## Caveats

1. **Sample size is very small (n = 8).** Correlation statistics have low power and should be interpreted with extreme caution. These results are exploratory, not confirmatory.

2. **Predicted expression and TPM measure different things.** AlphaGenome's REF_EXPR / ALT_EXPR are summed RNA-seq track values over a 1 MB window centred on the variant. This includes contributions from *all genes and regulatory elements* within that window, not just the target gene. TPM is a normalised per-gene measure. The two are not directly comparable in absolute terms — only relative trends are meaningful.

3. **All variants are Neutral.** Because AlphaGenome predicts essentially zero expression change for all 8 variants, we cannot assess whether the model's *predicted direction of change* agrees with reality. A meaningful validation of fold-change accuracy requires variants where the model predicts a non-trivial gain or loss.

4. **Single patient, single tissue.** All data comes from one TCGA LUAD patient. Results cannot be generalised without multi-patient analysis.

5. **No independent ground truth.** We are comparing model predictions against observational RNA-seq from the same patient. The RNA-seq TPM reflects the combined effect of *all* regulatory mechanisms (epigenetic, post-transcriptional, microenvironment), not just DNA-sequence-level effects. Agreement or disagreement with AlphaGenome does not confirm or refute the model's DNA-level predictions.

---

## Visualisation Notebook

All plots and the full influence/power analysis are in [`notebooks/prediction_vs_rnaseq.ipynb`](../notebooks/prediction_vs_rnaseq.ipynb). The notebook includes:

| Plot | Description |
|------|-------------|
| Scatter: ALT_EXPR vs TPM | Linear scale, coloured by vaccine priority, with regression line |
| Scatter: log-log | log₁₀ transform compresses the ERBB2 outlier |
| Scatter: ALT_EXPR vs raw counts | Avoids TPM normalisation artefacts |
| Bar chart: rescaled predicted vs observed | Min-max rescaled ALT_EXPR alongside TPM per gene |
| Heatmap: multi-metric summary | Z-scored log₂FC, VAF, TPM, NMD, priority per gene |
| ERBB2 influence analysis | Side-by-side scatter with/without ERBB2 (r drops 0.548 → 0.060) |
| Power curve | Shows n = 24 needed for 80% power at observed r ≈ 0.55 |

---

## Next Steps

- [ ] Expand to HIGH + MODERATE impact variants (46 total) for better statistical power
- [x] ~~Add GTEx lung baselines to distinguish tumour-specific from tissue-normal expression~~ → See GTEx section above.
- [ ] Consider per-gene normalisation of AlphaGenome output (divide by window gene count) for a fairer comparison with TPM
