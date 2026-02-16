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

| Gene | log₂FC | TPM | Raw Counts | VAF | NMD | Priority |
|------|--------|-----|------------|-----|-----|----------|
| FAM107A | −0.0001 | 10.14 | 608 | 0.236 | No | HIGH |
| LAMC3 | +0.0044 | 25.02 | 1,411 | 0.456 | No | HIGH |
| DOT1L | +0.0017 | 6.31 | 853 | 0.351 | No | HIGH |
| TTC7A | −0.0020 | 26.33 | 2,186 | 0.250 | Yes | LOW |
| TMTC1 | +0.0013 | 4.90 | 367 | 0.171 | No | MEDIUM |
| MMP25 | −0.0004 | 4.47 | 244 | 0.135 | No | MEDIUM |
| ERBB2 | −0.0023 | 222.22 | 11,791 | 0.360 | Yes | LOW |
| ELFN1-AS1 | −0.0000 | 0.44 | 38 | 0.415 | No | MEDIUM |

Full data: [`output/validation_table.csv`](../output/validation_table.csv)

---

## Correlation Analysis

We computed Pearson and Spearman correlations between AlphaGenome’s predicted fold-change (LOG2_FC, from `GeneMaskLFCScorer`) and the patient’s actual RNA-seq measurements. The question: **does the model’s predicted fold-change track real-world gene expression levels?**

### Results

| Comparison | Transform | n | Pearson r | p-value | Spearman ρ | p-value |
|------------|-----------|---|-----------|---------|------------|---------|
| LOG2_FC vs TPM | none | 8 | −0.472 | 0.237 | −0.357 | 0.385 |
| LOG2_FC vs TPM | log₁₀ | 8 | * | * | * | * |
| LOG2_FC vs raw counts | none | 8 | * | * | * | * |
| LOG2_FC vs TPM (expressed only) | none | 7 | * | * | * | * |

\* See `output/validation_correlations.csv` for full values.

Full data: [`output/validation_correlations.csv`](../output/validation_correlations.csv)

### Interpretation

**No correlations reach statistical significance** (all p > 0.05). However, the sample size (n = 8) gives very low statistical power — with only 8 data points, a correlation would need |r| ≈ 0.85 to reach significance at α = 0.05.

That said, several patterns are worth noting:

1. **Weak negative trend (r ≈ −0.47).** LOG2_FC shows a negative Pearson correlation with observed TPM. Genes with the most negative predicted fold-changes tend to have the highest TPM. This may reflect that highly-expressed genes have more room to decrease, or it may simply be noise.

2. **All fold-changes are near zero.** All 8 variants have |LOG2_FC| < 0.005, so there is essentially no predicted expression change for any variant. The correlation captures a directional tendency among very small values.

3. **ERBB2 dominates.** With TPM = 222 and LOG2_FC = −0.002, ERBB2 is a high-leverage point. Removing it changes the Pearson r from −0.472 to +0.138 — a dramatic shift. Any interpretation must acknowledge that a single outlier drives the correlation.

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

2. **Fold-change vs. absolute level.** AlphaGenome’s `GeneMaskLFCScorer` outputs a per-gene exon-masked log₂ fold-change (how much expression changes due to the variant), while TPM measures the absolute expression level. A gene can have high TPM but near-zero fold-change. Comparing fold-change to absolute level tests whether variant impact correlates with expression magnitude — a useful but indirect relationship.

3. **All variants are Neutral.** Because AlphaGenome predicts essentially zero expression change for all 8 variants, we cannot assess whether the model's *predicted direction of change* agrees with reality. A meaningful validation of fold-change accuracy requires variants where the model predicts a non-trivial gain or loss.

4. **Single patient, single tissue.** All data comes from one TCGA LUAD patient. Results cannot be generalised without multi-patient analysis.

5. **No independent ground truth.** We are comparing model predictions against observational RNA-seq from the same patient. The RNA-seq TPM reflects the combined effect of *all* regulatory mechanisms (epigenetic, post-transcriptional, microenvironment), not just DNA-sequence-level effects. Agreement or disagreement with AlphaGenome does not confirm or refute the model's DNA-level predictions.

---

## Visualisation Notebook

All plots and the full influence/power analysis are in [`notebooks/prediction_vs_rnaseq.ipynb`](../notebooks/prediction_vs_rnaseq.ipynb). The notebook includes:

| Plot | Description |
|------|-------------|
| Scatter: LOG2_FC vs TPM | Linear scale, coloured by vaccine priority, with regression line |
| Scatter: |LOG2_FC| vs log₁₀(TPM) | Magnitude of predicted fold-change vs log-transformed expression |
| Scatter: LOG2_FC vs raw counts | Avoids TPM normalisation artefacts |
| Correlation summary table | Pearson r, Spearman ρ, and p-values for all comparisons |
| Bar chart: LOG2_FC & TPM per gene | Dual-axis: predicted fold-change alongside observed TPM |
| Heatmap: multi-metric summary | Z-scored LOG2_FC, VAF, log₁₀(TPM), NMD, priority per gene |
| ERBB2 influence analysis | Side-by-side scatter with/without ERBB2 (r changes −0.472 → +0.138) |
| Power curve | Shows n = 33 needed for 80% power at observed |r| = 0.47 |
| GTEx paired bar chart | Tumour vs normal lung TPM per gene |
| GTEx divergence plot | log₂ tumour/GTEx ratio per gene |
| Dual scatter: LOG2_FC vs tumour & GTEx | Tumour r = −0.472, GTEx r = −0.072 |
| Three-source heatmap | Z-scored LOG2_FC, tumour TPM, GTEx TPM |

---

## Next Steps

- [ ] Expand to HIGH + MODERATE impact variants (~46 total) for better statistical power
- [x] ~~Add GTEx lung baselines to distinguish tumour-specific from tissue-normal expression~~ → See GTEx section above.
- [x] ~~Per-gene normalisation of AlphaGenome output~~ → Solved by switching to `GeneMaskLFCScorer` (exon-masked per-gene fold-change). See [GENE_DILUTION.md](GENE_DILUTION.md).
