# Pipeline Report: Evaluating AlphaGenome for Cancer Vaccine Target Selection

> **Audience:** This report is written for readers with an A-Level or early-undergraduate biology background. It assumes familiarity with DNA, RNA, gene expression, and mutations, but not with bioinformatics tools or deep learning models.

> **⚠️ Disclaimer:** All AlphaGenome outputs discussed in this report are **model predictions**, not experimentally validated results. They should be treated as computational hypotheses.

---

## 1. What Is This Project Trying to Achieve?

### The Problem

Cancer cells accumulate DNA mutations. Some of these mutations change the proteins a cancer cell produces, and those altered proteins — called **neoantigens** — can be recognised by the immune system. This is the basis of **personalised cancer vaccines**: if we know which mutant proteins a tumour is making, we can train the patient's immune system to attack cells displaying those proteins.

But there is a catch. Not every mutation produces a protein that the immune system can see. For a mutant protein to be a useful vaccine target, the gene carrying the mutation must actually be **expressed** — that is, the cell must be actively reading that gene and producing the corresponding RNA and protein. If a mutated gene is silenced (switched off), there is no mutant protein for the immune system to recognise, and a vaccine targeting it would be useless.

### The Question

This project asks: **Can a deep learning model called AlphaGenome predict whether a mutated gene will still be expressed?**

AlphaGenome, developed by Google DeepMind, is a DNA foundation model. You give it a stretch of DNA sequence, and it predicts what the resulting gene expression (RNA production) would look like. In principle, you could feed it the DNA sequence surrounding a tumour mutation and ask: "Does this mutation cause the gene to produce more RNA, less RNA, or about the same?"

If AlphaGenome's predictions are reliable, it could help prioritise which mutations are worth pursuing as vaccine targets — without needing to run expensive laboratory experiments for every single variant.

### The Approach — A Prototype Evaluation

This project is a **prototype** — a proof of concept to assess whether AlphaGenome's predictions are consistent with real-world gene expression data. We take data from a single lung cancer patient from The Cancer Genome Atlas (TCGA), run AlphaGenome on their tumour mutations, and then compare the predictions against what we actually observe in the patient's RNA sequencing (RNA-seq) data.

We are not yet claiming AlphaGenome works or does not work. We are building the computational infrastructure to ask that question rigorously, and showing what a preliminary evaluation looks like.

---

## 2. The Data: What Goes In

We use publicly available data from a single patient:

| Data Source | Description |
|-------------|-------------|
| **VCF file** (Variant Call Format) | A file listing every DNA mutation found in the patient's tumour, produced by comparing tumour DNA against normal DNA from the same patient. Contains 3,970 somatic (tumour-only) variants. |
| **RNA-seq CSV** | A file containing gene expression measurements from the patient's tumour tissue. For each gene, it reports the **TPM** (Transcripts Per Million) — a normalised measure of how actively that gene is being read. |
| **GTEx database** | A public database of gene expression in healthy tissues from hundreds of donors. We use it to look up what "normal lung" expression looks like for each gene. |

**Patient:** TCGA-05-4384, diagnosed with lung adenocarcinoma (LUAD).

**Variant calling:** The mutations were identified using a tool called **GATK MuTect2**, which is considered a gold-standard method for finding somatic mutations. It compares the tumour genome to the patient's own normal tissue to identify mutations that exist only in the cancer.

---

## 3. The Pipeline: Step by Step

The pipeline has six scripts, each performing one stage of analysis. Each stage reads from the output of the previous one, so they must be run in order.

### Step 1: Initial VCF Processing (`vcf_processing.py`)

**What it does:** Reads the raw VCF file containing all variant calls and extracts somatic variants — mutations present in the tumour but absent from normal tissue.

**Input:** The compressed VCF file from TCGA (approximately 31,000 raw variant calls).

**How it identifies somatic variants:** For each mutation, the VCF records the genotype (which alleles are present) for both the NORMAL and TUMOR samples. A somatic variant is one where:
- The **normal** sample has genotype 0/0 (two copies of the reference/normal allele — no mutation)
- The **tumour** sample carries at least one alternate allele (0/1 or 1/1)

**Output:** A tab-separated text file listing 3,970 somatic variants with their chromosome, position, reference allele, alternate allele, and genotypes.

**Why this matters:** Not all variants in a VCF file are somatic. Some are inherited (germline) variants that the patient has had since birth. We need to isolate the tumour-specific mutations because only those represent potential neoantigens.

---

### Step 2: VCF Filtering (`s2_vcf_filter.py`)

**What it does:** Applies three quality and relevance filters to narrow down from 3,970 somatic variants to a small, high-quality set.

**Input:** The raw VCF file + the patient's RNA-seq data.

**The three filters:**

| Filter | What It Checks | Why It Matters |
|--------|----------------|----------------|
| **PASS filter** | The variant passed MuTect2's internal quality checks (strand bias, mapping quality, contamination, etc.) | Removes technical artefacts — false mutations caused by sequencing errors rather than real biology. |
| **VEP Impact ≥ HIGH** | The variant is annotated by VEP (Variant Effect Predictor) as having HIGH functional impact: stop-gained (introduces a premature stop codon), frameshift (shifts the reading frame), or splice-site disruption. | LOW-impact variants (e.g., synonymous mutations that do not change the protein) are unlikely to produce novel neoantigens. HIGH-impact variants are the most likely to create detectably different proteins. |
| **Gene expressed in RNA-seq** | The gene affected by the variant appears in the patient's RNA-seq data. | If a gene does not appear in the RNA-seq at all, it may be on a chromosome region that was deleted in the tumour, or it may simply never be expressed in lung tissue. There is no point predicting expression for a gene that has no evidence of being active. |

**Key decision — why HIGH impact only?** This was a deliberate trade-off. HIGH-impact variants (frameshifts, stop-gained, splice disruptions) are the most functionally severe and most likely to produce recognisably different proteins. There are only 8 of them that pass all three filters, which gives us a small, manageable set for this prototype. The downside is that 8 variants is too few for robust statistics (more on this in Section 6). A future expansion to include MODERATE-impact variants (such as missense mutations, which change a single amino acid) would increase the set to approximately 46 variants.

**Output:** A filtered VCF file containing 8 HIGH-impact somatic variants in expressed genes.

---

### Step 3: AlphaGenome Expression Prediction (`s3_gene_expression_prediction.py`)

**What it does:** For each of the 8 filtered variants, sends the DNA sequence to AlphaGenome's API and records the model's predicted gene expression.

**Input:** The filtered VCF from Step 2, plus an API key for AlphaGenome.

**How AlphaGenome works (simplified):**
1. The script creates a **1 MB (1,048,576 base pair) window** of DNA sequence centred on each variant.
2. It asks AlphaGenome: "Given this 1 MB stretch of DNA, predict the RNA-seq signal across this region."
3. It does this twice: once with the **reference** (normal) allele, and once with the **alternate** (mutant) allele.
4. The predicted RNA-seq values across the entire 1 MB window are summed into a single number for each allele: **REF_EXPR** (predicted expression with normal DNA) and **ALT_EXPR** (predicted expression with mutant DNA).

**Why a 1 MB window?** AlphaGenome needs surrounding DNA context to make accurate predictions, because gene expression is controlled not just by the gene itself but by distant regulatory elements (enhancers, promoters, insulators) that can be hundreds of thousands of base pairs away. One megabase is the model's standard input size.

**Why sum over the window?** This is a simplification. AlphaGenome outputs a detailed predicted RNA-seq track (expression at every position in the window), but summing it gives a single number that can be easily compared across variants. The trade-off is that this sum captures expression from **all genes in the window**, not just the one we care about. For a gene sitting alone in a gene-sparse region, the sum mostly reflects that gene. For a gene in a crowded neighbourhood, other genes contribute to the total.

**Tissue specificity:** AlphaGenome allows you to specify which tissue type you want predictions for, using standardised ontology codes. We use `UBERON:0002048`, which corresponds to **Lung tissue**. Using the wrong tissue would produce misleading predictions — for example, a gene that is highly expressed in brain but silent in lung.

**Robustness:** The script includes automatic retry logic (if the API fails temporarily, it waits and tries again up to 3 times), checkpointing (if the script crashes partway through, it can resume from where it left off), and rate-limiting (a short pause between API calls to avoid overloading the server).

**Output:** A tab-separated file with 8 rows, one per variant, recording the chromosome, position, reference and alternate alleles, gene name, gene ID, REF_EXPR, and ALT_EXPR.

---

### Step 4: Biological Scoring (`s4_score_variants.py`)

**What it does:** Takes the raw AlphaGenome predictions from Step 3 and enriches each variant with additional biological information to produce a composite vaccine-suitability score.

**Input:** The raw predictions TSV, the filtered VCF, and the RNA-seq CSV.

**What it adds:**

| Metric | How It Is Calculated | What It Means Intuitively |
|--------|---------------------|--------------------------|
| **log₂ fold-change (LOG2_FC)** | $\log_2\!\left(\frac{\text{ALT\_EXPR}}{\text{REF\_EXPR}}\right)$ | How much expression changes when the mutation is introduced. A value of +1 means expression doubles; −1 means it halves; 0 means no change. |
| **Status** | Gain (log₂FC > +1), Loss (log₂FC < −1), or Neutral | A plain-English label for the direction and magnitude of change. |
| **VAF (Variant Allele Frequency)** | Fraction of tumour DNA reads carrying the mutation (from the VCF) | How prevalent the mutation is within the tumour. A VAF of 0.5 means roughly half the tumour cells carry it (a "clonal" mutation present in most cells). A VAF of 0.05 means only 5% of cells have it (a "subclonal" mutation). Clonal mutations make better vaccine targets because more tumour cells would be attacked. |
| **Observed TPM** | The gene's measured expression level in the patient's RNA-seq | Ground truth: is this gene actually active in this patient's tumour? |
| **NMD flag** | Whether VEP annotates the variant with Nonsense-Mediated Decay | NMD is a cellular quality-control mechanism that destroys RNA transcripts containing premature stop codons. If a mutation triggers NMD, the mutant RNA is degraded before it can be translated into protein — meaning no neoantigen is produced, making it a poor vaccine target. |
| **Vaccine Priority** | A composite HIGH / MEDIUM / LOW score | Combines all the above into an overall recommendation. |

**Vaccine priority logic:**
- **LOW** — automatically assigned if the variant triggers NMD (transcript is destroyed) or AlphaGenome predicts a large expression loss.
- **HIGH** — the gene is expressed (TPM ≥ 1), the variant is clonal (VAF ≥ 0.2), there is no NMD, and there is no predicted expression loss.
- **MEDIUM** — everything else.

**Why separate scoring from prediction?** A key architectural decision. The AlphaGenome API calls in Step 3 are slow and expensive. The scoring in Step 4 uses only local data (VCF fields, RNA-seq lookup) and runs in seconds. By separating them, you can adjust scoring thresholds, add new metrics, or re-run the analysis without re-querying the API. This saves both time and cost.

**Output:** A scored TSV file with 8 rows containing all original prediction data plus LOG2_FC, Status, VAF, TPM, NMD, and Vaccine Priority.

**Results at this stage:** All 8 variants were classified as **Neutral** (log₂FC between −0.005 and +0.005). The scoring layer differentiated them into 3 HIGH, 3 MEDIUM, and 2 LOW priority.

---

### Step 5: Validation Against RNA-seq (`s5_validate.py`)

**What it does:** Compares AlphaGenome's predicted expression values against the patient's actual measured gene expression, computing correlation statistics.

**Input:** The scored TSV from Step 4 and the RNA-seq CSV.

**The core question:** If AlphaGenome predicts that Gene X has a high expression sum and Gene Y has a low one, does that match what we see in the RNA-seq? In other words, does AlphaGenome at least get the **relative ordering** of genes right?

**What it computes:**

| Statistic | What It Measures |
|-----------|------------------|
| **Pearson r** | The strength of the linear relationship between predicted and observed expression. Ranges from −1 (perfect negative) to +1 (perfect positive). Values near 0 mean no linear relationship. |
| **Spearman ρ (rho)** | Like Pearson r, but based on ranks rather than raw values. Less sensitive to outliers. If the model gets the ordering right (highest predicted = highest observed), Spearman rho will be high even if the absolute values are on different scales. |
| **p-value** | The probability of seeing a correlation this strong (or stronger) by pure chance. A p-value below 0.05 is conventionally considered statistically significant. |

**Multiple comparisons are run:**
- ALT_EXPR vs. TPM (raw and log-transformed)
- REF_EXPR vs. TPM (raw and log-transformed)
- LOG2_FC vs. TPM (does the predicted change direction track expression level?)
- ALT_EXPR vs. raw read counts (an alternative to TPM that avoids normalisation assumptions)
- ALT_EXPR vs. TPM restricted to expressed genes only (TPM ≥ 1)

**Why both Pearson and Spearman?** Pearson r assumes a linear relationship and can be heavily influenced by outliers. In our data, the gene ERBB2 has a TPM of 222 while most others are below 30 — this single point has an outsized effect on Pearson r. Spearman ρ, which only looks at rank order, is more robust to such outliers.

**Why log-transform?** Gene expression data often spans several orders of magnitude (e.g., 0.44 TPM to 222 TPM). A log₁₀ transformation compresses this range, so that differences between lowly-expressed genes (4 vs. 10 TPM) are given similar visual and statistical weight as differences between highly-expressed genes (40 vs. 100 TPM). Without log-transformation, a single highly-expressed gene can dominate the entire analysis.

**Output:** Two CSV files:
1. A validation table joining predictions with RNA-seq data (per-variant).
2. A correlations table with Pearson r, Spearman ρ, and p-values for each comparison.

---

### Step 6: GTEx Normal Lung Baseline (`s6_gtex_baseline.py`)

**What it does:** For each gene in our set, looks up the typical expression level in **healthy lung tissue** from the GTEx database, then classifies whether the patient's tumour expression is normal, elevated, or reduced compared to healthy tissue.

**Input:** The validation table from Step 5.

**Why do we need this?** Imagine a gene has a TPM of 0.5 in the patient's tumour. Is this gene "silenced by the tumour," making it a poor vaccine target? Or is it just a gene that is normally not expressed in lung tissue, regardless of cancer? Without knowing the healthy baseline, we cannot tell the difference.

GTEx (the Genotype-Tissue Expression project) provides median gene expression values from hundreds of healthy donors across dozens of tissue types. By querying GTEx's Lung tissue data, we get a "normal" reference point for each gene.

**Classification logic:**

| Classification | Criteria | What It Means |
|----------------|----------|---------------|
| **Tissue-normal silence** | TPM < 1 in both tumour and GTEx | The gene is normally not expressed in lung. Its low expression in the tumour is not cancer-related. |
| **Tumour-specific silencing** | TPM < 1 in tumour, but ≥ 1 in GTEx | The gene is normally active in lung but has been switched off in this tumour. A potential loss-of-expression event. |
| **Tumour over-expression** | Tumour TPM ≥ 4× GTEx median | The tumour is producing much more of this gene's product than normal lung. Could indicate gene amplification or oncogenic activation. |
| **Comparable** | None of the above | Expression is broadly similar between tumour and normal lung. |

**Threshold choices:**
- **TPM ≥ 1 for "expressed"**: This is a widely used convention in the field. Below 1 TPM, a gene is considered to have minimal functional expression.
- **4× for "over-expression"**: A 4-fold increase over normal is a commonly used threshold for calling something meaningfully over-expressed, rather than just normal biological variation.

**Output:** A CSV file with each gene's tumour TPM, GTEx lung TPM, tumour/GTEx ratio, and silencing classification.

**Key findings:**
- **ELFN1-AS1** (the only gene with TPM < 1) was classified as "tissue-normal silence" — its GTEx median is just 0.08 TPM, meaning it is naturally inactive in healthy lung too.
- **ERBB2** showed "tumour over-expression" at 4.65× the GTEx median. ERBB2 (also known as HER2) is a well-known oncogene frequently amplified in cancers.
- **FAM107A** has a tumour/GTEx ratio of just 0.10 (10× under-expressed vs. normal lung), which may reflect tumour-specific down-regulation — FAM107A is a known tumour suppressor.

---

## 4. Results: What Did We Find?

### 4.1 AlphaGenome Predicts No Expression Changes

All 8 variants were classified as **Neutral**. The largest predicted fold-change was just 0.0044 (for LAMC3), far below the ±1.0 threshold needed for a Gain or Loss classification. In practical terms, AlphaGenome predicts that none of these HIGH-impact mutations — including frameshifts, stop-gained, and splice-site disruptions — meaningfully change gene expression at the DNA-sequence level.

This is not necessarily wrong. Many coding mutations alter protein function without affecting transcription. A frameshift mutation produces a truncated, non-functional protein, but the gene may still be transcribed at normal levels. The RNA is still made; it just encodes a broken protein. (Whether that RNA survives NMD is a separate question, handled by the NMD flag in Step 4.)

### 4.2 Scoring Layer Adds Biological Differentiation

Although AlphaGenome alone labels everything "Neutral," the biological scoring in Step 4 separates the 8 variants into three tiers:

| Priority | Variants | Reasoning |
|----------|----------|-----------|
| **HIGH** (3) | FAM107A, LAMC3, DOT1L | Clonal (VAF ≥ 0.2), expressed (TPM ≥ 1), no NMD, no predicted loss |
| **MEDIUM** (3) | TMTC1, MMP25, ELFN1-AS1 | Either subclonal (low VAF) or poorly expressed |
| **LOW** (2) | TTC7A, ERBB2 | Both flagged for NMD — the mutant transcript would likely be destroyed before producing protein |

This illustrates why the scoring layer matters. AlphaGenome models DNA-to-RNA effects but does not model post-transcriptional quality control like NMD. The scoring layer integrates information the model cannot access.

### 4.3 Weak Positive Correlation with Real Expression

The validation step found a **moderate positive Pearson correlation (r ≈ 0.55)** between AlphaGenome's predicted expression sum and the patient's observed TPM. Genes for which AlphaGenome predicted higher expression sums tended to have higher measured TPM values in the patient's tumour. However, **none of the correlations reached statistical significance** (all p-values > 0.05).

### 4.4 GTEx Baseline Context

The GTEx comparison revealed that 6 of the 8 genes had "comparable" expression between tumour and normal lung, one gene (ELFN1-AS1) was normally silent in both, and one gene (ERBB2) was dramatically over-expressed in the tumour.

---

## 5. Understanding the Notebook Charts

The Jupyter notebook ([notebooks/prediction_vs_rnaseq.ipynb](../notebooks/prediction_vs_rnaseq.ipynb)) contains 13 visualisations across two sections. Here is what each one shows, and what it would look like if AlphaGenome were performing well.

### Section 1–9: AlphaGenome vs. Patient RNA-seq

#### Chart 1 — Scatter: Predicted vs. Observed Expression (Linear Scale)

**What it shows:** Each dot is one gene. The x-axis is AlphaGenome's predicted expression (ALT_EXPR), and the y-axis is the patient's actual measured TPM. Dots are coloured by vaccine priority (red = HIGH, orange = MEDIUM, green = LOW). A grey dashed line shows the best-fit linear trend.

**What we see:** A mild upward trend — genes with higher predicted expression tend to have higher observed TPM. However, ERBB2 (TPM = 222) sits far from the other points and pulls the trend line towards it. The Pearson r of +0.548 is shown in the annotation box, but the p-value (0.160) indicates it is not statistically significant.

**What it would look like if AlphaGenome were useful:** Points would cluster tightly around the trend line (high r, low p-value), with no single point dominating. The correlation would hold even without ERBB2.

#### Chart 2 — Scatter: Log-Log Scale

**What it shows:** The same data as Chart 1, but with both axes log₁₀-transformed. This compresses the ERBB2 outlier so that the other 7 genes are easier to see.

**What we see:** The correlation weakens slightly on log scale (r ≈ 0.51), and the points are spread out more evenly. This tells us the apparent linear correlation in Chart 1 was partly driven by ERBB2's extreme position.

**What it would look like if AlphaGenome were useful:** Points would still follow a clear diagonal trend even after log transformation, meaning the model ranks low-, medium-, and high-expression genes correctly.

#### Chart 3 — Scatter: Predicted Expression vs. Raw Read Counts

**What it shows:** Same as Chart 1 but using raw sequencing read counts instead of TPM on the y-axis. Raw counts avoid potential artefacts introduced by TPM normalisation (which adjusts for gene length and library size).

**What we see:** A very similar pattern to Chart 1, with r ≈ 0.55. This reassures us that the trend is not an artefact of TPM normalisation.

**What it would look like if AlphaGenome were useful:** A strong positive correlation similar to or better than the TPM version.

#### Chart 4 — Correlation Summary Table

**What it shows:** A table listing Pearson r, Spearman ρ, and p-values for every comparison we computed. Significant p-values would be highlighted in green.

**What we see:** No green highlights — nothing is significant. Pearson r values range from about 0.43 to 0.55. Spearman ρ values are slightly lower (0.29 to 0.55).

**What it would look like if AlphaGenome were useful:** Multiple rows with green-highlighted p-values (below 0.05), especially for ALT_EXPR vs. TPM.

#### Chart 5 — Bar Chart: Rescaled Predicted vs. Observed Expression per Gene

**What it shows:** Paired bars for each gene: one for AlphaGenome's prediction (rescaled to fit the TPM range) and one for observed TPM. The original ALT_EXPR values are annotated above the predicted bars.

**What we see:** The bar heights do not consistently match. For example, ERBB2 has the highest TPM but a mid-range ALT_EXPR; DOT1L has the highest ALT_EXPR but only moderate TPM. This suggests the model does not accurately rank individual genes.

**What it would look like if AlphaGenome were useful:** The two bars for each gene would be similar heights, and the rank ordering would be preserved (tallest predicted bar = tallest observed bar).

#### Chart 6 — Heatmap: Multi-Metric Summary per Gene

**What it shows:** A colour-coded grid with genes as rows and scoring metrics as columns: log₂FC, VAF, log₁₀(TPM), NMD, and vaccine priority. Colours represent z-scores (how far each value is from the average across all 8 variants). Cell labels show the original values.

**What we see:** ERBB2 stands out with the highest TPM and the most negative log₂FC. TTC7A and ERBB2 are the only NMD-positive variants. The heatmap makes it easy to see which genes excel on which metrics.

#### Chart 7 — ERBB2 Influence Analysis

**What it shows:** Two side-by-side scatter plots. The left includes all 8 variants; the right excludes ERBB2. A table above shows how correlations change.

**What we see:** Removing ERBB2 drops the Pearson r from **0.548 to 0.060** — a dramatic collapse. This reveals that essentially the entire apparent correlation was driven by a single data point. The remaining 7 genes show almost no linear relationship between predicted and observed expression.

**Why this matters:** A result that depends entirely on one data point is unreliable. We need more variants to determine whether the correlation is genuine.

#### Chart 8 — Statistical Power Analysis

**What it shows:** Two panels:
- **Left:** A power curve — how statistical power increases with sample size, assuming the true correlation is r = 0.55. The red dashed line marks 80% power (the conventional threshold for "adequate" power). The green dashed line shows where this threshold is crossed (n = 24). The orange dashed line shows our current sample size (n = 8).
- **Right:** The minimum correlation we could detect at each sample size. At n = 8, we would need r ≈ 0.71 to detect a significant relationship — much higher than the r ≈ 0.55 we observe.

**What this tells us:** With only 8 variants, we simply do not have enough data to draw reliable conclusions. We would need at least **24 variants** to have an 80% chance of detecting a correlation as strong as r = 0.55.

### Section 11: GTEx Baseline Comparison

#### Chart 9 (11a) — Paired Bar Chart: Tumour vs. GTEx Expression

**What it shows:** For each gene, a pair of bars: tumour TPM (orange) and GTEx normal lung median TPM (green). A dashed line at TPM = 1 marks the expressed/silenced threshold. The y-axis uses a log scale to accommodate the wide range of values. Genes with non-"comparable" classifications are annotated.

**What we see:** Most genes have higher GTEx expression than tumour expression (the green bar is taller), suggesting the tumour has down-regulated many of these genes. ERBB2 is the exception — its tumour expression is much higher than normal lung. ELFN1-AS1 is below the TPM = 1 line in both conditions (tissue-normal silence).

#### Chart 10 (11b) — Divergence Plot: Tumour/GTEx Expression Ratio

**What it shows:** Horizontal bars showing the log₂ ratio of tumour TPM to GTEx TPM for each gene. The dashed vertical line at 0 represents equal expression. Bars extending to the right indicate tumour over-expression; bars to the left indicate under-expression. Colours correspond to the silencing classification. The raw ratio is annotated next to each bar.

**What we see:** ERBB2 extends to the right (4.65× over-expression). Most other genes extend to the left, with FAM107A showing the most extreme under-expression (0.10×, i.e., 10 times lower in the tumour than in normal lung). MMP25 and DOT1L are also notably under-expressed.

#### Chart 11 (11c) — Dual Scatter: AlphaGenome vs. Tumour TPM and GTEx TPM

**What it shows:** Two scatter plots side by side. The left plots AlphaGenome's predicted expression (ALT_EXPR) against the patient's tumour TPM. The right plots ALT_EXPR against GTEx normal lung TPM. Both include regression lines and correlation statistics.

**What we see:** AlphaGenome correlates better with tumour TPM (r = +0.548) than with GTEx normal lung TPM (r = +0.137). This suggests the model may be capturing something specific to this patient's tumour expression landscape, rather than simply reflecting normal lung biology. However, given the small sample size and the ERBB2 leverage effect, this difference is not statistically reliable.

**What it would look like if AlphaGenome were useful:** A strong correlation with tumour TPM (the actual data the model should be predicting) and a weaker correlation with GTEx TPM (which the model was not asked to predict).

#### Chart 12 (11d) — Three-Source Heatmap

**What it shows:** A heatmap with genes as rows and three columns: AlphaGenome's predicted expression, tumour TPM, and GTEx normal lung TPM. All values are log₁₀-transformed for comparable colour scaling. Cell labels show the original (untransformed) values.

**What we see:** The AlphaGenome column is uniformly dark red (high values) because the prediction sums are in the tens to hundreds of thousands, while TPM values are in the ones to hundreds. The relative ranking within each column is more informative than the absolute colours. The heatmap makes it easy to spot genes where tumour and GTEx expression diverge (e.g., ERBB2 and FAM107A).

---

## 6. What Would We Need to See to Trust AlphaGenome?

### The Core Problem: We Have Too Little Data

With only **8 variants**, we simply cannot draw confident conclusions about AlphaGenome's usefulness. Here is what we would need to see at a larger scale:

### 6.1 A Statistically Significant Correlation

Our power analysis shows that with the observed effect size (r ≈ 0.55), we would need at least **24 variants** to achieve 80% statistical power — that is, an 80% chance of correctly detecting a real correlation at the conventional significance threshold of p < 0.05.

The most immediate path to 24+ variants is to include **MODERATE-impact variants** (such as missense mutations, which change a single amino acid) alongside the current HIGH-impact set. This would increase the pool from 8 to approximately 46 variants.

### 6.2 A Correlation That Survives Outlier Removal

Currently, removing the single ERBB2 data point drops the Pearson r from 0.548 to 0.060. A reliable model should show consistent correlations that do not depend on any single gene. With 24+ variants, no individual gene should have such outsized influence.

### 6.3 Predicted Fold-Changes That Match Reality

All 8 current variants were predicted as Neutral (essentially zero expression change). To validate AlphaGenome's ability to predict **expression changes**, we need variants where the model predicts a meaningful gain or loss, and then check whether that matches the observed expression in the patient's RNA-seq.

The ideal result: variants predicted to lose expression should have lower TPM than variants predicted to maintain expression. If we could also access matched pre- and post-mutation expression data (which is not available in this single-patient design), we could directly validate the direction and magnitude of predicted fold-changes.

### 6.4 Multi-Patient Reproducibility

Our entire analysis comes from a single patient. Biological variability between patients is enormous — a correlation that holds in one patient might vanish in another. A convincing evaluation would run the same pipeline across **multiple TCGA patients** with the same cancer type and check whether the correlation is consistent.

### 6.5 A Stronger Correlation with Tumour Expression Than with Normal Tissue

Our preliminary GTEx comparison hints at this (r = 0.548 for tumour vs. r = 0.137 for GTEx), but the difference is not significant at n = 8. If AlphaGenome is truly modelling the variant's effect on expression in the tumour context, we would expect its predictions to correlate more strongly with tumour-specific expression than with healthy-tissue baselines.

---

## 7. Limitations and Honest Assessment

1. **Prototype, not production.** This pipeline demonstrates feasibility, not clinical utility. It is a framework for evaluation, not a finished tool.

2. **Scale mismatch.** AlphaGenome sums predicted expression over a 1 MB window containing multiple genes. TPM is a per-gene measure. Comparing the two is like comparing the total noise level of an entire street to the volume of one specific shop — they are related but not directly equivalent.

3. **All predictions were Neutral.** Because AlphaGenome predicts essentially zero expression change for all 8 variants, we can only assess whether the model captures expression *magnitude* (which gene is higher or lower), not expression *change* (which mutations increase or decrease expression). The more interesting and clinically useful validation — does the model correctly predict the direction and size of expression changes — remains untested.

4. **NMD is imperfect in cancer.** The NMD flag comes from VEP's annotation of the reference transcript. In reality, NMD efficiency varies between tissues and is often impaired in cancer cells. A variant flagged as NMD-triggering might actually produce a stable truncated transcript in a particular tumour.

5. **Single patient, single tissue.** Results from TCGA-05-4384 cannot be generalised without replication.

---

## 8. Conclusion

This prototype demonstrates a complete computational pipeline for evaluating whether AlphaGenome's DNA-based expression predictions could help identify cancer vaccine targets that are likely to be poorly expressed. The pipeline takes a patient's raw somatic mutation data, filters for functionally impactful variants, queries AlphaGenome for expression predictions, enriches with biological scoring, validates against real RNA-seq data, and contextualises against healthy-tissue baselines.

From 8 high-impact variants in a single lung cancer patient, we observe a moderate positive correlation (r ≈ 0.55) between AlphaGenome's predicted expression and the patient's actual RNA-seq measurements. However, this correlation is **not statistically significant**, is **driven by a single outlier (ERBB2)**, and is based on **too few data points** to draw conclusions.

**The infrastructure works; the verdict is pending.** The next step is to expand to ≥ 24 variants (by including MODERATE-impact mutations) and ideally to replicate across multiple patients. Only then can we assess whether AlphaGenome's predictions are reliable enough to inform vaccine target prioritisation in practice.
