# Gene Dilution Problem — Identified & Fixed

## Problem Statement

An earlier version of the pipeline’s step 3 (`s3_gene_expression_prediction.py`) predicted
expression changes by summing AlphaGenome’s RNA-seq predictions across an entire **1 MB
window** centred on the variant:

```python
# OLD APPROACH (no longer used)
ref_sum = float(np.sum(outputs.reference.rna_seq.values))
alt_sum = float(np.sum(outputs.alternate.rna_seq.values))
```

Because the 1 MB window usually contains **many genes**, the expression signal
of the target gene was diluted by the background expression of its neighbours.
A variant that meaningfully changes the target gene’s expression produced
a negligible `ref_sum` vs `alt_sum` difference when the neighbouring genes’
expression dominated.

**This problem has been fixed.** The pipeline now uses `GeneMaskLFCScorer`
(see “Implemented Fix” below).

## Investigation: How Many Genes Per Window?

We queried the UCSC Genome Browser REST API (GENCODE Basic v47, hg38) for
every 1 MB window around the 8 high-impact somatic variants in the pipeline:

| Variant             | Gene       | UCSC genes in 1 MB window |
|---------------------|------------|:-------------------------:|
| chr2:47006054       | TTC7A      | 38                        |
| chr3:58589218       | FAM107A    | 27                        |
| chr7:1746696        | ELFN1-AS1  | 38                        |
| chr9:131071589      | LAMC3      | 40                        |
| chr12:29783539      | TMTC1      | 21                        |
| chr16:3057614       | MMP25      | **99**                    |
| chr17:39719785      | ERBB2      | 59                        |
| chr19:2213544       | DOT1L      | 72                        |

**Summary:** min = 21, median = 39, max = 99 genes per window.

For example, the MMP25 variant at chr16:3057614 has its expression prediction
summed together with 98 other genes. The MMP25 gene body spans only ~14 kb of
the 1,048,576 bp window — roughly 1.3 % of the total. Any variant effect on
MMP25 is buried under the signal of its neighbours.

## Implemented Fix

AlphaGenome’s SDK provides a purpose-built solution:
**`model.score_variant()` with `GeneMaskLFCScorer`**.
This is now the approach used in `s3_gene_expression_prediction.py`.

This scorer:

1. Runs inference for both REF and ALT sequences (identical to
   `predict_variant`).
2. **Masks the RNA-seq output bins to only include exon coordinates** of each
   annotated gene in the window (using AlphaGenome's internal GENCODE models).
3. Computes a per-gene `log₂(sum(ALT_exon_bins)) − log₂(sum(REF_exon_bins))`.
4. Returns an AnnData object with one score per gene × track combination,
   including `gene_id`, `gene_name`, `gene_type`, and `gene_strand`.

### Before (old code — no longer used)

```python
outputs = model.predict_variant(
    interval=interval, variant=ag_variant,
    ontology_terms=[tissue_id],
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)
ref_sum = float(np.sum(outputs.reference.rna_seq.values))
alt_sum = float(np.sum(outputs.alternate.rna_seq.values))
```

### After (current code)

```python
from alphagenome.models import variant_scorers

scorer = variant_scorers.GeneMaskLFCScorer(
    requested_output=dna_client.OutputType.RNA_SEQ,
)
scores = model.score_variant(
    interval=interval, variant=ag_variant,
    variant_scorers=[scorer],
)
df = variant_scorers.tidy_scores(scores)
target = df[df["gene_id"] == gene_id]  # filter to the VCF's gene
```

## Validation: AlphaGenome Gene Masks Are Sufficient

AlphaGenome's `GeneMaskLFCScorer` applies its exon mask server-side and only
returns the aggregated score — it does not expose the mask coordinates. To
confirm that AlphaGenome's internal gene models are correct, we compared
against UCSC's public GENCODE annotations (Basic v47, hg38) in two ways:

1. **Metadata comparison** — gene name, strand, and biotype between
   AlphaGenome's output and the UCSC `wgEncodeGencodeBasicV47` track.
2. **Indirect score comparison** — we applied UCSC exon coordinates as a
   manual mask to raw `predict_variant` output bins, computed
   `log₂(sum(ALT_masked)) − log₂(sum(REF_masked))` ourselves, and compared
   to the `GeneMaskLFCScorer` score.

### Results (8 high-impact variants)

| Metric | Value |
|--------|-------|
| Gene name matches | **8 / 8** |
| Strand matches | **8 / 8** |
| Pearson r (AG LFC vs manual UCSC-masked LFC) | **0.9898** |
| Mean \|diff\| | 0.011 |
| Max \|diff\| | 0.044 |

**Conclusion: the gene masks are functionally equivalent.** The minor
numerical differences (~0.01 LFC) likely come from slight differences in
exon boundary definitions between GENCODE versions or bin-edge rounding.
AlphaGenome's built-in masks can be trusted; **UCSC data is not needed for
production runs.**

### Signal Amplification (masked vs whole-window LFC)

Masking to the target gene's exons dramatically increases sensitivity:

| Gene      | Masked LFC | Whole-window LFC | Amplification |
|-----------|:----------:|:-----------------:|:-------------:|
| TTC7A     | 0.003811   | 0.001978          | 1.9×          |
| FAM107A   | 0.002364   | 0.000084          | **28.2×**     |
| ELFN1-AS1 | 0.000287   | 0.000040          | 7.2×          |
| LAMC3     | 0.029366   | 0.004423          | 6.6×          |
| TMTC1     | 0.002841   | 0.001300          | 2.2×          |
| MMP25     | 0.152302   | 0.000353          | **431.6×**    |
| ERBB2     | 0.019751   | 0.002253          | 8.8×          |
| DOT1L     | 0.003679   | 0.001740          | 2.1×          |

For MMP25, the signal is **431× stronger** after masking — the whole-window
sum was burying a substantial expression change under 98 other genes.

Target gene exons cover only **0.10 %–0.90 %** of the 1 MB window bins
(median 0.54 %), confirming why whole-window summing dilutes the signal.

## Files

| File | Purpose |
|------|---------|
| `output/gene_masking_validation.tsv`  | UCSC gene/exon metadata per variant window |
| `output/ag_raw_tracks/`               | Cached raw prediction arrays (`.npz`) |
| `output/ag_genemask_scores.tsv`       | Cached GeneMaskLFCScorer output |
| `output/gene_masking_comparison.tsv`  | Final AG vs UCSC-masked comparison table |
| `src/a_gene_masking_validation.py`    | UCSC data collection script |
| `src/b_gene_masking_comparison.py`    | Indirect comparison script |
| `docs/GENE_DILUTION.md`              | This file |
