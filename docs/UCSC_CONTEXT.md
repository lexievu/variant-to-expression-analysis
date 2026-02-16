# Why We Don't Need UCSC Genome Browser Context Retrieval

The original execution plan called for fetching reference DNA context from the UCSC Genome Browser API (hg38), generating paired reference/mutated sequences, and passing them to AlphaGenome. After implementation, this step turned out to be **unnecessary**. This document explains why.

---

## How AlphaGenome handles sequence context

AlphaGenome's Python SDK provides two key abstractions:

```python
interval = genome.Interval("chr7", 55_000_000, 56_048_576)   # 1 MB window
variant  = genome.Variant("chr7", 55_191_822, "C", "T")      # somatic SNV

# Low-level API (used internally by score_variant):
outputs = model.predict_variant(
    interval=interval,
    variant=variant,
    ontology_terms=["UBERON:0002048"],              # Lung tissue
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)
```

> **Note:** The pipeline calls `model.score_variant()` with
> `GeneMaskLFCScorer` rather than `predict_variant` directly.
> The example above illustrates what happens *internally*.

Internally, `predict_variant`:

1. **Fetches the hg38 reference** for the given `Interval` — no external FASTA or API needed.
2. **Constructs the alternate sequence** by applying the `Variant` substitution.
3. **Runs inference twice** (reference + alternate) and returns both prediction tracks.

All sequence retrieval and mutation application happens server-side. The client only sends coordinates.

---

## What UCSC retrieval would have given us

| Use case | Benefit | Verdict |
|---|---|---|
| **REF allele validation** — confirm VCF `REF` matches hg38 | Catch genome-build mismatches | **Low priority.** MuTect2 already validates REF during variant calling. A spot-check of the VCF `##reference=` header line is sufficient. |
| **Gene/exon masking** — get exon coordinates for gene-level expression scoring | Focus predictions on target gene | **Not needed.** AlphaGenome's `GeneMaskLFCScorer` applies exon masks internally. We validated these against UCSC GENCODE (Pearson r = 0.9898, 8/8 gene name & strand matches). See `docs/GENE_DILUTION.md`. |
| **Alternative model input** — feed raw DNA strings to Enformer, Sei, or other sequence→expression models | Enable multi-model comparison | **Not needed now.** Only relevant if a second model is added. AlphaGenome takes coordinates, not strings. |
| **Sequence-context visualisation** — display ±50 bp flanking each variant for manual inspection | Useful for figures or reports | **Nice-to-have**, but a notebook one-liner with `pysam.FastaFile` against a local hg38 FASTA is simpler than HTTP calls to UCSC. |
| **Genome build mismatch detection** — verify VCF is hg38, not hg19 | Prevent silent coordinate errors | **Cheaper alternatives exist.** Check the VCF header or compare a handful of REF alleles against a local FASTA. |

---

## If we ever need it

If a future requirement does arise (e.g., adding Enformer as a second model), the recommended approach is:

1. **Local FASTA** — download `hg38.fa` from UCSC or NCBI and use `pysam.FastaFile` for O(1) lookups. This is faster, offline-capable, and more reliable than HTTP API calls.
2. **REF-allele QC only** — a lightweight check can be added to `s2_vcf_filter.py` without touching the prediction pipeline:
   ```python
   import pysam
   fasta = pysam.FastaFile("hg38.fa")
   expected_ref = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref))
   assert expected_ref == ref, f"REF mismatch at {chrom}:{pos}"
   ```

Neither of these requires the UCSC REST API.

---

## Decision

**UCSC Genome Browser data is not needed for production runs.** AlphaGenome
handles all sequence operations internally, and its built-in gene masks are
functionally equivalent to UCSC GENCODE annotations (validated with Pearson
r = 0.9898 across 8 variants; see `docs/GENE_DILUTION.md`). The UCSC REST
API was useful as a one-time validation tool but is not a runtime dependency.
This item is closed in the execution plan.
