DATA_PATH = "data/VCF_File/ff1df4a9-2318-4dba-8f34-cb69dde4360c/TCGA_LUAD.f368273c-bd2d-4b97-97ec-a04cb130af1e.wgs.GATK4_MuTect2_Pair.somatic_annotation.vcf.gz"
EXAMPLE_RNA_PATH = "data/Example_RNA.csv"
QUICK_10 = "output/high_confidence_10.vcf"
HIGH_IMPACT_VCF = "output/high_impact_variants.vcf"
RAW_PREDICTIONS = "output/raw_predictions.tsv"
SCORED_VARIANTS = "output/scored_variants.tsv"

# VEP impact levels to retain during filtering.
# HIGH  = stop_gained, frameshift, splice_donor/acceptor, etc.
# MODERATE = missense_variant, inframe_insertion/deletion, etc.
DEFAULT_IMPACT_LEVELS = {'HIGH'}