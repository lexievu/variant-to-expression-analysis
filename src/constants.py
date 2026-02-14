from pathlib import Path

# Project root — two levels up from this file (src/constants.py → project/).
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_PATH = str(PROJECT_ROOT / "data" / "VCF_File" / "ff1df4a9-2318-4dba-8f34-cb69dde4360c" / "TCGA_LUAD.f368273c-bd2d-4b97-97ec-a04cb130af1e.wgs.GATK4_MuTect2_Pair.somatic_annotation.vcf.gz")
EXAMPLE_RNA_PATH = str(PROJECT_ROOT / "data" / "Example_RNA.csv")
QUICK_10 = str(PROJECT_ROOT / "output" / "high_confidence_10.vcf")
HIGH_IMPACT_VCF = str(PROJECT_ROOT / "output" / "high_impact_variants.vcf")
RAW_PREDICTIONS = str(PROJECT_ROOT / "output" / "raw_predictions.tsv")
SCORED_VARIANTS = str(PROJECT_ROOT / "output" / "scored_variants.tsv")
VALIDATION_TABLE = str(PROJECT_ROOT / "output" / "validation_table.csv")
VALIDATION_CORRELATIONS = str(PROJECT_ROOT / "output" / "validation_correlations.csv")
GTEX_COMPARISON = str(PROJECT_ROOT / "output" / "gtex_comparison.csv")
DOTENV_PATH = str(PROJECT_ROOT / ".env")
LOG_DIR = PROJECT_ROOT / "log"
OUTPUT_DIR = PROJECT_ROOT / "output"

# VEP impact levels to retain during filtering.
# HIGH  = stop_gained, frameshift, splice_donor/acceptor, etc.
# MODERATE = missense_variant, inframe_insertion/deletion, etc.
DEFAULT_IMPACT_LEVELS = {'HIGH'}