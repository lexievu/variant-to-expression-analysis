import gzip
import csv
import logging
import pandas as pd
from global_variables import *
from cyvcf2 import VCF, Writer

# --- 1. CONFIGURATION ---
VCF_FILENAME = data_path 
RNA_FILENAME = example_rna_path
OUTPUT_FILENAME = 'output/high_impact_variants.vcf'
LOG_FILENAME = 'log/vcf_filter.log'
TUMOR_SAMPLE_ID = 'TCGA-05-4384-01A-01D-1751-02'

# VEP Consequence Levels to keep (Start with HIGH for "poorly expressed" targets)
TARGET_IMPACTS = {'HIGH'} # Add 'MODERATE' if you want missense variants later

# --- 2. LOGGING SETUP ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILENAME, mode='w'), # Writes to file
        logging.StreamHandler()                      # Optional: also prints to console
    ]
)

def load_rna_gene_ids(rna_file):
    """
    Loads RNA file and returns a set of Ensembl IDs (stripped of version).
    """
    logging.info(f"Loading RNA data from {rna_file}...")
    try:
        # Load CSV (assuming the file uploaded is CSV format based on extension)
        df = pd.read_csv(rna_file, comment='#')
        
        # Extract Gene IDs and strip version (e.g. ENSG00000000003.15 -> ENSG00000000003)
        # We assume the column name is 'gene_id' based on the file snippet
        valid_genes = set(df['gene_id'].apply(lambda x: str(x).split('.')[0]))
        
        logging.info(f"Loaded {len(valid_genes)} genes from RNA file.")
        return valid_genes
    except Exception as e:
        logging.error(f"Failed to load RNA file: {e}")
        return set()

def parse_csq_field(csq_string, valid_genes):
    """
    Parses the VEP CSQ string to find relevant transcripts.
    CSQ Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...
    Indices: Allele=0, Consequence=1, IMPACT=2, Gene=4
    """
    candidates = []
    # CSQ can have multiple transcripts separated by comma
    transcripts = csq_string.split(',')
    
    for tr in transcripts:
        fields = tr.split('|')
        if len(fields) < 5:
            continue
            
        consequence = fields[1]
        impact = fields[2]
        gene_id = fields[4]
        
        # Check if this transcript meets our criteria
        if impact in TARGET_IMPACTS and gene_id in valid_genes:
            candidates.append({
                'gene_id': gene_id,
                'impact': impact,
                'consequence': consequence
            })
            
    return candidates

def filter_vcf(vcf_file, valid_genes):
    logging.info(f"Starting VCF filtering on {vcf_file}...")
    count_kept = 0
    count_total = 0
    
    try:
        vcf = VCF(vcf_file)
        # 1. Setup Output Writer (copies header from input)
        w = Writer(OUTPUT_FILENAME, vcf)
        
        # 2. Find Tumor Sample Index
        tumor_idx = vcf.samples.index('TUMOR')

        for variant in vcf:
            count_total += 1
            # --- FILTER 1: PASS only ---
            # cyvcf2 'FILTER' is None if PASS (or "." in older VCFs), or the string "PASS"
            if variant.FILTER and variant.FILTER != 'PASS':
                continue

            # --- FILTER 2: Tumor Genotype Check ---
            # variant.genotypes is a list of [allele1, allele2, phased_bool] per sample
            # We check the tumor sample. If it's 0/0 (hom-ref) or -1/-1 (missing), skip.
            gt = variant.genotypes[tumor_idx]
            # gt[0] and gt[1] are the allele indices. 0=Ref, 1=Alt1, etc. -1=Missing.
            # We keep if at least one allele is > 0 (contains an ALT)
            if gt[0] <= 0 and gt[1] <= 0: 
                continue

            # --- FILTER 3: CSQ Parsing (VEP) ---
            csq_str = variant.INFO.get('CSQ')
            if not csq_str:
                continue

            # CSQ format: Allele|Consequence|IMPACT|SYMBOL|Gene|...
            # Indices: Consequence=1, IMPACT=2, Gene=4
            keep_variant = False
            
            for transcript in csq_str.split(','):
                fields = transcript.split('|')
                if len(fields) < 5: 
                    continue
                
                impact = fields[2]
                gene_id = fields[4]

                # Logic: Must be (HIGH Impact OR Splice Variant) AND (Gene in RNA file)
                is_high_impact = (impact in TARGET_IMPACTS)
                
                if is_high_impact and (gene_id in valid_genes):
                    keep_variant = True
                    break # Found one valid transcript, keep the whole variant line
            
            if keep_variant:
                w.write_record(variant)
                count_kept += 1
                if count_kept % 100 == 0:
                    logging.info(f"Kept {count_kept} variants so far...")

        w.close()
        vcf.close()
        logging.info(f"Finished. Scanned {count_total} variants. Wrote {count_kept} to {OUTPUT_FILENAME}.")

    except Exception as e:
        logging.error(f"Error processing VCF: {e}")

# --- 3. EXECUTION ---
if __name__ == "__main__":
    genes = load_rna_gene_ids(RNA_FILENAME)
    if genes:
        filter_vcf(VCF_FILENAME, genes)