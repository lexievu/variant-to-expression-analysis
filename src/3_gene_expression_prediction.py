import logging

import numpy as np
from cyvcf2 import VCF
from alphagenome.models import dna_client
from alphagenome.data import genome
import os
from dotenv import load_dotenv
from constants import DATA_PATH, EXAMPLE_RNA_PATH, QUICK_10
import utils

# Load variables from .env file into the environment
load_dotenv()

# --- 1. CONFIGURATION ---
LOG_FILENAME = "log/gene_expression_prediction.log"
API_KEY = os.getenv("ALPHAGENOME_API_KEY")
VCF_FILE = QUICK_10
OUTPUT_FILENAME = "output/alphagenome_hits.tsv"
TISSUE_ID = 'UBERON:0002048'  # Lung
DNA_SEQUENCE_LENGTH = 1_048_576

utils.setup_logging(LOG_FILENAME)

# --- 2. SETUP OUPUT FILE ---
# We open the file in 'w' (write) mode
outfile = open(OUTPUT_FILENAME, 'w')

# Write the Header Line
# \t means "tab", \n means "new line"
header = "CHROM\tPOS\tREF\tALT\tGENE\tGENE_ID\tREF_EXPR\tALT_EXPR\tLOG2_FC\tSTATUS\n"
outfile.write(header)

# --- 3. INITIALIZE ALPHAGENOME ---
logging.info(f"Connecting to AlphaGenome and writing to '{OUTPUT_FILENAME}'...")
model = dna_client.create(API_KEY)
vcf = VCF(VCF_FILE)

count_saved = 0


for variant in vcf:
    # (Optional) Skip variants that are not PASS if you haven't filtered yet
    if variant.FILTER is not None:
        continue

    # --- 4. PREPARE INPUTS ---
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alt = variant.ALT[0]
    
    gene_name = utils.get_gene_name(variant)
    gene_id = utils.get_gene_id(variant)  # Ensembl ID, version-stripped

    # Define Variant & Interval for AlphaGenome
    ag_variant = genome.Variant(chrom, pos, ref, alt)
    start = max(1, pos - DNA_SEQUENCE_LENGTH // 2)
    end = pos + DNA_SEQUENCE_LENGTH // 2
    interval = genome.Interval(chrom, start, end)

    try:
        # --- 5. PREDICT ---
        outputs = model.predict_variant(
            interval=interval,
            variant=ag_variant,
            ontology_terms=[TISSUE_ID],
            requested_outputs=[dna_client.OutputType.RNA_SEQ]
        )
        
        # --- 6. CALCULATE SCORES ---
        ref_track = outputs.reference.rna_seq.values
        alt_track = outputs.alternate.rna_seq.values
        
        ref_sum = np.sum(ref_track)
        alt_sum = np.sum(alt_track)
        
        # Calculate Fold Change (Log2)
        # Avoid division by zero
        if ref_sum == 0 and alt_sum == 0:
            log2_fc = 0.0
        else:
            log2_fc = np.log2((alt_sum + 1e-9) / (ref_sum + 1e-9))
        
        # --- 7. FILTER & SAVE "INTERESTING" HITS ---
        # Let's define "Interesting" as:
        # 1. Fold change > 1.0 (Doubling of expression)
        # 2. Fold change < -1.0 (Halving of expression)
        
        status = "Neutral"
        
        if log2_fc > 1.0:
            status = "Gain_of_Expression"
        elif log2_fc < -1.0:
            status = "Loss_of_Expression"

        # Format the line for the file
        # .4f means "4 decimal places"
        line = f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_name}\t{gene_id}\t{ref_sum:.2f}\t{alt_sum:.2f}\t{log2_fc:.4f}\t{status}\n"
        outfile.write(line)
        count_saved += 1
        
        # Optional: Print a progress indicator to terminal
        logging.info(f"Saved hit: {gene_name} ({status})")

    except Exception as e:
        logging.error(f"Error at {chrom}:{pos} -> {e}")

# --- 8. CLEANUP ---
outfile.close()
vcf.close()
logging.info(f"Done! Saved {count_saved} interesting variants to {OUTPUT_FILENAME}")