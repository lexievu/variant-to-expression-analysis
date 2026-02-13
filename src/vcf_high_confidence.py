from cyvcf2 import VCF, Writer
from global_variables import *

# 1. Load the VCF
vcf_input = VCF(data_path)

MIN_TLOD = 6.3  # Standard Mutect2 threshold
MIN_DEPTH = 20

count = 0

print("Filtering and writing to 'high_confidence.vcf'...")

with Writer('output/high_confidence_10.vcf', vcf_input) as vcf_output:
    for variant in vcf_input:
        
        # --- Filter Logic ---
        
        # 1. Must be PASS
        if variant.FILTER is not None:
            continue
            
        # 2. Must be Somatic (High TLOD)
        # Safely get TLOD, handling cases where it might be a list or None
        tlod = variant.INFO.get('TLOD')
        if isinstance(tlod, (list, tuple)):
            tlod = tlod[0]
            
        if tlod is None or tlod < MIN_TLOD:
            continue

        # 3. Must have good depth in Tumor
        # (Assuming Tumor is the 2nd sample, index 1. Adjust if needed.)
        tumor_idx = 1 
        tumor_depth = variant.format('DP')[tumor_idx][0]
        
        if tumor_depth < MIN_DEPTH:
            continue

        # --- Write the Variant ---
        vcf_output.write_record(variant)
        count += 1

        # Only save 10 variants here
        if (count == 10):
            break

# 4. Close the files to ensure data is flushed to disk
vcf_output.close()
vcf_input.close()

print(f"Success! Saved {count} variants.")