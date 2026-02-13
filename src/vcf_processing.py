from cyvcf2 import VCF
from global_variables import *
import utils

HIGH_EXPRESSION_GENE = ['SFTPB', 'CTNNB1', 'MUC5AC']
gene_names = set()

# 1. Setup Input 
vcf = VCF(data_path)

def write_to_text_file():
    # 2. Map sample names
    samples = vcf.samples
    normal_idx = samples.index('NORMAL')
    tumor_idx = samples.index('TUMOR')

    # 3. Write a Header Line
    # This makes it easy to read in Excel later
    output_file = open('output/somatic_variants.txt', 'w')
    header = "CHROM\tPOS\tREF\tALT\tNormal_GT\tTumor_GT\tType\tChange\n"
    output_file.write(header)

    count = 0

    for variant in vcf:
        # --- Logic from previous step ---
        gts = variant.genotypes
        normal_gt = gts[normal_idx]
        tumor_gt = gts[tumor_idx]
        
        # Normal is HomRef (0/0)
        is_normal_ref = (normal_gt[0] == 0) and (normal_gt[1] == 0)
        
        # Tumor is Variant (Not 0/0 and not Missing)
        is_tumor_variant = (tumor_gt[0] > 0) or (tumor_gt[1] > 0)
        
        if is_normal_ref and is_tumor_variant:
            # --- Format Data for Writing ---
            
            # Convert genotype lists [0,0,True] to strings "0/0" for readability
            # (This is optional, but makes it look standard)
            n_gt_str = f"{normal_gt[0]}/{normal_gt[1]}"
            t_gt_str = f"{tumor_gt[0]}/{tumor_gt[1]}"
            
            # Handle multiple ALTs safely by joining them (e.g., "A,T")
            alt_str = ",".join(variant.ALT)
            
            # Create the line
            line = f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{alt_str}\t{n_gt_str}\t{t_gt_str}\n"
            
            gene_name = utils.get_gene_name(variant)
            gene_names.add(gene_name)
            if (gene_name in HIGH_EXPRESSION_GENE):
                print(f"{gene_name}\t{line}")

            output_file.write(line)
            count += 1

    output_file.close()
    with open('output/gene_names.txt', 'w') as gene_name_file:
        gene_name_file.write('\n'.join(str(i) for i in gene_names))
    print(f"Done! Written {count} variants to 'output/somatic_variants.txt'")

write_to_text_file()

