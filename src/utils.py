def get_gene_name(variant):
    # Get Gene Name from VCF annotation (if available)
    # This assumes the CSQ format from your previous VCF header
    csq = variant.INFO.get('CSQ')
    gene_name = "."
    if csq:
        try:
            gene_name = csq.split('|')[3] # Index 3 is usually SYMBOL
        except IndexError:
            pass

    return gene_name