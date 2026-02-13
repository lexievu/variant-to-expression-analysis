import argparse
import logging
import pandas as pd
from cyvcf2 import VCF, Writer
import utils
from constants import DATA_PATH, EXAMPLE_RNA_PATH, DEFAULT_IMPACT_LEVELS

# --- 1. DEFAULTS ---
DEFAULT_OUTPUT = 'output/high_impact_variants.vcf'
LOG_FILENAME = 'log/vcf_filter.log'


def parse_args(argv=None):
    """Parse command-line arguments.

    Args:
        argv: Argument list (defaults to sys.argv[1:] when None).

    Returns:
        argparse.Namespace with vcf, rna, output, impact attributes.
    """
    parser = argparse.ArgumentParser(
        description='Filter a VCF for somatic, high-impact variants expressed in RNA-seq.'
    )
    parser.add_argument(
        '--vcf', default=DATA_PATH,
        help='Path to the input VCF file.'
    )
    parser.add_argument(
        '--rna', default=EXAMPLE_RNA_PATH,
        help='Path to the RNA-seq CSV with gene_id column.'
    )
    parser.add_argument(
        '--output', '-o', default=DEFAULT_OUTPUT,
        help=f'Output VCF path (default: {DEFAULT_OUTPUT}).'
    )
    parser.add_argument(
        '--impact', default=None,
        help=(
            'Comma-separated VEP impact levels to keep. '
            'E.g. "HIGH" or "HIGH,MODERATE". '
            f'Default: {",".join(sorted(DEFAULT_IMPACT_LEVELS))}.'
        )
    )
    args = parser.parse_args(argv)

    # Resolve impact set
    if args.impact is None:
        args.impact = DEFAULT_IMPACT_LEVELS
    else:
        args.impact = {s.strip().upper() for s in args.impact.split(',')}

    return args


def load_rna_gene_ids(rna_file):
    """Load RNA file and return a set of Ensembl IDs (stripped of version)."""
    logging.info(f"Loading RNA data from {rna_file}...")
    try:
        df = pd.read_csv(rna_file, comment='#')
        valid_genes = set(df['gene_id'].apply(lambda x: str(x).split('.')[0]))
        logging.info(f"Loaded {len(valid_genes)} genes from RNA file.")
        return valid_genes
    except Exception as e:
        logging.error(f"Failed to load RNA file: {e}")
        return set()


def parse_csq_field(variant, valid_genes, target_impacts):
    """Return CSQ transcripts matching *target_impacts* and *valid_genes*."""
    candidates = []
    for entry in utils.parse_csq(variant):
        if entry['impact'] in target_impacts and entry['gene_id_stripped'] in valid_genes:
            candidates.append({
                'gene_id': entry['gene_id_stripped'],
                'impact': entry['impact'],
                'consequence': entry['consequence'],
            })
    return candidates


def filter_vcf(vcf_file, valid_genes, output_file, target_impacts):
    """Filter *vcf_file* and write passing variants to *output_file*."""
    logging.info(f"Starting VCF filtering on {vcf_file}...")
    logging.info(f"Impact levels: {', '.join(sorted(target_impacts))}")
    count_kept = 0
    count_total = 0

    try:
        vcf = VCF(vcf_file)
        w = Writer(output_file, vcf)
        tumor_idx = vcf.samples.index('TUMOR')

        for variant in vcf:
            count_total += 1

            # --- FILTER 1: PASS only ---
            if variant.FILTER and variant.FILTER != 'PASS':
                continue

            # --- FILTER 2: Tumor carries at least one ALT allele ---
            gt = variant.genotypes[tumor_idx]
            if gt[0] <= 0 and gt[1] <= 0:
                continue

            # --- FILTER 3: CSQ impact + expressed gene ---
            hits = parse_csq_field(variant, valid_genes, target_impacts)
            if hits:
                w.write_record(variant)
                count_kept += 1
                if count_kept % 100 == 0:
                    logging.info(f"Kept {count_kept} variants so far...")

        w.close()
        vcf.close()
        logging.info(
            f"Finished. Scanned {count_total} variants. "
            f"Wrote {count_kept} to {output_file}."
        )

    except Exception as e:
        logging.error(f"Error processing VCF: {e}")


# --- ENTRY POINT ---
if __name__ == "__main__":
    args = parse_args()
    utils.setup_logging(LOG_FILENAME)
    logging.info(f"Config: impacts={sorted(args.impact)}, output={args.output}")
    genes = load_rna_gene_ids(args.rna)
    if genes:
        filter_vcf(args.vcf, genes, args.output, args.impact)