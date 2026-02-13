"""
Shared utilities for parsing VEP CSQ annotations from VCF variants.

CSQ field layout (pipe-delimited, comma-separated transcripts):
  Index 0: Allele
  Index 1: Consequence
  Index 2: IMPACT
  Index 3: SYMBOL  (gene name)
  Index 4: Gene    (Ensembl gene ID, e.g. ENSG00000000003.15)
"""

import logging


def setup_logging(log_file, level=logging.INFO):
    """Configure structured logging to *log_file* and the console.

    Safe to call multiple times — handlers are only added once.
    """
    logger = logging.getLogger()
    if not logger.handlers:
        logger.setLevel(level)
        fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh = logging.FileHandler(log_file, mode='w')
        fh.setFormatter(fmt)
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.addHandler(sh)


def get_gene_name(variant):
    """Return the gene symbol (SYMBOL, CSQ index 3) for a variant."""
    csq = variant.INFO.get('CSQ')
    if csq:
        try:
            return csq.split('|')[3]
        except IndexError:
            pass
    return "."


def get_gene_id(variant, strip_version=True):
    """Return the Ensembl gene ID (CSQ index 4) for a variant.

    Args:
        variant: A cyvcf2 Variant object with VEP CSQ annotations.
        strip_version: If True, strip the version suffix
                       (e.g. ENSG00000000003.15 → ENSG00000000003).

    Returns:
        The gene ID string, or "." if unavailable.
    """
    csq = variant.INFO.get('CSQ')
    if csq:
        try:
            gene_id = csq.split('|')[4]
            if strip_version and '.' in gene_id:
                gene_id = gene_id.split('.')[0]
            return gene_id
        except IndexError:
            pass
    return "."


def parse_csq(variant):
    """Parse all CSQ transcripts for a variant into a list of dicts.

    Returns a list of dicts, one per transcript annotation, with keys:
        allele, consequence, impact, gene_name, gene_id, gene_id_stripped
    Only transcripts with at least 5 fields are included.
    """
    csq = variant.INFO.get('CSQ')
    if not csq:
        return []

    results = []
    for transcript in csq.split(','):
        fields = transcript.split('|')
        if len(fields) < 5:
            continue
        gene_id_raw = fields[4]
        results.append({
            'allele': fields[0],
            'consequence': fields[1],
            'impact': fields[2],
            'gene_name': fields[3],
            'gene_id': gene_id_raw,
            'gene_id_stripped': gene_id_raw.split('.')[0] if '.' in gene_id_raw else gene_id_raw,
        })
    return results