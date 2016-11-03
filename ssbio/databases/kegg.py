import os.path as op
SEVEN_DAYS = 60 * 60 * 24 * 7
# import cachetools
from bioservices import KEGG
kegg = KEGG()
import io


def download_kegg_gene_metadata(organism_code, gene_id, outdir='', force_rerun=False):
    """Download the KEGG flatfile for a KEGG ID and return the path.

    Args:
        organism_code: KEGG organism three letter code
        gene_id: KEGG gene ID
        outdir: optional output directory of metadata

    Returns:
        Path to metadata file

    """
    outfile = op.join(outdir, '{}-{}.kegg'.format(organism_code, gene_id))
    if not op.exists(outfile) and not force_rerun:

        raw_text = kegg.get("{}:{}".format(organism_code, gene_id))
        if raw_text == 404:
            return

        with io.open(outfile, mode='wt', encoding='utf-8') as f:
            f.write(raw_text)

    return outfile


def download_kegg_aa_seq(organism_code, gene_id, outdir='', force_rerun=False):
    """Download a FASTA sequence of a protein from the KEGG database and return the path.

    Args:
        organism_code: the three letter KEGG code of your organism
        gene_id: the gene identifier
        outdir: optional path to output directory

    Returns:
        Path to FASTA file

    """
    outfile = op.join(outdir, '{}-{}.faa'.format(organism_code, gene_id))
    if not op.exists(outfile) and not force_rerun:

        raw_text = kegg.get("{}:{}".format(organism_code, gene_id), option='aaseq')
        if raw_text == 404:
            return

        with io.open(outfile, mode='wt', encoding='utf-8') as f:
           f.write(raw_text)

    return outfile


# @cachetools.func.ttl_cache(maxsize=800, ttl=SEVEN_DAYS)
def map_kegg_all_genes(organism_code, target_db):
    """Map all of an organism's gene IDs to the target database.

    This is faster than supplying a specific list of genes to map,
    plus there seems to be a limit on the number you can map with a manual REST query anyway.

    Args:
        organism_code: the three letter KEGG code of your organism
        target_db: ncbi-proteinid | ncbi-geneid | uniprot

    Returns:
        Dictionary of ID mapping

    """
    mapping = kegg.conv(target_db, organism_code)

    # strip the organism code from the keys and the identifier in the values
    new_mapping = {}
    for k,v in mapping.items():
        new_mapping[k.replace(organism_code + ':', '')] = str(v.split(':')[1])

    return new_mapping