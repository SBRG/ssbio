import io
import os.path as op
from collections import defaultdict
from bioservices import KEGG
import ssbio.utils
from ssbio.sequence.seqprop import SeqProp
import logging
from slugify import Slugify

log = logging.getLogger(__name__)
custom_slugify = Slugify(safe_chars='-_')
bs_kegg = KEGG()


class KEGGProp(SeqProp):
    def __init__(self, kegg_id, sequence_path=None, metadata_path=None):

        SeqProp.__init__(self, ident=kegg_id, sequence_path=sequence_path, metadata_path=metadata_path)

        if kegg_id:
            self.kegg = kegg_id

        if metadata_path:
            self.load_metadata_path(metadata_path)

    def load_metadata_path(self, metadata_path):
        SeqProp.load_metadata_path(self, metadata_path)
        self.update(parse_kegg_gene_metadata(metadata_path), overwrite=True)

    def download_seq_file(self, outdir, force_rerun=False):
        kegg_seq_file = download_kegg_aa_seq(gene_id=self.id,
                                             outdir=outdir,
                                             force_rerun=force_rerun)
        if kegg_seq_file:
            self.load_sequence_path(kegg_seq_file)
        else:
            log.warning('{}: no sequence file available'.format(self.id))

    def download_metadata_file(self, outdir, force_rerun=False):
        kegg_metadata_file = download_kegg_gene_metadata(gene_id=self.id,
                                                         outdir=outdir,
                                                         force_rerun=force_rerun)
        if kegg_metadata_file:
            self.load_metadata_path(kegg_metadata_file)
        else:
            log.warning('{}: no metadata file available'.format(self.id))


def download_kegg_gene_metadata(gene_id, outdir=None, force_rerun=False):
    """Download the KEGG flatfile for a KEGG ID and return the path.

    Args:
        gene_id: KEGG gene ID (with organism code), i.e. "eco:1244"
        outdir: optional output directory of metadata

    Returns:
        Path to metadata file

    """
    if not outdir:
        outdir = ''

    # Replace colon with dash in the KEGG gene ID
    outfile = op.join(outdir, '{}.kegg'.format(custom_slugify(gene_id)))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        raw_text = bs_kegg.get("{}".format(gene_id))
        if raw_text == 404:
            return

        with io.open(outfile, mode='wt', encoding='utf-8') as f:
            f.write(raw_text)

        log.debug('{}: downloaded KEGG metadata file'.format(outfile))
    else:
        log.debug('{}: KEGG metadata file already exists'.format(outfile))

    return outfile


def parse_kegg_gene_metadata(infile):
    """Parse the KEGG flatfile and return a dictionary of metadata.

    Dictionary keys are:
        refseq
        uniprot
        pdbs
        taxonomy

    Args:
        infile: Path to KEGG flatfile

    Returns:
        dict: Dictionary of metadata

    """
    metadata = defaultdict(str)

    with open(infile) as mf:
        kegg_parsed = bs_kegg.parse(mf.read())

    # TODO: additional fields can be parsed

    if 'DBLINKS' in kegg_parsed.keys():
        if 'UniProt' in kegg_parsed['DBLINKS']:
            unis = str(kegg_parsed['DBLINKS']['UniProt']).split(' ')
            # TODO: losing other uniprot ids by doing this
            if isinstance(unis, list):
                metadata['uniprot'] = unis[0]
            else:
                metadata['uniprot'] = unis
        if 'NCBI-ProteinID' in kegg_parsed['DBLINKS']:
            metadata['refseq'] = str(kegg_parsed['DBLINKS']['NCBI-ProteinID'])
    if 'STRUCTURE' in kegg_parsed.keys():
        metadata['pdbs'] = str(kegg_parsed['STRUCTURE']['PDB']).split(' ')
    else:
        metadata['pdbs'] = None
    if 'ORGANISM' in kegg_parsed.keys():
        metadata['taxonomy'] = str(kegg_parsed['ORGANISM'])

    return metadata


def download_kegg_aa_seq(gene_id, outdir=None, force_rerun=False):
    """Download a FASTA sequence of a protein from the KEGG database and return the path.

    Args:
        gene_id: the gene identifier
        outdir: optional path to output directory

    Returns:
        Path to FASTA file

    """
    if not outdir:
        outdir = ''

    outfile = op.join(outdir, '{}.faa'.format(custom_slugify(gene_id)))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        raw_text = bs_kegg.get("{}".format(gene_id), option='aaseq')
        if raw_text == 404:
            return

        with io.open(outfile, mode='wt', encoding='utf-8') as f:
           f.write(raw_text)

        log.debug('{}: downloaded KEGG FASTA file'.format(outfile))
    else:
        log.debug('{}: KEGG FASTA file already exists'.format(outfile))

    return outfile


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
    mapping = bs_kegg.conv(target_db, organism_code)

    # strip the organism code from the keys and the identifier in the values
    new_mapping = {}
    for k,v in mapping.items():
        new_mapping[k.replace(organism_code + ':', '')] = str(v.split(':')[1])

    return new_mapping