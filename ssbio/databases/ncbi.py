import os
import os.path as op

from Bio import Entrez
from ssbio import utils

date = utils.Date()

import logging
log = logging.getLogger(__name__)


def download_genome_sequence(genome_accession_or_id, seqtype, email, outdir='', outfile='', force_rerun=False):
    """Download the entire set of DNA or protein sequences from protein-encoding genes in a genome from NCBI.

    Saves a FASTA file in the optional directory specified.

    Args:
        genome_accession_or_id (str): RefSeq complete genome ID (e.g. "NC_000913") or GenBank ID (e.g. "U00096.3")
        seqtype (str): "nucl" or "prot" - if you want the coding sequences in amino acid or DNA format.
        email (str): mandatory email so NCBI knows who is accessing the data_dir
        outdir (str): optional output directory (default is the current directory)

    Returns:
        Path to downloaded FASTA file.

    """
    if seqtype == 'nucl':
        extension = 'fna'
        rettype = "fasta_cds_na"
    elif seqtype == 'prot':
        extension = 'faa'
        rettype = "fasta_cds_aa"
    else:
        raise ValueError('seqtype must be "nucl" or "prot"')

    # path and filename parsing
    if outfile:
        outfile = op.join(outdir, '{}.{}'.format(outfile, extension))
    else:
        # if no outfile is specified, default is "$GI.faa"
        outfile = op.join(outdir, '{}.{}'.format(genome_accession_or_id, extension))

    if not force_rerun:
        # return the path to the file if it was already downloaded
        if op.exists(outfile):
            log.debug('FASTA file already exists at {}'.format(outfile))
            return outfile

    Entrez.email = email

    # convert to NCBI GI
    convert_to_gi = Entrez.read(Entrez.esearch(db="nucleotide", term=genome_accession_or_id, retmode="xml"))
    gi = convert_to_gi['IdList']

    # download only the CDS protein FASTA
    # See: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
    record = Entrez.efetch(db="nucleotide", id=gi, rettype=rettype, retmode="fasta")

    # save the FASTA file
    with open(outfile, 'wt') as f:
        f.write(record.read())

    log.info('Saved FASTA file at {}'.format(outfile))
    return outfile