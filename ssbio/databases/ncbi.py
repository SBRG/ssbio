import os
import os.path as op

from Bio import Entrez
from ssbio import utils

date = utils.Date()

# def entrez_search

def download_all_protein_sequences(genome_accession_or_id, email, outdir='', outfile=''):
    """Download the entire set of amino acid sequences from protein-encoding genes in a genome from NCBI.

    Saves a FASTA file in the optional directory specified.

    Args:
        genome_accession_or_id: RefSeq complete genome ID (e.g. "NC_000913") or NCBI GI (e.g. "556503834")
        email: mandatory email so NCBI knows who is accessing the data
        outdir: optional output directory (default is the current directory)

    Returns:
        Path to downloaded FASTA file.

    """

    # path and filename parsing
    if outfile:
        outfile = op.join(outdir, '{}.faa'.format(outfile))
    else:
        # if no outfile is specified, default is "$DATE-$GI.faa"
        outfile = op.join(outdir, '{}-{}.faa'.format(date.short_date, genome_accession_or_id))

    # return the path to the file if it was already downloaded
    if op.exists(outfile):
        return outfile

    Entrez.email = email

    # convert to NCBI GI
    convert_to_gi = Entrez.read(Entrez.esearch(db="nucleotide", term=genome_accession_or_id, retmode="xml"))
    gi = convert_to_gi['IdList']

    # download only the CDS protein FASTA (https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/)
    record = Entrez.efetch(db="nucleotide", id=gi, rettype="fasta_cds_aa", retmode="fasta")

    # save the FASTA file
    with open(outfile, 'wt') as f:
        f.write(record.read())

    return outfile