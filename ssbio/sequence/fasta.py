from Bio.Alphabet import IUPAC
import os
import os.path as op
from Bio import SeqIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ssbio import utils
date = utils.Date()

# TODO: what if i want to write one fasta file for multiple sequences?


def write_fasta_file(seq_str, ident, description='',
                     extension='faa', outdir=None, overwrite=False, ignore_alphabet=False):
    '''
    This writes a fasta file for a single sequence string.
    It also checks if the file exists already and returns the filename.
    You can overwrite existing fasta files by setting overwrite=True

    Input:  seq_str (str) - amino acid string
            ident (str) - ID of the sequence
            outdir (str) - desired directory of file output
            overwrite (bool) - if you want to overwrite existing files
            ignore_alphabet (boolean): OPTIONAL check the alphabet to see if it contains valid amino acids.
    Output: Filename of fasta file
    '''

    if not description:
        description_final = '{}-ssbioSeq'.format(date.short_date)
    else:
        description_final = '{}_{}-ssbioSeq'.format(description, date.short_date)
    seq = load_seq_string(seq_str, ident=ident, desc=description_final, ignore_alphabet=ignore_alphabet)
    outfile = "{}.{}".format(ident, extension)

    if outdir:
        if not op.exists(outdir):
            os.mkdir(outdir)
        outfile = op.join(outdir, outfile)

    if op.isfile(outfile) and not overwrite:
        return outfile
    else:
        SeqIO.write(seq, outfile, "fasta")
        return outfile


def load_seq_string(seq_string, ident, name='', desc='', ignore_alphabet=False):
    """Load an amino acid sequence string.

    Args:
        seq_string (str): A protein amino acid sequence
        ident (str): Database identifier (ie. UniProt or PDB ID)
        name (str): OPTIONAL protein name (ie. gene name)
        desc (str): OPTIONAL description of this sequence (ie. catabolizes G6P)
        ignore_alphabet (boolean): OPTIONAL check the alphabet to see if it contains valid amino acids.

    Returns:
        my_seq (SeqRecord): Biopython SeqRecord object

    """
    my_seq = Seq(seq_string, Alphabet.IUPAC.extended_protein)

    # NOTE: using private _verify_alphabet method in Bio.Alphabet module
    if not Alphabet._verify_alphabet(my_seq) and not ignore_alphabet:
        raise ValueError('Sequence contains invalid characters')

    my_seq_record = SeqRecord(my_seq, id=ident, name=name, description=desc)

    return my_seq_record

def load_fasta_file(filename):
    """Load a FASTA file and returns the sequences as a list of SeqRecords

    Args:
        filename (str): Path to the FASTA file to load

    Returns:
        records (list): list of all sequences in the FASTA file as
        Biopython SeqRecord objects
    """
    if os.path.exists(filename):
        handle = open(filename, "r")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    else:
        raise IOError('File does not exist.')

def load_seq_str_as_temp_file(self, seq_str):
    """Load a sequence string and return a temporary path to the FASTA file.

    Args:
        seq_str:

    Returns:
        str: File path to FASTA file in temporary directory

    """
    pass


if __name__ == '__main__':
    print(write_fasta_file('ALALLALAL', ident='mypdb', outdir='/tmp/', overwrite=True))
