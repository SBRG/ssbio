from Bio.Alphabet import IUPAC
import os
from Bio import SeqIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ssbio import utils
date = utils.Date()

def write_fasta_file(seq_str, ident, description='{}-ssbioSeq'.format(date.short_date),
                     extension='faa', outpath=None, overwrite=False):
    '''
    This writes a fasta file for a SeqRecord object. It also checks if the file
    exists already and returns the filename.
    You can overwrite existing fasta files by setting overwrite=True

    Input:  seq_str (str) - amino acid string
            ident (str) - ID of the sequence
            outpath (str) - desired PATH of file output
            overwrite (bool) - if you want to overwrite existing files
    Output: Filename of fasta file
    '''
    seq = load_seq_string(seq_str, ident=ident, desc=description)
    outfile = "{}.{}".format(ident, extension)

    if outpath:
        outfile = os.path.join(outpath, outfile)

    if os.path.isfile(outfile) and not overwrite:
        return outfile
    else:
        SeqIO.write(seq, outfile, "fasta")
        return outfile


def load_seq_string(seq_string, ident, name='', desc=''):
    """Load an amino acid sequence string.

    Args:
        seq_string (str): A protein's amino acid sequence
        ident (str): Database identifier (ie. UniProt or PDB ID)
        name (str): OPTIONAL protein name (ie. gene name)
        desc (str): OPTIONAL description of this sequence (ie. catabolizes G6P)

    Returns:
        my_seq (SeqRecord): Biopython SeqRecord object

    """
    my_seq = Seq(seq_string, Alphabet.IUPAC.extended_protein)

    # NOTE: using private _verify_alphabet method in Bio.Alphabet module
    if not Alphabet._verify_alphabet(my_seq):
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




if __name__ == '__main__':
    print(write_fasta_file('ALALLALAL', ident='mypdb', outpath='/tmp/', overwrite=True))
