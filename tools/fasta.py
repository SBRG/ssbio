from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os

def write_fasta_file(seq_str, ident, outpath=None, overwrite=False):
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
    seq = SeqRecord(Seq(seq_str, IUPAC.protein), id=ident,
                    description='ProteinSequence')
    outfile = "%s.faa" % ident

    if outpath:
        outfile = os.path.join(outpath, outfile)

    if os.path.isfile(outfile) and not overwrite:
        return outfile
    else:
        SeqIO.write(seq, outfile, "fasta")
        return outfile

if __name__ == '__main__':
    print(write_fasta_file('ALALLALAL', ident='mypdb', outpath='/tmp/', overwrite=True))
