import os
from Bio import SeqIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SeqTools:
    """Various tools to load and write sequence files
    """

    def __init__(self):
        pass

    def load_seq_string(self, seq_string, ident, name='', desc=''):
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

    def load_fasta_file(self, filename):
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

    def write_fasta_file(self, seq_record, extension='faa',
                         custom_id='', outdir=None, overwrite=False):
        """Write a fasta file for a Biopython SeqRecord object.

        Also check if the file exists already and returns the filename.
        You can overwrite existing fasta files by setting overwrite=True

        Args:
            seq_record (SeqRecord): Biopython SeqRecord object
            extension (str): Desired extension of the file (default: .faa)
            custom_id (str): OPTIONAL user-provided ID (default: SeqRecord ID)
            outdir (str): OPTIONAL desired directory of file output (default: current wd)
            overwrite (bool): OPTIONAL if you want to overwrite existing files

        Returns:
            outfile (str): Filename of fasta file
        """
        if custom_id:
            ident = custom_id
        else:
            ident = seq_record.id

        outfile = "{}.{}".format(ident, extension)

        if outdir:
            outfile = os.path.join(outdir, outfile)

        if os.path.isfile(outfile) and not overwrite:
            return outfile
        else:
            SeqIO.write(seq_record, outfile, "fasta")
            return outfile


if __name__ == '__main__':
    s = SeqTools()
    pass

    #TODO: make this a executable script to
    #1) display the sequence of a fasta file in plain string format
    #2) write a string sequence argument to a file
