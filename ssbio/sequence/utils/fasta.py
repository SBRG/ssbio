import ssbio.utils
import tempfile
from Bio import SeqIO
import ssbio.sequence.utils


def write_fasta_file(seq_records, outname, outdir=None, outext='.faa', force_rerun=False):
    """Write a FASTA file for a SeqRecord or a list of SeqRecord objects.

    Args:
        seq_records (SeqRecord, list): SeqRecord or a list of SeqRecord objects
        outname: Name of the output file which will have outext appended to it
        outdir: Path to directory to output sequences to
        outext: Extension of FASTA file, default ".faa"
        force_rerun: If file should be overwritten if it exists

    Returns:
        str: Path to output FASTA file.

    """

    if not outdir:
        outdir = ''
    outfile = ssbio.utils.outfile_maker(inname='', outname=outname, outdir=outdir, outext=outext)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        SeqIO.write(seq_records, outfile, "fasta")

    return outfile


def write_fasta_file_from_dict(indict, outname, outdir=None, outext='.faa', force_rerun=False):
    """Write a FASTA file for a dictionary of IDs and their sequence strings.

    Args:
        indict: Input dictionary with keys as IDs and values as sequence strings
        outname: Name of the output file which will have outext appended to it
        outdir: Path to directory to output sequences to
        outext: Extension of FASTA file, default ".faa"
        force_rerun: If file should be overwritten if it exists

    Returns:
        str: Path to output FASTA file.

    """

    if not outdir:
        outdir = ''
    outfile = ssbio.utils.outfile_maker(inname='', outname=outname, outdir=outdir, outext=outext)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        seqs = []
        for i, s in indict.items():
            seq = ssbio.sequence.utils.cast_to_seq_record(s, id=i)
            seqs.append(seq)
        SeqIO.write(seqs, outfile, "fasta")

    return outfile


def write_seq_as_temp_fasta(seq):
    """Write a sequence as a temporary FASTA file

    Args:
        seq (str, Seq, SeqRecord): Sequence string, Biopython Seq or SeqRecord object

    Returns:
        str: Path to temporary FASTA file (located in system temporary files directory)

    """
    sr = ssbio.sequence.utils.cast_to_seq_record(seq, id='tempfasta')
    return write_fasta_file(seq_records=sr, outname='temp', outdir=tempfile.gettempdir(), force_rerun=True)


def load_fasta_file(filename):
    """Load a FASTA file and return the sequences as a list of SeqRecords

    Args:
        filename (str): Path to the FASTA file to load

    Returns:
        list: list of all sequences in the FASTA file as Biopython SeqRecord objects

    """

    with open(filename, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    return records


def load_fasta_file_as_dict_of_seqs(filename):
    """Load a FASTA file and return the sequences as a dict of {ID: sequence string}

    Args:
        filename (str): Path to the FASTA file to load

    Returns:
        dict: Dictionary of IDs to their sequence strings

    """

    results = {}
    records = load_fasta_file(filename)
    for r in records:
        results[r.id] = str(r.seq)

    return results


def fasta_files_equal(seq_file1, seq_file2):
    """Check equality of a FASTA file to another FASTA file

    Args:
        seq_file1: Path to a FASTA file
        seq_file2: Path to another FASTA file

    Returns:
        bool: If the sequences are the same

    """

    # Load already set representative sequence
    seq1 = SeqIO.read(open(seq_file1), 'fasta')

    # Load kegg sequence
    seq2 = SeqIO.read(open(seq_file2), 'fasta')

    # Test equality
    if str(seq1.seq) == str(seq2.seq):
        return True
    else:
        return False