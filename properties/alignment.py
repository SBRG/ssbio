from Bio.Emboss.Applications import NeedleCommandline
import os


def run_alignment(id_a, faa_a, id_b, faa_b):
    '''
    Runs the needle alignment program for two fasta files
    and writes the result to a file. Returns the filename.

    Input:  id_a - first id
            faa_a - first fasta file name
            id_b - second id
            faa_b - second fasta file name
    Output: alignment_file - file name of alignment
    '''

    alignment_file = "%s_%s_align.txt" % (id_a, id_b)

    if os.path.isfile(alignment_file):
        return alignment_file

    else:
        needle_cline = NeedleCommandline(
            asequence=faa_a, bsequence=faa_b, gapopen=10, gapextend=0.5, outfile=alignment_file)
        stdout, stderr = needle_cline()
        return alignment_file
