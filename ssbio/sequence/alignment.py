import os.path as op
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import AlignIO
from Bio.Emboss.Applications import NeedleCommandline


def run_needle_alignment_on_files(id_a, faa_a, id_b, faa_b, gapopen=10, gapextend=0.5,
                                  write_output=False, outdir='', outfile='', force_rerun=False):
    """Run the needle alignment program for two fasta files and return the raw alignment result.

    More info:
    EMBOSS needle: http://www.bioinformatics.nl/cgi-bin/emboss/help/needle
    Biopython wrapper: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84

    Args:
        id_a: ID of reference sequence
        faa_a: File path to reference sequence
        id_b: ID of sequence to be aligned
        faa_b: File path to sequence to be aligned
        gapopen: Gap open penalty is the score taken away when a gap is created
        gapextend: Gap extension penalty is added to the standard gap penalty for each base or residue in the gap
        write_output (bool): Default False, set to True if you want the alignment file saved
        outdir (str, optional): Path to output directory. Default is the current directory.
        outfile (str, optional): Name of output file. If not set, is {id_a}_{id_b}_align.txt
        force_rerun (bool): Default False, set to True if you want to rerun the alignment if outfile exists.

    Returns:
        str: Raw alignment result of the needle alignment in srspair format.

    """

    # If you don't want to save the output file, just run the alignment and return the raw results
    if not write_output:
        needle_cline = NeedleCommandline(asequence=faa_a, bsequence=faa_b,
                                         gapopen=gapopen, gapextend=gapextend,
                                         stdout=True, auto=True)
        raw_alignment_text, stderr = needle_cline()

    # If you do want to save an output file...
    else:
        # Make a default name if no outfile is set
        if not outfile:
            outfile = op.join(outdir, '{}_{}_align.txt'.format(id_a, id_b))
        else:
            outfile = op.join(outdir, outfile)

        # Check if the outfile already exists, and read it and return those results
        if op.exists(outfile) and not force_rerun:
            with open(outfile) as f:
                raw_alignment_text = f.read()

        # If it doesn't exist, or force_rerun=True, run the alignment
        else:
            needle_cline = NeedleCommandline(asequence=faa_a, bsequence=faa_b,
                                             gapopen=gapopen, gapextend=gapextend,
                                             outfile=outfile)
            stdout, stderr = needle_cline()
            with open(outfile) as f:
                raw_alignment_text = f.read()

    return raw_alignment_text


def run_needle_alignment_on_str(id_a, seq_a, id_b, seq_b, gapopen=10, gapextend=0.5,
                                write_output=False, outdir='', outfile='', force_rerun=False):
    """Run the needle alignment program for two strings and return the raw alignment result.

    More info:
    EMBOSS needle: http://www.bioinformatics.nl/cgi-bin/emboss/help/needle
    Biopython wrapper: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84
    Using strings as input: https://www.biostars.org/p/91124/

    Args:
        id_a: ID of reference sequence
        seq_a: String representation of reference sequence
        id_b: ID of sequence to be aligned
        seq_b: String representation of sequence to be aligned
        gapopen: Gap open penalty is the score taken away when a gap is created
        gapextend: Gap extension penalty is added to the standard gap penalty for each base or residue in the gap
        write_output (bool): Default False, set to True if you want the alignment file saved
        outdir (str, optional): Path to output directory. Default is the current directory.
        outfile (str, optional): Name of output file. If not set, is {id_a}_{id_b}_align.txt
        force_rerun (bool): Default False, set to True if you want to rerun the alignment if outfile exists.

    Returns:
        str: Raw alignment result of the needle alignment in srspair format.

    """
    # If you don't want to save the output file, just run the alignment and return the raw results
    if not write_output:
        needle_cline = NeedleCommandline(asequence="asis::"+seq_a, bsequence="asis::"+seq_b,
                                         gapopen=gapopen, gapextend=gapextend,
                                         stdout=True, auto=True)
        raw_alignment_text, stderr = needle_cline()

    # If you do want to save an output file...
    else:
        # Make a default name if no outfile is set
        if not outfile:
            outfile = op.join(outdir, '{}_{}_align.txt'.format(id_a, id_b))
        else:
            outfile = op.join(outdir, outfile)

        # Check if the outfile already exists, and read it and return those results
        if op.exists(outfile) and not force_rerun:
            with open(outfile) as f:
                raw_alignment_text = f.read()

        # If it doesn't exist, or force_rerun=True, run the alignment
        else:
            needle_cline = NeedleCommandline(asequence="asis::"+seq_a, bsequence="asis::"+seq_b,
                                             gapopen=gapopen, gapextend=gapextend,
                                             outfile=outfile)
            stdout, stderr = needle_cline()
            with open(outfile) as f:
                raw_alignment_text = f.read()

    return raw_alignment_text


def get_alignment_summary_df(alignment_file):
    """Create an alignment dataframe indicating where matches, mutations, insertions, and deletions are.
        Takes in a needle alignment file and returns a Pandas DataFrame.

    Args:
        alignment_file - needle alignment file between two amino acid sequences

    Retunrs:
        alignment_df - Pandas DataFrame of the alignment.

    """

    alignments = list(AlignIO.parse(alignment_file, "emboss"))
    alignment_df = pd.DataFrame(columns=['id_a', 'id_b', 'type', 'id_a_start', 'id_a_stop', 'count', 'id_a_aa', 'id_b_aa'])

    for alignment in alignments:
        a_seq_id = list(alignment)[0].id
        a_seq = str(list(alignment)[0].seq)
        b_seq_id = list(alignment)[1].id
        b_seq = str(list(alignment)[1].seq)

        idx_start = 1

        uniprot_shift = 0
        insertion_counter = 0
        for i, (a, b) in enumerate(zip(a_seq, b_seq)):
            if uniprot_shift != 0:
                new_i = i - uniprot_shift
            else:
                new_i = i

            if a == b and a != '-' and b != '-':
                aa_flag = 'match'
            if a != b and a == '-' and b != '-':
                aa_flag = 'insertion'
            if a != b and a != '-' and b == '-':
                aa_flag = 'deletion'
            if a != b and a != '-' and b == 'X':
                aa_flag = 'unresolved'
            elif a != b and a != '-' and b != '-':
                aa_flag = 'mutation'
            if i == 0:
                aa_flag_tracker = aa_flag

            if aa_flag == 'insertion':
                uniprot_shift += 1
                insertion_counter += 1

            if aa_flag != aa_flag_tracker:
                idx_end = new_i
                if aa_flag_tracker == 'match' or aa_flag_tracker == 'deletion' or aa_flag_tracker == 'unresolved':
                    appender = {}
                    appender['id_a'] = a_seq_id
                    appender['id_b'] = b_seq_id
                    appender['type'] = aa_flag_tracker
                    appender['id_a_start'] = idx_start
                    appender['id_a_stop'] = idx_end
                    appender['count'] = idx_end - idx_start + 1
                    alignment_df = alignment_df.append(appender, ignore_index=True)

                if aa_flag_tracker == 'insertion':
                    appender = {}
                    appender['id_a'] = a_seq_id
                    appender['id_b'] = b_seq_id
                    appender['type'] = aa_flag_tracker
                    appender['id_a_start'] = new_i
                    appender['id_a_stop'] = new_i
                    appender['count'] = insertion_counter
                    alignment_df = alignment_df.append(appender, ignore_index=True)
                    insertion_counter = 0
                idx_start = new_i + 1

            # TODO: test if this should be outside the for loop
            if aa_flag == 'mutation':
                appender = {}
                appender['id_a'] = a_seq_id
                appender['id_b'] = b_seq_id
                appender['type'] = aa_flag
                appender['id_a_start'] = new_i + 1
                appender['id_a_stop'] = new_i + 1
                appender['count'] = 1
                appender['id_a_aa'] = a
                appender['id_b_aa'] = b
                alignment_df = alignment_df.append(appender, ignore_index=True)
                idx_start = new_i + 1

            elif (i + 1) == len(a_seq) and aa_flag_tracker == aa_flag:
                idx_end = new_i + 1
                if aa_flag_tracker != 'mutation' and aa_flag_tracker != 'insertion':
                    appender = {}
                    appender['id_a'] = a_seq_id
                    appender['id_b'] = b_seq_id
                    appender['type'] = aa_flag_tracker
                    appender['id_a_start'] = idx_start
                    appender['id_a_stop'] = idx_end
                    appender['count'] = idx_end - idx_start + 1
                    alignment_df = alignment_df.append(appender, ignore_index=True)

            elif (i + 1) == len(a_seq) and aa_flag_tracker != aa_flag:
                idx_end = new_i + 1
                idx_start = new_i + 1
                if aa_flag_tracker != 'mutation' and aa_flag_tracker != 'insertion':
                    appender = {}
                    appender['id_a'] = a_seq_id
                    appender['id_b'] = b_seq_id
                    appender['type'] = aa_flag
                    appender['id_a_start'] = idx_start
                    appender['id_a_stop'] = idx_end
                    appender['count'] = idx_end - idx_start + 1
                    alignment_df = alignment_df.append(appender, ignore_index=True)

            aa_flag_tracker = aa_flag

    return alignment_df


def get_alignment_df(alignment_file, a_seq_id=None, b_seq_id=None):
    """Get a Pandas DataFrame of the Needle alignment results. Contains all positions of the sequences.

    Args:
        alignment_file:
        a_seq_id: Optional specification of the ID of the reference sequence
        b_seq_id: Optional specification of the ID of the aligned sequence

    Returns:
        Pandas DataFrame: all positions in the alignment

    """
    alignments = list(AlignIO.parse(alignment_file, "emboss"))

    appender = defaultdict(dict)
    idx = 0
    for alignment in alignments:
        if not a_seq_id:
            a_seq_id = list(alignment)[0].id
        a_seq = str(list(alignment)[0].seq)
        if not b_seq_id:
            b_seq_id = list(alignment)[1].id
        b_seq = str(list(alignment)[1].seq)

        a_idx = 1
        b_idx = 1

        for i, (a,b) in enumerate(zip(a_seq,b_seq)):
            if a == b and a != '-' and b != '-':
                aa_flag = 'match'
            if a != b and a == '-' and b != '-':
                aa_flag = 'insertion'
            if a != b and a != '-' and b == '-':
                aa_flag = 'deletion'
            if a != b and a != '-' and b == 'X':
                aa_flag = 'unresolved'
            if a != b and b != '-' and a == 'X':
                aa_flag = 'unresolved'
            elif a != b and a != '-' and b != '-':
                aa_flag = 'mutation'

            appender[idx]['id_a'] = a_seq_id
            appender[idx]['id_b'] = b_seq_id
            appender[idx]['type'] = aa_flag

            if aa_flag == 'match' or aa_flag == 'unresolved' or aa_flag == 'mutation':
                appender[idx]['id_a_aa'] = a
                appender[idx]['id_a_pos'] = a_idx
                appender[idx]['id_b_aa'] = b
                appender[idx]['id_b_pos'] = b_idx
                a_idx += 1
                b_idx += 1

            if aa_flag == 'deletion':
                appender[idx]['id_a_aa'] = a
                appender[idx]['id_a_pos'] = a_idx
                a_idx += 1

            if aa_flag == 'insertion':
                appender[idx]['id_b_aa'] = b
                appender[idx]['id_b_pos'] = b_idx
                b_idx += 1

            idx += 1

    alignment_df = pd.DataFrame.from_dict(appender, orient='index')
    alignment_df = alignment_df[['id_a', 'id_b', 'type', 'id_a_aa', 'id_a_pos', 'id_b_aa', 'id_b_pos']].fillna(value=np.nan)

    return alignment_df


def needle_statistics(infile):
    """Reads in a needle alignment file and spits out statistics of the alignment.

    Args:
        infile (str): Alignment file name

    Returns:
        dict: alignment_properties - a dictionary telling you the number of gaps, identity, etc.

    """

    alignments = list(AlignIO.parse(infile, "emboss"))

    f = open(infile)
    alignment_properties = defaultdict(dict)
    line = f.readline()

    for i in range(len(alignments)):
        while line.rstrip() != "#=======================================":
            line = f.readline()
            if not line:
                raise StopIteration

        while line[0] == "#":
            # Read in the rest of this alignment header,
            # try and discover the number of records expected and their length
            parts = line[1:].split(":", 1)
            key = parts[0].lower().strip()
            if key == '2':
                pdb_id = parts[1].strip()
            if key == 'identity':
                ident_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                ident_num = int(ident_parse[0].split('/')[0])
                ident_percent = float(ident_parse[1])
                alignment_properties[pdb_id]['identity'] = ident_num
                alignment_properties[pdb_id]['identity_percent'] = ident_percent
            if key == 'similarity':
                sim_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                sim_num = int(sim_parse[0].split('/')[0])
                sim_percent = float(sim_parse[1])
                alignment_properties[pdb_id]['similarity'] = sim_num
                alignment_properties[pdb_id]['similarity_percent'] = sim_percent
            if key == 'gaps':
                gap_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                gap_num = int(gap_parse[0].split('/')[0])
                gap_percent = float(gap_parse[1])
                alignment_properties[pdb_id]['gaps'] = gap_num
                alignment_properties[pdb_id]['gaps_percent'] = gap_percent
            if key == 'score':
                score = float(parts[1].strip())
                alignment_properties[pdb_id]['score'] = score

            # And read in another line...
            line = f.readline()

    return alignment_properties


# TODO: needleall has a bug
# def run_alignment_needleall(a_id, a_faa, b_id, b_faa):
#
#     from Bio.Emboss.Applications import NeedleallCommandline
#     import os.path
#
#     alignment_file = "%s-%s_align.txt" % (a_id, b_id)
#
#     if os.path.isfile(alignment_file):
#         # print 'Alignment %s file already exists' % alignment_file
#         return alignment_file
#
#     else:
#         print '**RUNNING ALIGNMENT FOR %kegg AND %kegg**' % (a_id, b_id)
#         needle_cline = NeedleallCommandline(asequence=a_faa, bsequence=b_faa, gapopen=10, gapextend=0.5, outfile=alignment_file)
#         stdout, stderr = needle_cline()
#         return alignment_file
