from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import os
import pandas as pd


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


def get_alignment_df(alignment_file):
    alignments = list(AlignIO.parse(alignment_file, "emboss"))
    alignment_df = pd.DataFrame(columns=[
                                'id_a', 'id_b', 'type', 'id_a_start', 'id_a_stop', 'count', 'id_a_aa', 'id_b_aa'])

    for alignment in alignments:
        #         if not switch:
        a_seq_id = list(alignment)[0].id
        a_seq = str(list(alignment)[0].seq)
        b_seq_id = list(alignment)[1].id
        b_seq = str(list(alignment)[1].seq)
    #         else:
    #             a_seq_id = list(alignment)[1].id
    #             a_seq = str(list(alignment)[1].seq)
    #             b_seq_id = list(alignment)[0].id
    #             b_seq = str(list(alignment)[0].seq)

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
# print '%s from %d to %d with %d members' % (appender['type'],
# appender['start'], appender['stop'], appender['count'])
                    alignment_df = alignment_df.append(
                        appender, ignore_index=True)

                if aa_flag_tracker == 'insertion':
                    appender = {}
                    appender['id_a'] = a_seq_id
                    appender['id_b'] = b_seq_id
                    appender['type'] = aa_flag_tracker
                    appender['id_a_start'] = new_i
                    appender['id_a_stop'] = new_i
                    appender['count'] = insertion_counter
# print '%s of length %d' % (aa_flag_tracker, insertion_counter)
                    alignment_df = alignment_df.append(
                        appender, ignore_index=True)
                    insertion_counter = 0
                idx_start = new_i + 1

            # HEY. THIS NEEDS TO BE OUTSIDE THE FOR LOOP
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
#                     print '%s at %d' % (aa_flag, new_i + 1)
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
# print '%s from %d to %d with %d members' % (appender['type'],
# appender['start'], appender['stop'], appender['count'])
                    alignment_df = alignment_df.append(
                        appender, ignore_index=True)

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
# print '%s from %d to %d with %d members' % (appender['type'],
# appender['start'], appender['stop'], appender['count'])
                    alignment_df = alignment_df.append(
                        appender, ignore_index=True)
# print new_i + 1, ':', a, b, 'last=%s, now=%s' %
# (aa_flag_tracker,aa_flag)

            aa_flag_tracker = aa_flag

    return alignment_df
