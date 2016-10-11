from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import os
import pandas as pd
import numpy as np
from collections import defaultdict

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

# SEE: https://www.biostars.org/p/91124/

def run_alignment2(a_id, a_seq, b_id, b_seq):
    '''
    Runs the needle alignment program and returns a raw text dump of the alignment

    Input:  a_id - sequence ID #1 (string)
            a_seq - sequence #1 (string)
            b_id - sequence ID #2 (string)
            b_seq - sequence #2 (string)
    Output: alignment_file - file name of alignment

    DEPENDENCIES:
    get_alignment_allpos_df
    '''

    from Bio.Emboss.Applications import NeedleCommandline

    alignment_file = "/tmp/%s_%s_align.txt" % (a_id, b_id)

    needle_cline = NeedleCommandline(asequence="asis::"+a_seq, bsequence="asis::"+b_seq, gapopen=10, gapextend=0.5, outfile=alignment_file)
    stdout, stderr = needle_cline()

    return get_alignment_allpos_df(alignment_file, a_id, b_id)


def get_alignment_df(alignment_file):
    '''
    Creates an alignment dataframe indicating where matches, mutations, insertions, and deletions are.
    Takes in a needle alignment file and returns a Pandas DataFrame.
    Input: alignment_file - needle alignment file between two amino acid sequences
    Output: alignment_df - Pandas DataFrame of the alignment.
    '''

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
# print new_i + 1, ':', a, b, 'last=%kegg, now=%kegg' %
# (aa_flag_tracker,aa_flag)

            aa_flag_tracker = aa_flag

    return alignment_df

def get_alignment_allpos_df(alignment_file, a_seq_id=None, b_seq_id=None):
    alignments = list(AlignIO.parse(alignment_file, "emboss"))

    appender = defaultdict(dict)
    idx = 0
    for alignment in alignments:
    #         if not switch:
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

#     alignment_df = alignment_df.join(alignment_df['id_a'].apply(lambda x: pd.Series(x.split('.')[1])))
#     alignment_df = alignment_df.rename(columns={0:'id_a_chain'})

#     alignment_df = alignment_df.join(alignment_df['id_b'].apply(lambda x: pd.Series(x.split('.')[1])))
#     alignment_df = alignment_df.rename(columns={0:'id_b_chain'})

    return alignment_df


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


def needle_reader(fl):
    '''
    Reads in a needle alignment file and spits out statistics of the alignment.

    Input: fl - alignment file name
    Output: alignment_properties - a dictionary telling you the number of gaps, identity, etc.
    '''

    alignments = list(AlignIO.parse(fl, "emboss"))

    f=open(fl)
    alignment_properties = defaultdict(dict)
    line = f.readline()

    for i in range(len(alignments)):
        #print alignments[i]
        while line.rstrip() != "#=======================================":
            line = f.readline()
            if not line:
                raise StopIteration

        while line[0] == "#":
            #Read in the rest of this alignment header,
            #try and discover the number of records expected
            #and their length
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

            #And read in another line...
            line = f.readline()

    return alignment_properties
