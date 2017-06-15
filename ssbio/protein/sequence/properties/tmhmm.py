from collections import defaultdict
import logging
log = logging.getLogger(__name__)

# Parsing TMHMM results
# There are two output formats: Long and short. Long output format
#
# For the long format (default), tmhmm gives some statistics and a list of the location of the predicted transmembrane helices and the predicted location of the intervening loop regions. Here is an example:
#
#     # COX2_BACSU Length: 278
#     # COX2_BACSU Number of predicted TMHs: 3
#     # COX2_BACSU Exp number of AAs in TMHs: 68.6888999999999
#     # COX2_BACSU Exp number, first 60 AAs: 39.8875
#     # COX2_BACSU Total prob of N-in: 0.99950
#     # COX2_BACSU POSSIBLE N-term signal sequence
#     COX2_BACSU TMHMM2.0 inside 1 6
#     COX2_BACSU TMHMM2.0 TMhelix 7 29
#     COX2_BACSU TMHMM2.0 outside 30 43
#     COX2_BACSU TMHMM2.0 TMhelix 44 66
#     COX2_BACSU TMHMM2.0 inside 67 86
#     COX2_BACSU TMHMM2.0 TMhelix 87 109
#     COX2_BACSU TMHMM2.0 outside 110 278
#
# If the whole sequence is labeled as inside or outside, the prediction is that it contains no membrane
# helices. It is probably not wise to interpret it as a prediction of location. The prediction gives the most probable location and orientation of transmembrane helices in the sequence. It is found by an algorithm called N-best (or 1-best in this case) that sums over all paths through the model with the same location and direction of the helices.
#
# The first few lines gives some statistics:
#
# Length: the length of the protein sequence.
#
# Number of predicted TMHs: The number of predicted transmembrane helices.
#
# Exp number of AAs in TMHs: The expected number of amino acids intransmembrane helices. If this number is larger than 18 it is very likely to be a transmembrane protein (OR have a signal peptide).
#
# Exp number, first 60 AAs: The expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein. If this number more than a few, you should be warned that a predicted transmembrane helix in the N-term could be a signal peptide.
#
# Total prob of N-in: The total probability that the N-term is on the cytoplasmic side of the membrane.
#
# POSSIBLE N-term signal sequence: a warning that is produced when "Exp number, first 60 AAs" is larger than 10.


def parse_tmhmm_long(tmhmm_results):
    with open(tmhmm_results) as f:
        lines = f.read().splitlines()

    infodict = defaultdict(dict)

    for l in lines:
        if 'Number of predicted TMHs:' in l:
            gene = l.split(' Number')[0].strip('# ')
            infodict[gene]['num_tm_helices'] = int(l.split(': ')[1])

        if 'WARNING' in l:
            log.warning('{}: no TMHMM predictions'.format(l))
            continue

        # TODO: POSSIBLE N-term signal sequence - parse this

        # Look for the lines without #, these are the TM predicted regions
        if '#' not in l:
            stuff = l.split()
            if stuff[1] == 'TMHMM2.0':

                gene = stuff[0]
                region = stuff[2]
                region_start = stuff[3]
                region_end = stuff[4]

                if 'sequence' in infodict[gene]:
                    tm_seq = infodict[gene]['sequence']
                else:
                    tm_seq = ''

                if region == 'outside':
                    info = 'O'
                elif region == 'inside':
                    info = 'I'
                elif region == 'TMhelix':
                    info = 'T'
                else:
                    log.error('{}: unknown region type'.format(info))
                    info = '-'

                for r in range(int(region_start), int(region_end) + 1):
                    tm_seq += info

                infodict[gene]['sequence'] = tm_seq

    return infodict


def label_TM_tmhmm_residue_numbers_and_leaflets(tmhmm_seq):
    """Determine the residue numbers of the TM-helix residues that cross the membrane and label them by leaflet.

    Args:
        tmhmm_seq: g.protein.representative_sequence.seq_record.letter_annotations['TM-tmhmm']

    Returns:
        leaflet_dict: a dictionary with leaflet_variable : [residue list] where the variable is inside or outside
        TM_boundary dict: outputs a dictionar with : TM helix number : [TM helix residue start , TM helix residue end]

    TODO:
        untested method!
    """

    TM_number_dict = {}
    T_index = []
    T_residue = []

    residue_count = 1
    for residue_label in tmhmm_seq:
        if residue_label == 'T':
            T_residue.append(residue_count)

        residue_count = residue_count + 1
    TM_number_dict.update({'T_residue': T_residue})

    # finding the TM boundaries
    T_residue_list = TM_number_dict['T_residue']

    count = 0
    max_count = len(T_residue_list) - 1
    TM_helix_count = 0
    TM_boundary_dict = {}

    while count <= max_count:
        # first residue = TM start
        if count == 0:
            TM_start = T_residue_list[count]
            count = count + 1
            continue
        # Last residue = TM end
        elif count == max_count:
            TM_end = T_residue_list[count]
            TM_helix_count = TM_helix_count + 1
            TM_boundary_dict.update({'TM_helix_' + str(TM_helix_count): [TM_start, TM_end]})
            break
        # middle residues need to be start or end
        elif T_residue_list[count] != T_residue_list[count + 1] - 1:
            TM_end = T_residue_list[count]
            TM_helix_count = TM_helix_count + 1
            TM_boundary_dict.update({'TM_helix_' + str(TM_helix_count): [TM_start, TM_end]})
            # new TM_start
            TM_start = T_residue_list[count + 1]
        count = count + 1
    # assign leaflet to proper TM residues O or I
    leaflet_dict = {}
    for leaflet in ['O', 'I']:
        leaflet_list = []
        for TM_helix, TM_residues in TM_boundary_dict.items():
            for residue_num in TM_residues:
                tmhmm_seq_index = residue_num - 1
                previous_residue = tmhmm_seq_index - 1
                next_residue = tmhmm_seq_index + 1
                # identify if the previous or next residue closest to the TM helix start/end is the proper leaflet
                if tmhmm_seq[previous_residue] == leaflet or tmhmm_seq[next_residue] == leaflet:
                    leaflet_list.append(residue_num)
        leaflet_dict.update({'tmhmm_leaflet_' + leaflet: leaflet_list})

    return TM_boundary_dict, leaflet_dict