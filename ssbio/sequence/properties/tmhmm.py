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

    infodict = defaultdict(str)
    for l in lines:
        # Look for the lines with tab separations, these are the TM predicted regions
        if '\t' not in l:
            continue

        stuff = l.split('\t')
        if stuff[1] == 'TMHMM2.0':
            gene = stuff[0]
            region = stuff[2:][0]
            region_resnums = stuff[3].split()
            region_start = region_resnums[0]
            region_end = region_resnums[1]

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
                infodict[gene] += info

    return infodict