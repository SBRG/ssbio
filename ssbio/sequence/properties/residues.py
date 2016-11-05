from Bio.SeqUtils.ProtParam import ProteinAnalysis


def sequence_properties(seq_str):
    analysed_seq = ProteinAnalysis(seq_str)
    return analysed_seq

#
# AAdict = {
#
#
# 'LYS': 'positive',
# 'ARG': 'positive',
# 'HIS': 'positive',
#
# 'ASP': 'negative',
# 'GLU': 'negative',
#
# 'LEU': 'nonpolar',
# 'TRP': 'nonpolar',
# 'VAL': 'nonpolar',
# 'PHE': 'nonpolar',
# 'PRO': 'nonpolar',
# 'ILE': 'nonpolar',
# 'GLY': 'nonpolar',
# 'ALA': 'nonpolar',
# 'MET': 'nonpolar',
#
# 'ASN': 'polar',
# 'THR': 'polar',
# 'TYR': 'polar',
# 'MSE': 'polar',
# 'SEC': 'polar',
# 'SER': 'polar',
# 'GLN': 'polar',
# 'CYS': 'polar',
# }
#
#
# def residue_props(seq_str):
#     """Return a dictionary of residue properties indicating the percentage of the respective property for a sequence.
#
#     Properties are: Polar, nonpolar, negative, positive.
#
#     Args:
#         pdb_file: PDB or MMCIF structure file
#
#     Returns:
#         dict: Dictonary of percentage (float) of properties
#     """
#
#     props = {}
#     polar = 0
#     nonpolar = 0
#     positive = 0
#     negative = 0
#     total = 0
#
#     for j in seq_str:
#         if j.resname in AAdict:
#             if AAdict[j.resname] == 'nonpolar':
#                 nonpolar = nonpolar + 1
#             elif AAdict[j.resname] == 'polar':
#                 polar = polar + 1
#             elif AAdict[j.resname] == 'positive':
#                 positive = positive + 1
#             elif AAdict[j.resname] == 'negative':
#                 negative = negative + 1
#             total = total + 1
#
#     props['ssb_per_NP'] = float(nonpolar) / float(total)
#     props['ssb_per_P'] = float(polar) / float(total)
#     props['ssb_per_pos'] = float(positive) / float(total)
#     props['ssb_per_neg'] = float(negative) / float(total)
#
#     return props