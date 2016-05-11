from ssbio.structure.pdbioext import PDBIOExt
from Bio import PDB

AAdict = {'CYS': 'polar',
          'ILE': 'nonpolar',
          'GLY': 'nonpolar',
          'SER': 'polar',
          'GLN': 'polar',
          'LYS': 'positive',
          'ASN': 'polar',
          'PRO': 'nonpolar',
          'ASP': 'negative',
          'THR': 'polar',
          'PHE': 'nonpolar',
          'ALA': 'nonpolar',
          'MET': 'nonpolar',
          'HIS': 'positive',
          'LEU': 'nonpolar',
          'ARG': 'positive',
          'TRP': 'nonpolar',
          'VAL': 'nonpolar',
          'GLU': 'negative',
          'TYR': 'polar',
          'MSE': 'polar',
          'SEC': 'polar'}


def residue_props(filename):
    """
    Returns a dictionary of residue properties indicating the
    percentage of the respective property.
    Properties are: Polar, nonpolar, negative, positive.
    Input: PDB or MMCIF structure file
    Output: Dictonary of percentage (float) of properties
    """

    my_structure = PDBIOExt(filename)
    model = my_structure.first_model

    props = {}
    polar = 0
    nonpolar = 0
    positive = 0
    negative = 0
    total = 0
    res_list = PDB.Selection.unfold_entities(model, 'R')

    for j in res_list:
        if j.resname in AAdict:
            if AAdict[j.resname] == 'nonpolar':
                nonpolar = nonpolar + 1
            elif AAdict[j.resname] == 'polar':
                polar = polar + 1
            elif AAdict[j.resname] == 'positive':
                positive = positive + 1
            elif AAdict[j.resname] == 'negative':
                negative = negative + 1
            total = total + 1

    props['ssb_per_NP'] = float(nonpolar) / float(total)
    props['ssb_per_P'] = float(polar) / float(total)
    props['ssb_per_pos'] = float(positive) / float(total)
    props['ssb_per_neg'] = float(negative) / float(total)

    return props

if __name__ == '__main__':
    import glob
    files = glob.glob('test_structures/*')
    for f in files:
        print(f)
        print(residue_props(f))
