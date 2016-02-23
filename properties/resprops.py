from Loader import Loader
from Bio import PDB
l = Loader()

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

    my_structure = l.structure_reader(filename)

    props = {}
    model = my_structure[0]
    polar = 0
    nonpolar = 0
    positive = 0
    negative = 0
    total = 0
    res_list = PDB.Selection.unfold_entities(my_structure, 'R')

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

    props['nonpolar'] = float(nonpolar) / float(total)
    props['polar'] = float(polar) / float(total)
    props['positive'] = float(positive) / float(total)
    props['negative'] = float(negative) / float(total)

    return props

if __name__ == '__main__':
    print(residue_props('test_structures/1u8f.pdb'))
