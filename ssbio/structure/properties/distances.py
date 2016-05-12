from ssbio.structure.pdbioext import PDBIOExt
from Bio import PDB
from Bio.PDB import Polypeptide
import ssbio.utils
import numpy as np
import cachetools

# calculate the magnitude of distance vector
def magni(a, b, c):
    return pow((pow(a, 2) + pow(b, 2) + pow(c, 2)), 1.0 / 2.0)

@cachetools.func.ttl_cache(maxsize=256)
def calculate_res_distance(res_1, res_2, pdb_file):
    """Calculate distance of one residue to another in a PDB file

    Args:
        res_1:
        res_2:
        pdb_file:

    Returns:

    """

    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model

    res_list = PDB.Selection.unfold_entities(model, 'R')

    ires_list = []
    res_chk_1 = ''
    res_chk_2 = ''
    for j in res_list:
        if j.id[1] in [res_1, res_2] and j.resname != 'HOH':
            ires_list.append(j)
            if res_chk_1 == '' and res_chk_2 == '':
                res_chk_1 = j.id[1]
            else:
                res_chk_2 = j.id[1]

    paired = ssbio.utils.combinations(ires_list, 2)
    # try:
    for k in paired:
        chainA = PDB.Selection.unfold_entities(k[0], 'C')[0]
        chainB = PDB.Selection.unfold_entities(k[1], 'C')[0]
        vec = list(
            np.array([x.get_coord() for x in k[0]]).mean(axis=0) - np.array([x.get_coord() for x in k[1]]).mean(
                axis=0))
        distance = magni(vec[0], vec[1], vec[2])

    return distance
    # except UnboundLocalError:
    #     return "Unknown interaction"