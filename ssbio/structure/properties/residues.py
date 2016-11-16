from ssbio.structure.pdbioext import PDBIOExt
from Bio import PDB
from Bio.PDB import Polypeptide
import ssbio.utils
import numpy as np
# import cachetools

from Bio.PDB import Polypeptide
from ssbio.structure.pdbioext import PDBIOExt
from Bio import PDB
# from Bio import Struct
from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB, ExposureCN

# TODO: network representation of structure ("neighborhood" stuff)?


# TODO: half sphere exposure
def hse_output(pdb_file):
    """
    The solvent exposure of an amino acid residue is important for analyzing,
    understanding and predicting aspects of protein structure and function [73].
    A residue's solvent exposure can be classified as four categories: exposed, partly exposed,
    buried and deeply buried residues. Hamelryck  et al. [73] established a new 2D measure that provides a
    different view of solvent exposure, i.e. half-sphere exposure (HSE). By conceptually dividing the sphere
    of a residue into two halves- HSE-up and HSE-down, HSE provides a more detailed description of an amino
    acid residue's spatial neighborhood. HSE is calculated by the hsexpo module implemented in the BioPython
    package [74] from a PDB file.

    http://onlinelibrary.wiley.com/doi/10.1002/prot.20379/abstract

    Args:
        pdb_file:

    Returns:

    """
    # Get the first model
    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model

    # Calculate HSEalpha
    exp_ca = HSExposureCA(model)
    # Calculate HSEbeta
    exp_cb = HSExposureCB(model)
    # Calculate classical coordination number
    exp_fs = ExposureCN(model)

    return


def get_pdb_seqs(pdb_file):
    """Get a dictionary of a PDB file's sequences.

    Special cases include:
        - Insertion codes. In the case of residue numbers like "15A", "15B", both residues are written out. Example: 9LPR
        - HETATMs. Currently written as an "X", or unknown amino acid.

    Args:
        pdb_file: Path to PDB file

    Returns:
        dict: Dictionary of:
        {chain_id: sequence}

    """
    # Get the first model
    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model

    structure_seqs = {}

    # Loop over each chain of the PDB
    for chain in model:
        chain_seq = ''
        tracker = 0

        # Loop over the residues
        for res in chain.get_residues():
            # NOTE: you can get the residue number too
            # res_num = res.id[1]

            # Double check if the residue name is a standard residue
            # If it is not a standard residue (ie. selenomethionine),
            # it will be filled in with an X on the next iteration)
            if Polypeptide.is_aa(res, standard=True):
                full_id = res.get_full_id()
                end_tracker = full_id[3][1]
                i_code = full_id[3][2]
                aa = Polypeptide.three_to_one(res.get_resname())

                # Tracker to fill in X's
                if end_tracker != (tracker + 1):
                    if i_code != ' ':
                        chain_seq += aa
                        tracker = end_tracker + 1
                        continue
                    else:
                        chain_seq += 'X' * (end_tracker - tracker - 1)

                chain_seq += aa
                tracker = end_tracker

            else:
                continue

        structure_seqs[chain.get_id()] = chain_seq

    return structure_seqs


# @cachetools.func.ttl_cache(maxsize=1000)
def get_pdb_res_starts(pdb_file):
    """Return a dictionary of the first residue number in each chain of a PDB file

    Args:
        pdb_file: path to PDB file

    Returns:
        start_residues: dictionary of {chainID: firstResNum, ...}

    """
    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model

    start_residues = {}
    for chain in model:
        residues = chain.get_residues()
        start_residues[chain.get_id()] = next(residues).get_id()[1]

    return start_residues

#
# def magni(a, b, c):
#     """Calculate the magnitude of distance vector
#     """
#     return pow((pow(a, 2) + pow(b, 2) + pow(c, 2)), 1.0 / 2.0)


# @cachetools.func.ttl_cache(maxsize=256)
# def calculate_res_distance(res_1, res_2, pdb_file):
#     """Calculate distance of one residue number to another in a PDB file
#
#     Args:
#         res_1: Residue number 1
#         res_2: Residue number 2
#         pdb_file: Path to PDB file
#
#     Returns:
#
#     """
#
#     my_structure = PDBIOExt(pdb_file)
#     model = my_structure.first_model
#
#     res_list = PDB.Selection.unfold_entities(model, 'R')
#
#     ires_list = []
#     res_chk_1 = ''
#     res_chk_2 = ''
#     for j in res_list:
#         if j.id[1] in [res_1, res_2] and j.resname != 'HOH':
#             ires_list.append(j)
#             if res_chk_1 == '' and res_chk_2 == '':
#                 res_chk_1 = j.id[1]
#             else:
#                 res_chk_2 = j.id[1]
#
#     paired = ssbio.utils.combinations(ires_list, 2)
#     try:
#         for k in paired:
#             chainA = PDB.Selection.unfold_entities(k[0], 'C')[0]
#             chainB = PDB.Selection.unfold_entities(k[1], 'C')[0]
#             vec = list(
#                 np.array([x.get_coord() for x in k[0]]).mean(axis=0) - np.array([x.get_coord() for x in k[1]]).mean(
#                     axis=0))
#             distance = magni(vec[0], vec[1], vec[2])
#
#         return distance
#     except UnboundLocalError:
#         log.error("Unknown interaction")
#         return None

# def count_ss_bond(filename, threshold_new=5):
#     """
#     Counts the number of sulfide bridges (S-S bonds formed by
#     cysteines in close proximity)
#     Input: PDB or mmCIF file (will be parsed into a Struct object)
#     Output: Int of number of SS bonds
#     """
#     s = Struct.read(filename)
#     my_structure = s.as_protein()
#     ss = my_structure.search_ss_bonds(threshold=threshold_new)
#
#     return len(list(ss))
#     return len(list(ss))