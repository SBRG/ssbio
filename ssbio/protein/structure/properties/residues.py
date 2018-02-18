"""
Structure Residues
==================
"""

import logging
from collections import defaultdict
from copy import deepcopy

from Bio.Alphabet import IUPAC
from Bio.PDB import Polypeptide
from Bio.PDB.HSExposure import ExposureCN, HSExposureCA, HSExposureCB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import ssbio.protein.sequence.utils
from ssbio.protein.structure.utils.structureio import StructureIO

log = logging.getLogger(__name__)


def search_ss_bonds(model, threshold=3.0):
    """ Searches S-S bonds based on distances
        between atoms in the structure (first model only).
        Average distance is 2.05A. Threshold is 3A default.
        Returns iterator with tuples of residues.

        ADAPTED FROM JOAO RODRIGUES' BIOPYTHON GSOC PROJECT (http://biopython.org/wiki/GSOC2010_Joao)
    """

    # Taken from http://docs.python.org/library/itertools.html
    # Python 2.4 does not include itertools.combinations

    def combinations(iterable, r):
        # combinations('ABCD', 2) --> AB AC AD BC BD CD
        # combinations(range(4), 3) --> 012 013 023 123
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = list(range(r))
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i + 1, r):
                indices[j] = indices[j - 1] + 1
            yield tuple(pool[i] for i in indices)

    cysteines = [r for r in model.get_residues() if r.get_resname() == 'CYS']

    pairs = combinations(cysteines, 2)  # Iterator with pairs

    bridges = []
    for cys_pair in pairs:
        try:
            if cys_pair[0]['SG'] - cys_pair[1]['SG'] < threshold:
                bridges.append(cys_pair)
        except KeyError:  # This will occur when a CYS residue is missing a SG atom for some reason
            continue

    infodict = {}
    if bridges:
        infodict = defaultdict(list)

        for disulfide_bridge in bridges:
            residue1 = disulfide_bridge[0]
            residue2 = disulfide_bridge[1]
            chain = residue1.get_parent().id
            infodict[chain].append((residue1.get_full_id()[3], residue2.get_full_id()[3]))

    return infodict


def resname_in_proximity(resname, model, chains, resnums, threshold=5):
    """Search within the proximity of a defined list of residue numbers and their chains for any specifed residue name.

    Args:
        resname (str): Residue name to search for in proximity of specified chains + resnums
        model: Biopython Model object
        chains (str, list): Chain ID or IDs to check
        resnums (int, list): Residue numbers within the chain to check
        threshold (float): Cutoff in Angstroms for returning True if a RESNAME is near

    Returns:
        bool: True if a RESNAME is within the threshold cutoff

    """
    residues = [r for r in model.get_residues() if r.get_resname() == resname]

    chains = ssbio.utils.force_list(chains)
    resnums = ssbio.utils.force_list(resnums)

    for chain in chains:
        for resnum in resnums:
            my_residue_last_atom = model[chain][resnum].child_list[-1]
            for rz in residues:
                distance = rz.child_list[-1] - my_residue_last_atom
                if distance < threshold:
                    # print(resnum, rz, distance)
                    return True

    return False


def get_structure_seqrecords(model):
    """Get a dictionary of a PDB file's sequences.

    Special cases include:
        - Insertion codes. In the case of residue numbers like "15A", "15B", both residues are written out. Example: 9LPR
        - HETATMs. Currently written as an "X", or unknown amino acid.

    Args:
        model: Biopython Model object of a Structure

    Returns:
        list: List of SeqRecords

    """

    structure_seq_records = []

    # Loop over each chain of the PDB
    for chain in model:
        tracker = 0
        chain_seq = ''
        chain_resnums = []

        # Loop over the residues
        for res in chain.get_residues():
            # NOTE: you can get the residue number too
            res_id = res.id
            res_num = res_id[1]
            res_icode = res_id[2]

            # Double check if the residue name is a standard residue
            # If it is not a standard residue (ie. selenomethionine),
            # it will be filled in with an X on the next iteration)
            if Polypeptide.is_aa(res, standard=True):
                end_tracker = res_num
                res_aa_one = Polypeptide.three_to_one(res.get_resname())

                # Tracker to fill in X's
                if end_tracker != (tracker + 1):
                    if res_icode != ' ':
                        chain_seq += res_aa_one
                        chain_resnums.append(res_num)
                        tracker = end_tracker + 1
                        continue
                    else:
                        multiplier = (end_tracker - tracker - 1)
                        chain_seq += 'X' * multiplier
                        # Residue numbers for unresolved or nonstandard residues are Infinite
                        chain_resnums.extend([float("Inf")] * multiplier)

                chain_seq += res_aa_one
                chain_resnums.append(res_num)
                tracker = end_tracker

            else:
                continue

        chain_seq_record = SeqRecord(Seq(chain_seq, IUPAC.protein), id=chain.get_id())
        chain_seq_record.letter_annotations['structure_resnums'] = chain_resnums
        structure_seq_records.append(chain_seq_record)

    return structure_seq_records


def get_structure_seqs(pdb_file, file_type):
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

    # TODO: Please check out capitalization of chain IDs in mmcif files. example: 5afi - chain "l" is present but
    # it seems like biopython capitalizes it to chain L

    # Get the first model
    my_structure = StructureIO(pdb_file)
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


def match_structure_sequence(orig_seq, new_seq, match='X', fill_with='X', ignore_excess=False):
    """Correct a sequence to match inserted X's in a structure sequence

    This is useful for mapping a sequence obtained from structural tools like MSMS or DSSP
        to the sequence obtained by the get_structure_seqs method.

    Examples:
        >>> structure_seq = 'XXXABCDEF'
        >>> prop_list = [4, 5, 6, 7, 8, 9]
        >>> match_structure_sequence(structure_seq, prop_list)
        ['X', 'X', 'X', 4, 5, 6, 7, 8, 9]

        >>> match_structure_sequence(structure_seq, prop_list, fill_with=float('Inf'))
        [inf, inf, inf, 4, 5, 6, 7, 8, 9]

        >>> structure_seq = '---ABCDEF---'
        >>> prop_list = ('H','H','H','C','C','C')
        >>> match_structure_sequence(structure_seq, prop_list, match='-', fill_with='-')
        ('-', '-', '-', 'H', 'H', 'H', 'C', 'C', 'C', '-', '-', '-')

        >>> structure_seq = 'ABCDEF---'
        >>> prop_list = 'HHHCCC'
        >>> match_structure_sequence(structure_seq, prop_list, match='-', fill_with='-')
        'HHHCCC---'

        >>> structure_seq = 'AXBXCXDXEXF'
        >>> prop_list = ['H', 'H', 'H', 'C', 'C', 'C']
        >>> match_structure_sequence(structure_seq, prop_list, match='X', fill_with='X')
        ['H', 'X', 'H', 'X', 'H', 'X', 'C', 'X', 'C', 'X', 'C']

    Args:
        orig_seq (str, Seq, SeqRecord): Sequence to match to
        new_seq (str, tuple, list): Sequence to fill in
        match (str): What to match
        fill_with: What to fill in when matches are found
        ignore_excess (bool): If excess sequence on the tail end of new_seq should be ignored

    Returns:
        str, tuple, list: new_seq which will match the length of orig_seq

    """
    if len(orig_seq) == len(new_seq):
        log.debug('Lengths already equal, nothing to fill in')
        return new_seq

    if not ignore_excess:
        if len(orig_seq) < len(new_seq):
            raise ValueError('Original sequence has a length less than the sequence provided to match to')
    else:
        log.debug('New sequence will be truncated to length of original sequence - information may be lost!')

    if not isinstance(new_seq, str) and not isinstance(new_seq, tuple) and not isinstance(new_seq, list):
        raise ValueError('Invalid sequence provided, must be string, tuple, or list')

    orig_seq = ssbio.protein.sequence.utils.cast_to_str(orig_seq)
    new_thing = deepcopy(new_seq)
    if isinstance(new_seq, tuple):
        new_thing = list(new_thing)

    for i, s in enumerate(orig_seq):
        if s == match:
            if isinstance(new_thing, str):
                new_thing = new_thing[:i] + fill_with + new_thing[i:]
            if isinstance(new_thing, list):
                new_thing.insert(i, fill_with)

    new_thing = new_thing[:len(orig_seq)]

    if isinstance(new_seq, tuple):
        new_thing = tuple(new_thing)

    return new_thing


def site_centroid(residues, model):
    """Get the XYZ coordinate of the center of a list of residues.

    Args:
        residues: List of residue numbers
        pdb_file: Path to PDB file

    Returns:
        tuple: (X, Y, Z) coordinate of centroid

    """
    pass


def distance_to_site(residue_of_interest, residues, model):
    """Calculate the distance between an amino acid and a group of amino acids.

    Args:
        residue_of_interest: Residue number you are interested in (ie. a mutation)
        residues: List of residue numbers

    Returns:
        float: Distance (in Angstroms) to the group of residues

    """
    centroid = site_centroid(residues, residue_of_interest)
    pass


# TODO: half sphere exposure
def hse_output(pdb_file, file_type):
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
    my_structure = StructureIO(pdb_file)
    model = my_structure.first_model

    # Calculate HSEalpha
    exp_ca = HSExposureCA(model)
    # Calculate HSEbeta
    exp_cb = HSExposureCB(model)
    # Calculate classical coordination number
    exp_fs = ExposureCN(model)

    return
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
#     my_structure = StructureIO(pdb_file)
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