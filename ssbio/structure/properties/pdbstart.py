import os
import cachetools
from ssbio.structure.pdbioext import PDBIOExt

@cachetools.func.ttl_cache(maxsize=1000)
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

if __name__ == '__main__':
    pass
