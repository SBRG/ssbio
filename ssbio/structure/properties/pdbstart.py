import os
from ssbio.structure.pdbioext import PDBIOExt

def get_pdb_res_starts(pdb_file):
    '''
    this uses biopython to get the first residue number in a pdb file.
    returns a list of tuples of pdb, chain, resnum.
    '''
    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model

    start_residues = []
    for chain in model:
        residues = chain.get_residues()
        start_residues.append((os.path.basename(pdb_file), chain.get_id(), next(residues).get_id()[1]))

    return start_residues

if __name__ == '__main__':
    import pandas as pd
    print(pd.DataFrame(get_pdb_res_starts('../test_files/structures/1u8f.pdb')))
