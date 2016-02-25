import os
from loader import Loader
l = Loader()

def get_pdb_res_starts(pdb_file):
    '''
    this uses biopython to get the first residue number in a pdb file.
    returns a list of tuples of pdb, chain, resnum.
    '''
    my_structure = l.structure_reader(pdb_file)

    start_residues = []
    for chain in my_structure[0]:
        residues = chain.get_residues()
        start_residues.append((os.path.basename(pdb_file), chain.get_id(), next(residues).get_id()[1]))

    return start_residues

if __name__ == '__main__':
    print(get_pdb_res_starts('../test_files/structures/1u8f.pdb'))
