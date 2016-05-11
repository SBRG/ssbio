from Bio import Struct

def count_ss_bond(filename, threshold_new=5):
    """
    Counts the number of sulfide bridges (S-S bonds formed by
    cysteines in close proximity)
    Input: PDB or mmCIF file (will be parsed into a Struct object)
    Output: Int of number of SS bonds
    """
    s = Struct.read(filename)
    my_structure = s.as_protein()
    ss = my_structure.search_ss_bonds(threshold=threshold_new)

    return len(list(ss))


if __name__ == '__main__':
    pass