import os
from tqdm import tqdm
import pandas as pd
from loader import Loader
l = Loader()


def count_ss_bond(filename, threshold_new=5):
    """
    Counts the number of sulfide bridges (S-S bonds formed by
    cysteines in close proximity)
    Input: PDB or mmCIF file (will be parsed into a Struct object)
    Output: Int of number of SS bonds
    """
    my_structure = l.structure_reader(filename)
    my_structure = my_structure.as_protein()
    ss = my_structure.search_ss_bonds(threshold=threshold_new)

    return len(list(ss))


if __name__ == '__main__':
    import glob
    files = glob.glob('../test_files/structures/*')
    print(files)
    for f in files:
        print(f)
        print(count_ss_bond(f))
