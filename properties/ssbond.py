import os
from tqdm import tqdm
import pandas as pd
from Loader import Loader
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


def count_ss_bond_batch(filenames, threshold_new=5):
    if isinstance(filenames, str):
        filenames = [filenames]

    results_list = []
    for f in tqdm(filenames):
        results = {}
        results['ssb_file'] = os.path.basename(f)
        results['ssb_cys_bridge'] = count_ss_bond(f)
        results_list.append(results)

    results_df = pd.DataFrame(results_list)

    return results_df

if __name__ == '__main__':
    import glob
    files = glob.glob('properties/test_structures/*')
    print(count_ss_bond_batch(files))
