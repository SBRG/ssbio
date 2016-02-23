import os
from Bio import PDB
from tqdm import tqdm
import pandas as pd

from Bio.PDB.DSSP import *
from Bio.PDB.Polypeptide import one_to_three
from Loader import Loader
l = Loader()


def calc_sasa(filename):
    """
    Calculation of SASA utilizing the DSSP program.

    DSSP must be installed for biopython to properly call it.
    Install using apt-get on Ubuntu
    or from: http://swift.cmbi.ru.nl/gv/dssp/

    Input: PDB or CIF structure file
    Output: SASA (integer) of structure
    """

    my_structure = l.structure_reader(filename)

    model = my_structure[0]
    try:
        dssp = PDB.DSSP(model, filename)
    except:
        return {}
    akeys = list(dssp)

    if len(akeys) == 0:
        akeys = [0] * 14

    df = pd.DataFrame(akeys)

    df.columns = ['dssp_index', 'aa', 'ss', 'relative_expo', 'phi', 'psi', 'NH_O_1_relidx', 'NH_O_1_energy',
                  'O_NH_1_relidx', 'O_NH_1_energy', 'NH_O_2_relidx', 'NH_O_2_energy', 'O_NH_2_relidx', 'O_NH_2_energy']

    droplist = []
    for j in df.index:
        if df.loc[j, 'relative_expo'] == 'NA':
            droplist.append(j)
    if len(droplist) > 0:
        dfad = df.drop(df.index[droplist])
    else:
        dfad = df

    dfad['aa_three'] = dfad['aa'].apply(one_to_three)
    dfad['max_acc'] = dfad['aa_three'].map(MAX_ACC.get)
    dfad[['relative_expo', 'max_acc']] = dfad[['relative_expo', 'max_acc']].astype(float)
    dfad['exposure_area'] = dfad['relative_expo'] * dfad['max_acc']

    infodict = {'ssb_sasa': dfad.exposure_area.sum(),
                'ssb_mean_rel_exposed': dfad.relative_expo.mean(),
                'ssb_size': len(akeys)}

    return infodict



if __name__ == '__main__':
    import glob
    files = glob.glob('test_structures/*')
    print(calc_sasa(files))
