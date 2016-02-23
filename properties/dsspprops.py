import os
from Bio import PDB
from tqdm import tqdm
import pandas as pd

from Bio.PDB.DSSP import *
from Bio.PDB.Polypeptide import one_to_three
from Loader import Loader
l = Loader()


def dssp_dataframe(filename):
    """
    Calculation of various properties utilizing the DSSP program.

    DSSP must be installed for biopython to properly call it.
    Install using apt-get on Ubuntu
    or from: http://swift.cmbi.ru.nl/gv/dssp/

    Input: PDB or CIF structure file
    Output: Dictionary of properties
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
    df.columns = ['dssp_index', 'aa', 'ss', 'relative_expo', 'phi', 'psi',
                  'NH_O_1_relidx', 'NH_O_1_energy', 'O_NH_1_relidx',
                  'O_NH_1_energy', 'NH_O_2_relidx', 'NH_O_2_energy',
                  'O_NH_2_relidx', 'O_NH_2_energy']

    return df


def calc_sasa(dssp_df):
    """
    Calculation of SASA utilizing the DSSP program.

    DSSP must be installed for biopython to properly call it.
    Install using apt-get on Ubuntu
    or from: http://swift.cmbi.ru.nl/gv/dssp/

    Input: PDB or CIF structure file
    Output: SASA (integer) of structure
    """
    droplist = []
    for j in dssp_df.index:
        if dssp_df.loc[j, 'relative_expo'] == 'NA':
            droplist.append(j)
    if len(droplist) > 0:
        dfad = dssp_df.drop(dssp_df.index[droplist])
    else:
        dfad = dssp_df

    dfad['aa_three'] = dfad['aa'].apply(one_to_three)
    dfad['max_acc'] = dfad['aa_three'].map(MAX_ACC.get)
    dfad[['relative_expo', 'max_acc']] = dfad[['relative_expo', 'max_acc']].astype(float)
    dfad['exposure_area'] = dfad['relative_expo'] * dfad['max_acc']

    infodict = {'ssb_sasa': dfad.exposure_area.sum(),
                'ssb_mean_rel_exposed': dfad.relative_expo.mean(),
                'ssb_size': len(dssp_df)}

    return infodict


def calc_sec_struct_composition(dssp_df):
    '''Calculates the percent of residues that are in a certain
    secondary structure. Returns a dictionary of this.
    '''
    expoinfo = {}
    Hn = 0
    Bn = 0
    En = 0
    Gn = 0
    In = 0
    Tn = 0
    Sn = 0
    Nn = 0
    Total = 0
    abh = dssp_df.ss.tolist()
    if len(abh) == 0:
        expoinfo = {}
    else:
        for i in abh:
            if i == '-':
                Nn = Nn + 1
            elif i == 'H':
                Hn = Hn + 1
            elif i == 'B':
                Bn = Bn + 1
            elif i == 'E':
                En = En + 1
            elif i == 'G':
                Gn = Gn + 1
            elif i == 'I':
                In = In + 1
            elif i == 'T':
                Tn = Tn + 1
            elif i == 'S':
                Sn = Sn + 1
            Total = Total + 1
        Nnp = float(Nn) / float(Total)
        Hnp = float(Hn) / float(Total)
        Bnp = float(Bn) / float(Total)
        Enp = float(En) / float(Total)
        Gnp = float(Gn) / float(Total)
        Inp = float(In) / float(Total)
        Tnp = float(Tn) / float(Total)
        Snp = float(Sn) / float(Total)
        expoinfo['ssb_per_irr'] = Nnp
        expoinfo['ssb_per_alpha'] = Hnp
        expoinfo['ssb_per_beta_bridge'] = Bnp
        expoinfo['ssb_per_ext_beta'] = Enp
        expoinfo['ssb_per_310_helix'] = Gnp
        expoinfo['ssb_per_5_helix'] = Inp
        expoinfo['ssb_per_hbond_turn'] = Tnp
        expoinfo['ssb_per_bent'] = Snp

    return expoinfo


if __name__ == '__main__':
    import glob
    files = glob.glob('test_structures/*')
    for f in files:
        t = dssp_dataframe(f)
        print(f)
        print(calc_sasa(t))
        print(calc_sec_struct_composition(t))
