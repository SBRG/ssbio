from Bio import PDB
import pandas as pd
import prody as pr

from Bio.PDB.DSSP import *
from Bio.PDB.Polypeptide import aa1
from Bio.PDB.Polypeptide import one_to_three
from ssbio.structure.pdbioext import PDBIOExt
import cachetools

AAdict = {'CYS': 'polar',
          'ILE': 'nonpolar',
          'GLY': 'nonpolar',
          'SER': 'polar',
          'GLN': 'polar',
          'LYS': 'positive',
          'ASN': 'polar',
          'PRO': 'nonpolar',
          'ASP': 'negative',
          'THR': 'polar',
          'PHE': 'nonpolar',
          'ALA': 'nonpolar',
          'MET': 'nonpolar',
          'HIS': 'positive',
          'LEU': 'nonpolar',
          'ARG': 'positive',
          'TRP': 'nonpolar',
          'VAL': 'nonpolar',
          'GLU': 'negative',
          'TYR': 'polar',
          'MSE': 'polar',
          'SEC': 'polar'}

@cachetools.func.ttl_cache(maxsize=300)
def dssp_dataframe(filename):
    """
    Calculation of various properties utilizing the DSSP program.

    DSSP must be installed for biopython to properly call it.
    Install using apt-get on Ubuntu
    or from: http://swift.cmbi.ru.nl/gv/dssp/

    Input: PDB or CIF structure file
    Output: Dictionary of properties
    """

    my_structure = PDBIOExt(filename)
    model = my_structure.first_model
    try:
        dssp = PDB.DSSP(model, filename)
    except:
        return pd.DataFrame()
    akeys = list(dssp)

    if len(akeys) == 0:
        akeys = [0] * 14

    df = pd.DataFrame(akeys)
    df.columns = ['dssp_index', 'aa', 'ss', 'relative_expo', 'phi', 'psi',
                  'NH_O_1_relidx', 'NH_O_1_energy', 'O_NH_1_relidx',
                  'O_NH_1_energy', 'NH_O_2_relidx', 'NH_O_2_energy',
                  'O_NH_2_relidx', 'O_NH_2_energy']

    df = df[df['aa'].isin(list(aa1))]
    df['aa_three'] = df['aa'].apply(one_to_three)
    df['max_acc'] = df['aa_three'].map(MAX_ACC.get)
    df[['relative_expo', 'max_acc']] = df[
        ['relative_expo', 'max_acc']].astype(float)
    df['exposure_area'] = df['relative_expo'] * df['max_acc']

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

    infodict = {'ssb_sasa': dssp_df.exposure_area.sum(),
                'ssb_mean_rel_exposed': dssp_df.relative_expo.mean(),
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


def calc_surface_buried(dssp_df):
    '''Calculates the percent of residues that are in the surface or buried,
    as well as if they are polar or nonpolar. Returns a dictionary of this.
    '''
    SN = 0
    BN = 0
    SP = 0
    SNP = 0
    SPo = 0
    SNe = 0
    BNP = 0
    BP = 0
    BPo = 0
    BNe = 0
    Total = 0

    sbinfo = {}

    df_min = dssp_df[['aa_three', 'exposure_area']]

    if len(df_min) == 0:
        return sbinfo
    else:
        for i, r in df_min.iterrows():
            res = r.aa_three
            area = r.exposure_area
            if res in AAdict:
                if AAdict[res] == 'nonpolar' and area > 3:
                    SNP = SNP + 1
                    SN = SN + 1
                elif AAdict[res] == 'polar' and area > 3:
                    SP = SP + 1
                    SN = SN + 1
                elif AAdict[res] == 'positive' and area > 3:
                    SPo = SPo + 1
                    SN = SN + 1
                elif AAdict[res] == 'negative' and area > 3:
                    SNe = SNe + 1
                    SN = SN + 1
                elif AAdict[res] == 'positive' and area <= 3:
                    BPo = BPo + 1
                    BN = BN + 1
                elif AAdict[res] == 'negative' and area <= 3:
                    BNe = BNe + 1
                    BN = BN + 1
                elif AAdict[res] == 'polar' and area <= 3:
                    BP = BP + 1
                    BN = BN + 1
                elif AAdict[res] == 'nonpolar' and area <= 3:
                    BNP = BNP + 1
                    BN = BN + 1
        Total = float(BN + SN)
        pSNP = float(SNP) / Total
        pSP = float(SP) / Total
        pSPo = float(SPo) / Total
        pSNe = float(SNe) / Total
        pBNP = float(BNP) / Total
        pBP = float(BP) / Total
        pBPo = float(BPo) / Total
        pBNe = float(BNe) / Total
        pBN = float(BN) / Total
        pSN = float(SN) / Total
        sbinfo['ssb_per_S_NP'] = pSNP
        sbinfo['ssb_per_S_P'] = pSP
        sbinfo['ssb_per_S_pos'] = pSPo
        sbinfo['ssb_per_S_neg'] = pSNe
        sbinfo['ssb_per_B_NP'] = pBNP
        sbinfo['ssb_per_B_P'] = pBP
        sbinfo['ssb_per_B_pos'] = pBPo
        sbinfo['ssb_per_B_neg'] = pBNe
        sbinfo['ssb_per_S'] = pSN
        sbinfo['ssb_per_B'] = pBN

        return sbinfo


def all_dssp_props(filename):
    '''Returns a large dictionary of SASA, secondary structure
    composition, and surface/buried composition. Values are computed using DSSP.
    Input: PDB or MMCIF filename
    Output: Dictionary of values obtained from dssp
    '''
    t = dssp_dataframe(filename)
    # print(t)
    sasa = calc_sasa(t)
    sstr = calc_sec_struct_composition(t)
    subu = calc_surface_buried(t)

    sasa.update(sstr)
    sasa.update(subu)

    return sasa

if __name__ == '__main__':
    import glob
    files = glob.glob('../test_files/structures/*')
    # print(files)
    for f in files:
        print(f)
        # print(all_dssp_props(f))


# TODO: convert these functions to use biopython?
def get_dssp_ss_content_multiplechains(prody_ag, chain):
    """Get secondary structure across chain

    Args:
        prody_ag: ProDy atomgroup object (parsed PDB file)
        chain (str): chain ID

    Returns:

    """
    a = prody_ag.getData('resnum')
    b = prody_ag.getData('name')
    c = prody_ag.getData('secondary')
    d = prody_ag.getData('chain')
    resid = []
    SSstr = []
    idx = [i for i, x in enumerate(d) if x == chain]
    for i in idx:
        if b[i] == 'CA':
            resid.append(i)
            SSstr.append(c[i])

    N = float(len(SSstr))
    helix_alpha = SSstr.count('H') / N
    helix_3_10 = SSstr.count('G') / N
    extended = SSstr.count('E') / N
    return helix_alpha, helix_3_10, extended


def get_ss_class(pdb_file, dssp_file, chain):
    """Define the secondary structure class of a PDB file at the specific chain

    Args:
        pdb_file:
        dssp_file:
        chain:

    Returns:

    """
    prag = pr.parsePDB(pdb_file)
    pr.parseDSSP(dssp_file, prag)
    alpha, threeTen, beta = get_dssp_ss_content_multiplechains(prag, chain)

    if alpha == 0 and beta > 0:
        classification = 'all-beta'
    elif beta == 0 and alpha > 0:
        classification = 'all-alpha'
    elif beta == 0 and alpha == 0:
        classification = 'mixed'
    elif float(alpha) / beta >= 20:
        classification = 'all-alpha'
    else:
        classification = 'mixed'

    return classification