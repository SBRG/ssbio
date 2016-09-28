import math
import cachetools
import scipy.constants

# author: Ke Chen

# R = scipy.constants.R # J/K.mol
# R is the molar gas constant in cal/K.mol
r_cal = scipy.constants.R / scipy.constants.calorie

# Oobatake dG
# dG and dCp from Table 8 in Oobatake paper.
# dG,dH in unit kcal/mol
# dCp,dS in unit cal/mol.K
oobatake_dictionary = {}
oobatake_dictionary['A'] = {'dG': -0.02, 'dCp': 14.22, 'dH': 0.51, 'dS': 1.82}
oobatake_dictionary['C'] = {'dG': 1.08, 'dCp': 9.41, 'dH': 5.21, 'dS': 13.85}
oobatake_dictionary['D'] = {'dG': -0.08, 'dCp': 2.73, 'dH': 0.18, 'dS': 0.86}
oobatake_dictionary['E'] = {'dG': -0.13, 'dCp': 3.17, 'dH': 0.05, 'dS': 0.65}
oobatake_dictionary['F'] = {'dG': 2.16, 'dCp': 39.06, 'dH': 6.82, 'dS': 15.64}
oobatake_dictionary['G'] = {'dG': 0.09, 'dCp': 4.88, 'dH': -0.23, 'dS': -1.07}
oobatake_dictionary['H'] = {'dG': 0.56, 'dCp': 20.05, 'dH': 0.79, 'dS': 0.75}
oobatake_dictionary['I'] = {'dG': -0.08, 'dCp': 41.98, 'dH': 0.19, 'dS': 0.89}
oobatake_dictionary['K'] = {'dG': -0.32, 'dCp': 17.68, 'dH': -1.45, 'dS': -3.78}
oobatake_dictionary['L'] = {'dG': -0.08, 'dCp': 38.26, 'dH': 0.17, 'dS': 0.81}
oobatake_dictionary['M'] = {'dG': 0.53, 'dCp': 31.67, 'dH': 2.89, 'dS': 7.93}
oobatake_dictionary['N'] = {'dG': -0.30, 'dCp': 3.91, 'dH': -2.03, 'dS': -5.83}
oobatake_dictionary['P'] = {'dG': -0.06, 'dCp': 23.69, 'dH': 0.02, 'dS': 0.30}
oobatake_dictionary['Q'] = {'dG': -0.23, 'dCp': 3.74, 'dH': -1.76, 'dS': -5.13}
oobatake_dictionary['R'] = {'dG': -0.71, 'dCp': 16.66, 'dH': -4.40, 'dS': -12.38}
oobatake_dictionary['S'] = {'dG': -0.40, 'dCp': 6.14, 'dH': -0.16, 'dS': 0.84}
oobatake_dictionary['T'] = {'dG': -0.24, 'dCp': 16.11, 'dH': 0.04, 'dS': 0.89}
oobatake_dictionary['V'] = {'dG': -0.06, 'dCp': 32.58, 'dH': 0.30, 'dS': 1.20}
oobatake_dictionary['W'] = {'dG': 1.78, 'dCp': 37.69, 'dH': 4.47, 'dS': 9.00}
oobatake_dictionary['Y'] = {'dG': 0.91, 'dCp': 30.54, 'dH': 3.73, 'dS': 9.46}
oobatake_dictionary['U'] = {'dG': 1.08, 'dCp': 9.41, 'dH': 5.21, 'dS': 13.85}  # assume U==C


def sum_of_dCp(seq):
    dCp_sum = 0
    for aa in seq:
        dCp_sum += oobatake_dictionary[aa]['dCp']
    return dCp_sum


def calculate_oobatake_dH(seq, temp):
    """Get dH in units cal/mol

    Args:
        seq: amino acid sequence
        temp: temperature in degrees C

    Returns:
        dH

    """
    dH = 0
    temp = temp + 273.15
    T0 = 298.15
    for aa in seq:
        H0 = oobatake_dictionary[aa]['dH'] * 1000
        dH += H0
    return dH + sum_of_dCp(seq) * (temp - T0)


def calculate_oobatake_dS(seq, temp):
    """Get dS in units cal/mol

    Args:
        seq: amino acid sequence
        temp: temperature in degrees C

    Returns:
        dS

    """
    dS = 0
    temp = temp + 273.15
    T0 = 298.15
    dCp_sum = sum_of_dCp(seq)
    for aa in seq:
        S0 = oobatake_dictionary[aa]['dS']
        dS += S0
    return dS + dCp_sum * math.log(temp / T0)


@cachetools.func.ttl_cache(maxsize=256)
def calculate_oobatake_dG(seq, temp):
    """Get dG in units cal/mol

    Args:
        seq: amino acid sequence
        temp: temperature in degrees C

    Returns:
        dG

    """
    T0 = 298.15
    dH = calculate_oobatake_dH(seq, temp)
    dS = calculate_oobatake_dS(seq, temp)
    dG = dH - (temp + 273.15) * dS

    # 563.552 - a correction for N- and C-terminal group (approximated from 7 examples in the paper)
    return dG - 563.552


@cachetools.func.ttl_cache(maxsize=128)
def calculate_dill_dG(seq_len, temp):
    """Get dG using Dill method in units J/mol

    Args:
        seq_len: length of amino acid sequence
        temp: temperature in degrees C

    Returns:
        dG in J/mol
    """
    Th = 373.5  # this quantity affects the up-and-down of the dG vs temperature curve (dG values)
    Ts = 385  # this quantity affects the left-and-right
    temp = temp + 273.15
    dH = (4.0 * seq_len + 143) * 1000
    dS = 13.27 * seq_len + 448
    dCp = (0.049 * seq_len + 0.85) * 1000
    dG = dH + dCp * (temp - Th) - temp * dS - temp * dCp * math.log(float(temp) / Ts)
    return dG

@cachetools.func.ttl_cache(maxsize=500)
def get_dG_at_T(seq, temp):
    """Predict dG at temperature T

    Args:
        seq: sequence string
        temp: temperature

    Returns:
        free energy of unfolding dG (cal/mol)
        equilibrium constant Keq (for ME2.0 coupling constraint)
        method used to calculate
    """

    oobatake = {}
    for t in range(20, 51):
        oobatake[t] = calculate_oobatake_dG(seq, t)

    stable = [i for i in oobatake.values() if i > 0]
    if len(stable) == 0:
        # if oobatake dG < 0 for all tempertures [20,50], use Dill dG
        # and convert the number from J/mol to cal/mol
        dG = 0.238846 * calculate_dill_dG(len(seq), temp)
        method='Dill'
    else:
        dG = oobatake[temp]
        method='Oobatake'

    keq = math.exp(-1 * dG / (r_cal * (temp + 273.15)))

    return dG, keq, method

if __name__ == '__main__':
    # load inputs from command line
    pass