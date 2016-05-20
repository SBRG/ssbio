import math

#R = 8.3144621 # J/K.mol
R = 1.9872041 # cal/K.mol

def get_dG_at_T(gene,T):
    """Predict dG at temperature T

    Args:
        gene: b number 
        T:    temperature

    Returns:
        free energy of unfolding dG (cal/mol)
        equilibrium constant Keq (for ME2.0 coupling constraint)

    """
    #seq = f(gene) # from ME2.0 or from structure

    oobatake = {}
    for t in range(20,51):
        oobatake[t] = calculate_Oobatake_dG(seq,t)

    stable=[i for i in oobatake.values() if i > 0]
    if len(stable) == 0:
        print "Use Dill!"
        # if oobatake dG < 0 for all tempertures [20,50], use Dill dG
        # and convert the number from J/mol to cal/mol
        dG = 0.238846 * calculate_Dill_dG(len(seq),T)
    else:
        dG = oobatake[T]

    keq = math.exp(-1*dG/(R*(T+273.15)))

    return dG,keq

## Dill dG 
def calculate_Dill_dG(N,T): #T in unit degreeC 
    # unit in J/mol
    Th = 373.5 # this quantity affects the up-and-down of the dG vs temperature curve (dG values)
    Ts = 385 # this quantity affects the left-and-right
    T = T + 273.15
    dH = (4.0*N+143)*1000
    dS = 13.27*N+448
    dCp = (0.049*N+0.85)*1000
    dG = dH + dCp*(T-Th) - T*dS - T*dCp*math.log(float(T)/Ts)
    return dG


## Oobatake dG
# dG and dCp from Table 8 in Oobatake paper.
# dG,dH in unit kcal/mol
# dCp,dS in unit cal/mol.K
Oobatake_dictionary = {} 
Oobatake_dictionary['A']={'dG':-0.02,'dCp':14.22,'dH':0.51,'dS':1.82}
Oobatake_dictionary['C']={'dG':1.08,'dCp':9.41,'dH':5.21,'dS':13.85}
Oobatake_dictionary['D']={'dG':-0.08,'dCp':2.73,'dH':0.18,'dS':0.86}
Oobatake_dictionary['E']={'dG':-0.13,'dCp':3.17,'dH':0.05,'dS':0.65}
Oobatake_dictionary['F']={'dG':2.16,'dCp':39.06,'dH':6.82,'dS':15.64}
Oobatake_dictionary['G']={'dG':0.09,'dCp':4.88,'dH':-0.23,'dS':-1.07}
Oobatake_dictionary['H']={'dG':0.56,'dCp':20.05,'dH':0.79,'dS':0.75}
Oobatake_dictionary['I']={'dG':-0.08,'dCp':41.98,'dH':0.19,'dS':0.89}
Oobatake_dictionary['K']={'dG':-0.32,'dCp':17.68,'dH':-1.45,'dS':-3.78}
Oobatake_dictionary['L']={'dG':-0.08,'dCp':38.26,'dH':0.17,'dS':0.81}
Oobatake_dictionary['M']={'dG':0.53,'dCp':31.67,'dH':2.89,'dS':7.93}
Oobatake_dictionary['N']={'dG':-0.30,'dCp':3.91,'dH':-2.03,'dS':-5.83}
Oobatake_dictionary['P']={'dG':-0.06,'dCp':23.69,'dH':0.02,'dS':0.30}
Oobatake_dictionary['Q']={'dG':-0.23,'dCp':3.74,'dH':-1.76,'dS':-5.13}
Oobatake_dictionary['R']={'dG':-0.71,'dCp':16.66,'dH':-4.40,'dS':-12.38}
Oobatake_dictionary['S']={'dG':-0.40,'dCp':6.14,'dH':-0.16,'dS':0.84}
Oobatake_dictionary['T']={'dG':-0.24,'dCp':16.11,'dH':0.04,'dS':0.89}
Oobatake_dictionary['V']={'dG':-0.06,'dCp':32.58,'dH':0.30,'dS':1.20}
Oobatake_dictionary['W']={'dG':1.78,'dCp':37.69,'dH':4.47,'dS':9.00}
Oobatake_dictionary['Y']={'dG':0.91,'dCp':30.54,'dH':3.73,'dS':9.46}
# assume U==C
Oobatake_dictionary['U']={'dG':1.08,'dCp':9.41,'dH':5.21,'dS':13.85}

def sum_of_dCp(seq):
    dCp_sum = 0
    for aa in seq:
        dCp_sum += Oobatake_dictionary[aa]['dCp']
    return dCp_sum

def calculate_Oobatake_dH(seq,T):
    # return dH in unit cal/mol
    dH = 0
    T = T + 273.15
    T0 = 298.15
    for aa in seq:
        H0 = Oobatake_dictionary[aa]['dH']*1000
        dH += H0
    return dH + sum_of_dCp(seq)*(T-T0)

def calculate_Oobatake_dS(seq,T):
    # return dS in unit cal/mol
    dS = 0
    T = T + 273.15
    T0 = 298.15
    dCp_sum = sum_of_dCp(seq)
    for aa in seq:
        S0 = Oobatake_dictionary[aa]['dS']
        dS += S0
    return dS + dCp_sum*math.log(T/T0)

def calculate_Oobatake_dG(seq,T): #T in unit degreeC 
    # return dG in unit cal/mol
    T0 = 298.15
    dH = calculate_Oobatake_dH(seq,T)
    dS = calculate_Oobatake_dS(seq,T)
    dG = dH - (T+273.15)*dS
    return dG - 563.552 # a correction for N- and C-terminal group (approximated from 7 examples in the paper)

