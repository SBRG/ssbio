import os
import glob
import pandas as pd
import numpy as np


def sequence_checker(structure_sequence, reference_sequence):



def p2f(x):
    return float(x.strip('%'))/100


def parse_procheck(quality_directory):
    '''
    Parses all PROCHECK files in a directory and returns a Pandas DataFrame of the results
    quality_directory: path to PROCHECK file output (.sum files)
    '''

    procheck_summaries = glob.glob(os.path.join(quality_directory, '*.sum'))

    if len(procheck_summaries) == 0:
        return pd.DataFrame()

    all_procheck = {}
    for summ in procheck_summaries:
        structure_id = os.path.basename(summ).split('.sum')[0]
        procheck_dict = {}
        with open(summ) as f_in:
            lines = (line.rstrip() for line in f_in) # All lines including the blank ones
            lines = (line for line in lines if line) # Non-blank lines
            for line in lines:
                if len(line.split()) > 1:
                    if line.split()[1] == 'Ramachandran':
                        procheck_dict['procheck_rama_favored'] = p2f(line.split()[3])
                        procheck_dict['procheck_rama_allowed'] = p2f(line.split()[5])
                        procheck_dict['procheck_rama_allowed_plus'] = p2f(line.split()[7])
                        procheck_dict['procheck_rama_disallowed'] = p2f(line.split()[9])
                    if line.split()[1] == 'G-factors':
                        procheck_dict['procheck_gfac_dihedrals'] = line.split()[3]
                        procheck_dict['procheck_gfac_covalent'] = line.split()[5]
                        procheck_dict['procheck_gfac_overall'] = line.split()[7]
        all_procheck[structure_id] = procheck_dict

    DF_PROCHECK = pd.DataFrame.from_dict(all_procheck, orient='index')

    return DF_PROCHECK


def parse_psqs(psqs_results_file):
    '''
    Parses a PSQS results file and returns a Pandas DataFrame of the results
    psqs_results_file: full path to psqs results file
    '''

    psqs_results = pd.read_csv(psqs_results_file, sep='\t', header=None)
    psqs_results['pdb_file'] = psqs_results[0].apply(lambda x: str(x).strip('./').strip('.pdb'))
    psqs_results = psqs_results.rename(columns = {1:'psqs_local', 2:'psqs_burial', 3:'psqs_contact', 4:'psqs_total'}).drop(0, axis=1)
    psqs_results['u_pdb'] = psqs_results['pdb_file'].apply(lambda x: x.upper() if len(x)==4 else np.nan)
    psqs_results['i_entry_name'] = psqs_results['pdb_file'].apply(lambda x: x.split('_model1')[0] if len(x)>4 else np.nan)
    psqs_results = psqs_results[pd.notnull(psqs_results.psqs_total)]

    return psqs_results
