import ssbio.utils
import os.path as op
from bs4 import BeautifulSoup
import itertools
from tqdm import tqdm
import numpy as np
import os
import pandas as pd

__author__ = "Anand Sastry"
__email__ = "avsastry@eng.ucsd.edu"


def run_fatcat(structure_path_1, structure_path_2, fatcat_sh, outdir='', silent=False, print_cmd=False, force_rerun=False):
    """Run FATCAT on two PDB files, and return the path of the XML result file.
    
    Args:
        structure_path_1 (str): Path to PDB file 
        structure_path_2 (str): Path to PDB file 
        fatcat_sh (str): Path to "runFATCAT.sh" executable script
        outdir (str): Path to where FATCAT XML output files will be saved
        silent (bool): If stdout should be silenced from showing up in Python console output
        print_cmd (bool): If command to run FATCAT should be printed to stdout
        force_rerun (bool): If FATCAT should be run even if XML output files already exist 

    Returns:
        str: Path to XML output file

    """
    filename1 = op.splitext(op.basename(structure_path_1))[0]
    filename2 = op.splitext(op.basename(structure_path_2))[0]

    if not op.exists(outdir):
        os.mkdir(outdir)
    outfile = op.join(outdir, filename1 + '__' + filename2 + '.xml')

    # Run FATCAT on the structures, print the XML of the result to stdout
    fatcat_cmd = '{} -file1 {} -file2 {} -outFile {}'.format(fatcat_sh, structure_path_1, structure_path_2, outfile)
    if print_cmd:
        print(fatcat_cmd)
    ssbio.utils.command_runner(fatcat_cmd, force_rerun_flag=force_rerun, outfile_checker=outfile, silent=silent)

    return outfile


def run_fatcat_all_by_all(list_of_structure_paths, fatcat_sh, outdir='', silent=True, force_rerun=False):
    """Run FATCAT on all pairs of structures given a list of structures.
    
    Args:
        list_of_structure_paths (list): List of PDB file paths 
        fatcat_sh (str): Path to "runFATCAT.sh" executable script
        outdir (str): Path to where FATCAT XML output files will be saved 
        silent (bool): If command to run FATCAT should be printed to stdout
        force_rerun (bool): If FATCAT should be run even if XML output files already exist 

    Returns:
        Pandas DataFrame: TM-scores (similarity) between all structures

    """
    structure_ids = {x: i for i, x in enumerate(list_of_structure_paths)}

    comps = itertools.combinations(list_of_structure_paths, 2)
    tm_score_matrix = np.eye(len(list_of_structure_paths))

    for pdb1, pdb2 in tqdm(comps):
        fatcat_file = run_fatcat(pdb1, pdb2, fatcat_sh, outdir=outdir, silent=silent, force_rerun=force_rerun)
        tm_score_matrix[structure_ids[pdb1], structure_ids[pdb2]] = parse_fatcat(fatcat_file)['tm_score']
        tm_score_matrix[structure_ids[pdb2], structure_ids[pdb1]] = parse_fatcat(fatcat_file)['tm_score']

    # Convert to dataframe with filenames
    filenames = [op.splitext(op.basename(x))[0] for x in list_of_structure_paths]
    tm_score_matrix_annotated = pd.DataFrame(data=tm_score_matrix, columns=filenames, index=filenames)

    return tm_score_matrix_annotated


def parse_fatcat(fatcat_xml):
    """Parse a FATCAT XML result file.
    
    Args:
        fatcat_xml (str): Path to FATCAT XML result file 

    Returns:
        dict: Parsed information from the output
        
    Todo:
        - Only returning TM-score at the moment

    """
    fatcat_results = {}

    # Parse output xml file
    with open(fatcat_xml, 'r') as f:
        soup = BeautifulSoup(f, 'lxml')

    # Find the tmScore of the alignment
    if soup.find('block'):
        fatcat_results['tm_score'] = float(soup.find('afpchain')['tmscore'])

    return fatcat_results
