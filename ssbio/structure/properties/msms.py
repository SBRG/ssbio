import json
import argparse
import logging
from Bio import PDB
from collections import defaultdict
from tqdm import tqdm
import pandas as pd

import ssbio.utils
from ssbio.structure.pdbioext import PDBIOExt

log = logging.getLogger(__name__)


def run_msms(pdb_file, outfile='', outdir='', outext='_msms.json', force_rerun=False):
    """Run MSMS (using Biopython) on a PDB file.

    Saves a JSON file for the PDB. The file contains:
        chain: chain ID
        resnum: residue number (PDB numbering)
        res_depth: average depth of all atoms in a residue
        ca_depth: depth of the alpha carbon atom

    It is formatted like so:
        {chain: {resnum1: [res_depth, ca_depth]},
                {resnum2: [res_depth, ca_depth]} }

    Depths are in units Angstroms. 1A = 10^-10 m = 1nm

    Args:
        pdb_file: Path to PDB file
        outfile: Optional name of output file (without extension)
        outdir: Optional output directory
        outext: Optional extension for the output file
        force_rerun: Rerun MSMS even if results exist already

    Returns:
        str: path to saved json file of residue and alpha-carbon depths

    """
    # Create the output file name
    outfile = ssbio.utils.outfile_name_maker(inname=pdb_file, outfile=outfile, outdir=outdir, outext=outext)

    if not ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        log.debug('{}: already ran MSMS and force_rerun={}'.format(outfile, force_rerun))
        return outfile

    else:
        # Load the structure
        my_structure = PDBIOExt(pdb_file)
        model = my_structure.first_model

        # Run MSMS with Biopython
        rd = PDB.ResidueDepth(model, pdb_file)

        # Create a json file of the results
        clean_rd = defaultdict(dict)
        for k, v in rd.property_dict.items():
            clean_rd[k[0]].update({k[1][1]: v})

        with open(outfile, 'w') as ff:
            json.dump(clean_rd, ff)

        return outfile


# TODO: What distance defines a surface/buried residue?
# TODO: function to take a list of residue numbers and return a dictionary of their location
def resnum_list_to_exposure(chain_id, residue_numbers, msms_json, exposed_cutoff=2.5):
    """Get a dictionary defining each residue number as surface or buried.

    Args:
        chain_id:
        residue_numbers:
        msms_json:
        exposed_cutoff:

    Returns:

    """
    results = {}

    with open(msms_json, 'r') as f:
        msms = json.load(f)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Run MSMS to calculate residue depth on a file, files, or a directory. Save a JSON file of results per structure.')
    p.add_argument('infile', help='PDB file, files, or folder', type=str, nargs='+')
    p.add_argument('--summary', action='store_true', help='Save a summary DataFrame of results.')

    args = p.parse_args()
    infiles = ssbio.utils.input_list_parser(args.infile)

    msms_errors = []

    for f in tqdm(infiles, leave=True):
        try:
            msms_stuff = run_msms(f)
        except:
            msms_errors.append(f)
            continue

    if args.summary:
        pass
        # TODO: what to save as a summary? Average residue depth of entire protein?
        log.info('Saved DF at: {}'.format('{}-df_msms_summary.csv'.format(ssbio.utils.Date().short_date)))

    if msms_errors:
        log.warning('Errors with: {}'.format(msms_errors))