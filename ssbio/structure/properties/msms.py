from Bio import PDB
from ssbio.structure.pdbioext import PDBIOExt
from collections import defaultdict
import os
import os.path as op
import pandas as pd
import ssbio.utils
import logging
log = logging.getLogger(__name__)
import json


def msms_output(pdb_file, outfile='', outdir='', outext='_msms.json', force_rerun=False):
    """Run MSMS (through Biopython) on a PDB file.

    Saves the result dataframe for the PDB file. The dataframe contains:
        chain: chain ID
        resnum: residue number (PDB numbering)
        res_depth: average depth of all atoms in a residue
        ca_depth: depth of the alpha carbon atom

    Depths are in units Angstroms. 1A = 10^-10 m = 1nm

    Args:
        pdb_file: path to PDB file

    Returns:
        str: path to saved Pandas DataFrame of residue and alpha-carbon depths.

    """
    # Create the output file name
    outfile = ssbio.utils.outfile_name_maker(inname=pdb_file, outfile=outfile, outdir=outdir, outext=outext)

    if not ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        log.debug('{}: already ran MSMS and force_rerun={}'.format(outfile, force_rerun))
        return outfile

    else:
        my_structure = PDBIOExt(pdb_file)
        model = my_structure.first_model
        rd = PDB.ResidueDepth(model, pdb_file)

        clean_rd = defaultdict(dict)
        for k,v in rd.property_dict.items():
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
    import argparse
    from tqdm import tqdm
    from ssbio import utils
    import pandas as pd

    p = argparse.ArgumentParser(description='Runs MSMS on a PDB file or folder')
    p.add_argument('infile', help='PDB file or folder', type=str, nargs='+')
    args = p.parse_args()

    print(args.infile)
    infiles = utils.input_list_parser(args.infile)

    msmsinfo = []
    redinfo = []

    msms_errors = []

    for f in tqdm(infiles, leave=True):
        try:
            msms_stuff = msms_output(f)
        except:
            msms_errors.append(f)
            continue

        msmsinfo.append([f, msms_stuff])

    # TODO: change behavior of dataframe stuff and printing errors
    DF_PROP_MSMS = pd.DataFrame(msmsinfo)
    DF_PROP_MSMS.columns = ['ssb_file', 'ssb_msms']
    DF_PROP_MSMS.to_csv('DF_PROP_MSMS.csv')

    print('Errors with: {}'.format(msms_errors))
    print('Saved DF at: {}'.format('DF_PROP_MSMS.csv'))