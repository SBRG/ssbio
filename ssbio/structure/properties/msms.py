from Bio import PDB
from ssbio.structure.pdbioext import PDBIOExt
import os
import os.path as op
import pandas as pd
import ssbio.utils
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
        return outfile

    else:
        my_structure = PDBIOExt(pdb_file)
        model = my_structure.first_model
        rd = PDB.ResidueDepth(model, pdb_file)

        clean_rd = {}
        for k,v in rd.property_dict.items():
            clean_rd[k[0]] = {k[1][1]: v}

        with open(outfile, 'w') as ff:
            json.dump(clean_rd, ff)

        return outfile

        # Old code to make dataframe of results
        # akeys = list(rd)
        # msms_pre_df = []
        # if len(akeys) == 0:
        #     raise AssertionError('Error running MSMS')
        # else:
        #     # Save the average and alpha-carbon residue depths
        #     for j in akeys:
                    # don't need to use, just check property_dict
        #         chain = [x.id for x in PDB.Selection.unfold_entities(j[0], 'C')][0]
        #         seqnum = j[0].id[1]
        #         re_depth = j[1][0]
        #         ca_depth = j[1][1]
        #
        #         if chain == ' ':
        #             chain = 'X'
        #
        #         msms_pre_df.append([chain, seqnum, re_depth, ca_depth])
        #
        # cols = ['chain', 'resnum', 'res_depth', 'ca_depth']
        # msms_df = pd.DataFrame.from_records(msms_pre_df)
        # msms_df = msms_df.rename(columns=cols)
        #
        # msms_df.to_csv(outfile)
        #
        # return outfile

def residue_depth(anslist):
    '''Computes the average residue and CA depth
    Returns a dictionary of "redepth" and "cadepth" (floats)
    '''
    depth_info = {}

    if anslist != ['NA', 'NA', 'NA', 'NA']:
        redepth = [x[2] for x in anslist]
        cadepth = [x[3] for x in anslist]
        redepth = [x for x in redepth if x is not None]
        cadepth = [x for x in cadepth if x is not None]
        depth_info['ssb_avg_res_depth'] = sum(redepth) / len(redepth)
        depth_info['ssb_ca_depth'] = sum(cadepth) / len(cadepth)
    else:
        depth_info['ssb_avg_res_depth'] = None
        depth_info['ssb_ca_depth'] = None

    return depth_info

if __name__ == '__main__':
    import argparse
    from tqdm import tqdm
    from ssbio import utils
    import pandas as pd

    p = argparse.ArgumentParser(description='Runs MSMS on a PDB file or folder')
    p.add_argument('infile', help='PDB file or folder', type=str)
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

        red_dict = residue_depth(msms_stuff)
        red_dict['ssb_file'] = f
        redinfo.append(red_dict)

    DF_PROP_MSMS = pd.DataFrame(msmsinfo)
    DF_PROP_MSMS.columns = ['ssb_file', 'ssb_msms']
    DF_PROP_MSMS.to_csv('DF_PROP_MSMS.csv')

    print('Errors with: {}'.format(msms_errors))
    print('Saved DF at: {}'.format('DF_PROP_MSMS.csv'))