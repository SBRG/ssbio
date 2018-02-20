import argparse
import logging
import os.path as op
import pandas as pd
from Bio import PDB
from tqdm import tqdm

import ssbio.utils
from ssbio.protein.structure.utils.structureio import StructureIO

log = logging.getLogger(__name__)


def get_msms_df(model, pdb_id, outfile=None, outdir=None, outext='_msms.df', force_rerun=False):
    """Run MSMS (using Biopython) on a Biopython Structure Model.

    Depths are in units Angstroms. 1A = 10^-10 m = 1nm. Returns a dictionary of::

        {
            chain_id:{
                        resnum1_id: (res_depth, ca_depth),
                        resnum2_id: (res_depth, ca_depth)
                     }
        }

    Args:
        model: Biopython Structure Model

    Returns:
        Pandas DataFrame: ResidueDepth property_dict, reformatted

    """
    # XTODO: need to deal with temporary surface/vertex files in tmp directory when running on a large scale --
    # XTODO: will run into inode limits! Also, some valuable information is in these MSMS output files that we should save.

    # Create the output file name
    outfile = ssbio.utils.outfile_maker(inname=pdb_id, outname=outfile, outdir=outdir, outext=outext)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        # Run MSMS with Biopython
        try:
            rd = PDB.ResidueDepth(model)
        except AssertionError:
            log.error('{}: unable to run MSMS'.format(pdb_id))
            return pd.DataFrame()

        # Reorganize the results into a csv file
        appender = []
        for k in rd.property_keys:
            x = rd.property_dict[k]
            chain = k[0]
            residue = k[1]
            het = residue[0]
            resnum = residue[1]
            icode = residue[2]
            resdepth = x[0]
            cadepth = x[1]
            appender.append((chain, resnum, icode, resdepth, cadepth))

        df = pd.DataFrame.from_records(appender, columns=['chain', 'resnum', 'icode', 'res_depth', 'ca_depth'])
        df.to_csv(outfile)
    else:
        log.debug('{}: already ran MSMS and force_rerun={}, loading results'.format(outfile, force_rerun))
        df = pd.read_csv(outfile, index_col=0)

    return df


def get_msms_df_on_file(pdb_file, outfile=None, outdir=None, outext='_msms.df', force_rerun=False):
    """Run MSMS (using Biopython) on a PDB file.

    Saves a CSV file of:
        chain: chain ID
        resnum: residue number (PDB numbering)
        icode: residue insertion code
        res_depth: average depth of all atoms in a residue
        ca_depth: depth of the alpha carbon atom

    Depths are in units Angstroms. 1A = 10^-10 m = 1nm

    Args:
        pdb_file: Path to PDB file
        outfile: Optional name of output file (without extension)
        outdir: Optional output directory
        outext: Optional extension for the output file
        outext: Suffix appended to json results file
        force_rerun: Rerun MSMS even if results exist already

    Returns:
        Pandas DataFrame: ResidueDepth property_dict, reformatted

    """
    # Create the output file name
    outfile = ssbio.utils.outfile_maker(inname=pdb_file, outname=outfile, outdir=outdir, outext=outext)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        # Load the structure
        my_structure = StructureIO(pdb_file)
        model = my_structure.first_model
        df = get_msms_df(model, pdb_id=op.splitext(op.basename(pdb_file))[0],
                         outfile=outfile, outdir=outdir, outext=outext, force_rerun=force_rerun)
    else:
        log.debug('{}: already ran MSMS and force_rerun={}, loading results'.format(outfile, force_rerun))
        df = pd.read_csv(outfile, index_col=0)

    return df


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Run MSMS to calculate residue depth on a file, files, or a directory. Save a JSON file of results per structure.')
    p.add_argument('infile', help='PDB file, files, or folder', type=str, nargs='+')
    p.add_argument('--summary', action='store_true', help='Save a summary DataFrame of results.')

    args = p.parse_args()
    infiles = ssbio.utils.input_list_parser(args.infile)

    msms_errors = []

    for f in tqdm(infiles):
        msms_stuff = get_msms_df_on_file(f)
        if msms_stuff.empty:
            msms_errors.append(f)

    if args.summary:
        pass
        # TODO: what to save as a summary? Average residue depth of entire protein?
        log.info('Saved DF at: {}'.format('{}-df_msms_summary.csv'.format(ssbio.utils.todays_short_date())))

    if msms_errors:
        log.warning('Errors with: {}'.format(msms_errors))