from Bio import PDB
from ssbio.structure.pdbioext import PDBIOExt
import os
import os.path as op
import pandas as pd


def msms_output(pdb_file, outdir=None):
    """Run MSMS (through Biopython) on a PDB file.

    Saves the result dataframe

    Args:
        pdb_file:

    Returns:

    """
    basename = op.splitext(op.basename(pdb_file))[0]

    my_structure = PDBIOExt(pdb_file)
    model = my_structure.first_model
    rd = PDB.ResidueDepth(model, pdb_file)

    akeys = list(rd)
    if len(akeys) == 0:
        akeys = [('<Residue UNKNOWN het=  resseq=0 icode= >', (0, 0))]

    anslist = []
    if akeys == [('NA', ('NA', 'NA'))]:
        anslist = ["NA", "NA", "NA", "NA"]
        # print "warning: error at index:"+str(len(msmslist))

    else:
        for j in akeys:
            chain = [x.id for x in PDB.Selection.unfold_entities(j[0], 'C')][0]
            seqnum = j[0].id[1]
            re_depth = j[1][0]
            ca_depth = j[1][1]
            if chain == ' ':
                chain = 'X'
            anslist.append([chain, seqnum, re_depth, ca_depth])

    msms_df = pd.DataFrame(anslist)
    msms_df = msms_df.rename(columns={0:'chain',1:'resnum',2:'res_depth',3:'ca_depth'})
    if outdir:
        outfile = op.join(outdir, basename + '_msms.df')
    else:
        outfile = basename + '_msms.df'

    msms_df.to_csv(outfile)

    return outfile

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