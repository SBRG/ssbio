import os.path as op
import pandas as pd
import shutil
import os
import logging
log = logging.getLogger(__name__)
import time

def organize_itasser_models(raw_dir, copy_to_dir, rename_model_to=''):
    """Reorganize the raw information from I-TASSER modeling.

    - Copies the model1.pdb file and COACH results (if they exist) to specified folder
    - Parses modeling information such as C-scores, template, etc
    - Also will parse binding site (COACH) information if it exists

    Returns:
        dict: summary of I-TASSER results for a sequence
    """

    hom_dict = {}

    # Best homology model is "model1.pdb"
    model = 'model1.pdb'
    old_model_path = op.join(raw_dir, model)
    new_model_path = op.join(copy_to_dir, model)
    if rename_model_to:
        new_model_path = op.join(copy_to_dir, '{}.pdb'.format(rename_model_to))
    if op.exists(old_model_path):
        if not op.exists(new_model_path):
            shutil.copy2(old_model_path, new_model_path)
        hom_dict['model_file'] = op.basename(new_model_path)
        hom_dict['model_date'] = time.strftime('%Y-%m-%d', time.gmtime(os.path.getmtime(old_model_path)))
    else:
        log.debug('{}: No homology model exists'.format(raw_dir))
        return hom_dict

    ### Other modeling results
    dest_itasser_dir = op.join(copy_to_dir, 'itasser')
    if rename_model_to:
        dest_itasser_dir = op.join(copy_to_dir, '{}_itasser'.format(rename_model_to))
    if not op.exists(dest_itasser_dir) and op.exists(old_model_path):
        os.mkdir(dest_itasser_dir)

    # init.dat
    # File size is usually pretty big, don't copy for now, just parse
    init_dat = 'init.dat'
    old_init_dat = op.join(raw_dir, init_dat)
    new_init_dat = op.join(dest_itasser_dir, init_dat)
    if op.exists(old_init_dat):
        hom_dict.update(parse_init_dat(old_init_dat))
        # if not op.exists(new_init_dat):
        #     shutil.copy2(old_init_dat, new_init_dat)

    # seq.dat
    seq_dat = 'seq.dat'
    old_seq_dat = op.join(raw_dir, seq_dat)
    new_seq_dat = op.join(dest_itasser_dir, seq_dat)
    if op.exists(old_seq_dat):
        if not op.exists(new_seq_dat):
            shutil.copy2(old_seq_dat, new_seq_dat)
            # TODO: parse seq.dat

    # cscore
    cscore = 'cscore'
    old_cscore = op.join(raw_dir, cscore)
    new_cscore = op.join(dest_itasser_dir, cscore)
    if op.exists(old_cscore):
        hom_dict.update(parse_cscore(old_cscore))
        if not op.exists(new_cscore):
            shutil.copy2(old_cscore, new_cscore)

    ### COACH (binding site prediction) results
    model1_coach_folder = op.join(raw_dir, 'model1', 'coach')
    dest_coach_dir = op.join(dest_itasser_dir, 'coach')
    if rename_model_to:
        dest_coach_dir = op.join(dest_itasser_dir, '{}_coach'.format(rename_model_to))
    if not op.exists(dest_coach_dir) and op.exists(model1_coach_folder):
        os.mkdir(dest_coach_dir)

    # Bsites.inf
    bsites_inf = 'Bsites.inf'
    old_bsites_inf = op.join(model1_coach_folder, bsites_inf)
    new_bsites_inf = op.join(dest_coach_dir, bsites_inf)
    if op.exists(old_bsites_inf):
        if not op.exists(new_bsites_inf):
            shutil.copy2(old_bsites_inf, new_bsites_inf)
        # TODO: parse Bsites.inf

    # EC.dat
    ec = 'EC.dat'
    old_ec = op.join(model1_coach_folder, ec)
    new_ec = op.join(dest_coach_dir, ec)
    if op.exists(old_ec):
        if not op.exists(new_ec):
            shutil.copy2(old_ec, new_ec)
        # TODO: parse EC

    go_mf = 'GO_MF.dat'
    old_go_mf = op.join(model1_coach_folder, go_mf)
    new_go_mf = op.join(dest_coach_dir, go_mf)
    if op.exists(old_go_mf):
        if not op.exists(new_go_mf):
            shutil.copy2(old_go_mf, new_go_mf)
        # TODO: parse GO_MF

    go_bp = 'GO_BP.dat'
    old_go_bp = op.join(model1_coach_folder, go_bp)
    new_go_bp = op.join(dest_coach_dir, go_bp)
    if op.exists(old_go_bp):
        if not op.exists(new_go_bp):
            shutil.copy2(old_go_bp, new_go_bp)
        # TODO: parse GO_BP

    go_cc = 'GO_CC.dat'
    old_go_cc = op.join(model1_coach_folder, go_cc)
    new_go_cc = op.join(dest_coach_dir, go_cc)
    if op.exists(old_go_cc):
        if not op.exists(new_go_cc):
            shutil.copy2(old_go_cc, new_go_cc)
        # TODO: parse GO_CC

    return hom_dict

def parse_init_dat(infile):
    """Parse the main init.dat file which contains the modeling results

    The first line of the file init.dat contains stuff like:
    "120 easy  40   8"
    # TODO: what is this?

    The other lines look like this:
    "     161   11.051   1  1guqA MUSTER"
    and getting the first 10 gives you the top 10 templates used in modeling

    Args:
        infile:

    Returns:

    """
    init_dict = {}

    with open(infile, 'r') as f:
        # Get first 2 lines of file
        head = [next(f).strip() for x in range(2)]

    summary = head[0].split()
    difficulty = summary[1]

    top_template_info = head[1].split()
    top_template_pdbchain = top_template_info[3]
    top_template_pdb = top_template_pdbchain[:4]
    top_template_chain = top_template_pdbchain[4:]

    init_dict['difficulty'] = difficulty
    init_dict['top_template_pdb'] = top_template_pdb
    init_dict['top_template_chain'] = top_template_chain

    return init_dict

def parse_seq_dat(infile):
    """Parse the secondary structure predictions in seq.dat

    Args:
        infile: Path to seq.dat

    Returns:
        dict: Dictionary of residue number to secondary structure

    """
    pass

def parse_cscore(infile):
    """Parse the cscore file to return a dictionary of scores.

    Args:
        infile: Path to cscore

    Returns:
        dict: Dictionary of scores

    """
    cscore_dict = {}
    with open(infile, 'r') as f:
        lines = f.readlines()
        for ll in lines:
            # Look for the first line that starts with model1
            if ll.lower().startswith('model1'):
                l = ll.split()

                cscore = l[1]
                tmscore_full = l[2].split('+-')
                tmscore = tmscore_full[0]
                tmscore_err = tmscore_full[1]
                rmsd_full = l[3].split('+-')
                rmsd = rmsd_full[0]
                rmsd_err = rmsd_full[1]

                cscore_dict['c_score'] = float(cscore)
                cscore_dict['tm_score'] = float(tmscore)
                cscore_dict['tm_score_err'] = float(tmscore_err)
                cscore_dict['rmsd'] = float(rmsd)
                cscore_dict['rmsd_err'] = float(rmsd_err)

    return cscore_dict

def parse_bsites_inf(infile):
    """Parse the Bsites.inf output file of COACH and returns a dataframe

    Bsites.inf contains the summary of COACH clustering results after all
        other prediction algorithms have finished
    For each site (cluster), there are three lines:
    Line 1: site number, c-score of coach prediction, cluster size
    Line 2: algorithm, PDB ID, ligand ID, center of binding site (cartesian coordinates),
        c-score of the algorithm's prediction, binding residues from single template
    Line 3: Statistics of ligands in the cluster

    Args:
        infile (str): Path to Bsites.inf

    Returns:
        bsites_inf_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            site_num: cluster which is the consensus binding site
            c_score: confidence score of the cluster prediction
            cluster_size: number of predictions within this cluster
            algorithm: main? algorithm used to make the prediction
            pdb_template_id: PDB ID of the template used to make the prediction
            pdb_template_chain: chain of the PDB which has the ligand
            pdb_ligand: predicted ligand to bind
            binding_location_coords: centroid of the predicted ligand position in the homology model
            c_score_method: confidence score for the main algorithm
            binding_residues: predicted residues to bind the ligand
            ligand_cluster_counts: number of predictions per ligand

    """

    bsites_inf_df = pd.DataFrame()

    with open(infile) as pp:
        lines = list(filter(None, (line.rstrip() for line in pp)))

        for i in range(len(lines) // 3):
            bsites_site_dict = {}

            line1 = lines[i * 3].split('\t')
            line2 = lines[i * 3 + 1].split('\t')
            line3 = lines[i * 3 + 2]

            bsites_site_dict['site_num'] = line1[0]
            bsites_site_dict['c_score'] = float(line1[1])
            bsites_site_dict['cluster_size'] = line1[2]

            bsites_site_dict['algorithm'] = line2[0]
            bsites_site_dict['pdb_template_id'] = line2[1][:4]
            bsites_site_dict['pdb_template_chain'] = line2[1][4]
            bsites_site_dict['pdb_ligand'] = line2[2]
            bsites_site_dict['binding_location_coords'] = tuple(float(x) for x in line2[3].split())

            # TODO: what's the difference between this c-score and the cluster's c-score?
            # how is the cluster's c-score computed? it's not the average c-score of all methods
            # also why are some COFACTOR c-scores >1?
            # 160411 - seems like the COFACTOR "BS-score" is being reported here, not its c-score...
            tmp_split = line2[4].split(' :')
            bsites_site_dict['c_score_method'] = tmp_split[0]
            bsites_site_dict['binding_residues'] = tmp_split[1]

            bsites_site_dict['ligand_cluster_counts'] = line3

            bsites_inf_df = bsites_inf_df.append(bsites_site_dict, ignore_index=True)

    if len(bsites_inf_df) == 0:
        return bsites_inf_df

    bsites_inf_df = bsites_inf_df[['site_num', 'c_score', 'cluster_size', 'algorithm',
                                   'pdb_template_id', 'pdb_template_chain', 'pdb_ligand',
                                   'binding_location_coords', 'c_score_method', 'binding_residues',
                                   'ligand_cluster_counts']]
    return bsites_inf_df  # .set_index('site_num')


def parse_ec(infile):
    """Parse the EC.dat output file of COACH and returns a dataframe

    EC.dat contains the predicted EC number and active residues.
        The columns are: PDB_ID, TM-score, RMSD, Sequence identity,
        Coverage, Confidence score, EC number, and Active site residues

    Args:
        infile (str): Path to EC.dat

    Returns:
        ec_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            pdb_template_id: PDB ID of the template used to make the prediction
            pdb_template_chain: chain of the PDB which has the ligand
            tm_score: TM-score of the template to the model (similarity score)
            rmsd: RMSD of the template to the model (also a measure of similarity)
            seq_ident: percent sequence identity
            seq_coverage: percent sequence coverage
            c_score: confidence score of the EC prediction
            ec_number: predicted EC number
            binding_residues: predicted residues to bind the ligand

    """

    ec_df = pd.read_table(infile, delim_whitespace=True,
                          names=['pdb_template', 'tm_score', 'rmsd', 'seq_ident', 'seq_coverage',
                                 'c_score', 'ec_number', 'binding_residues'])
    ec_df['pdb_template_id'] = ec_df['pdb_template'].apply(lambda x: x[:4])
    ec_df['pdb_template_chain'] = ec_df['pdb_template'].apply(lambda x: x[4])

    ec_df = ec_df[['pdb_template_id', 'pdb_template_chain', 'tm_score', 'rmsd',
                   'seq_ident', 'seq_coverage', 'c_score', 'ec_number', 'binding_residues']]
    ec_df['c_score'] = pd.to_numeric(ec_df.c_score, errors='coerce')
    return ec_df  # .set_index('pdb_template_id')


def parse_go(infile):
    """Parse a GO output file from COACH and returns a dataframe

    The columns in all files are: GO terms, Confidence score, Name of GO terms.
        GO_MF.dat - GO terms in 'molecular function'
        GO_BP.dat - GO terms in 'biological process'
        GO_CC.dat - GO terms in 'cellular component'

    Args:
        infile (str): Path to any GO file

    Returns:
        go_df (Pandas DataFrame): Organized dataframe of results, columns defined below
            go_id: GO term ID
            go_term: GO term text
            c_score: confidence score of the GO prediction

    """

    go_df = pd.DataFrame()

    with open(infile) as go_file:
        for line in go_file.readlines():
            go_dict = {}

            go_split = line.split()
            go_dict['go_id'] = go_split[0]
            go_dict['c_score'] = go_split[1]
            go_dict['go_term'] = ' '.join(go_split[2:])

            go_df = go_df.append(go_dict, ignore_index=True)

    if len(go_df) == 0:
        return go_df

    go_df = go_df[['go_id', 'go_term', 'c_score']]
    go_df['c_score'] = pd.to_numeric(go_df.c_score, errors='coerce')
    return go_df  # .set_index('go_id')

class ITASSERParse():
    '''
    Parse I-TASSER & COACH results
    '''

    def __init__(self):
        pass
