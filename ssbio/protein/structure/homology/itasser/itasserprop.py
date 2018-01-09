import logging
import os
import os.path as op
import shutil
import time

import pandas as pd
from ssbio.protein.structure.utils.cleanpdb import CleanPDB
from ssbio.protein.structure.utils.structureio import StructureIO

import ssbio.utils
from ssbio.protein.structure.structprop import StructProp

log = logging.getLogger(__name__)


class ITASSERProp(StructProp):

    """Parse all available information for a local I-TASSER modeling run.

    Initializes a class to collect I-TASSER modeling information and optionally copy results to a new directory.
    SEE: https://zhanglab.ccmb.med.umich.edu/papers/2015_1.pdf for detailed information.

    Args:
        ident (str): ID of I-TASSER modeling run
        original_results_path (str): Path to I-TASSER modeling folder
        coach_results_folder (str): Path to original COACH results
        model_to_use (str): Which I-TASSER model to use. Default is "model1"

    """

    # TODO: parse 1) predicted SS (seq.dat) 2) predicted solvent acc (exp.dat) 3) B-factor (BFP.dat?)
    # Other? Top 10 proteins with similar structure in PDB: $datadir/model1/cofactor/similarpdb_model1.lst

    def __init__(self, ident, original_results_path, coach_results_folder='model1/coach', model_to_use='model1'):

        if not op.exists(original_results_path):
            raise OSError('{}: folder does not exist'.format(original_results_path))

        StructProp.__init__(self, ident, is_experimental=False)

        self._attrs_to_copy = []

        self.original_results_path = original_results_path
        self.model_to_use = model_to_use

        ### MODEL
        self.model_date = None
        original_model_path = op.join(original_results_path, '{}.pdb'.format(model_to_use))
        if not op.exists(original_model_path):
            raise IOError('{}: no homology model generated'.format(original_model_path))
        else:
            self.load_structure_path(structure_path=original_model_path, file_type='pdb')

        ### MODELING RESULTS
        self.difficulty = None
        self.top_template_pdb = None
        self.top_template_chain = None

        # Parse init.dat
        self._init_dat_path = op.join(original_results_path, 'init.dat')
        if op.exists(self._init_dat_path):
            self.update(parse_init_dat(self._init_dat_path))

        # Parse seq.dat
        self._seq_dat_path = op.join(original_results_path, 'seq.dat')
        if op.exists(self._seq_dat_path):
            self._attrs_to_copy.append('_seq_dat_path')
            # TODO: parse seq.dat and store in modeling_results
            self.update(parse_seq_dat(self._seq_dat_path))

        # Parse exp.dat
        self._exp_dat_path = op.join(original_results_path, 'exp.dat')
        if op.exists(self._exp_dat_path):
            self._attrs_to_copy.append('_exp_dat_path')
            # TODO: parse exp.dat and store in modeling_results
            self.update(parse_exp_dat(self._exp_dat_path))

        # Parse BFP.dat
        self._bfp_dat_path = op.join(original_results_path, 'BFP.dat')
        if op.exists(self._bfp_dat_path):
            self._attrs_to_copy.append('_bfp_dat_path')
            # TODO: parse BFP.dat and store in modeling_results
            self.update(parse_bfp_dat(self._bfp_dat_path))

        # Parse cscore
        self.c_score = None
        self.tm_score = None
        self.tm_score_err = None
        self.rmsd = None
        self.rmsd_err = None
        self._cscore_path = op.join(original_results_path, 'cscore')
        if op.exists(self._cscore_path):
            self._attrs_to_copy.append('_cscore_path')
            self.update(parse_cscore(self._cscore_path))


        ### COACH RESULTS
        self.coach_bsites = []
        self.coach_ec = []
        self.coach_go_mf = []
        self.coach_go_bp = []
        self.coach_go_cc = []

        coach_folder = op.join(original_results_path, coach_results_folder)
        self._coach_bsites_inf_path = op.join(coach_folder, 'Bsites.inf')
        self._coach_ec_dat_path = op.join(coach_folder, 'EC.dat')
        self._coach_go_mf_dat_path = op.join(coach_folder, 'GO_MF.dat')
        self._coach_go_bp_dat_path = op.join(coach_folder, 'GO_BP.dat')
        self._coach_go_cc_dat_path = op.join(coach_folder, 'GO_CC.dat')

        if op.isdir(coach_folder):
            # Parse Bsites.inf
            if op.exists(self._coach_bsites_inf_path):
                parsed_bsites = parse_coach_bsites_inf(self._coach_bsites_inf_path)
                if parsed_bsites:
                    self._attrs_to_copy.append('_coach_bsites_inf_path')
                    self.coach_bsites = parsed_bsites

            # Parse EC.dat
            if op.exists(self._coach_ec_dat_path):
                parsed_ec = parse_coach_ec(self._coach_ec_dat_path)
                if parsed_ec:
                    self._attrs_to_copy.append('_coach_ec_dat_path')
                    self.coach_ec = parsed_ec

            # Parse GO_MF.dat
            if op.exists(self._coach_go_mf_dat_path):
                parsed_go_mf = parse_coach_go(self._coach_go_mf_dat_path)
                if parsed_go_mf:
                    self._attrs_to_copy.append('_coach_go_mf_dat_path')
                    self.coach_go_mf = parsed_go_mf

            # Parse GO_BP.dat
            if op.exists(self._coach_go_bp_dat_path):
                parsed_go_bp = parse_coach_go(self._coach_go_bp_dat_path)
                if parsed_go_bp:
                    self._attrs_to_copy.append('_coach_go_bp_dat_path')
                    self.coach_go_bp = parsed_go_bp

            # Parse GO_CC.dat
            if op.exists(self._coach_go_cc_dat_path):
                parsed_go_cc = parse_coach_go(self._coach_go_cc_dat_path)
                if parsed_go_cc:
                    self._attrs_to_copy.append('_coach_go_cc_dat_path')
                    self.coach_go_cc = parsed_go_cc

    def load_structure_path(self, structure_path, file_type='pdb'):
        StructProp.load_structure_path(self, structure_path=structure_path, file_type=file_type)
        # Additionally save model creation date
        self.model_date = time.strftime('%Y-%m-%d', time.gmtime(os.path.getmtime(self.structure_path)))

    def copy_results(self, copy_to_dir, rename_model_to=None, force_rerun=False):
        """Copy the raw information from I-TASSER modeling to a new folder.

        Copies all files in the list _attrs_to_copy.

        Args:
            copy_to_dir (str): Directory to copy the minimal set of results per sequence.
            rename_model_to (str): New file name (without extension)
            force_rerun (bool): If existing models and results should be overwritten.

        """
        # Save path to the structure and copy it if specified
        if not rename_model_to:
            rename_model_to = self.model_to_use

        new_model_path = op.join(copy_to_dir, '{}.pdb'.format(rename_model_to))

        if self.structure_path:
            if ssbio.utils.force_rerun(flag=force_rerun, outfile=new_model_path):
                # Clean and save it
                custom_clean = CleanPDB()
                my_pdb = StructureIO(self.structure_path)
                new_model_path = my_pdb.write_pdb(custom_selection=custom_clean,
                                                  custom_name=rename_model_to,
                                                  out_dir=copy_to_dir,
                                                  force_rerun=force_rerun)

            # Update the structure_path to be the new clean file
            self.load_structure_path(structure_path=new_model_path, file_type='pdb')

            # Other modeling results - store in a new folder
            dest_itasser_dir = op.join(copy_to_dir, '{}_itasser'.format(rename_model_to))
            if not op.exists(dest_itasser_dir):
                os.mkdir(dest_itasser_dir)

            for attr in self._attrs_to_copy:
                old_file_path = getattr(self, attr)
                new_file_path = op.join(dest_itasser_dir, op.basename(old_file_path))
                if ssbio.utils.force_rerun(flag=force_rerun, outfile=new_file_path):
                    shutil.copy2(old_file_path, new_file_path)
                    log.debug('{}: copied from {}'.format(new_file_path, old_file_path))
                else:
                    log.debug('{}: file already exists'.format(new_file_path))
                setattr(self, attr, new_file_path)

    @property
    def df_coach_bsites(self):
        df_cols = ['site_num', 'c_score', 'cluster_size', 'algorithm',
                   'pdb_template_id', 'pdb_template_chain', 'pdb_ligand',
                   'binding_location_coords', 'c_score_method', 'binding_residues',
                   'ligand_cluster_counts']

        bsites_inf_df = pd.DataFrame.from_records(self.coach_bsites, columns=df_cols).drop_duplicates().reset_index(drop=True)

        if bsites_inf_df.empty:
            log.warning('Empty dataframe')
            return bsites_inf_df
        else:
            bsites_inf_df['c_score'] = pd.to_numeric(bsites_inf_df.c_score, errors='coerce')
            bsites_inf_df['cluster_size'] = pd.to_numeric(bsites_inf_df.cluster_size, errors='coerce')
            return ssbio.utils.clean_df(bsites_inf_df)

    @property
    def df_coach_ec(self):
        if self._coach_ec_dat_path:
            return parse_coach_ec_df(self._coach_ec_dat_path)
        else:
            log.warning('Empty dataframe')
            return pd.DataFrame()

    @property
    def df_coach_go(self):
        cols = ['go_id', 'go_term', 'c_score']

        go_all_df = pd.DataFrame()

        for go_list in [self.coach_go_mf, self.coach_go_cc, self.coach_go_bp]:
            go_df = pd.DataFrame.from_records(go_list, columns=cols).drop_duplicates().reset_index(drop=True)
            go_df['c_score'] = pd.to_numeric(go_df.c_score, errors='coerce')

            if go_all_df.empty:
                go_all_df = go_df
            else:
                go_all_df.append(go_df)

        return go_all_df

    def get_dict(self, only_attributes=None, exclude_attributes=None, df_format=False):
        """Summarize the I-TASSER run in a dictionary containing modeling results and top predictions from COACH

        Args:
            only_attributes (str, list): Attributes that should be returned. If not provided, all are returned.
            exclude_attributes (str, list): Attributes that should be excluded.
            df_format (bool): If dictionary values should be formatted for a dataframe
                (everything possible is transformed into strings, int, or float -
                if something can't be transformed it is excluded)

        Returns:
            dict: Dictionary of attributes

        """

        to_exclude = ['coach_bsites', 'coach_ec', 'coach_go_mf', 'coach_go_bp', 'coach_go_cc']
        if not exclude_attributes:
            excluder = to_exclude
        else:
            excluder = ssbio.utils.force_list(exclude_attributes)
            excluder.extend(to_exclude)

        summary_dict = StructProp.get_dict(self, only_attributes=only_attributes,
                                           exclude_attributes=excluder,
                                           df_format=df_format)

        if self.coach_bsites:
            tmp = {'top_bsite_' + k:v for k, v in self.coach_bsites[0].items()}
            summary_dict.update(tmp)

        if self.coach_ec:
            tmp = {'top_ec_' + k: v for k, v in self.coach_ec[0].items()}
            summary_dict.update(tmp)

        if self.coach_go_mf:
            tmp = {'top_go_mf_' + k: v for k, v in self.coach_go_mf[0].items()}
            summary_dict.update(tmp)

        if self.coach_go_bp:
            tmp = {'top_go_bp_' + k: v for k, v in self.coach_go_bp[0].items()}
            summary_dict.update(tmp)

        if self.coach_go_cc:
            tmp = {'top_go_cc_' + k: v for k, v in self.coach_go_cc[0].items()}
            summary_dict.update(tmp)

        return summary_dict


def parse_init_dat(infile):
    """Parse the main init.dat file which contains the modeling results

    The first line of the file init.dat contains stuff like::
    
        "120 easy  40   8"
    
    The other lines look like this::
    
        "     161   11.051   1  1guqA MUSTER"
    
    and getting the first 10 gives you the top 10 templates used in modeling

    Args:
        infile (stt): Path to init.dat

    Returns:
        dict: Dictionary of parsed information

    """

    # TODO: would be nice to get top 10 templates instead of just the top
    init_dict = {}

    log.debug('{}: reading file...'.format(infile))
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
        infile (str): Path to seq.dat

    Returns:
        list: List of secondary structure predictions for all residues

    """

    # TODO: parser for seq.dat
    seq_dict = {}
    return seq_dict


def parse_exp_dat(infile):
    """Parse the solvent accessibility predictions in exp.dat

    Args:
        infile (str): Path to exp.dat

    Returns:
        list: List of solvent accessibility predictions for all residues

    """

    # TODO: parser for exp.dat
    exp_dict = {}
    return exp_dict


def parse_bfp_dat(infile):
    """Parse the B-factor predictions in BFP.dat

    Args:
        infile (str): Path to BFP.dat

    Returns:
        list: List of B-factor predictions for all residues

    """

    # TODO: parser for BFP.dat
    bfp_dict = {}
    return bfp_dict


def parse_cscore(infile):
    """Parse the cscore file to return a dictionary of scores.

    Args:
        infile (str): Path to cscore

    Returns:
        dict: Dictionary of scores

    """

    cscore_dict = {}

    with open(infile, 'r') as f:
        for ll in f.readlines():
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


def parse_coach_bsites_inf(infile):
    """Parse the Bsites.inf output file of COACH and return a list of rank-ordered binding site predictions

    Bsites.inf contains the summary of COACH clustering results after all other prediction algorithms have finished
    For each site (cluster), there are three lines:
    
        - Line 1: site number, c-score of coach prediction, cluster size
        - Line 2: algorithm, PDB ID, ligand ID, center of binding site (cartesian coordinates), 
          c-score of the algorithm's prediction, binding residues from single template
        - Line 3: Statistics of ligands in the cluster

    C-score information:

        - "In our training data, a prediction with C-score>0.35 has average false positive and false negative rates below
          0.16 and 0.13, respectively." (https://zhanglab.ccmb.med.umich.edu/COACH/COACH.pdf)

    Args:
        infile (str): Path to Bsites.inf

    Returns:
        list: Ranked list of dictionaries, keys defined below
        
            - ``site_num``: cluster which is the consensus binding site
            - ``c_score``: confidence score of the cluster prediction
            - ``cluster_size``: number of predictions within this cluster
            - ``algorithm``: main? algorithm used to make the prediction
            - ``pdb_template_id``: PDB ID of the template used to make the prediction
            - ``pdb_template_chain``: chain of the PDB which has the ligand
            - ``pdb_ligand``: predicted ligand to bind
            - ``binding_location_coords``: centroid of the predicted ligand position in the homology model
            - ``c_score_method``: confidence score for the main algorithm
            - ``binding_residues``: predicted residues to bind the ligand
            - ``ligand_cluster_counts``: number of predictions per ligand

    """

    bsites_results = []

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

        bsites_results.append(bsites_site_dict)

    return bsites_results


def parse_coach_ec_df(infile):
    """Parse the EC.dat output file of COACH and return a dataframe of results

    EC.dat contains the predicted EC number and active residues.
    The columns are: PDB_ID, TM-score, RMSD, Sequence identity,
    Coverage, Confidence score, EC number, and Active site residues

    Args:
        infile (str): Path to EC.dat

    Returns:
        DataFrame: Pandas DataFrame summarizing EC number predictions

    """

    ec_df = pd.read_table(infile, delim_whitespace=True,
                          names=['pdb_template', 'tm_score', 'rmsd', 'seq_ident', 'seq_coverage',
                                 'c_score', 'ec_number', 'binding_residues'])

    ec_df['pdb_template_id'] = ec_df['pdb_template'].apply(lambda x: x[:4])
    ec_df['pdb_template_chain'] = ec_df['pdb_template'].apply(lambda x: x[4])

    ec_df = ec_df[['pdb_template_id', 'pdb_template_chain', 'tm_score', 'rmsd',
                   'seq_ident', 'seq_coverage', 'c_score', 'ec_number', 'binding_residues']]
    ec_df['c_score'] = pd.to_numeric(ec_df.c_score, errors='coerce')

    return ec_df


def parse_coach_ec(infile):
    """Parse the EC.dat output file of COACH and return a list of rank-ordered EC number predictions

    EC.dat contains the predicted EC number and active residues.
    The columns are: PDB_ID, TM-score, RMSD, Sequence identity,
    Coverage, Confidence score, EC number, and Active site residues

    Args:
        infile (str): Path to EC.dat

    Returns:
         list: Ranked list of dictionaries, keys defined below

            - ``pdb_template_id``: PDB ID of the template used to make the prediction
            - ``pdb_template_chain``: chain of the PDB which has the ligand
            - ``tm_score``: TM-score of the template to the model (similarity score)
            - ``rmsd``: RMSD of the template to the model (also a measure of similarity)
            - ``seq_ident``: percent sequence identity
            - ``seq_coverage``: percent sequence coverage
            - ``c_score``: confidence score of the EC prediction
            - ``ec_number``: predicted EC number
            - ``binding_residues``: predicted residues to bind the ligand

    """
    return list(parse_coach_ec_df(infile).to_dict(orient='index').values())


def parse_coach_go(infile):
    """Parse a GO output file from COACH and return a rank-ordered list of GO term predictions

    The columns in all files are: GO terms, Confidence score, Name of GO terms. The files are:
        
        - GO_MF.dat - GO terms in 'molecular function'
        - GO_BP.dat - GO terms in 'biological process'
        - GO_CC.dat - GO terms in 'cellular component'

    Args:
        infile (str): Path to any COACH GO prediction file

    Returns:
        Pandas DataFrame: Organized dataframe of results, columns defined below
            
            - ``go_id``: GO term ID
            - ``go_term``: GO term text
            - ``c_score``: confidence score of the GO prediction

    """
    go_list = []

    with open(infile) as go_file:
        for line in go_file.readlines():
            go_dict = {}

            go_split = line.split()
            go_dict['go_id'] = go_split[0]
            go_dict['c_score'] = go_split[1]
            go_dict['go_term'] = ' '.join(go_split[2:])

            go_list.append(go_dict)

    return go_list