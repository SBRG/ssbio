import os.path as op
import pandas as pd
import shutil
import os
import logging
log = logging.getLogger(__name__)
import time

from ssbio.structure.pdbioext import PDBIOExt
from ssbio.structure.cleanpdb import CleanPDB
import ssbio.utils

class ITASSERParse():
    """Parse all available information for an I-TASSER modeling run.
    """

    _main_files_to_copy = ['seq.dat', 'cscore']
    _coach_files_to_copy = ['Bsites.inf', 'EC.dat', 'GO_MF.dat', 'GO_BP.dat', 'GO_CC.dat']

    def __init__(self, original_results_path, create_dfs=False,
                 coach_results_folder='model1/coach', model_to_use='model1'):
        """Initialize a class to collect I-TASSER modeling information and optionally copy results to a new directory.

        Args:
            original_results_path: Path to I-TASSER modeling folder
            create_dfs: If data frames should be created for COACH results
            coach_results_folder: Path to original COACH results
            model_to_use: Which I-TASSER model to use. Default is "model1"
        """
        self.results_path = original_results_path
        self.model_to_use = model_to_use
        self.create_dfs = create_dfs

        # For summarizing main modeling results
        self.modeling_results = {}

        ### MODEL1.PDB THINGS
        old_model_path = op.join(original_results_path, '{}.pdb'.format(model_to_use))
        if not op.exists(old_model_path):
            log.debug('{}: No homology model available'.format(old_model_path))
            self.structure_path = None
        else:
            self.structure_path = old_model_path
            # Save model creation date
            self.modeling_results['model_date'] = time.strftime('%Y-%m-%d', time.gmtime(os.path.getmtime(old_model_path)))

        ### MODELING RESULTS
        # Parse init.dat
        init_dat_path = op.join(original_results_path, 'init.dat')
        if op.exists(init_dat_path):
            self.modeling_results.update(self.parse_init_dat(init_dat_path))

        # Parse seq.dat
        seq_dat_path = op.join(original_results_path, 'seq.dat')
        if op.exists(seq_dat_path):
            pass
            # TODO: parse seq.dat and store in modeling_results

        # Parse cscore
        cscore_path = op.join(original_results_path, 'cscore')
        if op.exists(cscore_path):
            self.modeling_results.update(self.parse_cscore(cscore_path))

        # Parse COACH results if available
        coach_folder = op.join(original_results_path, coach_results_folder)
        if op.isdir(coach_folder):

            # Parse Bsites.inf
            bsites_inf_path = op.join(coach_folder, 'Bsites.inf')
            if op.exists(bsites_inf_path):
                parsed_bsites = self.parse_bsites_inf(infile=bsites_inf_path, create_df=create_dfs)
                if parsed_bsites:
                    tmp = {'top_bsite_'+k: v for k, v in parsed_bsites[0].items()}
                    self.modeling_results.update(tmp)

            # Parse EC.dat
            ec_dat_path = op.join(coach_folder, 'EC.dat')
            if op.exists(ec_dat_path):
                parsed_ec = self.parse_ec(infile=ec_dat_path, create_df=create_dfs)
                if parsed_ec:
                    tmp = {'top_ec_' + k: v for k, v in parsed_ec[0].items()}
                    self.modeling_results.update(tmp)

            # Parse GO_MF.dat
            go_mf_dat_path = op.join(coach_folder, 'GO_MF.dat')
            if op.exists(go_mf_dat_path):
                parsed_go_mf = self.parse_go(infile=go_mf_dat_path, create_df=create_dfs)
                if parsed_go_mf:
                    tmp = {'top_go_mf_' + k: v for k, v in parsed_go_mf[0].items()}
                    self.modeling_results.update(tmp)

            # Parse GO_BP.dat
            go_bp_dat_path = op.join(coach_folder, 'GO_BP.dat')
            if op.exists(go_bp_dat_path):
                parsed_go_bp = self.parse_go(infile=go_bp_dat_path, create_df=create_dfs)
                if parsed_go_bp:
                    tmp = {'top_go_bp_' + k: v for k, v in parsed_go_bp[0].items()}
                    self.modeling_results.update(tmp)

            # Parse GO_CC.dat
            go_cc_dat_path = op.join(coach_folder, 'GO_CC.dat')
            if op.exists(go_cc_dat_path):
                parsed_go_cc = self.parse_go(infile=go_cc_dat_path, create_df=create_dfs)
                if parsed_go_cc:
                    tmp = {'top_go_cc_' + k: v for k, v in parsed_go_cc[0].items()}
                    self.modeling_results.update(tmp)

    def copy_results(self, copy_to_dir, rename_model_to=None, force_rerun=False):
        """Copy the raw information from I-TASSER modeling to a new folder.

        Copies files in the variables _main_files_to_copy and _coach_files_to_copy.

        Args:
            copy_to_dir: Directory to copy the minimal set of results per sequence.
            rename_model_to: New file name (without extension)
            force_rerun: If existing models and results should be overwritten.

        """
        # Save path to the structure and copy it if specified
        if not rename_model_to:
            rename_model_to = self.model_to_use

        new_model_path = op.join(copy_to_dir, '{}.pdb'.format(rename_model_to))

        if self.structure_path:
            if ssbio.utils.force_rerun(flag=force_rerun, outfile=new_model_path):
                # Clean and save it
                custom_clean = CleanPDB()
                my_pdb = PDBIOExt(self.structure_path, file_type='pdb')
                new_model_path = my_pdb.write_pdb(custom_selection=custom_clean,
                                                  custom_name=rename_model_to,
                                                  out_dir=copy_to_dir)
            # Update the structure_path to be the copied, clean file
            self.structure_path = new_model_path

            # Other modeling results - store in a new folder
            dest_itasser_dir = op.join(copy_to_dir, '{}_itasser'.format(rename_model_to))
            if not op.exists(dest_itasser_dir):
                os.mkdir(dest_itasser_dir)

            for f in self._main_files_to_copy:
                old_file_path = op.join(self.results_path, f)
                new_file_path = op.join(dest_itasser_dir, f)
                if op.exists(old_file_path):
                    if ssbio.utils.force_rerun(flag=force_rerun, outfile=new_file_path):
                        shutil.copy2(old_file_path, new_file_path)
                        log.debug('{}: copied from {}'.format(new_file_path, old_file_path))
                    else:
                        log.debug('{}: file already exists'.format(new_file_path))

            for ff in self._coach_files_to_copy:
                old_file_path = op.join(self.results_path, 'model1/coach/', ff)
                new_file_path = op.join(dest_itasser_dir, ff)
                if op.exists(old_file_path):
                    if ssbio.utils.force_rerun(flag=force_rerun, outfile=new_file_path):
                        shutil.copy2(old_file_path, new_file_path)
                        log.debug('{}: copied from {}'.format(new_file_path, old_file_path))
                    else:
                        log.debug('{}: file already exists'.format(new_file_path))

    def save_dataframes(self, outdir):
        """Save all attributes that start with "df_" into a specified directory.

        Args:
            outdir: Path to output directory.

        Returns:

        """
        if self.create_dfs:
            # Get list of attributes that start with "df_"
            dfs = list(filter(lambda x: x.startswith('df_'), dir(self)))
            for df in dfs:
                outpath = ssbio.utils.outfile_maker(inname=df, outext='.csv', outdir=outdir)
                my_df = getattr(self, df)
                my_df.to_csv(outpath)
                log.debug('{}: saved dataframe'.format(outpath))
        else:
            raise NameError('Nothing to save! Create dataframes by specifying create_dfs=True.')

        log.debug('Saved {} dataframes at {}'.format(len(dfs), outdir))

    def parse_init_dat(self, infile):
        """Parse the main init.dat file which contains the modeling results

        The first line of the file init.dat contains stuff like:
        "120 easy  40   8"
        The other lines look like this:
        "     161   11.051   1  1guqA MUSTER"
        and getting the first 10 gives you the top 10 templates used in modeling

        Args:
            infile: Path to init.dat
            copy_to_folder: Optional folder to copy init.dat to

        Returns:
            dict: Dictionary of parsed information

        """
        init_dict = {}

        log.debug('{}: Reading file...')
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

    def parse_seq_dat(self, infile):
        """Parse the secondary structure predictions in seq.dat

        Args:
            infile: Path to seq.dat

        Returns:
            dict: Dictionary of residue number to secondary structure

        """
        pass

    def parse_cscore(self, infile):
        """Parse the cscore file to return a dictionary of scores.

        Args:
            infile: Path to cscore

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

    def parse_bsites_inf(self, infile, create_df=False):
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
            create_df (bool): If a Pandas DataFrame should be saved to the object under the attribute "df_bsites"

        Returns:
            list: Ranked list of dictionaries, keys defined below
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

        pre_bsites_inf_df = []

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

            pre_bsites_inf_df.append(bsites_site_dict)

        if create_df:
            cols = ['site_num', 'c_score', 'cluster_size', 'algorithm',
                    'pdb_template_id', 'pdb_template_chain', 'pdb_ligand',
                    'binding_location_coords', 'c_score_method', 'binding_residues',
                    'ligand_cluster_counts']
            bsites_inf_df = pd.DataFrame.from_records(pre_bsites_inf_df, columns=cols).drop_duplicates().reset_index(drop=True)
            bsites_inf_df['c_score'] = pd.to_numeric(bsites_inf_df.c_score, errors='coerce')
            bsites_inf_df['cluster_size'] = pd.to_numeric(bsites_inf_df.cluster_size, errors='coerce')

            self.df_bsites = bsites_inf_df
            log.debug('Created binding site prediction DataFrame in attribute "df_bsites')

        return pre_bsites_inf_df


    def parse_ec(self, infile, create_df=False):
        """Parse the EC.dat output file of COACH and returns a dataframe

        EC.dat contains the predicted EC number and active residues.
            The columns are: PDB_ID, TM-score, RMSD, Sequence identity,
            Coverage, Confidence score, EC number, and Active site residues

        Args:
            infile (str): Path to EC.dat

        Returns:
             list: Ranked list of dictionaries, keys defined below
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

        if create_df:
            self.df_ec = ec_df
            log.debug('Created EC prediction DataFrame in attribute "df_ec')

        return list(ec_df.to_dict(orient='index').values())


    def parse_go(self, infile, create_df=False):
        """Parse a GO output file from COACH and returns a dataframe

        The columns in all files are: GO terms, Confidence score, Name of GO terms.
            GO_MF.dat - GO terms in 'molecular function'
            GO_BP.dat - GO terms in 'biological process'
            GO_CC.dat - GO terms in 'cellular component'

        Args:
            infile (str): Path to any GO file

        Returns:
            Pandas DataFrame: Organized dataframe of results, columns defined below
                go_id: GO term ID
                go_term: GO term text
                c_score: confidence score of the GO prediction

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

        if create_df:
            cols = ['go_id', 'go_term', 'c_score']
            go_df = pd.DataFrame.from_records(go_list, columns=cols).drop_duplicates().reset_index(drop=True)
            go_df['c_score'] = pd.to_numeric(go_df.c_score, errors='coerce')

            if hasattr(self, 'df_go'):
                self.df_go.append(go_df)
            else:
                self.df_go = go_df
            log.debug('Created GO prediction DataFrame in attribute "df_go"')

        return go_list

