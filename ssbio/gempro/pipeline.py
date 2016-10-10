import argparse
import os
import os.path as op
import shutil
import warnings
import glob

from cobra.io import load_matlab_model
from cobra.io.sbml import create_cobra_model_from_sbml_file

from Bio import SeqIO

from bioservices.uniprot import UniProt

bsup = UniProt()
from bidict import bidict

from ssbio import utils
from ssbio.cobra.utils import is_spontaneous, true_num_reactions, true_num_genes

date = utils.Date()

import ssbio.databases.kegg
import ssbio.databases.uniprot
import ssbio.databases.pdb
import pandas as pd

from tqdm import tqdm
import requests

from ssbio.databases.pdb import top_pdb_blast_hit

from bioservices import KEGG

bs_kegg = KEGG()


class GEMPRO():
    """Generic class to represent all information of a GEM-PRO for a GEM

    Main steps are:
    1. Preparation of folders
    2. Automatic ID mapping
    3. QC/QA
    4. Homology modeling
    5. Representative structures

    Each step may generate a report and also request additional files if mappings are incomplete
    """

    def __init__(self, gem_name, root_dir):
        self.root_dir = root_dir
        self.model_dir = op.join(root_dir, gem_name)
        # model_files - directory where original gems and gem-related files are stored
        self.model_files = op.join(self.model_dir, 'model_files')

    def prep_folders(self):
        """Prepare all folders for a GEM-PRO project
        """
        # data - directory where all data will be stored
        self.data = op.join(self.model_dir, 'data')

        # notebooks - directory where ipython notebooks will be stored for manual analyses
        self.notebooks = op.join(self.model_dir, 'notebooks')

        # # figures - directory where all figures will be stored
        self.figures = op.join(self.model_dir, 'figures')

        # missing - directory where missing reports go
        self.missing = op.join(self.model_dir, 'missing')

        # struct_files - directory where structure related files will be downloaded/are located
        self.struct_files = op.join(self.model_dir, 'structure_files')
        self.struct_single_chain = op.join(self.struct_files, 'by_gene')

        # seq_files - sequence related files are stored here
        self.seq_files = op.join(self.model_dir, 'sequence_files')

        for directory in [self.model_dir, self.data, self.notebooks, self.figures, self.missing,
                          self.model_files, self.struct_files, self.struct_single_chain, self.seq_files]:
            if not op.exists(directory):
                os.makedirs(directory)
                print('[PREP] Created directory: {}'.format(directory))
            else:
                print('[PREP] Directory already exists: {}'.format(directory))

    def load_model(self, gem_file, file_type):
        """Load the GEM using COBRApy. Accept SBML or MAT files as input

        Args:
            file_type: if your model is in "SBML" (or "XML"), or "Matlab" formats

        Returns: None - loads the model in the GEMPRO object instance

        """
        extension = op.splitext(gem_file)[1]
        if file_type.replace('.', '').lower() == 'xml' or file_type.replace('.', '').lower() == 'sbml':
            self.model = create_cobra_model_from_sbml_file(gem_file, print_time=False)
            print('[PREP] Loaded model: {}'.format(gem_file))
        elif extension.replace('.', '').lower() == 'mat':
            self.model = load_matlab_model(gem_file)
            print('[PREP] Loaded model: {}'.format(gem_file))

        # place a copy of the current used model in model_files
        if not op.exists(op.join(self.model_files, op.basename(gem_file))):
            shutil.copy(gem_file, self.model_files)
            print('[PREP] Copied model file to model_files')
        self.gem_file = op.join(self.model_files, op.basename(gem_file))

        # obtain list of all gene ids, excluding spontaneous things
        # TODO: account for isoforms
        self.genes_unfiltered = [x.id for x in self.model.genes]
        self.genes = [y for y in self.genes_unfiltered if not is_spontaneous(y)]

        # print information on the number of things
        print('[MODEL] Number of reactions: {}'.format(len(self.model.reactions)))
        print('[MODEL] Number of reactions linked to a gene: {}'.format(true_num_reactions(self.model)))
        print('[MODEL] Number of genes (excluding spontaneous): {}'.format(true_num_genes(self.model)))
        print('[MODEL] Number of metabolites: {}'.format(len(self.model.metabolites)))

    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the KEGG service. Also download all metadata and sequences

        Args:
            kegg_organism_code (str): The three letter KEGG code of your organism
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model genes.
            force_rerun (bool): Flag to indicate if function should be rerun (True) or if
                cached files that match the format *_kegg_df.csv should be used instead (default, False).

        Returns: path to metadata/mapping dataframe

        """

        # all data will be stored in a dataframe
        kegg_pre_df = []

        # first map all of the organism's kegg genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        self.kegg_missing_genes = []

        kegg_df_outfile = op.join(self.data, '{}-kegg_df.csv'.format(date.short_date))
        err_file = op.join(self.missing, '{}-kegg_mapping.err'.format(date.short_date))

        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-kegg_df.csv', self.data)
            saved_err = utils.rank_dated_files('*-kegg_mapping.err', self.missing)

            if saved_err:
                to_load_err = saved_err[0]
                err_date = op.basename(to_load_err).split('-')[0]

                with open(to_load_err, mode='rt', encoding='utf-8') as error_file:
                    self.kegg_missing_genes = error_file.read().splitlines()
                    print('[INFO] Loaded KEGG mapping errors dated {}'.format(err_date))

            if saved_dfs:
                to_load_df = saved_dfs[0]
                df_date = op.basename(to_load_df).split('-')[0]
                self.kegg_df = pd.read_csv(to_load_df, index_col=0)
                print('[INFO] Loaded KEGG mapping dataframe dated {}'.format(df_date))
                return

        for g in tqdm(self.genes):

            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[g]
            else:
                kegg_g = g

            kegg_dict = {}
            kegg_dict['m_gene'] = g
            kegg_dict['k_kegg_id'] = kegg_g

            # make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # download kegg metadata
            metadata_file = ssbio.databases.kegg.download_kegg_gene_metadata(organism_code=kegg_organism_code,
                                                                             gene_id=kegg_g,
                                                                             out_dir=gene_folder)

            # download kegg sequence
            sequence_file = ssbio.databases.kegg.download_kegg_aa_seq(organism_code=kegg_organism_code,
                                                                      gene_id=kegg_g,
                                                                      out_dir=gene_folder)

            # if there is no FASTA for this gene, consider it missing and do not add it to the dataframe
            if not sequence_file:
                self.kegg_missing_genes.append(g)
                continue
            else:
                kegg_dict['k_seq_file'] = op.basename(sequence_file)
                kegg_dict['k_seq_len'] = len(SeqIO.read(open(sequence_file), "fasta"))

            if metadata_file:
                kegg_dict['k_metadata_file'] = op.basename(metadata_file)

                # parse PDB mapping from metadata file
                # in some cases, KEGG IDs do not map to a UniProt ID - this is to ensure we get the PDB mapping
                with open(metadata_file) as mf:
                    kegg_parsed = bs_kegg.parse(mf.read())
                if 'STRUCTURE' in kegg_parsed.keys():
                    kegg_dict['k_pdb_id'] = kegg_parsed['STRUCTURE']
                    # TODO: there is a lot more you can get from the KEGG metadata file, examples below (consider saving it in the DF)
                    # 'DBLINKS': {'NCBI-GeneID': '100763844', 'NCBI-ProteinID': 'XP_003514445'}
                    # 'ORTHOLOGY': {'K00473': 'procollagen-lysine,2-oxoglutarate 5-dioxygenase 1 [EC:1.14.11.4]'},
                    # 'PATHWAY': {'cge00310': 'Lysine degradation'},

            if kegg_g in kegg_to_uniprot.keys():
                kegg_dict['k_uniprot_acc'] = kegg_to_uniprot[kegg_g]

            kegg_pre_df.append(kegg_dict)

        # save a log file for the genes that could not be mapped
        if self.kegg_missing_genes:
            warnings.warn(
                '{} genes were not able to be mapped to a KEGG entry. Please see {} for a list of these genes'.format(
                    len(self.kegg_missing_genes),
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.kegg_missing_genes))

        # save a dataframe of the file mapping info
        cols = ['m_gene', 'k_kegg_id', 'k_uniprot_acc', 'k_pdb_id', 'k_metadata_file', 'k_seq_file',
                'k_seq_len']
        self.kegg_df = pd.DataFrame(kegg_pre_df)[cols].set_index('m_gene')
        self.kegg_df.to_csv(kegg_df_outfile)
        print('[INFO] Saved KEGG information at {}'.format(kegg_df_outfile))

        return kegg_df_outfile

    def uniprot_mapping_and_metadata(self, model_gene_source, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the UniProt mapping service. Also download all metadata and sequences.

        Args:
            model_gene_source (str): the database source of your model gene IDs, see http://www.uniprot.org/help/programmatic_access
                Common model gene sources are:
                Ensembl Genomes - ENSEMBLGENOME_ID
                Entrez Gene (GeneID) - P_ENTREZGENEID
                RefSeq Protein - P_REFSEQ_AC
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model genes.
            force_rerun (bool): Flag to indicate if function should be rerun (True) or if
                cached files that match the format *_uniprot_df.csv should be used instead (default, False).

        Returns: path to metadata/mapping dataframe

        """

        # allow model gene --> custom ID mapping
        if custom_gene_mapping:
            # create a reverse mapping dictionary to map from custom ID -> model gene later
            rev_map = bidict(custom_gene_mapping)
            genes_to_map = list(custom_gene_mapping.values())
        else:
            genes_to_map = self.genes

        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-uniprot_df.csv', self.data)
            saved_err = utils.rank_dated_files('*-uniprot_mapping.err', self.missing)

            if saved_err:
                to_load_err = saved_err[0]
                err_date = op.basename(to_load_err).split('-')[0]

                with open(to_load_err, mode='rt', encoding='utf-8') as error_file:
                    self.uniprot_missing_genes = error_file.read().splitlines()
                    print('[INFO] Loaded UniProt mapping errors dated {}'.format(err_date))

            if saved_dfs:
                to_load_df = saved_dfs[0]
                df_date = op.basename(to_load_df).split('-')[0]
                self.uniprot_df = pd.read_csv(to_load_df, index_col=0)
                print('[INFO] Loaded UniProt mapping dataframe dated {}'.format(df_date))
                return

        genes_to_uniprots = bsup.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)
        self.uniprot_missing_genes = list(set(genes_to_map).difference(genes_to_uniprots.keys()))

        # save a log file for the genes that could not be mapped
        if self.uniprot_missing_genes:
            err_file = op.join(self.missing, '{}-uniprot_mapping.err'.format(date.short_date))
            warnings.warn(
                '[INFO] {} genes were not able to be mapped to a UniProt entry. Please see {} for a list of these genes'.format(
                    len(self.uniprot_missing_genes),
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.uniprot_missing_genes))

        # all data will be stored in a dataframe
        uniprot_pre_df = []
        for g, v in tqdm(genes_to_uniprots.items()):
            # TODO: reverse mapping to custom ids could probably be done in a better way
            if custom_gene_mapping:
                original_m_gene = rev_map[g]
            else:
                original_m_gene = g

            for mapped_uniprot in v:
                uniprot_dict = {}
                uniprot_dict['m_gene'] = original_m_gene
                uniprot_dict['u_uniprot_acc'] = mapped_uniprot

                # make the gene specific folder under the sequence_files directory
                gene_folder = op.join(self.seq_files, g)
                if not op.exists(gene_folder):
                    os.mkdir(gene_folder)

                # download uniprot metadata
                metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot, filetype='txt',
                                                                              outdir=gene_folder)

                # download uniprot sequence
                sequence_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot,
                                                                              filetype='fasta',
                                                                              outdir=gene_folder)

                uniprot_dict['u_seq_file'] = op.basename(sequence_file)
                uniprot_dict['u_metadata_file'] = op.basename(metadata_file)

                # adding additional uniprot metadata
                metadata = ssbio.databases.uniprot.parse_uniprot_txt_file(metadata_file)
                uniprot_dict.update(metadata)

                uniprot_pre_df.append(uniprot_dict)

        # save a dataframe of the file mapping info
        cols = ['m_gene', 'u_uniprot_acc', 'u_gene_name', 'u_reviewed', 'u_seq_len', 'u_kegg_id', 'u_refseq',
                'u_ec_number', 'u_pfam', 'u_description', 'u_entry_version', 'u_seq_version', 'u_metadata_file',
                'u_seq_file']
        self.uniprot_df = pd.DataFrame(uniprot_pre_df)[cols].set_index('m_gene')
        uniprot_df_outfile = op.join(self.data, '{}-uniprot_df.csv'.format(date.short_date))
        self.uniprot_df.to_csv(uniprot_df_outfile)
        print('[INFO] Saved UniProt information at {}'.format(uniprot_df_outfile))

        return uniprot_df_outfile

    def cosonlidate_mappings(self):
        """Combine information from KEGG, UniProt, and manual mappings.

        Returns:

        """
        # compare gene -> uniprot mappings

        # compare sequence lengths

        # if different, compare sequences?

        # sometimes KEGG won't map to a uniprot ID, but will map to a PDB - this needs to be saved somewhere


        pass

    def blast_seqs_to_pdb(self, custom_df, evalue=0.001, link=False):
        """BLAST each gene sequence to the PDB. Results are saved as dataframes in the structure_files folder.

        Returns: Path to summary results dataframe which details best hit for each gene.

        """

        top_blast_pre_df = []

        for i, r in tqdm(custom_df.iterrows()):
            g = r['m_gene']
            seq = r['u_seq']

            # make the gene specific folder under the structure_files directory
            gene_folder = op.join(self.struct_single_chain, str(g))
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            blast_df = ssbio.databases.pdb.blast_pdb(seq, evalue=evalue, link=link)
            if blast_df.empty:
                continue
            else:
                # save the blast_df
                blast_df.to_csv(op.join(gene_folder, str(g) + '_blast_df.csv'))
                # save the top blast hit
                top_hit = blast_df.loc[0].to_dict()
                top_hit['m_gene'] = g
                top_blast_pre_df.append(top_hit)

        top_blast_df = pd.DataFrame(top_blast_pre_df)
        reorg = ['m_gene', 'hit_pdb', 'hit_pdb_chain', 'hit_evalue', 'hit_score', 'hit_num_ident', 'hit_percent_ident']
        self.top_blast_df = top_blast_df[reorg]
        top_blast_df_outfile = op.join(self.data, '{}-pdb_blast_df.csv'.format(date.short_date))
        self.top_blast_df.to_csv(top_blast_df_outfile)

        return top_blast_df_outfile

    def organize_itasser_models(self, raw_dir):
        """Reorganize the raw information from I-TASSER modeling.

        - Copies the model1.pdb file to the matching gene folder in structure_files/by_gene/
        - Parses modeling information such as C-scores, template, etc
        - Also will parse binding site (COACH) information if it exists

        TODO: what else to parse?
        TODO: finish summary df stuff

        Returns: path to summary results dataframe
        """
        hom_pre_df = []

        mode11_counter = 0
        model1_bsites_counter = 0

        for g in tqdm(self.genes):
            hom_dict = {}
            g_files = glob.glob(op.join(raw_dir, g + '*'))

            ### homology model results
            model_folder = op.join(raw_dir, g)
            # best homology model is "model1.pdb"
            model1_file = op.join(model_folder, 'model1.pdb')
            init_dat = 'init.dat'
            init_dat_path = op.join(model_folder, init_dat)

            ### coach (binding site prediction) results
            model1_coach_folder = op.join(model_folder, 'model1', 'coach')
            bsites_inf = 'Bsites.inf'
            bsites_inf_path = op.join(model1_coach_folder, bsites_inf)
            ec = 'EC.dat'
            ec_path = op.join(model1_coach_folder, ec)
            go_mf = 'GO_MF.dat'
            go_mf_path = op.join(model1_coach_folder, go_mf)
            go_bp = 'GO_BP.dat'
            go_bp_path = op.join(model1_coach_folder, go_bp)
            go_cc = 'GO_CC.dat'
            go_cc_path = op.join(model1_coach_folder, go_cc)

            # if the homology model has been completed
            if op.exists(model1_file):
                # make the destination gene folder
                dest_gene_dir = op.join(self.struct_single_chain, g)
                if not op.exists(dest_gene_dir):
                    os.mkdir(dest_gene_dir)

                # the new model file will be named by the gene
                dest = op.join(dest_gene_dir, '{}.pdb'.format(g))
                if not op.exists(dest):
                    shutil.copy2(model1_file, dest)
                hom_dict['m_gene'] = g
                hom_dict['i_model'] = '{}.pdb'.format(g)
                hom_pre_df.append(hom_dict)
                mode11_counter += 1

            # if the binding site predictions have been completed
            if op.exists(model1_coach_folder):
                # make the destination folder
                dest_coach_dir = op.join(self.struct_single_chain, g, 'coach')
                if not op.exists(dest_coach_dir):
                    os.mkdir(dest_coach_dir)

                new_bsites_inf = op.join(dest_coach_dir, bsites_inf)
                if op.exists(bsites_inf_path):
                    if not op.exists(new_bsites_inf):
                        shutil.copy2(bsites_inf_path, new_bsites_inf)
                model1_bsites_counter += 1

                new_ec = op.join(dest_coach_dir, ec)
                if op.exists(ec_path) and not op.exists(new_ec):
                    shutil.copy2(ec_path, new_ec)

                new_go_mf = op.join(dest_coach_dir, go_mf)
                if op.exists(go_mf_path) and not op.exists(new_go_mf):
                    shutil.copy2(go_mf_path, new_go_mf)

                new_go_bp = op.join(dest_coach_dir, go_bp)
                if op.exists(go_bp_path) and not op.exists(new_go_bp):
                    shutil.copy2(go_bp_path, new_go_bp)

                new_go_cc = op.join(dest_coach_dir, go_cc)
                if op.exists(go_cc_path) and not op.exists(new_go_cc):
                    shutil.copy2(go_cc_path, new_go_cc)

        print('[INFO] {} homology models completed and copied to {}.'.format(mode11_counter, self.struct_single_chain))
        if model1_bsites_counter:
            print('[INFO] {} binding site predictions completed and copied.'.format(model1_bsites_counter))

        self.homology_df = pd.DataFrame(hom_pre_df)[['m_gene', 'i_model']]
        self.homology_df.to_csv(op.join(self.data, '{}-homology_models_df.csv'.format(date.short_date)))
        return op.join(self.data, '{}-homology_models_df.csv'.format(date.short_date))

    def run_pipeline(self):
        """Run the entire GEM-PRO pipel ine.

        Options include:
        ...

        Returns:

        """
        self.prep_folders()
        self.load_model()
        # map gene IDs
        # check for missing maps, request manual mapping file?
        pass


if __name__ == '__main__':
    # run the GEM-PRO pipeline!

    # parse arguments
    p = argparse.ArgumentParser(description='Runs the GEM-PRO pipeline')
    p.add_argument('gemfile', help='Path to the GEM file')
    p.add_argument('gemname', help='Name you would like to use to refer to this GEM')
    p.add_argument('rootdir', help='Directory where GEM-PRO files should be stored')

    args = p.parse_args()

    my_gem_pro = GEMPRO(args.gemfile, args.gemname, args.rootdir)
    my_gem_pro.run_pipeline()
