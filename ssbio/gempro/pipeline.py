import argparse
import os
import os.path as op
import shutil
import warnings

from cobra.io import load_matlab_model
from cobra.io.sbml import create_cobra_model_from_sbml_file

from Bio import SeqIO

from bioservices.uniprot import UniProt

bsup = UniProt()

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
        # data - directory where all data frames will be stored (all stages)
        self.data = op.join(self.model_dir, 'data')

        # notebooks - directory where ipython notebooks will be stored for manual analyses
        self.notebooks = op.join(self.model_dir, 'notebooks')

        # # figures - directory where all figures will be stored (all stages)
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

    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_list=None):
        """Map all genes in the model to UniProt IDs using the KEGG service. Also download all metadata and sequences

        Args:
            kegg_organism_code: the three letter KEGG code of your organism

        Returns: path to metadata/mapping dataframe

        """

        # all data will be stored in a dataframe
        kegg_pre_df = []

        # first map all of the organism's genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        self.kegg_missing_genes = []

        if custom_gene_list:
            genes_to_map = custom_gene_list
        else:
            genes_to_map = self.genes

        for g in tqdm(genes_to_map):
            kegg_dict = {}
            kegg_dict['k_kegg_id'] = g

            # make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # download kegg metadata
            metadata_file = ssbio.databases.kegg.download_kegg_gene_metadata(organism_code=kegg_organism_code,
                                                                             gene_id=g,
                                                                             out_dir=gene_folder)

            # download kegg sequence
            sequence_file = ssbio.databases.kegg.download_kegg_aa_seq(organism_code=kegg_organism_code,
                                                                      gene_id=g,
                                                                      out_dir=gene_folder)

            # if there is no FASTA for this gene, consider it missing
            if sequence_file:
                kegg_dict['k_seq_file'] = op.basename(sequence_file)
                kegg_dict['k_seq_len'] = len(SeqIO.read(open(sequence_file), "fasta"))
            else:
                self.kegg_missing_genes.append(g)

            if metadata_file:
                kegg_dict['k_metadata_file'] = op.basename(metadata_file)

            if g in kegg_to_uniprot.keys():
                kegg_dict['k_uniprot_acc'] = kegg_to_uniprot[g]

            kegg_pre_df.append(kegg_dict)

        # save a log file for the genes that could not be mapped
        if self.kegg_missing_genes:
            err_file = op.join(self.missing, '{}_kegg_mapping.err'.format(date.short_date))
            warnings.warn(
                '{} genes were not able to be mapped to a Kegg entry. Please see {} for a list of these genes'.format(
                    len(self.kegg_missing_genes),
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.kegg_missing_genes))

        # save a dataframe of the file mapping info
        # TODO: add more info from kegg metadata parsing
        self.kegg_df = pd.DataFrame(kegg_pre_df)[['k_kegg_id','k_uniprot_acc','k_metadata_file','k_seq_file','k_seq_len']]
        kegg_df_outfile = op.join(self.data, '{}-kegg_df.csv'.format(date.short_date))
        self.kegg_df.to_csv(kegg_df_outfile)
        print('[INFO] Saved KEGG information at {}'.format(kegg_df_outfile))

        return kegg_df_outfile

    def uniprot_mapping_and_metadata(self, model_gene_source, custom_gene_list=None):
        """Map all genes in the model to UniProt IDs using the UniProt mapping service. Also download all metadata and sequences.

        Common model gene sources are:
            Ensembl Genomes - ENSEMBLGENOME_ID
            Entrez Gene (GeneID) - P_ENTREZGENEID
            RefSeq Protein - P_REFSEQ_AC

        Args:
            model_gene_source: the database source of your model gene IDs, see http://www.uniprot.org/help/programmatic_access

        Returns: path to metadata/mapping dataframe

        """

        if custom_gene_list:
            genes_to_map = custom_gene_list
        else:
            genes_to_map = self.genes

        genes_to_uniprots = bsup.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)
        self.uniprot_missing_genes = list(set(genes_to_map).difference(genes_to_uniprots.keys()))

        # save a log file for the genes that could not be mapped
        if self.uniprot_missing_genes:
            err_file = op.join(self.missing, '{}_uniprot_mapping.err'.format(date.short_date))
            warnings.warn(
                '[INFO] {} genes were not able to be mapped to a UniProt entry. Please see {} for a list of these genes'.format(
                    len(self.uniprot_missing_genes),
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.uniprot_missing_genes))

        # all data will be stored in a dataframe
        uniprot_pre_df = []

        for g, v in tqdm(genes_to_uniprots.items()):
            for mapped_uniprot in v:
                uniprot_dict = {}
                uniprot_dict['m_gene'] = g
                uniprot_dict['u_uniprot_acc'] = mapped_uniprot

                # make the gene specific folder under the sequence_files directory
                gene_folder = op.join(self.seq_files, g)
                if not op.exists(gene_folder):
                    os.mkdir(gene_folder)

                # download uniprot metadata
                metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot, filetype='txt',
                                                                              outdir=gene_folder)

                # download uniprot sequence
                sequence_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot, filetype='fasta',
                                                                              outdir=gene_folder)

                uniprot_dict['u_seq_file'] = op.basename(sequence_file)
                uniprot_dict['u_seq_len'] = len(SeqIO.read(open(sequence_file), "fasta"))
                uniprot_dict['u_metadata_file'] = op.basename(metadata_file)

                uniprot_pre_df.append(uniprot_dict)

        # save a dataframe of the file mapping info
        # TODO: add more info from uniprot metadata parsing
        self.uniprot_df = pd.DataFrame(uniprot_pre_df)[
            ['m_gene', 'u_uniprot_acc', 'u_metadata_file', 'u_seq_file', 'u_seq_len']]
        uniprot_df_outfile = op.join(self.data, '{}-uniprot_df.csv'.format(date.short_date))
        self.uniprot_df.to_csv(uniprot_df_outfile)
        print('[INFO] Saved KEGG information at {}'.format(uniprot_df_outfile))

        return uniprot_df_outfile


    def combine_kegg_uniprot_manual(self):
        """Combine information from KEGG, UniProt, and manual mappings.

        Returns:

        """
        pass

    def blast_seqs_to_pdb(self, custom_df, evalue=0.001, link=False):
        """BLAST each gene sequence to the PDB. Results are saved as dataframes in the structure_files folder.

        Returns: Path to summary results dataframe which details best hit for each gene.

        """

        top_blast_pre_df = []

        for i,r in tqdm(custom_df.iterrows()):
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
        reorg = ['m_gene','hit_pdb', 'hit_pdb_chain', 'hit_evalue', 'hit_score', 'hit_num_ident', 'hit_percent_ident']
        self.top_blast_df = top_blast_df[reorg]
        top_blast_df_outfile = op.join(self.data, '{}-pdb_blast_df.csv'.format(date.short_date))
        self.top_blast_df.to_csv(top_blast_df_outfile)

        return top_blast_df_outfile



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
