import argparse
import os
import os.path as op
import shutil
import warnings
import glob
import numpy as np

from Bio import SeqIO

from bioservices.uniprot import UniProt

bsup = UniProt()
from bidict import bidict

import ssbio.cobra.utils
from ssbio import utils

date = utils.Date()

import ssbio.databases.kegg
import ssbio.databases.uniprot
import ssbio.databases.pdb
import ssbio.sequence.fasta
import pandas as pd

from tqdm import tqdm
import requests

import ssbio.databases.pdb

import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(__name__)

from bioservices import KEGG
bs_kegg = KEGG()


class GEMPRO(object):
    """Generic class to represent all information of a GEM-PRO for a GEM.

    Main steps are:
    1. Preparation of folders
    2. Automatic ID mapping
    3. QC/QA
    4. Homology modeling
    5. Representative structures

    Each step may generate a report and also allow input of manual mappings if applicable.
    """

    def __init__(self, gem_name, root_dir, gem_file_path=None, gem_file_type=None, genes=None):
        """Initialize the GEM-PRO project with a GEM or a list of genes.

        Specify the name of your project, along with the root directory where a folder with that name will be created.

        Args:
            gem_name:
            root_dir:
            gem_file_path:
            gem_file_type:
            genes:
        """
        self.root_dir = root_dir
        self.model_dir = op.join(root_dir, gem_name)

        # model_files - directory where original gems and gem-related files are stored
        self.model_files = op.join(self.model_dir, 'model_files')

        # data - directory where all data will be stored
        self.data = op.join(self.model_dir, 'data')

        # notebooks_dir - directory where ipython notebooks_dir will be stored for manual analyses
        self.notebooks_dir = op.join(self.model_dir, 'notebooks')

        # # figures - directory where all figures will be stored
        self.figures = op.join(self.model_dir, 'figures')

        # struct_files - directory where structure related files will be downloaded/are located
        self.struct_files = op.join(self.model_dir, 'structure_files')
        self.struct_single_chain = op.join(self.struct_files, 'by_gene')

        # seq_files - sequence related files are stored here
        self.seq_files = op.join(self.model_dir, 'sequence_files')

        # Create directory tree
        for directory in [self.model_dir, self.data, self.notebooks_dir, self.figures,
                          self.model_files, self.struct_files, self.struct_single_chain, self.seq_files]:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('Created directory: {}'.format(directory))
            else:
                log.info('Directory already exists: {}'.format(directory))

        if gem_file_path and gem_file_type:
            self.model = ssbio.cobra.utils.model_loader(gem_file_path, gem_file_type)
            log.info('Loaded model: {}'.format(gem_file_path))

            # Place a copy of the current used model in model_files
            if not op.exists(op.join(self.model_files, op.basename(gem_file_path))):
                shutil.copy(gem_file_path, self.model_files)
                log.debug('Copied model file to model_files')
            self.gem_file = op.join(self.model_files, op.basename(gem_file_path))

            # Obtain list of all gene ids, excluding spontaneous things
            # TODO: account for isoforms
            genes_unfiltered = [x.id for x in self.model.genes]
            genes = [y for y in genes_unfiltered if not ssbio.cobra.utils.is_spontaneous(y)]
            self.genes = list(set(genes))

            # Log information on the number of things
            log.info('Number of reactions: {}'.format(len(self.model.reactions)))
            log.info(
                'Number of reactions linked to a gene: {}'.format(ssbio.cobra.utils.true_num_reactions(self.model)))
            log.info(
                'Number of genes (excluding spontaneous): {}'.format(ssbio.cobra.utils.true_num_genes(self.model)))
            log.info('Number of metabolites: {}'.format(len(self.model.metabolites)))

        if genes:
            self.genes = genes
            log.info('Number of genes: {}'.format(len(genes)))

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

        # all data_dir will be stored in a dataframe
        kegg_pre_df = []

        # first map all of the organism's kegg genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        kegg_df_outfile = op.join(self.data, '{}-kegg_df.csv'.format(date.short_date))

        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-kegg_df.csv', self.data)

            if saved_dfs:
                kegg_df_outfile = saved_dfs[0]
                df_date = op.basename(kegg_df_outfile).split('-')[0]
                self.kegg_df = pd.read_csv(kegg_df_outfile, index_col=0)
                self.kegg_missing = self.kegg_df[pd.isnull(self.kegg_df.k_kegg_id)].m_gene.unique().tolist()
                log.info('Loaded KEGG mapping dataframe dated {}'.format(df_date))
                return

        for g in tqdm(list(set(self.genes))):

            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[g]
            else:
                kegg_g = g

            kegg_dict = {}
            kegg_dict['m_gene'] = g

            # always make the gene specific folder under the sequence_files directory
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

            # if there is no FASTA for this gene, consider it missing but still add empty info to the DF
            if sequence_file:
                kegg_dict['k_kegg_id'] = kegg_g
                kegg_dict['k_seq_file'] = op.basename(sequence_file)
                kegg_dict['k_seq_len'] = len(SeqIO.read(open(sequence_file), "fasta"))

            if metadata_file:
                kegg_dict['k_metadata_file'] = op.basename(metadata_file)

                # parse PDB mapping from metadata file
                # in some cases, KEGG IDs do not map to a UniProt ID - this is to ensure we get the PDB mapping
                with open(metadata_file) as mf:
                    kegg_parsed = bs_kegg.parse(mf.read())
                if 'STRUCTURE' in kegg_parsed.keys():
                    kegg_dict['k_pdb_id'] = ','.join(kegg_parsed['STRUCTURE']['PDB'].split(' '))
                    # TODO: there is a lot more you can get from the KEGG metadata file, examples below (consider saving it in the DF)
                    # 'DBLINKS': {'NCBI-GeneID': '100763844', 'NCBI-ProteinID': 'XP_003514445'}
                    # 'ORTHOLOGY': {'K00473': 'procollagen-lysine,2-oxoglutarate 5-dioxygenase 1 [EC:1.14.11.4]'},
                    # 'PATHWAY': {'cge00310': 'Lysine degradation'},
                else:
                    kegg_dict['k_pdb_id'] = np.nan

            if kegg_g in kegg_to_uniprot.keys():
                kegg_dict['k_uniprot_acc'] = kegg_to_uniprot[kegg_g]

            kegg_pre_df.append(kegg_dict)

        # save a dataframe of the file mapping info
        cols = ['m_gene', 'k_kegg_id', 'k_uniprot_acc', 'k_seq_len', 'k_pdb_id', 'k_metadata_file', 'k_seq_file']
        self.kegg_df = pd.DataFrame(kegg_pre_df)[cols]
        self.kegg_df.to_csv(kegg_df_outfile)

        log.info('Saved KEGG information at {}'.format(kegg_df_outfile))

        # info on genes that could not be mapped
        self.kegg_missing = self.kegg_df[pd.isnull(self.kegg_df.k_kegg_id)].m_gene.unique().tolist()
        log.warning('{} gene(s) could not be mapped. Inspect the "kegg_missing" attribute.'.format(
                len(self.kegg_missing)))

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

        Returns:
            Path to metadata/mapping dataframe

        """

        # allow model gene --> custom ID mapping
        if custom_gene_mapping:
            # create a reverse mapping dictionary to map from custom ID -> model gene later
            genes_to_map = list(custom_gene_mapping.values())
            rev_map = bidict(custom_gene_mapping).inv
        else:
            genes_to_map = list(set(self.genes))

        # TODO: allow loading of previous dataframe + adding of new genes (instead of all force rerun)
        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-uniprot_df.csv', self.data)

            if saved_dfs:
                to_load_df = saved_dfs[0]
                df_date = op.basename(to_load_df).split('-')[0]
                self.uniprot_df = pd.read_csv(to_load_df, index_col=0)
                self.uniprot_missing = self.uniprot_df[
                    pd.isnull(self.uniprot_df.u_uniprot_acc)].m_gene.unique().tolist()
                print('[INFO] Loaded UniProt mapping dataframe dated {}'.format(df_date))
                return

        genes_to_uniprots = bsup.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)
        self.uniprot_missing = list(set(genes_to_map).difference(genes_to_uniprots.keys()))


        # all data_dir will be stored in a dataframe
        uniprot_pre_df = []
        for g in tqdm(genes_to_map):
            if custom_gene_mapping:
                original_m_gene = rev_map[g]
            else:
                original_m_gene = g

            # always make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.seq_files, original_m_gene)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            if g not in list(genes_to_uniprots.keys()):
                # append empty information for a gene that cannot be mapped
                uniprot_dict = {}
                uniprot_dict['m_gene'] = original_m_gene
                uniprot_pre_df.append(uniprot_dict)
            else:
                for mapped_uniprot in genes_to_uniprots[g]:
                    uniprot_dict = {}
                    uniprot_dict['m_gene'] = original_m_gene
                    uniprot_dict['u_uniprot_acc'] = mapped_uniprot

                    # download uniprot metadata
                    metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot,
                                                                                  filetype='txt',
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
        self.uniprot_df = pd.DataFrame(uniprot_pre_df)[cols]
        uniprot_df_outfile = op.join(self.data, '{}-uniprot_df.csv'.format(date.short_date))
        self.uniprot_df.to_csv(uniprot_df_outfile)
        print('[INFO] Saved UniProt information at {}'.format(uniprot_df_outfile))

        # info on genes that could not be mapped
        print('[WARN] {} gene(s) could not be mapped. Inspect the "uniprot_missing" attribute.'.format(
                len(self.uniprot_missing)))

        return uniprot_df_outfile

    def manual_uniprot_mapping(self, infile):
        """Read a manual input file of model gene IDs --> UniProt IDs.

        This allows for mapping of the missing genes, or overriding of automatic mappings.

        Args:
            infile: Path to .csv file with the headers:
                m_gene - model gene ID
                u_uniprot_acc - manually determined UniProt ID

        """
        # read the manual mapping file
        # TODO: what happens if you use this function and uniprot_mapping has not been done ?

        man_df = pd.read_csv(infile)

        uniprot_pre_df = []

        # iterate over the rows
        for i, r in man_df.iterrows():
            g = r['m_gene']
            u = r['u_uniprot_acc']

            uniprot_dict = {}
            uniprot_dict['m_gene'] = g
            uniprot_dict['u_uniprot_acc'] = u

            # make the gene folder
            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # download uniprot metadata
            metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=u, filetype='txt',
                                                                          outdir=gene_folder)

            # download uniprot sequence
            sequence_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=u,
                                                                          filetype='fasta',
                                                                          outdir=gene_folder)

            uniprot_dict['u_seq_file'] = op.basename(sequence_file)
            uniprot_dict['u_metadata_file'] = op.basename(metadata_file)

            # adding additional uniprot metadata
            metadata = ssbio.databases.uniprot.parse_uniprot_txt_file(metadata_file)
            uniprot_dict.update(metadata)

            uniprot_pre_df.append(uniprot_dict)
            # remove the entry from the missing genes list and also the DF
            # if it is in the DF but not marked as missing, we do not remove the information already mapped
            if g in self.uniprot_missing:
                self.uniprot_missing.remove(g)
            if g in self.uniprot_df.m_gene.tolist():
                get_index = self.uniprot_df[self.uniprot_df.m_gene == g].index
                for i in get_index:
                    self.uniprot_df = self.uniprot_df.drop(i)

        # add all new entries to the uniprot_df
        self.uniprot_df = self.uniprot_df.append(uniprot_pre_df, ignore_index=True)

        print('[INFO] Updated UniProt dataframe.')

    def manual_seq_mapping(self, infile):
        """Read a manual input file of model gene IDs --> protein sequences.

        Args:
            infile: Path to .csv file with the headers:
                m_gene - model gene ID
                m_sequence - manually specified protein sequence

        Returns:
            Path to mapping dataframe

        """
        # read the manual mapping file
        man_df = pd.read_csv(infile)

        manual_pre_df = []

        # save the sequence information in individual FASTA files
        for i, r in man_df.iterrows():
            g = r['m_gene']
            s = r['m_sequence']

            # make the a manual_mapping_df with added info seq_len & seq_file
            manual_dict = {}
            manual_dict['m_gene'] = g
            manual_dict['m_seq_len'] = len(s)

            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            seq_file = ssbio.sequence.fasta.write_fasta_file(seq_str=s, ident=g, outdir=gene_folder)
            manual_dict['m_seq_file'] = op.basename(seq_file)

            manual_pre_df.append(manual_dict)

        cols = ['m_gene', 'm_seq_len', 'm_seq_file']
        self.manual_df = pd.DataFrame(manual_pre_df)[cols]
        manual_df_outfile = op.join(self.data, '{}-manual_df.csv'.format(date.short_date))
        self.manual_df.to_csv(manual_df_outfile)
        print('[INFO] Saved manually defined sequence information at {}'.format(manual_df_outfile))

        return manual_df_outfile

    def consolidate_mappings(self):
        """Combine information from KEGG, UniProt, and manual mappings.

        Manual mappings override all existing mappings.
        UniProt mappings override KEGG mappings.

        Returns:

        """
        # TODO: add more information printouts
        if hasattr(self, 'manual_df') and not hasattr(self, 'kegg_df') and not hasattr(self, 'uniprot_df'):
            self.mapping_df = self.manual_df
            cols = self.mapping_df.columns
            print('[INFO] Mapping information is solely the manually provided sequences.')

        elif hasattr(self, 'manual_df') and not hasattr(self, 'kegg_df'):
            # merge uniprot and manual info
            filter_uniprot = self.uniprot_df[~self.uniprot_df.m_gene.isin(self.uniprot_missing)]
            filter_uniprot = filter_uniprot[~filter_uniprot.m_gene.isin(self.manual_df.m_gene.tolist())]
            filter_uniprot = filter_uniprot.rename(
                    columns={'u_uniprot_acc': 'm_uniprot_acc', 'u_seq_len': 'm_seq_len', 'u_seq_file': 'm_seq_file'})
            mapping_df = filter_uniprot.merge(self.manual_df, how='outer')
            cols = ['m_gene', 'm_uniprot_acc', 'm_seq_len', 'm_seq_file', 'u_gene_name', 'u_reviewed',
                    'u_kegg_id', 'u_refseq', 'u_ec_number', 'u_pfam', 'u_description',
                    'u_entry_version', 'u_seq_version', 'u_metadata_file']

        elif hasattr(self, 'manual_df') and not hasattr(self, 'uniprot_df'):
            # merge kegg and manual info
            filter_kegg = self.kegg_df[~self.kegg_df.m_gene.isin(self.kegg_missing)]
            filter_kegg = filter_kegg[~filter_kegg.m_gene.isin(self.manual_df.m_gene.tolist())]
            filter_kegg = filter_kegg.rename(
                    columns={'k_uniprot_acc': 'm_uniprot_acc', 'k_seq_len': 'm_seq_len', 'k_seq_file': 'm_seq_file'})
            mapping_df = filter_kegg.merge(self.manual_df, how='outer')
            cols = ['m_gene', 'm_uniprot_acc', 'm_seq_len', 'm_seq_file', 'k_kegg_id', 'k_pdb_id', 'k_metadata_file']

        else:
            # take out missing entries - we do not want these to overwrite any existing mappings
            filter_kegg = self.kegg_df[~self.kegg_df.m_gene.isin(self.kegg_missing)]
            filter_uniprot = self.uniprot_df[~self.uniprot_df.m_gene.isin(self.uniprot_missing)]

            if hasattr(self, 'manual_df'):
                # take out entries that will be replaced by manual mapping
                filter_kegg = filter_kegg[~filter_kegg.m_gene.isin(self.manual_df.m_gene.tolist())]
                filter_uniprot = filter_uniprot[~filter_uniprot.m_gene.isin(self.manual_df.m_gene.tolist())]

            filter_kegg = filter_kegg.rename(columns={'k_uniprot_acc': 'm_uniprot_acc', 'k_seq_len': 'm_seq_len'})
            filter_uniprot = filter_uniprot.rename(columns={'u_uniprot_acc': 'm_uniprot_acc', 'u_seq_len': 'm_seq_len'})

            # these are unique kegg gene->uniprot mappings AND ones that contain a gene->pdb mapping
            # NOTE: this will override any manual uniprot mappings that have a gene->pdb mapping
            a = filter_kegg[(pd.notnull(filter_kegg.k_pdb_id)) & (
                ~filter_kegg.m_uniprot_acc.isin(filter_uniprot.m_uniprot_acc.tolist()))]

            # these are unique kegg gene mappings
            b = filter_kegg[~filter_kegg.m_gene.isin(filter_uniprot.m_gene)]

            # this is all the kegg mappings to add
            one = a.merge(b, how='outer')

            # take out uniprot mappings to be replaced by kegg mappings
            two = filter_uniprot[~filter_uniprot.m_gene.isin(a.m_gene.tolist())]

            # merge uniprot mappings with kegg info
            three = two.merge(filter_kegg.drop('k_seq_file', axis=1), how='left')

            # merge with kegg mappings to add
            mapping_df = three.rename(columns={'u_seq_file': 'm_seq_file'}).merge(
                    one.rename(columns={'k_seq_file': 'm_seq_file'}), how='outer')

            if hasattr(self, 'manual_df'):
                # last, merge manual
                mapping_df = mapping_df.merge(self.manual_df, how='outer')

            cols = ['m_gene', 'm_uniprot_acc', 'm_seq_len', 'm_seq_file', 'u_gene_name', 'u_reviewed',
                    'u_kegg_id', 'u_refseq', 'u_ec_number', 'u_pfam', 'u_description',
                    'u_entry_version', 'u_seq_version', 'u_metadata_file', 'k_kegg_id', 'k_pdb_id', 'k_metadata_file']

        self.mapping_df = pd.DataFrame(mapping_df)[cols]
        mapping_df_outfile = op.join(self.data, '{}-mapping_df.csv'.format(date.short_date))
        self.mapping_df.to_csv(mapping_df_outfile)

        print('[INFO] Merged ID mapping information at {}'.format(mapping_df_outfile))
        return mapping_df_outfile

    def map_uniprot_to_pdb(self, seq_ident_cutoff=0.95, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Returns:
            Path to PDB mapping dataframe

        """
        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-ranked_pdbs.csv', self.data)

            if saved_dfs:
                to_load_df = saved_dfs[0]
                df_date = op.basename(to_load_df).split('-')[0]
                self.ranked_pdbs_df = pd.read_csv(to_load_df, index_col=0)
                print('[INFO] Loaded UniProt -> PDB mapping dataframe dated {}'.format(df_date))
                return

        best_structures_pre_df = []

        for i, r in tqdm(self.mapping_df[pd.notnull(self.mapping_df.m_uniprot_acc)].iterrows()):
            g = r['m_gene']
            u = r['m_uniprot_acc']
            best_structures = ssbio.databases.pdb.best_structures(u)
            if best_structures:
                rank = 1
                for best_structure in best_structures:
                    best_structure_dict = {}
                    best_structure_dict['m_gene'] = g
                    best_structure_dict['m_uniprot_acc'] = u
                    best_structure_dict['pdb_id'] = best_structure['pdb_id']
                    best_structure_dict['pdb_chain_id'] = best_structure['chain_id']
                    best_structure_dict['experimental_method'] = best_structure['experimental_method']
                    best_structure_dict['pdb_resolution'] = best_structure['resolution']
                    best_structure_dict['pdb_start'] = best_structure['start']
                    best_structure_dict['pdb_end'] = best_structure['end']
                    best_structure_dict['seq_coverage'] = best_structure['coverage']
                    best_structure_dict['unp_start'] = best_structure['unp_start']
                    best_structure_dict['unp_end'] = best_structure['unp_end']
                    best_structure_dict['tax_id'] = best_structure['tax_id']
                    best_structure_dict['rank'] = rank
                    best_structures_pre_df.append(best_structure_dict)
                    rank += 1

        cols = ['m_gene', 'm_uniprot_acc', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'pdb_resolution',
                'pdb_start', 'pdb_end', 'seq_coverage', 'unp_start', 'unp_end', 'tax_id', 'rank']
        self.ranked_pdbs_df = pd.DataFrame(best_structures_pre_df)[cols]
        # filter for high % seq ident proteins only
        # TODO: double check the seq_coverage number, does it provide the right %?
        self.ranked_pdbs_df = self.ranked_pdbs_df[self.ranked_pdbs_df.seq_coverage >= seq_ident_cutoff]
        ranked_pdbs_df_outfile = op.join(self.data, '{}-ranked_pdbs.csv'.format(date.short_date))
        self.ranked_pdbs_df.to_csv(ranked_pdbs_df_outfile)

        print('[INFO] Saved UniProt -> PDB mapping at {}'.format(ranked_pdbs_df_outfile))
        return ranked_pdbs_df_outfile

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0.90, all_genes=False, force_rerun=False, evalue=0.001, display_link=False):
        """BLAST each gene sequence to the PDB. Results are saved as dataframes in the structure_files folder.

        Returns: Path to summary results dataframe which details best hit for each gene.

        """
        # TODO: what if force_rerun=False but all_genes=True? and the last run was all_genes=False?
        # answer -- use logging functions...
        if not force_rerun:
            # check for already saved mappings
            # load the latest mapping if there are multiple
            saved_dfs = utils.rank_dated_files('*-top_pdb_blast_df.csv', self.data)

            if saved_dfs:
                to_load_df = saved_dfs[0]
                df_date = op.basename(to_load_df).split('-')[0]
                self.top_blast_df = pd.read_csv(to_load_df, index_col=0)
                print('[INFO] Loaded PDB BLAST top hits dated {}'.format(df_date))
                return

        # save the best blast hits per gene
        top_blast_pre_df = []

        # if all_genes=False, BLAST only genes without a uniprot->pdb mapping
        # requires ranked_pdbs_df to have been made
        if hasattr(self, 'ranked_pdbs_df') and all_genes == False:
            genes_with_pdbs = self.ranked_pdbs_df.m_gene.unique().tolist()
            all_genes = self.mapping_df.m_gene.unique().tolist()
            genes_without_pdbs = [x for x in all_genes if x not in genes_with_pdbs]
            seqs_to_map = self.mapping_df[self.mapping_df.m_gene.isin(genes_without_pdbs)]
        else:
            # if all_genes=True, BLAST all sequences we have
            seqs_to_map = self.mapping_df[pd.notnull(self.mapping_df.m_seq_file)]

        for i, r in tqdm(seqs_to_map.iterrows()):
            g = r['m_gene']
            seq_file = r['m_seq_file']
            seq_file_path = op.join(self.seq_files, g, seq_file)

            # make the gene specific folder under the structure_files directory
            gene_folder = op.join(self.struct_single_chain, str(g))
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # read the sequence
            # TODO: replace with biopython sequence loading..
            seq = ssbio.sequence.fasta.load_fasta_file(seq_file_path)
            seq_str = str(seq[0].seq)

            # BLAST the sequence to the PDB
            # but do not run the request if the output file already exists and force_rerun is False
            blast_outfile = op.join(gene_folder, str(g) + '_blast_df.csv')
            if op.exists(blast_outfile) and force_rerun == False:
                blast_df = pd.read_csv(blast_outfile, index_col=0)
            else:
                blast_df = ssbio.databases.pdb.blast_pdb(seq_str, evalue=evalue, link=display_link)

                # if no blast results are returned, move on to the next sequence
                if blast_df.empty:
                    continue

                # save the blast_df to the structure_files/by_gene folder
                blast_df.to_csv(blast_outfile)

            # retrieve the top hits
            # if the top BLAST hit % sequence identity is less than the cutoff, do not save it as a top hit
            top_hits_df = blast_df[blast_df.hit_percent_ident >= seq_ident_cutoff]

            for i,r in top_hits_df.iterrows():
                hit = r.to_dict()
                hit['m_gene'] = g
                hit['m_seq_file'] = seq_file
                top_blast_pre_df.append(hit)

        top_blast_df = pd.DataFrame(top_blast_pre_df)
        reorg = ['m_gene', 'm_seq_file', 'hit_pdb', 'hit_pdb_chain', 'hit_evalue', 'hit_score', 'hit_num_ident',
                 'hit_percent_ident', 'hit_num_similar', 'hit_percent_similar', 'hit_num_gaps', 'hit_percent_gaps']
        self.top_blast_df = top_blast_df[reorg]
        top_blast_df_outfile = op.join(self.data, '{}-top_pdb_blast_df.csv'.format(date.short_date))
        self.top_blast_df.to_csv(top_blast_df_outfile)

        print('[INFO] Saved PDB BLAST top hits at {}'.format(top_blast_df_outfile))
        return top_blast_df_outfile

    # def download_structures(self):
    #     """Download all structures which have been mapped to our genes.
    #
    #     Returns:
    #
    #     """
    #     genes_plus_structure_df = self.ranked_pdbs_df[['m_gene','pdb_id','pdb_chain_id']]
    #     if hasattr(self, 'top_blast_df'):
    #         blasted_genes = self.top_blast_df.m_gene.tolist()
    #         genes_with_structure.extend(blasted_genes)
    #
    #     for g in self.genes:
    #         if g in genes_with_structure:
    #             pdbs_for_this_gene =


    def organize_itasser_models(self, raw_dir, custom_name_mapping=None):
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
        """Run the entire GEM-PRO pipeline.

        Options include:
        ...

        Returns:

        """
        self.prep_folders()
        self.load_model()
        # map gene IDs
        # check for missing maps, request manual mapping file?
        pass

#
# if __name__ == '__main__':
#     # run the GEM-PRO pipeline!
#
#     # parse arguments
#     p = argparse.ArgumentParser(description='Runs the GEM-PRO pipeline')
#     p.add_argument('gemfile', help='Path to the GEM file')
#     p.add_argument('gemname', help='Name you would like to use to refer to this GEM')
#     p.add_argument('rootdir', help='Directory where GEM-PRO files should be stored')
#
#     args = p.parse_args()
#
#     my_gem_pro = GEMPRO(args.gemfile, args.gemname, args.rootdir)
#     my_gem_pro.run_pipeline()
