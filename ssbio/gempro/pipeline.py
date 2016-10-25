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
from cobra.core import Gene
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(__name__)

from bioservices import KEGG
bs_kegg = KEGG()

from cobra.core import DictList



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

    def __init__(self, gem_name, root_dir, gem_file_path=None, gem_file_type=None, genes_list=None):
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

        # Main attributes
        cols = ['gene', 'uniprot_acc', 'kegg_id', 'seq_len', 'pdbs', 'seq_file', 'metadata_file']
        self.df_sequence_mapping = pd.DataFrame(columns=cols)

        # Load the model
        if gem_file_path and gem_file_type:
            self.model = ssbio.cobra.utils.model_loader(gem_file_path, gem_file_type)
            log.info('Loaded model: {}'.format(gem_file_path))

            # Place a copy of the current used model in model_files
            if not op.exists(op.join(self.model_files, op.basename(gem_file_path))):
                shutil.copy(gem_file_path, self.model_files)
                log.debug('Copied model file to model_files')
            self.gem_file = op.join(self.model_files, op.basename(gem_file_path))

            # Obtain list of all gene ids
            # TODO: only will work for prokaryotic genomes for now
            self.genes = self.model.genes

            # Log information on the number of things
            log.info('Number of reactions: {}'.format(len(self.model.reactions)))
            log.info(
                'Number of reactions linked to a gene: {}'.format(ssbio.cobra.utils.true_num_reactions(self.model)))
            log.info(
                'Number of genes (excluding spontaneous): {}'.format(ssbio.cobra.utils.true_num_genes(self.model)))
            log.info('Number of metabolites: {}'.format(len(self.model.metabolites)))

        # Or, load a list of gene IDs
        elif genes_list:
            self.genes = genes_list
            log.info('Number of genes: {}'.format(len(self._genes)))

        # If neither a model or a list of gene IDs is input, you can still add_genes_by_id later
        else:
            self.genes = []
            log.warning('No model or list of genes input.')

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, genes_list):
        """Set the genes attribute to be a DictList of COBRApy Gene objects.

        Extra "annotation" fields will be added to the objects.

        Args:
            genes_list: DictList of COBRApy Gene objects, or list of gene IDs

        """

        if not isinstance(genes_list, DictList):
            tmp_list = []
            for x in list(set(genes_list)):
                new_gene = Gene(id=x)

                # TODO: clean up this dict
                new_gene.annotation['sequence'] = {'kegg'          : {'uniprot_acc'  : None,
                                                                      'kegg_id'      : None,
                                                                      'seq_len'      : 0,
                                                                      'pdbs'         : [],
                                                                      'seq_file'     : None,
                                                                      'metadata_file': None},
                                                   'uniprot'       : {},
                                                   'representative': {'uniprot_acc'  : None,
                                                                      'kegg_id'      : None,
                                                                      'seq_len'      : 0,
                                                                      'pdbs'         : [],
                                                                      'seq_file'     : None,
                                                                      'metadata_file': None}}
                new_gene.annotation['structure'] = {'ranked_pdbs'   : [],
                                                    'blast_pdbs'    : [],
                                                    'representative': {}}
                tmp_list.append(new_gene)
            self._genes = DictList(tmp_list)
        else:
            for x in genes_list:
                x.annotation['sequence'] = {'kegg'          : {'uniprot_acc'  : None,
                                                               'kegg_id'      : None,
                                                               'seq_len'      : 0,
                                                               'pdbs'         : [],
                                                               'seq_file'     : None,
                                                               'metadata_file': None},
                                            'uniprot'       : {},
                                            'representative': {'uniprot_acc'  : None,
                                                               'kegg_id'      : None,
                                                               'seq_len'      : 0,
                                                               'pdbs'         : [],
                                                               'seq_file'     : None,
                                                               'metadata_file': None}}
                x.annotation['structure'] = {'ranked_pdbs'   : [],
                                             'blast_pdbs'    : [],
                                             'representative': {}}
            self._genes = genes_list

    def add_genes_by_id(self, genes_list):
        """Add gene IDs manually into our GEM-PRO project

        Args:
            genes_list (list): List of gene IDs as strings.

        """
        new_genes = []
        for x in list(set(genes_list)):
            new_gene = Gene(id=x)
            new_gene.annotation['sequence'] = {'kegg'          : {'uniprot_acc'  : None,
                                                                  'kegg_id'      : None,
                                                                  'seq_len'      : 0,
                                                                  'pdbs'         : [],
                                                                  'seq_file'     : None,
                                                                  'metadata_file': None},
                                               'uniprot'       : {},
                                               'representative': {'uniprot_acc'  : None,
                                                                  'kegg_id'      : None,
                                                                  'seq_len'      : 0,
                                                                  'pdbs'         : [],
                                                                  'seq_file'     : None,
                                                                  'metadata_file': None}}
            new_gene.annotation['structure'] = {'ranked_pdbs'   : [],
                                                'blast_pdbs'    : [],
                                                'representative': {}}
            new_genes.append(new_gene)

        # Add unique genes only
        self.genes.union(new_genes)

    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the KEGG service. Also download all metadata and sequences.

        Saves the following:
            - df_kegg_mapping attribute
            - df_kegg_mapping.csv in data folder
            - KEGG metadata in each Gene object

        Adds to the df_sequence_mapping attribute, overwriting any existing genes.

        Args:
            kegg_organism_code (str): The three letter KEGG code of your organism
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model genes.

        Returns:
            Path to metadata/mapping dataframe

        """

        # all data will be stored in a dataframe
        kegg_pre_df = []

        # first map all of the organism's kegg genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        kegg_df_outfile = op.join(self.data, 'df_kegg_mapping.csv')

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[gene_id]
            else:
                kegg_g = gene_id

            # Always make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.seq_files, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # For saving the KEGG dataframe
            kegg_dict = {}

            # Download kegg metadata
            metadata_file = ssbio.databases.kegg.download_kegg_gene_metadata(organism_code=kegg_organism_code,
                                                                             gene_id=kegg_g,
                                                                             outdir=gene_folder,
                                                                             force_rerun=force_rerun)
            if metadata_file:
                kegg_dict['metadata_file'] = op.basename(metadata_file)

                # parse PDB mapping from metadata file
                # in some cases, KEGG IDs do not map to a UniProt ID - this is to ensure we get the PDB mapping
                with open(metadata_file) as mf:
                    kegg_parsed = bs_kegg.parse(mf.read())
                    if 'STRUCTURE' in kegg_parsed.keys():
                        kegg_dict['pdbs'] = kegg_parsed['STRUCTURE']['PDB'].split(' ')
                        # TODO: there is a lot more you can get from the KEGG metadata file, examples below
                        # (consider saving it in the DF and in gene)
                        # 'DBLINKS': {'NCBI-GeneID': '100763844', 'NCBI-ProteinID': 'XP_003514445'}
                        # 'ORTHOLOGY': {'K00473': 'procollagen-lysine,2-oxoglutarate 5-dioxygenase 1 [EC:1.14.11.4]'},
                        # 'PATHWAY': {'cge00310': 'Lysine degradation'},

            # Download kegg sequence
            sequence_file = ssbio.databases.kegg.download_kegg_aa_seq(organism_code=kegg_organism_code,
                                                                      gene_id=kegg_g,
                                                                      outdir=gene_folder,
                                                                      force_rerun=force_rerun)
            if sequence_file:
                kegg_dict['kegg_id'] = kegg_organism_code + ':' + kegg_g
                # kegg_dict['seq'] = str(SeqIO.read(open(sequence_file), 'fasta').seq)
                kegg_dict['seq_file'] = op.basename(sequence_file)
                kegg_dict['seq_len'] = len(SeqIO.read(open(sequence_file), "fasta"))

            if kegg_g in kegg_to_uniprot.keys():
                kegg_dict['uniprot_acc'] = kegg_to_uniprot[kegg_g]

            # Save in Gene
            g.annotation['sequence']['kegg'].update(kegg_dict)

            # Save in dataframe
            kegg_dict['gene'] = gene_id
            kegg_pre_df.append(kegg_dict)

            log.debug('Saved KEGG information for gene {}'.format(gene_id))

        # Save a dataframe of the file mapping info
        cols = ['gene', 'uniprot_acc', 'kegg_id', 'seq_len', 'pdbs', 'seq_file', 'metadata_file']
        self.df_kegg_mapping = pd.DataFrame.from_records(kegg_pre_df, columns=cols)
        self.df_kegg_mapping.to_csv(kegg_df_outfile)

        log.info('Saved KEGG dataframe at {}'.format(kegg_df_outfile))

        # Info on genes that could not be mapped
        self.missing_kegg_mapping = self.df_kegg_mapping[pd.isnull(self.df_kegg_mapping.kegg_id)].gene.unique().tolist()
        if len(self.missing_kegg_mapping) > 0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_kegg_mapping" attribute.'.format(len(self.missing_kegg_mapping)))

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

        Returns:
            Path to metadata/mapping dataframe

        """

        # Allow model gene --> custom ID mapping ({'TM_1012':'TM1012'})
        if custom_gene_mapping:
            genes_to_map = list(custom_gene_mapping.values())
        else:
            genes_to_map = [x.id for x in self.genes]

        # Map all IDs first to available UniProts
        genes_to_uniprots = bsup.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)

        uniprot_pre_df = []
        for g in tqdm(self.genes):
            gene_id = str(g.id)
            if custom_gene_mapping and gene_id in custom_gene_mapping.keys():
                uniprot_gene = custom_gene_mapping[gene_id]
            else:
                uniprot_gene = gene_id

            # Always make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.seq_files, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            uniprot_dict = {}

            if uniprot_gene not in list(genes_to_uniprots.keys()):
                # Append empty information for a gene that cannot be mapped
                uniprot_dict['gene'] = gene_id
                uniprot_pre_df.append(uniprot_dict)
            else:
                for mapped_uniprot in genes_to_uniprots[uniprot_gene]:

                    uniprot_dict['uniprot_acc'] = mapped_uniprot

                    # Download uniprot metadata
                    metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot,
                                                                                  filetype='txt',
                                                                                  outdir=gene_folder,
                                                                                  force_rerun=force_rerun)

                    # Download uniprot sequence
                    sequence_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=mapped_uniprot,
                                                                                  filetype='fasta',
                                                                                  outdir=gene_folder,
                                                                                  force_rerun=force_rerun)

                    uniprot_dict['seq_file'] = op.basename(sequence_file)
                    uniprot_dict['metadata_file'] = op.basename(metadata_file)

                    # Adding additional uniprot metadata
                    metadata = ssbio.databases.uniprot.parse_uniprot_txt_file(metadata_file)
                    uniprot_dict.update(metadata)

                    # Save in Gene
                    g.annotation['sequence']['uniprot'][mapped_uniprot] = uniprot_dict

                    # Save in dataframe
                    uniprot_dict['gene'] = gene_id
                    uniprot_pre_df.append(uniprot_dict)

        # Save a dataframe of the UniProt metadata
        if hasattr(self, 'df_uniprot_mapping'):
            self.df_uniprot_mapping = self.df_uniprot_mapping.append(uniprot_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            cols = ['gene', 'uniprot_acc', 'seq_len', 'seq_file', 'pdbs', 'gene_name', 'reviewed', 'kegg_id', 'refseq',
                    'ec_number', 'pfam', 'description', 'entry_version', 'seq_version', 'metadata_file']
            self.df_uniprot_mapping = pd.DataFrame.from_records(uniprot_pre_df, columns=cols)

        uniprot_df_outfile = op.join(self.data, 'df_uniprot_mapping.csv')
        self.df_uniprot_mapping.to_csv(uniprot_df_outfile)
        log.info('Saved UniProt information at {}'.format(uniprot_df_outfile))

        self.missing_uniprot_mapping = self.df_uniprot_mapping[pd.isnull(self.df_uniprot_mapping.kegg_id)].gene.unique().tolist()
        # Info on genes that could not be mapped
        if len(self.missing_uniprot_mapping)>0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_uniprot_mapping" attribute.'.format(
                    len(self.missing_uniprot_mapping)))

        return uniprot_df_outfile

    def manual_uniprot_mapping(self, gene_to_uniprot_dict):
        """Read a manual dictionary of model gene IDs --> UniProt IDs.

        This allows for mapping of the missing genes, or overriding of automatic mappings.

        Args:
            gene_to_uniprot_dict: Dictionary of mappings with key:value pairs:
                <gene_id1>:<uniprot_id1>,
                <gene_id2>:<uniprot_id2>,
                ...
        """
        uniprot_pre_df = []

        for g,u in gene_to_uniprot_dict.items():
            gene = self.genes.get_by_id(g)

            uniprot_dict = {}
            uniprot_dict['uniprot_acc'] = u

            # Make the gene folder
            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # Download uniprot metadata
            metadata_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=u, filetype='txt',
                                                                          outdir=gene_folder)

            # Download uniprot sequence
            sequence_file = ssbio.databases.uniprot.download_uniprot_file(uniprot_id=u,
                                                                          filetype='fasta',
                                                                          outdir=gene_folder)

            uniprot_dict['seq_file'] = op.basename(sequence_file)
            uniprot_dict['metadata_file'] = op.basename(metadata_file)

            # Adding additional uniprot metadata
            metadata = ssbio.databases.uniprot.parse_uniprot_txt_file(metadata_file)
            uniprot_dict.update(metadata)

            # Add info to Gene
            your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
            for_saving = {your_key: uniprot_dict[your_key] for your_key in your_keys}
            gene.annotation['sequence']['representative'].update(for_saving)

            uniprot_dict['gene'] = g
            uniprot_pre_df.append(uniprot_dict)

            # Remove the entry from the missing genes list and also the DF
            # If it is in the DF but not marked as missing, we do not remove the information already mapped
            if hasattr(self, 'df_uniprot_mapping'):
                if g in self.missing_uniprot_mapping:
                    self.missing_uniprot_mapping.remove(g)
                if g in self.df_uniprot_mapping.gene.tolist():
                    get_index = self.df_uniprot_mapping[self.df_uniprot_mapping.gene == g].index
                    for i in get_index:
                        self.df_uniprot_mapping = self.df_uniprot_mapping.drop(i)

        # Save a dataframe of the UniProt metadata
        # Add all new entries to the df_uniprot_mapping
        if hasattr(self, 'df_uniprot_mapping'):
            self.df_uniprot_mapping = self.df_uniprot_mapping.append(uniprot_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            cols = ['gene', 'uniprot_acc', 'seq_len', 'seq_file', 'pdbs', 'gene_name', 'reviewed',
                    'kegg_id', 'refseq', 'ec_number', 'pfam', 'description', 'entry_version', 'seq_version',
                    'metadata_file']
            self.df_uniprot_mapping = pd.DataFrame.from_records(uniprot_pre_df, columns=cols)

        uniprot_df_outfile = op.join(self.data, 'df_uniprot_mapping.csv')
        self.df_uniprot_mapping.to_csv(uniprot_df_outfile)
        log.info('Saved UniProt information at {}'.format(uniprot_df_outfile))

        return uniprot_df_outfile

    def manual_seq_mapping(self, gene_to_seq_dict):
        """Read a manual input dictionary of model gene IDs --> protein sequences.

        Save the sequence in the Gene object.

        Args:
            gene_to_seq_dict: Dictionary of mappings with key:value pairs:
                <gene_id1>:<protein_seq1>,
                <gene_id2>:<protein_seq2>,
                ...

        """

        # save the sequence information in individual FASTA files
        for g, s in gene_to_seq_dict.items():
            gene = self.genes.get_by_id(g)

            gene.annotation['sequence']['representative']['seq_len'] = len(s)

            gene_folder = op.join(self.seq_files, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            seq_file = ssbio.sequence.fasta.write_fasta_file(seq_str=s, ident=g, outdir=gene_folder)
            gene.annotation['sequence']['representative']['seq_file'] = op.basename(seq_file)

            log.info('Loaded manually defined sequence information for {}.'.format(g))

    def set_representative_sequences(self):
        """Combine information from KEGG, UniProt, and manual mappings.

        Manual mappings override all existing mappings.
        UniProt mappings override KEGG mappings except when KEGG mappings have PDBs associated with them and UniProt doesn't.

        Returns:
            Path to df_sequence_mapping dataframe.

        """
        # TODO: clean up this code!
        seq_mapping_pre_df = []

        for gene in self.genes:
            g = gene.id

            genedict = {}
            genedict['gene'] = g
            genedict['pdbs'] = []

            gene_folder = op.join(self.seq_files, g)

            seq_prop = gene.annotation['sequence']

            # If a representative sequence has already been set, nothing needs to be done
            if seq_prop['representative']['seq_len'] > 0:
                genedict['metadata_file'] = seq_prop['representative']['metadata_file']
                genedict['pdbs'] = seq_prop['representative']['pdbs']
                genedict['uniprot_acc'] = seq_prop['representative']['uniprot_acc']
                if isinstance(seq_prop['representative']['kegg_id'], list):
                    genedict['kegg_id'] = ';'.join(seq_prop['representative']['kegg_id'])
                else:
                    genedict['kegg_id'] = seq_prop['representative']['kegg_id']
                genedict['seq_len'] = seq_prop['representative']['seq_len']
                genedict['seq_file'] = seq_prop['representative']['seq_file']
                log.debug('Representative sequence already set for {}'.format(g))

            # If there is a KEGG annotation and no UniProt annotations
            elif seq_prop['kegg']['seq_len'] > 0 and len(seq_prop['uniprot']) == 0:
                kegg_prop = seq_prop['kegg']
                seq_prop['representative'].update(kegg_prop)
                genedict.update(kegg_prop)
                log.debug('Representative sequence set from KEGG for {}'.format(g))

            # If there are UniProt annotations and no KEGG annotations
            elif seq_prop['kegg']['seq_len'] == 0 and len(seq_prop['uniprot']) > 0:

                # If there are multiple uniprots rank them by the sum of reviewed + num_pdbs
                uniprots = seq_prop['uniprot'].keys()
                u_ranker = []
                for u in uniprots:
                    u_ranker.append((u, seq_prop['uniprot'][u]['reviewed'] + len(seq_prop['uniprot'][u]['pdbs'])))
                sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
                best_u = sorted_by_second[0][0]

                uni_prop = seq_prop['uniprot'][best_u]
                your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
                for_saving = { your_key: uni_prop[your_key] for your_key in your_keys }
                seq_prop['representative'].update(for_saving)
                genedict.update(for_saving)
                genedict['kegg_id'] = ';'.join(genedict['kegg_id'])
                log.debug('Representative sequence set from UniProt for {} using {}'.format(g, best_u))

            # If there are both UniProt and KEGG annotations..
            elif seq_prop['kegg']['seq_len'] > 0 and len(seq_prop['uniprot']) > 0:
                # Use KEGG if the mapped UniProt is unique, and it has PDBs
                kegg_prop = seq_prop['kegg']
                if len(kegg_prop['pdbs']) > 0 and kegg_prop['uniprot_acc'] not in seq_prop['uniprot'].keys():
                    seq_prop['representative'].update(kegg_prop)
                    genedict.update(kegg_prop)
                    log.debug('Representative sequence set from KEGG for {}'.format(g))
                else:
                    # If there are multiple uniprots rank them by the sum of reviewed + num_pdbs
                    uniprots = seq_prop['uniprot'].keys()
                    u_ranker = []
                    for u in uniprots:
                        u_ranker.append((u, seq_prop['uniprot'][u]['reviewed'] + len(seq_prop['uniprot'][u]['pdbs'])))
                    sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
                    best_u = sorted_by_second[0][0]

                    uni_prop = seq_prop['uniprot'][best_u]
                    your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
                    for_saving = {your_key: uni_prop[your_key] for your_key in your_keys}
                    seq_prop['representative'].update(for_saving)
                    genedict.update(for_saving)
                    genedict['kegg_id'] = ';'.join(genedict['kegg_id'])
                    log.debug('Representative sequence set from UniProt for {} using {}'.format(g, best_u))

            if genedict['pdbs']:
                genedict['pdbs'] = ';'.join(genedict['pdbs'])
            else:
                genedict['pdbs'] = None
            seq_mapping_pre_df.append(genedict)

        cols = ['gene', 'uniprot_acc', 'kegg_id', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
        tmp = pd.DataFrame.from_records(seq_mapping_pre_df, columns=cols)
        self.missing_final = tmp[pd.isnull(tmp.seq_file)].gene.unique().tolist()
        # Info on genes that could not be mapped
        if len(self.missing_final) > 0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_final" attribute.'.format(
                    len(self.missing_final)))

        self.df_sequence_mapping = tmp[pd.notnull(tmp.seq_file)].reset_index(drop=True)
        self.df_sequence_mapping.fillna(value=np.nan, inplace=True)
        mapping_df_outfile = op.join(self.data, 'df_sequence_mapping.csv')
        self.df_sequence_mapping.to_csv(mapping_df_outfile)
        log.info('Merged ID mapping information at {}'.format(mapping_df_outfile))

        return mapping_df_outfile

    def map_uniprot_to_pdb(self, seq_ident_cutoff=0, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Returns:
            Path to PDB mapping dataframe

        """
        best_structures_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            uniprot_id = g.annotation['sequence']['representative']['uniprot_acc']

            if uniprot_id:
                best_structures = ssbio.databases.pdb.best_structures(uniprot_id,
                                                                      outfile='{}_best_structures.json'.format(uniprot_id),
                                                                      outdir=op.join(self.seq_files, gene_id),
                                                                      force_rerun=force_rerun)

                if best_structures:
                    ranked_pdbs = []
                    rank = 1

                    for best_structure in best_structures:
                        best_structure_dict = {}
                        best_structure_dict['gene'] = gene_id
                        best_structure_dict['uniprot_acc'] = uniprot_id
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

                        # Allow filter for high % seq coverage proteins only
                        # seq_coverage provides the num_ident_in_pdb / total_uniprot_seq_length
                        if best_structure['coverage'] >= seq_ident_cutoff:
                            ranked_pdbs.append(best_structure_dict)

                        # Always add all mapped structures to the dataframe though
                        best_structures_pre_df.append(best_structure_dict)

                        rank += 1

                    g.annotation['structure']['ranked_pdbs'] = ranked_pdbs
                    log.debug('Saved PDB ranking for {}'.format(gene_id))

        cols = ['gene', 'uniprot_acc', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'pdb_resolution',
                'pdb_start', 'pdb_end', 'seq_coverage', 'unp_start', 'unp_end', 'tax_id', 'rank']
        self.df_pdb_ranking = pd.DataFrame.from_records(best_structures_pre_df, columns=cols)
        df_pdb_ranking_outfile = op.join(self.data, 'df_pdb_ranking.csv')
        self.df_pdb_ranking.to_csv(df_pdb_ranking_outfile)

        log.info('Saved UniProt -> PDB mapping at {}'.format(df_pdb_ranking_outfile))
        return df_pdb_ranking_outfile

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, all_genes=False, force_rerun=False, display_link=False):
        """BLAST each gene sequence to the PDB. Results are saved as dataframes in the structure_files folder.

        Returns: Path to summary results dataframe which details best hit for each gene.

        """

        # Save the best blast hits per gene
        top_blast_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(str(g.id))
            seq_file = g.annotation['sequence']['representative']['seq_file']
            seq_dir, seq_name, seq_ext = utils.split_folder_and_path(seq_file)

            # Check if a representative sequence was set
            if not seq_file:
                log.warning('No sequence set for {}'.format(gene_id))
                continue

            # If all_genes=False, BLAST only genes without a uniprot->pdb mapping
            already_ranked_pdbs = g.annotation['structure']['ranked_pdbs']
            if already_ranked_pdbs and not all_genes:
                log.debug('Skipping BLAST for {}, structures already mapped and all_genes flag is False'.format(gene_id))
                continue

            # Read the sequence
            seq_file_path = op.join(self.seq_files, gene_id, seq_file)
            seq_record = SeqIO.read(open(seq_file_path), "fasta")
            seq_str = str(seq_record.seq)

            # Make the gene specific folder under the structure_files directory
            gene_folder = op.join(self.struct_single_chain, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # BLAST the sequence to the PDB
            blast_results, blast_df = ssbio.databases.pdb.blast_pdb_df(seq_str,
                                                                       xml_outfile='{}_blast_pdb.xml'.format(seq_name),
                                                                       xml_outdir=op.join(self.seq_files, gene_id),
                                                                       force_rerun=force_rerun,
                                                                       evalue=evalue,
                                                                       link=display_link)

            # If no blast results are returned, move on to the next sequence
            if not blast_results:
                log.error('No BLAST results for {}'.format(gene_id))
                continue

            # Save blast hits in gene annotation
            g.annotation['structure']['blast_pdbs'] = blast_results

            # Save top blast hits in summary dataframe
            self.df_blast_pdb = blast_df[blast_df.hit_percent_ident >= seq_ident_cutoff]

    # def download_structures(self):
    #     """Download all structures which have been mapped to our genes.
    #
    #     Returns:
    #
    #     """
    #     genes_plus_structure_df = self.df_pdb_ranking[['m_gene','pdb_id','pdb_chain_id']]
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

        log.info('{} homology models completed and copied to {}.'.format(mode11_counter, self.struct_single_chain))
        if model1_bsites_counter:
            log.info('{} binding site predictions completed and copied.'.format(model1_bsites_counter))

        self.homology_df = pd.DataFrame(hom_pre_df)[['m_gene', 'i_model']]
        self.homology_df.to_csv(op.join(self.data, '{}-homology_models_df.csv'.format(date.short_date)))
        return op.join(self.data, '{}-homology_models_df.csv'.format(date.short_date))

    def run_pipeline(self):
        """Run the entire GEM-PRO pipeline.

        Options include:
        ...

        Returns:

        """
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
