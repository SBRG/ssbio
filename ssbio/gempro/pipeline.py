import argparse
import os
import os.path as op
import shutil
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

import ssbio.databases.pdb
from cobra.core import Gene
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)
log = logging.getLogger(__name__)

from bioservices import KEGG
bs_kegg = KEGG()

from cobra.core import DictList
from more_itertools import unique_everseen
# from dotmap import DotMap

"""
notes

force_rerun
if possible, files that require a web request should be downloaded and used as a cache
- metadata files
- sequence files
- structure files (pdb, sifts)
- blast to the pdb
- best_structures api
also, files that require runtime of a program
- dssp
- msms
- cleaned pdbs
- minimized pdbs
- homology models

each function should have a force_rerun flag to download the newest files
also there should be a global setting?

if dataframes are made, it should always return the actual dataframe
it should be left up to the user if they want to save that information
only important dfs should be saved, like the gempro final
everything else should be saved in gene annotations
"""


# TODO using these classes in the annotation field will work
class StructureProp(object):
    def __init__(self, ranked_pdbs=None, blast_pdbs=None, pdb=None, representative=None):
        if not ranked_pdbs:
            ranked_pdbs = []
        if not blast_pdbs:
            blast_pdbs = []
        if not pdb:
            pdb = {}
        if not representative:
            representative = {'pdb_id'             : None,
                              'resolution'         : float('inf'),
                              'original_pdb_file'  : None,
                              'original_mmcif_file': None,
                              'clean_pdb_file'     : None}
        self.ranked_pdbs = ranked_pdbs
        self.blast_pdbs = blast_pdbs
        self.pdb = pdb
        self.representative = representative


class SequenceProp(object):
    def __init__(self, kegg=None, uniprot=None, representative=None):
        if not kegg:
            kegg = {'uniprot_acc'  : None,
                    'kegg_id'      : None,
                    'seq_len'      : 0,
                    'pdbs'         : [],
                    'seq_file'     : None,
                    'metadata_file': None}
        if not uniprot:
            uniprot = {}
        if not representative:
            representative = {'uniprot_acc'  : None,
                              'kegg_id'      : None,
                              'seq_len'      : 0,
                              'pdbs'         : [],
                              'seq_file'     : None,
                              'metadata_file': None}
        self.kegg = kegg
        self.uniprot = uniprot
        self.representative = representative


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
        self.base_dir = op.join(root_dir, gem_name)

        # model_dir - directory where original gems and gem-related files are stored
        self.model_dir = op.join(self.base_dir, 'model')

        # data_dir - directory where all data (dataframes, etc) will be stored
        self.data_dir = op.join(self.base_dir, 'data')

        # notebooks_dir - directory where ipython notebooks_dir will be stored for manual analyses
        self.notebooks_dir = op.join(self.base_dir, 'notebooks')

        # # figure_dir - directory where all figure_dir will be stored
        self.figure_dir = op.join(self.base_dir, 'figures')

        # structure_dir - directory where structure related files will be downloaded/are located
        self.structure_dir = op.join(self.base_dir, 'structures')
        self.structure_single_chain_dir = op.join(self.structure_dir, 'by_gene')

        # sequence_dir - sequence related files are stored here
        self.sequence_dir = op.join(self.base_dir, 'sequences')

        list_of_dirs = [self.base_dir,
                        self.model_dir,
                        self.data_dir,
                        self.notebooks_dir,
                        self.figure_dir,
                        self.structure_dir,
                        self.structure_single_chain_dir,
                        self.sequence_dir]

        # Create directory tree
        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('Created directory: {}'.format(directory))
            else:
                log.debug('Directory already exists: {}'.format(directory))

        # Main attributes
        cols = ['gene', 'uniprot_acc', 'kegg_id', 'seq_len', 'pdbs', 'seq_file', 'metadata_file']
        self.df_sequence_mapping = pd.DataFrame(columns=cols)
        self.missing_mapping = []

        # Load the model
        if gem_file_path and gem_file_type:
            self.model = ssbio.cobra.utils.model_loader(gem_file_path, gem_file_type)
            log.info('Loaded model: {}'.format(gem_file_path))

            # Place a copy of the current used model in model_dir
            if not op.exists(op.join(self.model_dir, op.basename(gem_file_path))):
                shutil.copy(gem_file_path, self.model_dir)
                log.debug('Copied model file to model directory')
            self.gem_file = op.join(self.model_dir, op.basename(gem_file_path))

            # Obtain list of all gene ids
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

                # TODO: if we don't define dictionaries like this we risk pointing to the same object for some reason
                # is there a better way to do this?
                if 'sequence' not in new_gene.annotation.keys():
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
                if 'structure' not in new_gene.annotation.keys():
                    new_gene.annotation['structure'] = {'ranked_pdbs'   : [],
                                                        'blast_pdbs'    : [],
                                                        'pdb'           : {},
                                                        'representative': {'pdb_id'             : None,
                                                                           'resolution'         : float('inf'),
                                                                           'original_pdb_file'  : None,
                                                                           'original_mmcif_file': None,
                                                                           'clean_pdb_file'     : None}}
                tmp_list.append(new_gene)
            self._genes = DictList(tmp_list)
        else:
            for x in genes_list:
                if 'sequence' not in x.annotation.keys():
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
                if 'structure' not in x.annotation.keys():
                    x.annotation['structure'] = {'ranked_pdbs'   : [],
                                                 'blast_pdbs'    : [],
                                                 'pdb'           : {},
                                                 'representative': {'pdb_id'             : None,
                                                                    'resolution'         : float('inf'),
                                                                    'original_pdb_file'  : None,
                                                                    'original_mmcif_file': None,
                                                                    'clean_pdb_file'     : None}}
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
                                                'pdb'           : {},
                                                'representative': {'pdb_id'             : None,
                                                                   'resolution'         : float('inf'),
                                                                   'original_pdb_file'  : None,
                                                                   'original_mmcif_file': None,
                                                                   'clean_pdb_file'     : None}}
            new_genes.append(new_gene)

        # Add unique genes only
        self.genes.union(new_genes)
        # TODO: report union length?
        # log.info('Added {} genes to the list of genes'.format(len(new_genes)))

    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the KEGG service.

        This function does these things:
            1. Download all metadata and sequence files in the sequences directory
            2. Saves KEGG metadata in each Gene object under the "sequence" key
            3. Returns a Pandas DataFrame of mapping results

        Args:
            kegg_organism_code (str): The three letter KEGG code of your organism
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model gene IDs.
            force_rerun (bool): If you want to overwrite any existing mappings and files

        Returns:
            Pandas DataFrame: of KEGG mapping results

        """

        # all data_dir will be stored in a dataframe
        kegg_pre_df = []

        # first map all of the organism's kegg genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[gene_id]
            else:
                kegg_g = gene_id

            # Always make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.sequence_dir, gene_id)
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

            log.debug('{}: Loaded KEGG information for gene'.format(gene_id))

        # Save a dataframe of the file mapping info
        cols = ['gene', 'uniprot_acc', 'kegg_id', 'seq_len', 'pdbs', 'seq_file', 'metadata_file']
        self.df_kegg_metadata = pd.DataFrame.from_records(kegg_pre_df, columns=cols)
        log.info('Created KEGG metadata dataframe.')

        # Info on genes that could not be mapped
        self.missing_kegg_mapping = self.df_kegg_metadata[pd.isnull(self.df_kegg_metadata.kegg_id)].gene.unique().tolist()
        if len(self.missing_kegg_mapping) > 0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_kegg_mapping" attribute.'.format(len(self.missing_kegg_mapping)))

    def uniprot_mapping_and_metadata(self, model_gene_source, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the UniProt mapping service. Also download all metadata and sequences.

        Args:
            model_gene_source (str): the database source of your model gene IDs, see http://www.uniprot.org/help/programmatic_access
                Common model gene sources are:
                Ensembl Genomes - ENSEMBLGENOME_ID (i.e. E. coli b-numbers)
                Entrez Gene (GeneID) - P_ENTREZGENEID
                RefSeq Protein - P_REFSEQ_AC
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model genes.
            force_rerun (bool): If you want to overwrite any existing mappings and files

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
            gene_folder = op.join(self.sequence_dir, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            uniprot_dict = {}

            if uniprot_gene not in list(genes_to_uniprots.keys()):
                # Append empty information for a gene that cannot be mapped
                uniprot_dict['gene'] = gene_id
                uniprot_pre_df.append(uniprot_dict)
                log.debug('{}: Unable to map to UniProt'.format(gene_id))
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
                    if 'pdbs' not in uniprot_dict:
                        uniprot_dict['pdbs'] = []
                    g.annotation['sequence']['uniprot'][mapped_uniprot] = uniprot_dict

                    # Add info to dataframe
                    # TODO: empty pdb lists should be NaN in the dataframe
                    uniprot_dict['gene'] = gene_id
                    uniprot_pre_df.append(uniprot_dict)

        # Save a dataframe of the UniProt metadata
        if hasattr(self, 'df_uniprot_metadata'):
            self.df_uniprot_metadata = self.df_uniprot_metadata.append(uniprot_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            cols = ['gene', 'uniprot_acc', 'seq_len', 'seq_file', 'pdbs', 'gene_name', 'reviewed', 'kegg_id', 'refseq',
                    'ec_number', 'pfam', 'description', 'entry_version', 'seq_version', 'metadata_file']
            self.df_uniprot_metadata = pd.DataFrame.from_records(uniprot_pre_df, columns=cols)
            log.info('Created UniProt metadata dataframe.')

        self.missing_uniprot_mapping = self.df_uniprot_metadata[pd.isnull(self.df_uniprot_metadata.kegg_id)].gene.unique().tolist()
        # Info on genes that could not be mapped
        if len(self.missing_uniprot_mapping) > 0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_uniprot_mapping" attribute.'.format(
                    len(self.missing_uniprot_mapping)))

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
            gene_folder = op.join(self.sequence_dir, g)
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
            # TODO: Setting as representative for now, but should also save in uniprot key
            if 'pdbs' not in uniprot_dict:
                uniprot_dict['pdbs'] = []
            your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
            for_saving = {your_key: uniprot_dict[your_key] for your_key in your_keys}
            gene.annotation['sequence']['representative'].update(for_saving)

            # Add info to dataframe
            # TODO: empty pdb lists should be NaN in the dataframe
            uniprot_dict['gene'] = g
            uniprot_pre_df.append(uniprot_dict)

            # Remove the entry from the missing genes list and also the DF
            # If it is in the DF but not marked as missing, we do not remove the information already mapped
            if hasattr(self, 'df_uniprot_metadata'):
                if g in self.missing_uniprot_mapping:
                    self.missing_uniprot_mapping.remove(g)
                if g in self.df_uniprot_metadata.gene.tolist():
                    get_index = self.df_uniprot_metadata[self.df_uniprot_metadata.gene == g].index
                    for i in get_index:
                        self.df_uniprot_metadata = self.df_uniprot_metadata.drop(i)

        # Save a dataframe of the UniProt metadata
        # Add all new entries to the df_uniprot_metadata
        if hasattr(self, 'df_uniprot_metadata'):
            self.df_uniprot_metadata = self.df_uniprot_metadata.append(uniprot_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            cols = ['gene', 'uniprot_acc', 'seq_len', 'seq_file', 'pdbs', 'gene_name', 'reviewed',
                    'kegg_id', 'refseq', 'ec_number', 'pfam', 'description', 'entry_version', 'seq_version',
                    'metadata_file']
            self.df_uniprot_metadata = pd.DataFrame.from_records(uniprot_pre_df, columns=cols)
            log.info('Created UniProt metadata dataframe.')

    def manual_seq_mapping(self, gene_to_seq_dict):
        """Read a manual input dictionary of model gene IDs --> protein sequences.

        Save the sequence in the Gene object.

        Args:
            gene_to_seq_dict: Dictionary of mappings with key:value pairs:
                {<gene_id1>:<protein_seq1>,
                 <gene_id2>:<protein_seq2>,
                 ...}

        """

        # save the sequence information in individual FASTA files
        for g, s in gene_to_seq_dict.items():
            gene = self.genes.get_by_id(g)

            gene.annotation['sequence']['representative']['seq_len'] = len(s)

            gene_folder = op.join(self.sequence_dir, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            seq_file = ssbio.sequence.fasta.write_fasta_file(seq_str=s, ident=g, outdir=gene_folder)
            gene.annotation['sequence']['representative']['seq_file'] = op.basename(seq_file)

            log.info('{}: Loaded manually defined sequence information'.format(g))

    def set_representative_sequence(self):
        """Combine information from KEGG, UniProt, and manual mappings. Saves a DataFrame of results.

        Manual mappings override all existing mappings.
        UniProt mappings override KEGG mappings except when KEGG mappings have PDBs associated with them and UniProt doesn't.

        Returns:
            Pandas DataFrame: dataframe of ID mapping

        """
        # TODO: clean up this code!
        seq_mapping_pre_df = []

        for gene in self.genes:
            g = gene.id

            genedict = {'gene': g,
                        'pdbs': []}

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

            # If there is a KEGG annotation and no UniProt annotations, set KEGG as representative
            elif seq_prop['kegg']['seq_len'] > 0 and len(seq_prop['uniprot']) == 0:
                kegg_prop = seq_prop['kegg']
                seq_prop['representative'].update(kegg_prop)
                genedict.update(kegg_prop)
                log.debug('{}: Representative sequence set from KEGG'.format(g))

            # If there are UniProt annotations and no KEGG annotations, set UniProt as representative
            elif seq_prop['kegg']['seq_len'] == 0 and len(seq_prop['uniprot']) > 0:

                # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
                # This way, UniProts with PDBs get ranked to the top, or if no PDBs, reviewed entries
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
                log.debug('{}: Representative sequence set from UniProt using {}'.format(g, best_u))

            # If there are both UniProt and KEGG annotations...
            elif seq_prop['kegg']['seq_len'] > 0 and len(seq_prop['uniprot']) > 0:
                # Use KEGG if the mapped UniProt is unique, and it has PDBs
                kegg_prop = seq_prop['kegg']
                if len(kegg_prop['pdbs']) > 0 and kegg_prop['uniprot_acc'] not in seq_prop['uniprot'].keys():
                    seq_prop['representative'].update(kegg_prop)
                    genedict.update(kegg_prop)
                    log.debug('{}: Representative sequence set from KEGG'.format(g))
                else:
                    # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
                    uniprots = seq_prop['uniprot'].keys()
                    u_ranker = []
                    for u in uniprots:
                        u_ranker.append((u, seq_prop['uniprot'][u]['reviewed'] + len(seq_prop['uniprot'][u]['pdbs'])))
                    sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
                    best_u = sorted_by_second[0][0]

                    uni_prop = seq_prop['uniprot'][best_u]
                    your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
                    for_saving = {your_key: uni_prop[your_key] for your_key in your_keys if your_key in uni_prop}
                    seq_prop['representative'].update(for_saving)
                    genedict.update(for_saving)

                    # For saving in dataframe, save as string
                    if 'kegg_id' in genedict:
                        genedict['kegg_id'] = ';'.join(genedict['kegg_id'])

                    log.debug('{}: Representative sequence set from UniProt using {}'.format(g, best_u))

            # For saving in dataframe, save as string
            if genedict['pdbs']:
                genedict['pdbs'] = ';'.join(genedict['pdbs'])
            else:
                genedict['pdbs'] = None
            seq_mapping_pre_df.append(genedict)

        cols = ['gene', 'uniprot_acc', 'kegg_id', 'pdbs', 'seq_len', 'seq_file', 'metadata_file']
        tmp = pd.DataFrame.from_records(seq_mapping_pre_df, columns=cols)

        # Info on genes that could not be mapped
        self.missing_mapping = tmp[pd.isnull(tmp.seq_file)].gene.unique().tolist()
        if len(self.missing_mapping) > 0:
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_mapping" attribute.'.format(
                    len(self.missing_mapping)))

        self.df_sequence_mapping = tmp[pd.notnull(tmp.seq_file)].reset_index(drop=True)
        self.df_sequence_mapping.fillna(value=np.nan, inplace=True)
        mapping_df_outfile = op.join(self.data_dir, 'df_sequence_mapping.csv')
        self.df_sequence_mapping.to_csv(mapping_df_outfile)
        log.info('{}: Saved ID mapping information. Inspect the df_sequence_mapping attribute for more info.'.format(mapping_df_outfile))

    def map_uniprot_to_pdb(self, seq_ident_cutoff=0, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Returns:
            Pandas DataFrame: PDB mapping dataframe

        """
        best_structures_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            uniprot_id = g.annotation['sequence']['representative']['uniprot_acc']

            if uniprot_id:
                best_structures = ssbio.databases.pdb.best_structures(uniprot_id,
                                                                      outfile='{}_best_structures.json'.format(uniprot_id),
                                                                      outdir=op.join(self.sequence_dir, gene_id),
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

                    # TODO: gene is a key, not really needed in gene annotation
                    g.annotation['structure']['ranked_pdbs'] = ranked_pdbs
                    log.debug('{}: Loaded PDB ranking'.format(gene_id))
                    if len(best_structures) != len(ranked_pdbs):
                        log.info('{}: {} PDBs filtered out based on sequence identity cutoff'.format(gene_id, len(best_structures)-len(ranked_pdbs)))
                    log.info('{}: {} PDBs mapped'.format(gene_id, len(ranked_pdbs)))
                else:
                    log.info('{}: No PDBs mapped'.format(gene_id))

        cols = ['gene', 'uniprot_acc', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'pdb_resolution',
                'pdb_start', 'pdb_end', 'seq_coverage', 'unp_start', 'unp_end', 'tax_id', 'rank']
        self.df_pdb_ranking = pd.DataFrame.from_records(best_structures_pre_df, columns=cols)

        log.info('Completed UniProt -> best PDB mapping')

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, all_genes=False, force_rerun=False, display_link=False):
        """BLAST each gene sequence to the PDB. Raw BLAST results (XML files) are saved per gene in the "structures" folder.

        Returns:

        """

        for g in tqdm(self.genes):
            gene_id = str(str(g.id))
            seq_file = g.annotation['sequence']['representative']['seq_file']

            # Check if a representative sequence was set
            if not seq_file:
                log.warning('{}: No sequence set'.format(gene_id))
                continue

            seq_dir, seq_name, seq_ext = utils.split_folder_and_path(seq_file)

            # If all_genes=False, BLAST only genes without a uniprot->pdb mapping
            already_ranked_pdbs = g.annotation['structure']['ranked_pdbs']
            if already_ranked_pdbs and not all_genes:
                log.debug('Skipping BLAST for {}, structures already mapped and all_genes flag is False'.format(gene_id))
                continue

            # Read the sequence
            seq_file_path = op.join(self.sequence_dir, gene_id, seq_file)
            seq_record = SeqIO.read(open(seq_file_path), "fasta")
            seq_str = str(seq_record.seq)

            # Make the gene specific folder under the structure_files directory
            gene_folder = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # BLAST the sequence to the PDB
            blast_results = ssbio.databases.pdb.blast_pdb(seq_str,
                                                          outfile='{}_blast_pdb.xml'.format(seq_name),
                                                          outdir=op.join(self.sequence_dir, gene_id),
                                                          force_rerun=force_rerun,
                                                          evalue=evalue,
                                                          seq_ident_cutoff=seq_ident_cutoff,
                                                          link=display_link)

            if not blast_results:
                log.debug('No BLAST results for {}'.format(gene_id))
            else:
                log.debug('{}: {} PDBs BLASTed'.format(gene_id, len(blast_results)))

            # Save blast hits in gene annotation
            g.annotation['structure']['blast_pdbs'] = blast_results

    def set_representative_structure(self):
        pass

    """
    now we may have 2 things : 'blast_pdbs' and 'ranked_pdbs'
    blast_pdbs has a list of stuff, look for 'hit_pdb'
    ranked_pdbs has it by 'pdb_id'

    Take these two lists and make a big list of

    potential_structures = list of tuples (pdb_id, chain, coverage?, resolution?)

    """

    def pdb_downloader_and_metadata(self, all_pdbs=False, force_rerun=False):
        """Download structures which have been mapped to our genes. Gets both PDB and MMCIF files.

        Args:
            all_pdbs (bool): Default False, if True, all PDBs in the ranked and blasted pdb fields are downloaded.
                If False, only the representative PDB is downloaded.

        """
        pdb_pre_df = []

        for g in tqdm(self.genes):
            gene_id = g.id

            # Make the gene directory for structures
            gene_struct_dir = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(gene_struct_dir):
                os.mkdir(gene_struct_dir)

            if all_pdbs:
                # Check if we have any PDBs
                if not g.annotation['structure']['blast_pdbs'] and not g.annotation['structure']['ranked_pdbs']:
                    log.debug('{}: No structures available - no structures will be downloaded'.format(gene_id))
                    continue

                # TODO: keep chain information?
                # Get list of BLASTed PDBs
                blasted_pdbs = [x['hit_pdb'].lower() for x in g.annotation['structure']['blast_pdbs']]

                # Get list of PDBs from best_structures
                all_pdbs = [c['pdb_id'].lower() for c in g.annotation['structure']['ranked_pdbs']]
                all_pdbs.extend(blasted_pdbs)

                # Make sure theoretical or obsolete pdbs are filtered out (obsolete is replaced)
                pdbs_to_download = ssbio.databases.pdb.update_pdb_list(all_pdbs)
                log.debug('{}: Gene is mapped to these PDBs - {}'.format(gene_id, pdbs_to_download))

            else:
                # Check if we have a representative structure
                if not g.annotation['structure']['representative']['pdb_id']:
                    log.debug('{}: No representative structure available - no structure will be downloaded'.format(gene_id))
                    continue

                # Download representative structure only!
                pdbs_to_download = [g.annotation['structure']['representative']['pdb_id'].lower()]

            # Download the PDBs
            for p in pdbs_to_download:
                info_dict = {}

                log.debug('{}: Downloading PDB and mmCIF'.format(p))
                pdb_file = ssbio.databases.pdb.download_structure(pdb_id=p, file_type='pdb', outdir=gene_struct_dir, force_rerun=force_rerun)
                cif_file = ssbio.databases.pdb.download_structure(pdb_id=p, file_type='cif', outdir=gene_struct_dir, force_rerun=force_rerun)

                # Parse the mmCIF header
                cif_dict = ssbio.databases.pdb.parse_mmcif_header(cif_file)

                info_dict.update(cif_dict)

                # Save annotation info
                info_dict['pdb_file'] = op.basename(pdb_file)
                info_dict['mmcif_file'] = op.basename(cif_file)
                g.annotation['structure']['pdb'][p] = info_dict

                info_dict['gene'] = gene_id
                pdb_pre_df.append(info_dict)

        # Save a dataframe of the UniProt metadata
        if hasattr(self, 'df_pdb_metadata'):
            self.df_pdb_metadata = self.df_pdb_metadata.append(pdb_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing PDB dataframe.')
        else:
            cols = ['gene', 'organism', 'experiment', 'resolution', 'chemicals', 'pdb_file', 'mmcif_file']
            self.df_pdb_metadata = pd.DataFrame.from_records(pdb_pre_df, columns=cols)
            log.info('Created PDB metadata dataframe.')

    def pdb_ranking(self):

        for g in tqdm(self.genes):
            gene_id = g.id

            # Get list of BLASTed PDBs and their chains
            blasted_pdbs_c = [(x['hit_pdb'].lower(), x['hit_pdb_chains']) for x in g.annotation['structure']['blast_pdbs']]
            # Split the chain up
            blasted_pdbs = []
            for b in blasted_pdbs_c:
                for c in b[1]:
                    blasted_pdbs.append((b[0], c))

            # Get list of PDBs from best_structures
            ranked_pdbs = [(c['pdb_id'].lower(), c['pdb_chain_id'].upper()) for c in
                           g.annotation['structure']['ranked_pdbs']]

            # Preserve the rankings from best_structures
            total_pdbs = list(unique_everseen(ranked_pdbs))

            # Add additional BLASTed PDBs to the list
            for x in blasted_pdbs:
                if x not in total_pdbs:
                    total_pdbs.append(x)












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
                dest_gene_dir = op.join(self.structure_single_chain_dir, g)
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
                dest_coach_dir = op.join(self.structure_single_chain_dir, g, 'coach')
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

        log.info('{} homology models completed and copied to {}.'.format(mode11_counter, self.structure_single_chain_dir))
        if model1_bsites_counter:
            log.info('{} binding site predictions completed and copied.'.format(model1_bsites_counter))

        self.homology_df = pd.DataFrame(hom_pre_df)[['m_gene', 'i_model']]
        self.homology_df.to_csv(op.join(self.data_dir, '{}-homology_models_df.csv'.format(date.short_date)))
        return op.join(self.data_dir, '{}-homology_models_df.csv'.format(date.short_date))

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
