import os
import copy
import os.path as op
import numpy as np
import shutil
import pandas as pd
try:
    from IPython.display import clear_output
    have_ipython = True
    from tqdm import tqdm_notebook as tqdm
except ImportError:
    have_ipython = False
    from tqdm import tqdm

from collections import OrderedDict

from bioservices.uniprot import UniProt
from bioservices import KEGG

from Bio import SeqIO
import requests
from ssbio import utils
import ssbio.cobra.utils
import ssbio.databases.kegg
import ssbio.databases.uniprot
import ssbio.databases.pdb
import ssbio.sequence.fasta
from ssbio.itasser.itasserparse import ITASSERParse
import ssbio.structure.properties.residues
import ssbio.structure.properties.quality
from ssbio.structure.cleanpdb import CleanPDB
from ssbio.structure.pdbioext import PDBIOExt
import ssbio.sequence.properties.residues
import ssbio.structure.properties.msms

from cobra.core import Gene
from cobra.core import DictList

import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)
log = logging.getLogger(__name__)

date = utils.Date()
bs_unip = UniProt()
bs_kegg = KEGG()


class RepresentativeSequence(object):
    pass


class RepresentativeStructure(object):
    pass


# TODO using these classes in the annotation field will work
class StructureProp(object):
    def __init__(self, homology=None, pdb=None, representative=None):
        if not homology:
            homology = {}
        if not pdb:
            pdb = OrderedDict()
        if not representative:
            representative = {'pdb_id'             : None,
                              'resolution'         : float('inf'),
                              'original_pdb_file'  : None,
                              'original_mmcif_file': None,
                              'clean_pdb_file'     : None,
                              'properties': {}}
        self.homology = homology
        self.pdb = pdb
        self.representative = representative


class SequenceProp(object):
    def __init__(self, kegg_id=None, kegg_sequence_file=None, kegg_metadata_file=None,
                 uniprot_ids=None, representative=None):

        self.kegg = KEGGProp(kegg_id, kegg_sequence_file, kegg_metadata_file)

        if not uniprot_ids:
            uniprot_ids = []
        uniprots = {}
        for x in uniprot_ids:
            uniprots[x] = UniProtProp(uniprot_acc=x)
        self.uniprot = uniprots


class SeqProp(object):
    def __init__(self, ident=None, sequence_file=None, metadata_file=None):
        if not ident:
            ident = ''
        if not sequence_file:
            sequence_file = ''
        if not metadata_file:
            metadata_file = ''

        self.ident = ident
        self.kegg = None
        self.refseq = None
        self.uniprot = None
        self.pdbs = []
        self.seq_len = 0
        self.seq_file = op.basename(sequence_file)
        self.metadata_file = op.basename(metadata_file)

        if sequence_file:
            self.seq_len = len(SeqIO.read(open(sequence_file), "fasta"))


class UniProtProp(SeqProp):
    def __init__(self, uniprot_acc=None, sequence_file=None, metadata_file=None):

        SeqProp.__init__(self, sequence_file=sequence_file, metadata_file=metadata_file)

        if uniprot_acc:
            self.uniprot = uniprot_acc

        if metadata_file:
            metadata = ssbio.databases.uniprot.parse_uniprot_txt_file(metadata_file)


class KEGGProp(SeqProp):
    def __init__(self, kegg_id=None, sequence_file=None, metadata_file=None):

        SeqProp.__init__(self, sequence_file=sequence_file, metadata_file=metadata_file)

        if kegg_id:
            self.kegg = kegg_id

    def get_dict(self):
        return self.__dict__.copy()


class GEMPRO(object):
    """Generic class to represent all information of a GEM-PRO for a GEM.

    Main steps are:
    1. Mapping of sequence IDs
        a. With KEGG mapper
        b. With UniProt mapper
        c. Using both
        d. Allowing manual gene ID --> protein sequence entry
        e. Allowing manual gene ID --> UniProt ID
    2. Consolidating sequence IDs and setting a representative sequence
    3. Mapping of representative sequence --> structures
        a. With UniProt --> ranking of PDB structures
        b. BLAST representative sequence --> PDB database
    3. QC/QA
    4. Homology modeling
    5. Representative structures

    Each step may generate a report and also allow input of manual mappings if applicable.
    """

    def __init__(self, gem_name, root_dir, gem_file_path=None, gem_file_type=None, genes_list=None, genes_and_sequences=None):
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

        # Or, load a dictionary of genes and their sequences
        elif genes_and_sequences:
            self.genes = list(genes_and_sequences.keys())
            self.manual_seq_mapping(genes_and_sequences)

        # If neither a model or genes are input, you can still add_genes_by_id later
        else:
            self.genes = []
            log.warning('No model or genes input')

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
                x = str(x)
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
                                                                          'properties'   : {},
                                                                          'seq_file'     : None,
                                                                          'metadata_file': None}}
                if 'structure' not in new_gene.annotation.keys():
                    new_gene.annotation['structure'] = {'homology'      : {},
                                                        'pdb'           : OrderedDict(),
                                                        'representative': {'structure_id'     : None,
                                                                           'seq_coverage'     : 0,
                                                                           'original_pdb_file': None,
                                                                           'clean_pdb_file'   : None,
                                                                           'properties': {}}}
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
                                                                   'properties'   : {},
                                                                   'seq_file'     : None,
                                                                   'metadata_file': None}}
                if 'structure' not in x.annotation.keys():
                    x.annotation['structure'] = {'homology'      : {},
                                                 'pdb'           : OrderedDict(),
                                                 'representative': {'structure_id'     : None,
                                                                    'seq_coverage'     : 0,
                                                                    'original_pdb_file': None,
                                                                    'clean_pdb_file'   : None,
                                                                    'properties': {}}}
            self._genes = genes_list

    def add_genes_by_id(self, genes_list):
        """Add gene IDs manually into our GEM-PRO project.

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
                                                                  'properties'   : {},
                                                                  'seq_file'     : None,
                                                                  'metadata_file': None}}
            new_gene.annotation['structure'] = {'homology'      : {},
                                                'pdb'           : OrderedDict(),
                                                'representative': {'structure_id'     : None,
                                                                   'seq_coverage'     : 0,
                                                                   'original_pdb_file': None,
                                                                   'clean_pdb_file'   : None,
                                                                   'properties': {}}}
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

            # Download kegg metadata
            metadata_file = ssbio.databases.kegg.download_kegg_gene_metadata(organism_code=kegg_organism_code,
                                                                             gene_id=kegg_g,
                                                                             outdir=gene_folder,
                                                                             force_rerun=force_rerun)

            # Download kegg sequence
            sequence_file = ssbio.databases.kegg.download_kegg_aa_seq(organism_code=kegg_organism_code,
                                                                      gene_id=kegg_g,
                                                                      outdir=gene_folder,
                                                                      force_rerun=force_rerun)
            kegg_prop = KEGGProp(kegg_id=kegg_organism_code + ':' + kegg_g,
                                 sequence_file=sequence_file,
                                 metadata_file=metadata_file)
            kegg_dict = kegg_prop.get_dict()

            # Update potentially old UniProt ID
            if kegg_g in kegg_to_uniprot.keys():
                kegg_prop.uniprot_acc = kegg_to_uniprot[kegg_g]

            # Save in Gene
            g.annotation['sequence']['kegg'] = kegg_dict

            if sequence_file:
                # Also check if KEGG sequence matches a potentially set representative sequence
                # Do not add any info if a UniProt ID was already mapped though, we want to use that
                if g.annotation['sequence']['representative']['seq_len'] > 0 and not g.annotation['sequence']['representative']['uniprot_acc']:
                    ### Check if sequences are the same

                    # Load already set representative sequence
                    rep_seq_file = g.annotation['sequence']['representative']['seq_file']
                    rep_seq_path = op.join(self.sequence_dir, gene_id, rep_seq_file)
                    rep_seq = SeqIO.read(open(rep_seq_path), 'fasta')

                    # Load kegg sequence
                    kegg_seq = SeqIO.read(open(sequence_file), 'fasta')

                    # Test equality
                    if str(rep_seq.seq) == str(kegg_seq.seq):
                        # If equal, fill in representative sequence fields with kegg metadata but still use
                        # the already set sequence file
                        your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'metadata_file']
                        for_saving = {your_key: kegg_dict[your_key] for your_key in your_keys if
                                      your_key in kegg_dict}
                        g.annotation['sequence']['representative'].update(for_saving)
                    else:
                        log.warning('{}: manually added sequence does not match mapped KEGG sequence.'.format(gene_id))

            # Save in dataframe
            kegg_dict['gene'] = gene_id
            kegg_pre_df.append(kegg_dict)

            log.debug('{}: Loaded KEGG information for gene'.format(gene_id))

        # Save a dataframe of the file mapping info
        cols = ['gene', 'uniprot_acc', 'kegg_id', 'seq_len', 'pdbs', 'seq_file', 'metadata_file']
        self.df_kegg_metadata = pd.DataFrame.from_records(kegg_pre_df, columns=cols)
        log.info('Created KEGG metadata dataframe. See the "df_kegg_metadata" attribute.')

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

        """

        # Allow model gene --> custom ID mapping ({'TM_1012':'TM1012'})
        if custom_gene_mapping:
            genes_to_map = list(custom_gene_mapping.values())
        else:
            genes_to_map = [x.id for x in self.genes]

        # Map all IDs first to available UniProts
        genes_to_uniprots = bs_unip.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)

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

                    uniprot_dict['uniprot_acc'] = str(mapped_uniprot)

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
                    # TODO: fix default representation for pdbs
                    if 'pdbs' not in uniprot_dict:
                        uniprot_dict['pdbs'] = []
                    g.annotation['sequence']['uniprot'][str(mapped_uniprot)] = uniprot_dict

                    # Also check if UniProt sequence matches a potentially set representative sequence
                    if g.annotation['sequence']['representative']['seq_len'] > 0:
                        ### Check if sequences are the same

                        # Load already set representative sequence
                        rep_seq_file = g.annotation['sequence']['representative']['seq_file']
                        rep_seq_path = op.join(self.sequence_dir, gene_id, rep_seq_file)
                        rep_seq = SeqIO.read(open(rep_seq_path), 'fasta')

                        # Load uniprot sequence
                        unip_seq = SeqIO.read(open(sequence_file), 'fasta')

                        # Test equality
                        if str(rep_seq.seq) == str(unip_seq.seq):
                            # If equal, fill in representative sequence fields with uniprot metadata but still use
                            # the already set sequence file
                            your_keys = ['kegg_id', 'uniprot_acc', 'pdbs', 'metadata_file']
                            for_saving = {your_key: uniprot_dict[your_key] for your_key in your_keys if
                                          your_key in uniprot_dict}
                            g.annotation['sequence']['representative'].update(for_saving)
                        else:
                            log.warning('{}: manually added sequence does not match mapped UniProt sequence.'.format(gene_id))

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
            log.info('Created UniProt metadata dataframe. See the "df_uniprot_metadata" attribute.')

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
            for_saving = {your_key: uniprot_dict[your_key] for your_key in your_keys if your_key in uniprot_dict}
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

    # TODO: should also have a seq --> uniprot id function (has to be 100% match) (what about organism?)
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
            g = str(g)
            gene = self.genes.get_by_id(g)

            gene.annotation['sequence']['representative']['seq_len'] = len(s)

            gene_folder = op.join(self.sequence_dir, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            seq_file = ssbio.sequence.fasta.write_fasta_file(seq_str=s, ident=g, outdir=gene_folder)
            gene.annotation['sequence']['representative']['seq_file'] = op.basename(seq_file)

            log.debug('{}: Loaded manually defined sequence information'.format(g))

        log.info('Loaded in {} sequences'.format(len(gene_to_seq_dict)))

    def set_representative_sequence(self):
        """Combine information from KEGG, UniProt, and manual mappings. Saves a DataFrame of results.

        Manual mappings override all existing mappings.
        UniProt mappings override KEGG mappings except when KEGG mappings have PDBs associated with them and UniProt doesn't.

        """
        # TODO: clean up this code!
        # TODO: check erol's notebook to see why uniprot_accs are not being set
        seq_mapping_pre_df = []

        for gene in self.genes:
            g = gene.id

            genedict = {'pdbs': []}

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
                for_saving = { your_key: uni_prop[your_key] for your_key in your_keys if your_key in uni_prop}
                seq_prop['representative'].update(for_saving)
                genedict.update(for_saving)
                if 'kegg_id' in genedict:
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
            genedict['gene'] = g
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
        # mapping_df_outfile = op.join(self.data_dir, 'df_sequence_mapping.csv')
        # self.df_sequence_mapping.to_csv(mapping_df_outfile)
        log.info('Created sequence mapping dataframe. See the "df_sequence_mapping" attribute.')

    def map_uniprot_to_pdb(self, seq_ident_cutoff=0, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Creates a summary dataframe accessible by the attribute "df_pdb_ranking".

        """
        best_structures_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            uniprot_id = g.annotation['sequence']['representative']['uniprot_acc']

            if not uniprot_id:
                # Check if a representative sequence was set
                log.warning('{}: No representative UniProt ID set, cannot use best structures API'.format(gene_id))
                continue
            else:
                best_structures = ssbio.databases.pdb.best_structures(uniprot_id,
                                                                      outfile='{}_best_structures.json'.format(uniprot_id),
                                                                      outdir=op.join(self.sequence_dir, gene_id),
                                                                      seq_ident_cutoff=seq_ident_cutoff,
                                                                      force_rerun=force_rerun)

                if best_structures:
                    rank = 1
                    to_add_to_annotation = OrderedDict()

                    for best_structure in best_structures:
                        currpdb = str(best_structure['pdb_id'].lower())
                        currchain = str(best_structure['chain_id'])

                        # TODO: fix caching
                        # pdb_rel = ssbio.databases.pdb.get_release_date(currpdb)

                        best_structure_dict = {}
                        best_structure_dict['pdb_id'] = currpdb
                        best_structure_dict['pdb_chain_id'] = currchain
                        best_structure_dict['uniprot_acc'] = uniprot_id
                        best_structure_dict['experimental_method'] = best_structure['experimental_method']
                        best_structure_dict['resolution'] = best_structure['resolution']
                        best_structure_dict['seq_coverage'] = best_structure['coverage']
                        # best_structure_dict['release_date'] = pdb_rel
                        best_structure_dict['taxonomy_id'] = best_structure['tax_id']
                        best_structure_dict['pdb_start'] = best_structure['start']
                        best_structure_dict['pdb_end'] = best_structure['end']
                        best_structure_dict['unp_start'] = best_structure['unp_start']
                        best_structure_dict['unp_end'] = best_structure['unp_end']
                        best_structure_dict['rank'] = rank

                        # For saving in the Gene annotation
                        to_add_to_annotation[str(currpdb) + ';' + str(currchain)] = best_structure_dict.copy()

                        # For saving in the summary dataframe
                        best_structure_dict['gene'] = gene_id
                        best_structures_pre_df.append(best_structure_dict)

                        rank += 1

                    # If structure annotation exists already, remove existing
                    # (pdb,chain) keys and use the best_structure annotation instead
                    # NOTE: sometimes, the (pdb,chain) key is not unique - some chains have a fused protein with
                    # perhaps another organism's protein fused. Currently not considering these as unique, but
                    # TODO: when cleaning, these should be removed
                    if g.annotation['structure']['pdb']:
                        current_annotation = g.annotation['structure']['pdb']
                        temp_annotation = OrderedDict([k, v] for k, v in current_annotation.items() if k not in to_add_to_annotation)
                        to_add_to_annotation.update(temp_annotation)

                    g.annotation['structure']['pdb'] = to_add_to_annotation

                    log.debug('{}: {} PDB/chain pairs mapped'.format(gene_id, to_add_to_annotation))
                else:
                    log.debug('{}: No PDB/chain pairs mapped'.format(gene_id))

        cols = ['gene', 'uniprot_acc', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'resolution', 'seq_coverage',
                'taxonomy_id', 'pdb_start', 'pdb_end', 'unp_start', 'unp_end', 'rank'] #'release_date',
        self.df_pdb_ranking = pd.DataFrame.from_records(best_structures_pre_df, columns=cols)

        # TODO: also report genes with no PDB?

        log.info('Completed UniProt -> best PDB mapping. See the "df_pdb_ranking" attribute.')

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, all_genes=False, force_rerun=False, display_link=False):
        """BLAST each gene sequence to the PDB. Raw BLAST results (XML files) are saved per gene in the "structures" folder.

        Returns:

        """
        blast_results_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)
            seq_file = g.annotation['sequence']['representative']['seq_file']

            # Check if a representative sequence was set
            if not seq_file:
                log.warning('{}: No sequence set'.format(gene_id))
                continue

            seq_dir, seq_name, seq_ext = utils.split_folder_and_path(seq_file)

            # If all_genes=False, BLAST only genes without a uniprot->pdb mapping
            already_has_pdbs = g.annotation['structure']['pdb']
            if already_has_pdbs and not all_genes:
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

            if blast_results:
                to_add_to_annotation = OrderedDict()

                for blast_result in blast_results:
                    pdb = str(blast_result['hit_pdb'].lower())
                    chains = blast_result['hit_pdb_chains']

                    # pdb_rez = ssbio.databases.pdb.get_resolution(pdb)
                    # pdb_rel = ssbio.databases.pdb.get_release_date(pdb)

                    for chain in chains:
                        chain = str(chain)
                        blast_dict = {}
                        blast_dict['pdb_id'] = pdb
                        blast_dict['pdb_chain_id'] = chain
                        # blast_dict['resolution'] = pdb_rez
                        # blast_dict['release_date'] = pdb_rel
                        blast_dict['blast_score'] = blast_result['hit_score']
                        blast_dict['blast_evalue'] = blast_result['hit_evalue']
                        blast_dict['seq_coverage'] = blast_result['hit_percent_ident']
                        blast_dict['seq_similar'] = blast_result['hit_percent_similar']
                        blast_dict['seq_num_coverage'] = blast_result['hit_num_ident']
                        blast_dict['seq_num_similar'] = blast_result['hit_num_similar']

                        # For saving in Gene annotation
                        to_add_to_annotation[pdb + ';' + chain] = blast_dict.copy()

                        # For saving in summary dataframe
                        blast_dict['gene'] = gene_id

                        blast_results_pre_df.append(blast_dict)

                # If structure annotation exists already, remove existing
                # (pdb,chain) keys from BLAST results and append rest to the end
                if g.annotation['structure']['pdb']:
                    to_add_to_annotation = OrderedDict([k, v] for k, v in to_add_to_annotation.items() if k not in g.annotation['structure']['pdb'])
                if to_add_to_annotation:
                    log.info('{}: Adding {} PDBs from BLAST results.'.format(gene_id, len(to_add_to_annotation)))
                g.annotation['structure']['pdb'].update(to_add_to_annotation)

                log.debug('{}: {} PDBs BLASTed'.format(gene_id, len(blast_results)))
            else:
                log.debug('No BLAST results for {}'.format(gene_id))

        cols = ['gene', 'pdb_id', 'pdb_chain_id', 'blast_score', 'blast_evalue', # 'resolution', 'release_date',
                'seq_coverage', 'seq_similar', 'seq_num_coverage', 'seq_num_similar']
        self.df_pdb_blast = pd.DataFrame.from_records(blast_results_pre_df, columns=cols)

        log.info('Completed sequence --> PDB BLAST. See the "df_pdb_blast" attribute.')
        # TODO: log.info for counts - num pdbs with no blast hits, number with (instead of in the for loop)

    def manual_homology_models(self, input_dict):
        """Copy homology models and manually defined information per model to the GEM-PRO project.

        Args:
            input_dict: Dictionary of dictionaries of gene names to homology model IDs and information. Input a dict of:
                {model_gene: {homology_model_id1: {'model_file': '/path/to/homology/model',
                                                  'other_info': 'other_info_here',
                                                  ...},
                              homology_model_id2: {'model_file': '/path/to/homology/model',
                                                  'other_info': 'other_info_here',
                                                  ...}}}

        """
        counter = 0
        for g in tqdm(self.genes):
            gene_id = g.id

            if gene_id not in input_dict:
                continue

            for hid, hdict in input_dict[gene_id].items():
                if 'model_file' not in hdict:
                    raise KeyError('"model_file" must be a key in the manual input dictionary.')

                # Make the destination structure folder
                dest_gene_dir = op.join(self.structure_single_chain_dir, gene_id)
                if not op.exists(dest_gene_dir):
                    os.mkdir(dest_gene_dir)

                # Just copy the file to the structure directory and store the file name
                if not op.exists(op.join(dest_gene_dir, op.basename(hdict['model_file']))):
                    shutil.copy2(hdict['model_file'], dest_gene_dir)

                g.annotation['structure']['homology'][hid] = hdict
                log.debug('{}: updated homology model information and copied model file.'.format(gene_id))
            counter += 1

        log.info('Updated homology model information for {} genes.'.format(counter))

    # TODO: code to prep homology modeling
    def prep_itasser_models(self, outdir, itasser_installation, itlib_location, runtype, print_exec=False, **kwargs):
        pass
        # import ssbio.itasser.itasserprep as ip
        # # if running on ssb1-3
        # runtype = 'torque'
        # itasser_loc = '~/software/I-TASSER4.4'
        # itlib_loc = '~/software/ITLIB'
        # local_dir = '/home/nathan/projects_unsynced/thermatoga_me/structure_files/homology_models'
        # run_from_dir = '~/projects/thermatoga_me/structure_files/homology_models'
        #
        # for ident, seq in genes_and_sequences:
        #     ip.ITASSERPrep(ident, seq, local_dir,
        #                    itasser_loc=itasser_loc, itlib_loc=itlib_loc,
        #                    runtype=runtype, print_exec=True, data_dir=run_from_dir)

        # # if running on NERSC
        # runtype = 'slurm'
        # itasser_loc = '~/projects/software/I-TASSER4.4'
        # itlib_loc = '~/projects/software/ITLIB'
        # email = 'nmih@ucsd.edu'
        # username = 'nmih'
        # project = 'm1244'
        # data_dir = '~/projects/thermatoga_me/structure_files/homology_models'
        #
        # for ident, seq in genes_and_sequences:
        #     ip.ITASSERPrep(ident, seq, '/home/nathan/projects_unsynced/thermatoga_me/structure_files/homology_models',
        #                    itasser_loc=itasser_loc, itlib_loc=itlib_loc,
        #                    runtype=runtype, print_exec=True, data_dir=data_dir, project_id=project,
        #                    slurm_email=email, slurm_username=username)

    def get_itasser_models(self, homology_raw_dir, custom_itasser_name_mapping=None, force_rerun=False):
        """Copy generated homology models from a directory to the GEM-PRO directory.

        Args:
            homology_raw_dir: Root directory of I-TASSER folders.
            custom_itasser_name_mapping: Use this if your I-TASSER folder names differ from your model gene names.
                Input a dict of {model_gene: ITASSER_folder}.

        """
        itasser_pre_df = []

        for g in tqdm(self.genes):
            gene_id = g.id
            itasser_name = gene_id + '_model1'

            # Make the destination structure folder and other results folder
            dest_gene_dir = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(dest_gene_dir):
                os.mkdir(dest_gene_dir)

            dest_itasser_extra_dir = op.join(dest_gene_dir, '{}_itasser'.format(itasser_name))
            if not op.exists(dest_itasser_extra_dir):
                os.mkdir(dest_itasser_extra_dir)

            if custom_itasser_name_mapping and gene_id in custom_itasser_name_mapping:
                hom_id = custom_itasser_name_mapping[gene_id]
            else:
                hom_id = gene_id
            orig_itasser_dir = op.join(homology_raw_dir, hom_id)

            itasser_parse = ITASSERParse(original_results_path=orig_itasser_dir, create_dfs=True)
            itasser_parse.copy_results(copy_to_dir=dest_gene_dir, rename_model_to=itasser_name, force_rerun=force_rerun)
            itasser_parse.save_dataframes(outdir=dest_itasser_extra_dir)

            if itasser_parse.modeling_results:
                # Always set sequence coverage to 100% for an ITASSER model
                # TODO: should probably just do the alignment to make sure (representative sequences can change)
                copied = copy.deepcopy(itasser_parse.modeling_results)
                copied['seq_coverage'] = 1
                copied['model_file'] = op.basename(itasser_parse.structure_path)
                g.annotation['structure']['homology'][gene_id] = copied

                copied['gene'] = gene_id
                itasser_pre_df.append(copied)
            else:
                log.debug('{}: No homology model available.'.format(gene_id))

        cols = ['gene', 'model_file', 'model_date', 'difficulty',
                'top_template_pdb', 'top_template_chain', 'c_score',
                'tm_score', 'tm_score_err', 'rmsd', 'rmsd_err',
                'top_bsite_site_num', 'top_bsite_c_score', 'top_bsite_cluster_size', 'top_bsite_algorithm',
                'top_bsite_pdb_template_id', 'top_bsite_pdb_template_chain', 'top_bsite_pdb_ligand',
                'top_bsite_binding_location_coords', 'top_bsite_c_score_method', 'top_bsite_binding_residues',
                'top_bsite_ligand_cluster_counts',
                'top_ec_pdb_template_id', 'top_ec_pdb_template_chain', 'top_ec_tm_score', 'top_ec_rmsd',
                'top_ec_seq_ident', 'top_ec_seq_coverage', 'top_ec_c_score', 'top_ec_ec_number',
                'top_ec_binding_residues',
                'top_go_mf_go_id', 'top_go_mf_go_term', 'top_go_mf_c_score', 'top_go_bp_go_id', 'top_go_bp_go_term',
                'top_go_bp_c_score', 'top_go_cc_go_id', 'top_go_cc_go_term', 'top_go_cc_c_score']
        self.df_itasser = pd.DataFrame.from_records(itasser_pre_df, columns=cols)

        log.info('Completed copying of {} I-TASSER models to GEM-PRO directory. See the "df_itasser" attribute.'.format(len(self.df_itasser)))

    def set_representative_structure(self, always_use_homology=True, sort_homology_by='seq_coverage',
                                     allow_missing_on_termini=0.1, allow_mutants=True, allow_deletions=False,
                                     allow_insertions=False, allow_unresolved=True, force_rerun=False):
        """Set the representative structure for a gene.

        Each gene can have a combination of the following:
        - Homology model(s)
        - Ranked PDBs
        - BLASTed PDBs

        If the always_use_homology flag is true, homology models are always set as representative when they exist.
            If there are multiple homology models, we rank by default by the seq_coverage key. Other parameters:
            - c_score
            - model_date
            - tm_score

        """

        for g in tqdm(self.genes):
            gene_id = str(g.id)

            has_homology = False
            has_pdb = False
            use_homology = False
            use_pdb = False

            if len(g.annotation['structure']['homology']) > 0:
                has_homology = True
            if len(g.annotation['structure']['pdb']) > 0:
                has_pdb = True

            # If there are no structures at all, move on
            if not has_pdb and not has_homology:
                log.debug('{}: No structures available - no representative structure will be set.'.format(gene_id))
                continue

            # If we mark to always use homology, use it if it exists
            if always_use_homology:
                if has_homology:
                    use_homology = True
                elif has_pdb:
                    use_pdb = True
            # If we don't always want to use homology, use PDB if it exists
            else:
                if has_homology and has_pdb:
                    use_pdb = True
                elif has_homology and not has_pdb:
                    use_homology = True
                elif has_pdb and not has_homology:
                    use_pdb = True

            structure_set = False
            gene_seq_dir = op.join(self.sequence_dir, gene_id)
            gene_struct_dir = op.join(self.structure_single_chain_dir, gene_id)

            # Make the destination structure folder
            if not op.exists(gene_struct_dir):
                os.mkdir(gene_struct_dir)

            if use_pdb:
                try:
                    # Get the representative sequence
                    seq_file = g.annotation['sequence']['representative']['seq_file']
                    seq_file_path = op.join(gene_seq_dir, seq_file)
                    seq_record = SeqIO.read(open(seq_file_path), "fasta")
                    ref_seq = str(seq_record.seq)

                    # Put PDBs through QC/QA
                    all_pdbs_and_chains = list(g.annotation['structure']['pdb'].keys())
                    convert_to_dict = utils.DefaultOrderedDict(list)
                    for x in all_pdbs_and_chains:
                        convert_to_dict[x.split(';')[0]].append(x.split(';')[1])

                    for pdb, chains in convert_to_dict.items():
                        # Download the PDB
                        try:
                            file_type = 'pdb'
                            pdb_file = ssbio.databases.pdb.download_structure(pdb_id=pdb, file_type=file_type, header=False,
                                                                              outdir=gene_struct_dir,
                                                                              force_rerun=force_rerun)

                            log.debug('{}: downloaded PDB file'.format(pdb))
                        except requests.exceptions.HTTPError:
                            file_type = 'cif'
                            pdb_file = ssbio.databases.pdb.download_structure(pdb_id=pdb, file_type=file_type, header=False,
                                                                              outdir=gene_struct_dir,
                                                                              force_rerun=force_rerun)
                            log.debug('{}: download mmCIF file'.format(pdb))
                        else:
                            log.warning('{}: PDB file not available'.format(pdb))
                            continue

                        # Get the sequences of the chains
                        chain_to_seq = ssbio.structure.properties.residues.get_pdb_seqs(pdb_file, file_type=file_type)

                        for chain in chains:
                            chain_seq = chain_to_seq[chain]

                            # Compare representative sequence to structure sequence
                            found_good_pdb = ssbio.structure.properties.quality.sequence_checker(reference_id=gene_id,
                                                                                                 reference_sequence=ref_seq,
                                                                                                 structure_id=pdb+'_'+chain,
                                                                                                 structure_sequence=chain_seq,
                                                                                                 allow_missing_on_termini=allow_missing_on_termini,
                                                                                                 allow_mutants=allow_mutants,
                                                                                                 allow_deletions=allow_deletions,
                                                                                                 allow_insertions=allow_insertions,
                                                                                                 allow_unresolved=allow_unresolved,
                                                                                                 outfile='{}_{}.needle'.format(gene_id, pdb+'_'+chain),
                                                                                                 outdir=gene_seq_dir,
                                                                                                 force_rerun=force_rerun)

                            # If found_good_pdb = True, set as representative
                            # If not, move on to the next potential PDB
                            if found_good_pdb:
                                orig_pdb_data = g.annotation['structure']['pdb'][str(pdb) + ';' + str(chain)]
                                g.annotation['structure']['representative']['structure_id'] = str(pdb) + ';' + str(chain)
                                g.annotation['structure']['representative']['seq_coverage'] = orig_pdb_data['seq_coverage']
                                g.annotation['structure']['representative']['original_pdb_file'] = op.basename(pdb_file)

                                # Clean it
                                custom_clean = CleanPDB(keep_chains=chain)
                                my_pdb = PDBIOExt(pdb_file, file_type=file_type)
                                default_cleaned_pdb = my_pdb.write_pdb(custom_selection=custom_clean,
                                                                       out_suffix='_{}_clean'.format(chain),
                                                                       out_dir=gene_struct_dir)
                                default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

                                g.annotation['structure']['representative']['clean_pdb_file'] = default_cleaned_pdb_basename

                                # TODO: download SIFTS file
                                # TODO: mutate it based on the alignment
                                # TODO: run xleap commands too?

                                structure_set = True
                                raise StopIteration
                except StopIteration:
                    log.debug('{}: Found representative PDB ({}, {})'.format(gene_id, pdb, chain))
                    continue
                else:
                    if has_homology:
                        use_homology = True

            # If we are to use homology, save its information in the representative structure field
            if use_homology and not structure_set:
                hm = g.annotation['structure']['homology']
                # Sort the available homology models by the specified field
                sorted_homology_ids = sorted(hm, key=lambda x: hm[x][sort_homology_by], reverse=True)

                top_homology = sorted_homology_ids[0]
                original_pdb_file = g.annotation['structure']['homology'][top_homology]['model_file']
                seq_coverage = g.annotation['structure']['homology'][top_homology]['seq_coverage']

                g.annotation['structure']['representative']['structure_id'] = top_homology
                g.annotation['structure']['representative']['seq_coverage'] = seq_coverage
                g.annotation['structure']['representative']['original_pdb_file'] = original_pdb_file

                # Clean it
                custom_clean = CleanPDB()
                my_pdb = PDBIOExt(op.join(gene_struct_dir, original_pdb_file), file_type='pdb')
                default_cleaned_pdb = my_pdb.write_pdb(custom_selection=custom_clean,
                                                       out_suffix='_clean',
                                                       out_dir=gene_struct_dir)
                default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

                g.annotation['structure']['representative']['clean_pdb_file'] = default_cleaned_pdb_basename

                structure_set = True
            else:
                log.debug('{}: No representative PDB'.format(gene_id))

        # TODO: summary of results

    def pdb_downloader_and_metadata(self, force_rerun=False):
        """Download ALL structures which have been mapped to our genes. Gets PDB file and mmCIF header and
            creates a metadata table.

        Args:
            force_rerun (bool): If you want to overwrite any existing mappings and files

        """
        pdb_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)

            # Make the gene directory for structures
            gene_struct_dir = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(gene_struct_dir):
                os.mkdir(gene_struct_dir)

            # Check if we have any PDBs
            if len(g.annotation['structure']['pdb']) == 0:
                log.debug('{}: No structures available - no structures will be downloaded'.format(gene_id))
                continue

            # Download the PDBs
            for k, v in g.annotation['structure']['pdb'].items():
                p = v['pdb_id']

                log.debug('{}: Downloading PDB and mmCIF'.format(p))
                pdb_file = ssbio.databases.pdb.download_structure(pdb_id=p, file_type='pdb', header=False, outdir=gene_struct_dir, force_rerun=force_rerun)
                cif_file = ssbio.databases.pdb.download_structure(pdb_id=p, file_type='cif', header=True, outdir=gene_struct_dir, force_rerun=force_rerun)

                # Parse the mmCIF header
                cif_dict = ssbio.databases.pdb.parse_mmcif_header(cif_file)

                # Save annotation info
                cif_dict['pdb_file'] = op.basename(pdb_file)
                cif_dict['mmcif_header'] = op.basename(cif_file)
                g.annotation['structure']['pdb'][k].update(cif_dict)

                adder = g.annotation['structure']['pdb'][k].copy()
                adder['chemicals'] = ';'.join(adder['chemicals'])
                # TODO: error when taxonomy name does not exist for some PDBs
                if isinstance(adder['taxonomy_name'], list):
                    adder['taxonomy_name'] = ';'.join(adder['taxonomy_name'])
                adder['gene'] = gene_id
                pdb_pre_df.append(adder)

        # Save a dataframe of the PDB metadata
        if hasattr(self, 'df_pdb_metadata'):
            self.df_pdb_metadata = self.df_pdb_metadata.append(pdb_pre_df, ignore_index=True).drop_duplicates().reset_index(drop=True)
            log.info('Updated existing PDB dataframe.')
        else:
            cols = ['gene', 'pdb_id', 'pdb_chain_id', 'taxonomy_name', 'experimental_method',
                    'resolution', 'seq_coverage', 'chemicals', 'rank', 'release_date', 'blast_evalue', 'blast_score',
                    'seq_similar', 'seq_num_coverage', 'seq_num_similar',
                    'pdb_file', 'mmcif_header']
            self.df_pdb_metadata = pd.DataFrame.from_records(pdb_pre_df, columns=cols).drop_duplicates().reset_index(drop=True)
            log.info('Created PDB metadata dataframe.')

    # TODO: prepare genes for homology modeling
    def homology_modeling_prep(self):
        pass

    def calculate_sequence_properties(self, force_rerun=False):
        """Run EMBOSS pepstats and biopython ProteinAnalysis on all representative sequences

        """
        seq_prop_pre_df = []
        for g in tqdm(self.genes):
            gene_id = str(g.id)
            gene_folder = op.join(self.sequence_dir, gene_id)
            repseq = g.annotation['sequence']['representative']['seq_file']

            if not repseq:
                log.warning('{}: No representative sequence available'.format(gene_id))
                continue
            else:
                repseq_path = op.join(gene_folder, repseq)
                pepstats_file = ssbio.sequence.properties.residues.emboss_pepstats_on_fasta(infile=repseq_path,
                                                                                            outdir=gene_folder,
                                                                                            force_rerun=force_rerun)
                pepstats_parsed = ssbio.sequence.properties.residues.emboss_pepstats_parser(pepstats_file)

                # TODO: also run biopython ProteinAnalysis
                g.annotation['sequence']['representative']['properties'].update(pepstats_parsed.copy())

                pepstats_parsed['gene'] = gene_id
                seq_prop_pre_df.append(pepstats_parsed)

        cols = ['gene', 'percent_acidic', 'percent_aliphatic', 'percent_aromatic', 'percent_basic', 'percent_charged',
                'percent_non-polar', 'percent_polar', 'percent_small', 'percent_tiny']

        self.df_sequence_properties = pd.DataFrame.from_records(seq_prop_pre_df, columns=cols).drop_duplicates().reset_index(drop=True)
        log.info('Created sequence property dataframe. See the "df_sequence_properties" attribute.')

    def calculate_residue_depth(self, force_rerun=False):
        """Run MSMS on all representative structures, save a json file per structure that contains
            the average residue depth per residue, and the alpha-carbon depth. Units are in Angstroms.

        Args:
            force_rerun (bool): If you want to overwrite any existing files

        """
        self.missing_msms = {}

        for g in tqdm(self.genes):
            gene_id = g.id
            gene_structure_dir = op.join(self.structure_single_chain_dir, gene_id)
            clean_pdb = g.annotation['structure']['representative']['clean_pdb_file']

            if not clean_pdb:
                log.warning('{}: No representative structure available'.format(gene_id))
                continue
            else:
                infile = op.join(gene_structure_dir, clean_pdb)

                try:
                    outfile = ssbio.structure.properties.msms.run_msms(pdb_file=infile,
                                                                       outdir=gene_structure_dir,
                                                                       force_rerun=force_rerun)
                # Log if MSMS was not run on a file
                except AssertionError:
                    self.missing_msms[gene_id] = infile
                    log.warning('{}: MSMS failed to run on {}'.format(gene_id, clean_pdb))
                    continue

                g.annotation['structure']['representative']['properties']['msms_file'] = op.basename(outfile)

        if len(self.missing_msms) > 0:
            log.warning('{} gene(s) failed calculation. Inspect the "missing_msms" attribute.'.format(len(self.missing_msms)))

        log.info('Completed calculations of residue depth')

    def run_pipeline(self, sequence_mapping_engine='', structure_mapping_engine='', **kwargs):
        """Run the entire GEM-PRO pipeline.

        Options include:
        ...

        Returns:

        """
        current_sequence_mapping_engines = ['kegg', 'uniprot', 'all']
        if sequence_mapping_engine not in current_sequence_mapping_engines:
            raise ValueError('Sequence mapping engine not available')

        if sequence_mapping_engine == 'kegg' or sequence_mapping_engine == 'all':
            if not kwargs['kegg_organism_code']:
                raise TypeError('kegg_organism_code needed')

            print('Running KEGG mapping...')
            self.kegg_mapping_and_metadata(kegg_organism_code=kwargs['kegg_organism_code'])

        if sequence_mapping_engine == 'uniprot' or sequence_mapping_engine == 'all':
            if not kwargs['model_gene_source']:
                raise TypeError('UniProt model_gene_source needed')

            print('Running UniProt mapping...')
            self.uniprot_mapping_and_metadata(model_gene_source=kwargs['model_gene_source'])

        print('Setting representative sequences...')
        self.set_representative_sequence()

        current_structure_mapping_engines = ['uniprot', 'itasser', 'all'] # 'blast', 'manual_homology'
        if structure_mapping_engine not in current_structure_mapping_engines:
            raise ValueError('Structure mapping engine not available')

        if structure_mapping_engine == 'uniprot' or structure_mapping_engine == 'all':
            print('Mapping UniProt IDs to the PDB...')
            self.map_uniprot_to_pdb()

        if structure_mapping_engine == 'itasser' or structure_mapping_engine == 'all':
            print('Copying I-TASSER models...')
            self.get_itasser_models(homology_raw_dir=kwargs['homology_raw_dir'],
                                    custom_itasser_name_mapping=kwargs['custom_itasser_name_mapping'])

        print('Setting representative structures...')
        self.set_representative_structure(always_use_homology=kwargs['always_use_homology'])


if __name__ == '__main__':
    pass
    # Run the GEM-PRO pipeline!

    # Parse arguments
#     p = argparse.ArgumentParser(description='Runs the GEM-PRO pipeline')
#     p.add_argument('gemfile', help='Path to the GEM file')
#     p.add_argument('gemname', help='Name you would like to use to refer to this GEM')
#     p.add_argument('rootdir', help='Directory where GEM-PRO files should be stored')
#
#     args = p.parse_args()
#
#     my_gem_pro = GEMPRO(args.gemfile, args.gemname, args.rootdir)
#     my_gem_pro.run_pipeline()
