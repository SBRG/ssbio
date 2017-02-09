import os
import os.path as op
import shutil
import numpy as np
import pandas as pd
from Bio import SeqIO
from bioservices import KEGG
from bioservices import UniProt
from cobra.core import DictList
from cobra.core import Gene
from slugify import slugify
from ssbio.core.genepro import GenePro
import ssbio.cobra.utils
import ssbio.databases.kegg
import ssbio.databases.pdb
import ssbio.databases.uniprot
import ssbio.sequence.properties.residues
import ssbio.sequence.utils.fasta
import ssbio.structure.properties.msms
import ssbio.structure.properties.quality
import ssbio.structure.properties.residues
from ssbio import utils
from ssbio.sequence.seqprop import SeqProp
from ssbio.structure.homology.itasser.itasserprep import ITASSERPrep
if utils.is_ipynb():
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm
import sys
import logging


logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger("requests").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.ERROR)
log = logging.getLogger(__name__)
date = utils.Date()
bs_unip = UniProt()
bs_kegg = KEGG()


def __new__(cls, *args, **kwargs):
    """Casting Gene objects into GenePro objects
        This replaces any new instance of Genes with GenePros. Even when you load a model later
        using COBRApy methods. Use with caution!
    See http://stackoverflow.com/questions/3464061/cast-base-class-to-derived-class-python-or-more-pythonic-way-of-extending-class

    Returns:
        GenePro: a Gene object with a .protein attribute
    """
    if cls == Gene:
        return object.__new__(GenePro)
    return object.__new__(cls)
Gene.__new__ = staticmethod(__new__)


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
    4. Homology modeling
        a. Mapping to existing models
        b. Preparation for running I-TASSER
        c. Parsing I-TASSER runs
    5. Running QC/QA on structures and setting a representative structure

    Each step may generate reports as a Pandas DataFrame and various logging messages.

    """

    def __init__(self, gem_name, root_dir=None, gem_file_path=None, gem_file_type=None,
                 genes_list=None, genes_and_sequences=None, genome_path=None):
        """Initialize the GEM-PRO project with a GEM (genome-scale model) or a list of genes.

        Specify the name of your project, along with the root directory where a folder with that name will be created.

        Args:
            gem_name (str): The name of your GEM or just your project in general.
                This will be the name of the main folder that is created in root_dir.
            root_dir (str): Path to where the folder named gem_name will be created.
                Default is current working directory.
            gem_file_path (str): Path to GEM file
            gem_file_type (str): GEM model type - 'sbml' (or 'xml'), 'mat', or 'json' format
            genes_list (list): List of gene IDs that you want to map
            genes_and_sequences (dict): Dictionary of gene IDs and their amino acid sequence strings
            genome_path (str): Simple reference link to the genome FASTA file (CDS)

        """

        list_of_dirs = []

        # Create directory paths
        if not root_dir:
            root_dir = os.getcwd()
        self.root_dir = root_dir
        self.base_dir = op.join(root_dir, gem_name)
        list_of_dirs.append(self.base_dir)

        # model_dir - directory where original GEMs and GEM-related files are stored
        self.model_dir = op.join(self.base_dir, 'model')
        list_of_dirs.append(self.model_dir)

        # data_dir - directory where all data (dataframes and more) will be stored
        self.data_dir = op.join(self.base_dir, 'data')
        list_of_dirs.append(self.data_dir)

        # sequence_dir - sequence related files are stored here
        self.sequence_dir = op.join(self.base_dir, 'sequences')
        list_of_dirs.append(self.sequence_dir)

        # structure_dir - directory where structure related files will be downloaded/are located
        self.structure_dir = op.join(self.base_dir, 'structures')
        self.structure_single_chain_dir = op.join(self.structure_dir, 'by_gene')
        list_of_dirs.append(self.structure_dir)
        list_of_dirs.append(self.structure_single_chain_dir)

        # Make the directories
        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('{}: created directory'.format(directory))
            else:
                log.debug('{}: directory already exists'.format(directory))

        # Load the model
        if gem_file_path and gem_file_type:
            self.model = ssbio.cobra.utils.model_loader(gem_file_path, gem_file_type)
            log.info('{}: loaded model'.format(gem_file_path))
            log.warning('IMPORTANT: All Gene objects have been transformed into GenePro objects, and will be for any new ones')

            # Place a copy of the current used model in model_dir
            if not op.exists(op.join(self.model_dir, op.basename(gem_file_path))):
                shutil.copy(gem_file_path, self.model_dir)
                log.debug('Copied model file to model directory')
            self.gem_file = op.join(self.model_dir, op.basename(gem_file_path))

            # Obtain list of all gene ids
            self.genes = self.model.genes

            # Log information on the number of things
            log.info('{}: number of reactions'.format(len(self.model.reactions)))
            log.info('{}: number of reactions linked to a gene'.format(ssbio.cobra.utils.true_num_reactions(self.model)))
            log.info('{}: number of genes (excluding spontaneous)'.format(ssbio.cobra.utils.true_num_genes(self.model)))
            log.info('{}: number of metabolites'.format(len(self.model.metabolites)))

        # Or, load a list of gene IDs
        elif genes_list:
            self.genes = genes_list
            log.info('{}: number of genes'.format(len(self.genes)))

        # Or, load a dictionary of genes and their sequences
        elif genes_and_sequences:
            self.genes = list(genes_and_sequences.keys())
            self.manual_seq_mapping(genes_and_sequences)
            log.info('{}: number of genes'.format(len(self.genes)))

        # If neither a model or genes are input, you can still add_genes_by_id later
        else:
            self.genes = []
            log.warning('No model or genes input')

        # Other attributes and dataframes
        self.genome_path = genome_path
        self.df_kegg_metadata = pd.DataFrame()
        self.df_uniprot_metadata = pd.DataFrame()
        self.df_sequence_properties = pd.DataFrame()

        self.df_itasser = pd.DataFrame()
        self.df_pdb_ranking = pd.DataFrame()
        self.df_pdb_blast = pd.DataFrame()
        self.df_pdb_metadata = pd.DataFrame()

    @property
    def genes_with_structures(self):
        return DictList(x for x in self.genes if x.protein.num_structures > 0)

    @property
    def genes_with_experimental_structures(self):
        """Return a DictList of all genes that have at least one experimental structure"""
        return DictList(x for x in self.genes_with_structures if x.protein.num_structures_experimental > 0)

    @property
    def genes_with_homology_models(self):
        """Return a DictList of all homology models in self.structures"""
        return DictList(x for x in self.genes_with_structures if x.protein.num_structures_homology > 0)

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, genes_list):
        """Set the genes attribute to be a DictList of GenePro objects.

        A "protein" attribute will be added to the Genes.

        Args:
            genes_list: DictList of COBRApy Gene objects, or list of gene IDs

        """

        if not isinstance(genes_list, DictList):
            tmp_list = []
            for x in list(set(genes_list)):
                x = str(x)
                new_gene = GenePro(id=x)
                tmp_list.append(new_gene)
            self._genes = DictList(tmp_list)
        else:
            self._genes = genes_list

    def add_genes_by_id(self, genes_list):
        """Add gene IDs manually into our GEM-PRO project.

        Args:
            genes_list (list): List of gene IDs as strings.

        """

        new_genes = []
        for x in list(set(genes_list)):
            new_gene = GenePro(id=x)
            new_genes.append(new_gene)

        # Add unique genes only
        self.genes.union(new_genes)
        # TODO: report union length?
        # log.info('Added {} genes to the list of genes'.format(len(new_genes)))

    ####################################################################################################################
    ### SEQUENCE RELATED METHODS ###
    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_mapping=None, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the KEGG service.

        This function does these things:
            1. Download all metadata and sequence files in the sequences directory
            2. Creates a KEGGProp object in the protein.sequences attribute
            3. Returns a Pandas DataFrame of mapping results

        Args:
            kegg_organism_code (str): The three letter KEGG code of your organism
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model gene IDs.
            force_rerun (bool): If you want to overwrite any existing mappings and files

        """

        # all data will be stored in a dataframe
        kegg_pre_df = []
        # TODO: update dataframe columns (there might be more now)
        df_cols = ['gene', 'kegg', 'refseq', 'uniprot', 'num_pdbs', 'pdbs', 'seq_len', 'sequence_file', 'metadata_file']

        # First map all of the organism's KEGG genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        counter = 0
        kegg_missing = []

        for g in tqdm(self.genes):
            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[g.id]
            else:
                kegg_g = g.id

            # Always make the gene specific folder under the sequence_files directory
            gene_folder = op.join(self.sequence_dir, g.id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            kegg_prop = g.protein.load_kegg(kegg_id=g.id, kegg_organism_code=kegg_organism_code,
                                            download=True, outdir=gene_folder, force_rerun=force_rerun)

            # Update potentially old UniProt ID
            if kegg_g in kegg_to_uniprot.keys():
                kegg_prop.uniprot = kegg_to_uniprot[kegg_g]
                if g.protein.representative_sequence:
                    if g.protein.representative_sequence.kegg == kegg_prop.kegg:
                        g.protein.representative_sequence.uniprot = kegg_to_uniprot[kegg_g]

            # Keep track of missing mappings - missing means no sequence is available
            if not kegg_prop.sequence_file:
                kegg_missing.append(g.id)
            else:
                counter += 1

            # Get a (deep)copy of the KEGG properties
            kegg_dict = kegg_prop.get_dict(df_format=True, only_keys=df_cols)

            # Save in dataframe
            kegg_dict['gene'] = g.id
            kegg_pre_df.append(kegg_dict)

            log.debug('{}: loaded KEGG information for gene'.format(g.id))

        log.info('{} genes mapped to KEGG'.format(counter))

        # Save a dataframe of the file mapping info
        self.df_kegg_metadata = pd.DataFrame.from_records(kegg_pre_df, columns=df_cols)
        self.df_kegg_metadata.fillna(value=np.nan, inplace=True)
        log.info('Created KEGG metadata dataframe. See the "df_kegg_metadata" attribute.')

        # Info on genes that could not be mapped
        if len(kegg_missing) > 0:
            self.missing_kegg_mapping = list(set(kegg_missing))
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
        df_cols = ['gene', 'uniprot', 'reviewed', 'gene_name', 'kegg', 'refseq', 'num_pdbs', 'pdbs', 'ec_number', 'pfam',
                   'seq_len', 'description', 'entry_version', 'seq_version', 'sequence_file', 'metadata_file']

        counter = 0
        uniprot_missing = []
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

            if uniprot_gene not in list(genes_to_uniprots.keys()):
                # Append empty information for a gene that cannot be mapped
                uniprot_prop = SeqProp(ident=gene_id)
                uniprot_dict = uniprot_prop.get_dict(df_format=True, only_keys=df_cols)
                uniprot_dict['gene'] = gene_id
                uniprot_pre_df.append(uniprot_dict)
                uniprot_missing.append(gene_id)
                log.debug('{}: unable to map to UniProt'.format(gene_id))
            else:
                for mapped_uniprot in genes_to_uniprots[uniprot_gene]:
                    uniprot_prop = g.protein.load_uniprot(uniprot_id=mapped_uniprot, download=True,
                                                          outdir=gene_folder, force_rerun=force_rerun)

                    # Keep track of missing ones
                    if not uniprot_prop.sequence_file:
                        uniprot_missing.append(gene_id)
                    else:
                        counter += 1

                    uniprot_dict = uniprot_prop.get_dict(df_format=True, only_keys=df_cols)

                    # Add info to dataframe
                    uniprot_dict['gene'] = gene_id
                    uniprot_pre_df.append(uniprot_dict)

        # Save a dataframe of the UniProt metadata
        if not self.df_uniprot_metadata.empty:
            self.df_uniprot_metadata = self.df_uniprot_metadata.append(uniprot_pre_df, ignore_index=True).drop_duplicates().reset_index(drop=True)
            self.df_uniprot_metadata.fillna(value=np.nan, inplace=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            self.df_uniprot_metadata = pd.DataFrame.from_records(uniprot_pre_df, columns=df_cols)
            self.df_uniprot_metadata.fillna(value=np.nan, inplace=True)
            log.info('Created UniProt metadata dataframe. See the "df_uniprot_metadata" attribute.')

        # Info on genes that could not be mapped
        if len(uniprot_missing) > 0:
            self.missing_uniprot_mapping = list(set(uniprot_missing))
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
        df_cols = ['gene', 'uniprot', 'seq_len', 'sequence_file', 'num_pdbs', 'pdbs', 'gene_name', 'reviewed',
                    'kegg', 'refseq', 'ec_number', 'pfam', 'description', 'entry_version', 'seq_version',
                    'metadata_file']

        for g,u in gene_to_uniprot_dict.items():
            gene = self.genes.get_by_id(g)

            # Make the gene folder
            gene_folder = op.join(self.sequence_dir, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            uniprot_prop = gene.protein.load_uniprot(uniprot_id=u,
                                                     outdir=gene_folder, download=True,
                                                     set_as_representative=True)

            uniprot_dict = uniprot_prop.get_dict(df_format=True, only_keys=df_cols)

            # Add info to dataframe
            uniprot_dict['gene'] = g
            uniprot_pre_df.append(uniprot_dict)

            # Remove the entry from the missing genes list and also the DF
            # If it is in the DF but not marked as missing, we do not remove the information already mapped
            if not self.df_uniprot_metadata.empty:
                if g in self.missing_uniprot_mapping:
                    self.missing_uniprot_mapping.remove(g)
                if g in self.df_uniprot_metadata.gene.tolist():
                    get_index = self.df_uniprot_metadata[self.df_uniprot_metadata.gene == g].index
                    for i in get_index:
                        self.df_uniprot_metadata = self.df_uniprot_metadata.drop(i)

        # Save a dataframe of the UniProt metadata
        # Add all new entries to the df_uniprot_metadata
        if not self.df_uniprot_metadata.empty:
            self.df_uniprot_metadata = self.df_uniprot_metadata.append(uniprot_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing UniProt dataframe.')
        else:
            self.df_uniprot_metadata = pd.DataFrame.from_records(uniprot_pre_df, columns=df_cols)
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
        # Save the sequence information in individual FASTA files
        for g, s in gene_to_seq_dict.items():
            g = str(g)
            gene = self.genes.get_by_id(g)

            gene_folder = op.join(self.sequence_dir, g)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            manual_info = gene.protein.load_manual_sequence_str(ident=g, seq_str=s, outdir=gene_folder, set_as_representative=True)
            log.debug('{}: loaded manually defined sequence information'.format(g))

        log.info('Loaded in {} sequences'.format(len(gene_to_seq_dict)))

    def set_representative_sequence(self):
        """Combine information from KEGG, UniProt, and manual mappings. Saves a DataFrame of results.

        Manual mappings override all existing mappings.
        UniProt mappings override KEGG mappings except when KEGG mappings have PDBs associated with them and UniProt doesn't.

        """
        seq_mapping_pre_df = []
        df_cols = ['gene', 'uniprot', 'kegg', 'num_pdbs', 'pdbs', 'seq_len', 'sequence_file', 'metadata_file']

        sequence_missing = []
        counter = 0
        for gene in self.genes:
            g = gene.id
            repseq = gene.protein.set_representative_sequence()

            if not repseq:
                sequence_missing.append(g)
            elif not repseq.sequence_path:
                sequence_missing.append(g)
            else:
                counter += 1

            # For saving the dataframe
            gene_dict = gene.protein.representative_sequence.get_dict(df_format=True, only_keys=df_cols)
            gene_dict['gene'] = g
            seq_mapping_pre_df.append(gene_dict)

        tmp = pd.DataFrame.from_records(seq_mapping_pre_df, columns=df_cols)

        # TODO: info on genes that could be mapped!

        # Info on genes that could not be mapped
        if len(sequence_missing) > 0:
            self.missing_repseq = list(set(sequence_missing))
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_repseq" attribute.'.format(
                    len(self.missing_repseq)))

        self.df_representative_sequences = tmp[pd.notnull(tmp.sequence_file)].reset_index(drop=True)
        self.df_representative_sequences.fillna(value=np.nan, inplace=True)
        log.info('Created sequence mapping dataframe. See the "df_representative_sequences" attribute.')

    def write_representative_sequences_file(self, outname, outdir=None):
        """Write all the model's sequences as a single FASTA file

        Args:
            outname: Name of FASTA file without the extension
            outdir: Path to output file
            use_original_ids: If the original sequence IDs in the FASTA file should be used, or the model IDs

        """
        if not outdir:
            outdir = self.sequence_dir

        # TODO: not sure if i like this outname / outdir / outext thing
        outfile = op.join(outdir, outname + '.faa')

        tmp = []
        for x in self.genes:
            repseq = x.protein.representative_sequence.seq_record
            if repseq:
                repseq.id = x.id
                tmp.append(repseq)

        SeqIO.write(tmp, outfile, "fasta")

        log.info('{}: wrote all representative sequences to file'.format(outfile))
        return outfile
    ### END SEQUENCE RELATED METHODS ###
    ####################################################################################################################

    ####################################################################################################################
    ### STRUCTURE RELATED METHODS ###
    def map_uniprot_to_pdb(self, seq_ident_cutoff=0, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Creates a summary dataframe accessible by the attribute "df_pdb_ranking".

        """
        best_structures_pre_df = []

        # First get all UniProt IDs and check if they have PDBs
        all_representative_uniprots = []
        for g in self.genes:
            if g.protein.representative_sequence:
                uniprot_id = g.protein.representative_sequence.uniprot
                if uniprot_id:
                    all_representative_uniprots.append(uniprot_id)
        uniprots_to_pdbs = bs_unip.mapping(fr='ACC', to='PDB_ID', query=all_representative_uniprots)

        for g in tqdm(self.genes):
            if not g.protein.representative_sequence:
                # Check if a representative sequence was set
                log.warning('{}: no representative sequence set, cannot use best structures API'.format(g.id))
                continue

            uniprot_id = g.protein.representative_sequence.uniprot
            if not uniprot_id:
                log.warning('{}: no representative UniProt ID set, cannot use best structures API'.format(g.id))
                continue
            else:
                if uniprot_id in uniprots_to_pdbs:
                    best_structures = ssbio.databases.pdb.best_structures(uniprot_id,
                                                                          outname='{}_best_structures'.format(slugify(uniprot_id)),
                                                                          outdir=op.join(self.sequence_dir, g.id),
                                                                          seq_ident_cutoff=seq_ident_cutoff,
                                                                          force_rerun=force_rerun)

                    if best_structures:
                        rank = 1
                        for best_structure in best_structures:
                            currpdb = str(best_structure['pdb_id'].lower())
                            currchain = str(best_structure['chain_id'])

                            # load_pdb will append this protein to the list
                            new_pdb = g.protein.load_pdb(pdb_id=currpdb, mapped_chains=currchain)
                            # Also add this chain to the chains attribute so we can save the info we get from best_structures
                            new_pdb.add_chain_ids(currchain)

                            pdb_specific_keys = ['experimental_method','resolution']
                            chain_specific_keys = ['coverage','start','end','unp_start','unp_end']

                            new_pdb.update(best_structure, only_keys=pdb_specific_keys)
                            new_chain = new_pdb.chains.get_by_id(currchain)
                            new_chain.update(best_structure, only_keys=chain_specific_keys)

                            # For saving in the summary dataframe
                            pdb_dict = new_pdb.get_dict_with_chain(chain=currchain, df_format=True,
                                                                   chain_keys=chain_specific_keys)
                            pdb_dict['gene'] = g.id
                            pdb_dict['uniprot'] = uniprot_id
                            pdb_dict['pdb_id'] = currpdb
                            pdb_dict['pdb_chain_id'] = currchain
                            pdb_dict['rank'] = rank
                            best_structures_pre_df.append(pdb_dict)

                            rank += 1

                        log.debug('{}, {}: {} PDB/chain pairs mapped'.format(g.id, uniprot_id, len(best_structures)))
                    else:
                        log.debug('{}, {}: no PDB/chain pairs mapped'.format(g.id, uniprot_id))
                else:
                    log.debug('{}, {}: no PDBs available'.format(g.id, uniprot_id))

        cols = ['gene', 'uniprot', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'resolution', 'coverage',
                'taxonomy_name', 'start', 'end', 'unp_start', 'unp_end', 'rank']
        self.df_pdb_ranking = pd.DataFrame.from_records(best_structures_pre_df, columns=cols)

        log.info('Completed UniProt -> best PDB mapping. See the "df_pdb_ranking" attribute.')
        log.info('{}: number of genes with at least one structure'.format(len(self.genes_with_experimental_structures)))
        log.info('{}: number of genes with no structures'.format(len(self.genes)-len(self.genes_with_experimental_structures)))

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, all_genes=False, force_rerun=False, display_link=False):
        """BLAST each gene sequence to the PDB. Raw BLAST results (XML files) are saved per gene in the "structures" folder
        """
        blast_results_pre_df = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)

            # Check if a representative sequence was set
            if not g.protein.representative_sequence:
                log.warning('{}: no representative sequence set, cannot BLAST'.format(gene_id))
                continue

            seq_name = g.protein.representative_sequence.id
            seq_str = g.protein.representative_sequence.seq_str

            if not seq_str:
                log.warning('{}: no representative sequence set, cannot BLAST'.format(gene_id))
                continue

            # If all_genes=False, BLAST only genes without a uniprot->pdb mapping

            if g.protein.num_structures > 0 and not all_genes:
                log.debug('{}: skipping BLAST, {} structures already mapped and all_genes flag is False'.format(g.protein.num_structures,
                                                                                                                gene_id))
                continue

            # Make the gene specific folder under the structure_files directory
            gene_folder = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(gene_folder):
                os.mkdir(gene_folder)

            # BLAST the sequence to the PDB
            blast_results = ssbio.databases.pdb.blast_pdb(seq_str,
                                                          outfile='{}_blast_pdb.xml'.format(slugify(seq_name)),
                                                          outdir=op.join(self.sequence_dir, gene_id),
                                                          force_rerun=force_rerun,
                                                          evalue=evalue,
                                                          seq_ident_cutoff=seq_ident_cutoff,
                                                          link=display_link)

            if blast_results:
                # Filter for new BLASTed PDBs
                pdbs = [x['hit_pdb'].lower() for x in blast_results]
                new_pdbs = [y for y in pdbs if not g.protein.structures.has_id(y)]
                if new_pdbs:
                    log.info('{}: adding {} PDBs from BLAST results'.format(gene_id, len(new_pdbs)))
                blast_results = [z for z in blast_results if z['hit_pdb'].lower() in new_pdbs]

                for blast_result in blast_results:
                    pdb = blast_result['hit_pdb'].lower()
                    chains = blast_result['hit_pdb_chains']

                    for chain in chains:
                        # load_pdb will append this protein to the list
                        new_pdb = g.protein.load_pdb(pdb_id=pdb, mapped_chains=chain)
                        new_pdb.add_chain_ids(chain)
                        new_chain = new_pdb.chains.get_by_id(chain)
                        new_chain.update(blast_result, only_keys=['hit_score', 'hit_evalue', 'hit_percent_similar',
                                                                  'hit_percent_ident', 'hit_num_ident', 'hit_num_similar'])
                        # For saving in the summary dataframe
                        pdb_dict = new_pdb.get_dict_with_chain(chain=chain, df_format=True)
                        pdb_dict['gene'] = gene_id
                        pdb_dict['pdb_id'] = pdb
                        pdb_dict['pdb_chain_id'] = chain
                        blast_results_pre_df.append(pdb_dict)

                log.debug('{}: {} PDBs BLASTed'.format(gene_id, len(blast_results)))
            else:
                log.debug('{}: no BLAST results'.format(gene_id))

        cols = ['gene', 'pdb_id', 'pdb_chain_id', 'hit_score', 'hit_evalue', 'hit_percent_similar',
                'hit_percent_ident', 'hit_num_ident', 'hit_num_similar']
        self.df_pdb_blast = pd.DataFrame.from_records(blast_results_pre_df, columns=cols)

        log.info('Completed sequence --> PDB BLAST. See the "df_pdb_blast" attribute.')
        # TODO: log.info for counts - num pdbs with no blast hits, number with (instead of in the for loop)

    def manual_homology_models(self, input_dict, clean=True, force_rerun=False):
        """Copy homology models to the GEM-PRO project.

        Args:
            input_dict: Dictionary of dictionaries of gene names to homology model IDs and information. Input a dict of:
                {model_gene: {homology_model_id1: {'model_file': '/path/to/homology/model'},
                              homology_model_id2: {'model_file': '/path/to/homology/model'}
                             },
                }

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

                new_homology = g.protein.load_generic_structure(hid, hdict['model_file'])

                if clean:
                    new_homology.structure_path = new_homology.clean_structure(outdir=dest_gene_dir, force_rerun=force_rerun)
                else:
                    if ssbio.utils.force_rerun(force_rerun, op.basename(hdict['model_file'])):
                        # Just copy the file to the structure directory and store the file name
                        shutil.copy2(hdict['model_file'], dest_gene_dir)

                # TODO: how to handle other info?
                new_homology.update(hdict)

                log.debug('{}: updated homology model information and copied model file.'.format(gene_id))
            counter += 1

        # TODO: dataframe for this stuff

        log.info('Updated homology model information for {} genes.'.format(counter))

    def get_itasser_models(self, homology_raw_dir, custom_itasser_name_mapping=None, force_rerun=False):
        """Copy generated homology models from a directory to the GEM-PRO directory.

        Args:
            homology_raw_dir: Root directory of I-TASSER folders.
            custom_itasser_name_mapping: Use this if your I-TASSER folder names differ from your model gene names.
                Input a dict of {model_gene: ITASSER_folder}.

        """
        itasser_pre_df = []
        counter = 0

        df_cols = ['gene', 'model_file', 'model_date', 'difficulty',
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

        for g in tqdm(self.genes):
            gene_id = g.id

            if custom_itasser_name_mapping and gene_id in custom_itasser_name_mapping:
                hom_id = custom_itasser_name_mapping[gene_id]
            else:
                hom_id = gene_id

            # The name of the actual pdb file will be $GENEID_model1.pdb
            new_itasser_name = hom_id + '_model1'
            orig_itasser_dir = op.join(homology_raw_dir, hom_id)

            try:
                itasser_prop = g.protein.load_itasser_folder(ident=hom_id, itasser_folder=orig_itasser_dir,
                                                             create_dfs=True, force_rerun=force_rerun)
            except OSError:
                log.debug('{}: no homology model available.'.format(gene_id))
                continue

            if itasser_prop.model_file:
                counter += 1

                # Make the destination structure folder and other results folder
                dest_gene_dir = op.join(self.structure_single_chain_dir, gene_id)
                if not op.exists(dest_gene_dir):
                    os.mkdir(dest_gene_dir)

                dest_itasser_extra_dir = op.join(dest_gene_dir, '{}_itasser'.format(new_itasser_name))
                if not op.exists(dest_itasser_extra_dir):
                    os.mkdir(dest_itasser_extra_dir)

                itasser_prop.copy_results(copy_to_dir=dest_gene_dir, rename_model_to=new_itasser_name, force_rerun=force_rerun)
                itasser_prop.save_dataframes(outdir=dest_itasser_extra_dir)

                new_homology_dict = itasser_prop.get_dict(df_format=True, only_keys=df_cols)

                # TODO: should do the alignment to make calc seq_coverage (representative sequences can change)
                # copied['seq_coverage'] = 1

                new_homology_dict['gene'] = gene_id
                itasser_pre_df.append(new_homology_dict)
            else:
                log.debug('{}: homology modeling unfinished'.format(gene_id))

        # TODO: write a "dataframe_update" method to 1) create a new df if empty, 2) or append to an existing one but
        # checks for existing gene IDs (has to be an empty row value tho?)
        if not self.df_itasser.empty:
            self.df_itasser = self.df_itasser.append(itasser_pre_df, ignore_index=True).reset_index(drop=True)
            log.info('Updated existing I-TASSER dataframe.')
        else:
            self.df_itasser = pd.DataFrame.from_records(itasser_pre_df, columns=df_cols)
        # Drop columns that don't have anything in them
        self.df_itasser.dropna(axis=1, how='all', inplace=True)

        log.info('Completed copying of {} I-TASSER models to GEM-PRO directory'.format(counter))
        log.info('{} I-TASSER models total'.format(len(self.df_itasser)))
        log.info('See the "df_itasser" attribute for a summary dataframe')

    def set_representative_structure(self, seq_ident_cutoff=0.5, always_use_homology=False,
                                     allow_missing_on_termini=0.2, allow_mutants=True, allow_deletions=False,
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
        df_cols = ['gene', 'id', 'is_experimental', 'reference_seq', 'reference_seq_top_coverage']

        rep_struct_pre_df = []
        structure_missing = []

        for g in tqdm(self.genes):
            gene_id = str(g.id)

            gene_seq_dir = op.join(self.sequence_dir, gene_id)
            gene_struct_dir = op.join(self.structure_single_chain_dir, gene_id)
            if not op.exists(gene_struct_dir):
                os.mkdir(gene_struct_dir)

            repstruct = g.protein.set_representative_structure(seq_outdir=gene_seq_dir,
                                                               struct_outdir=gene_struct_dir,
                                                               seq_ident_cutoff=seq_ident_cutoff,
                                                               always_use_homology=always_use_homology,
                                                               allow_missing_on_termini=allow_missing_on_termini,
                                                               allow_mutants=allow_mutants,
                                                               allow_deletions=allow_deletions,
                                                               allow_insertions=allow_insertions,
                                                               allow_unresolved=allow_unresolved,
                                                               force_rerun=force_rerun)

            if not repstruct:
                structure_missing.append(gene_id)
                continue

            repdict = repstruct.get_dict(df_format=True)
            repdict['gene'] = gene_id
            rep_struct_pre_df.append(repdict)

        # Info on genes that could not be mapped
        if len(structure_missing) > 0:
            self.missing_repstruct = list(set(structure_missing))
            log.warning('{} gene(s) could not be mapped. Inspect the "missing_repstruct" attribute.'.format(
                    len(self.missing_repstruct)))

        self.df_representative_structures = pd.DataFrame.from_records(rep_struct_pre_df, columns=df_cols)
        self.df_representative_structures.fillna(value=np.nan, inplace=True)
        log.info('Created representative structures dataframe. See the "df_representative_structures" attribute.')

    def prep_itasser_models(self, itasser_installation, itlib_folder, runtype, create_in_dir=None,
                            execute_from_dir=None, all_genes=False, print_exec=False, **kwargs):
        # TODO: kwargs for slurm options

        if not create_in_dir:
            self.structure_homology_models = op.join(self.structure_dir, 'homology_models')
        else:
            self.structure_homology_models = create_in_dir
        if not op.exists(self.structure_homology_models):
            os.makedirs(self.structure_homology_models)
            log.info('Created directory: {}'.format(self.structure_homology_models))

        if not execute_from_dir:
            execute_from_dir = self.structure_homology_models

        counter = 0
        for g in self.genes:
            gene_id = g.id

            repseq = g.protein.representative_sequence
            if not repseq:
                log.warning('{}: No representative sequence available'.format(gene_id))
                continue

            # Skip homology modeling if there is a structure set
            repstruct = g.protein.representative_structure
            if repstruct and not all_genes:
                log.debug('{}: Representative structure available, skipping homology modeling'.format(gene_id))
                continue

            ITASSERPrep(ident=gene_id, seq_str=repseq.seq_str, root_dir=self.structure_homology_models,
                        itasser_path=itasser_installation, itlib_path=itlib_folder,
                        runtype=runtype, print_exec=print_exec, execute_dir=execute_from_dir, walltime=kwargs['walltime'])
            counter += 1

        log.info('Prepared I-TASSER modeling folders for {} genes'.format(counter))

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
            if g.protein.num_structures_experimental == 0:
                log.debug('{}: No structures available - no structures will be downloaded'.format(gene_id))
                continue

            # Download the PDBs
            for s in g.protein.get_experimental_structures():

                log.debug('{}: Downloading PDB or mmCIF file'.format(s.id))
                s.download_structure_file(outdir=gene_struct_dir, force_rerun=force_rerun)
                s.download_cif_header_file(outdir=gene_struct_dir, force_rerun=force_rerun)

                infodict = s.get_dict(df_format=True)
                infodict['gene'] = gene_id
                infodict['pdb_id'] = s.id
                pdb_pre_df.append(infodict)

        # Save a dataframe of the PDB metadata
        if not self.df_pdb_metadata.empty:
            self.df_pdb_metadata = self.df_pdb_metadata.append(pdb_pre_df, ignore_index=True).drop_duplicates().reset_index(drop=True)
            log.info('Updated existing PDB dataframe.')
        else:
            cols = ['gene', 'pdb_id', 'pdb_title','description', 'experimental_method',
                    'resolution', 'chemicals', 'date', 'taxonomy_name',
                    'structure_file']
            self.df_pdb_metadata = pd.DataFrame.from_records(pdb_pre_df, columns=cols)
            log.info('Created PDB metadata dataframe.')

    # def run_pipeline(self, sequence_mapping_engine='', structure_mapping_engine='', **kwargs):
    #     """Run the entire GEM-PRO pipeline.
    #
    #     Options include:
    #     ...
    #
    #     Returns:
    #
    #     """
    #     current_sequence_mapping_engines = ['kegg', 'uniprot', 'all']
    #     if sequence_mapping_engine not in current_sequence_mapping_engines:
    #         raise ValueError('Sequence mapping engine not available')
    #
    #     if sequence_mapping_engine == 'kegg' or sequence_mapping_engine == 'all':
    #         if 'kegg_organism_code' not in kwargs:
    #             raise TypeError('kegg_organism_code needed')
    #
    #     if sequence_mapping_engine == 'uniprot' or sequence_mapping_engine == 'all':
    #         if 'model_gene_source' not in kwargs:
    #             raise TypeError('UniProt model_gene_source needed')
    #
    #     print('Running KEGG mapping...')
    #     # TODO: better way of passing optional arguments to these things - yaml? There are too many if statements to check for things
    #     if 'custom_gene_mapping' in kwargs:
    #         self.kegg_mapping_and_metadata(kegg_organism_code=kwargs['kegg_organism_code'], custom_gene_mapping=kwargs['custom_gene_mapping'])
    #     else:
    #         self.kegg_mapping_and_metadata(kegg_organism_code=kwargs['kegg_organism_code'])
    #
    #     print('Running UniProt mapping...')
    #     self.uniprot_mapping_and_metadata(model_gene_source=kwargs['model_gene_source'])
    #
    #     print('Setting representative sequences...')
    #     self.set_representative_sequence()
    #
    #     current_structure_mapping_engines = ['uniprot', 'itasser', 'all', 'none'] # 'blast', 'manual_homology'
    #     if structure_mapping_engine not in current_structure_mapping_engines:
    #         raise ValueError('Structure mapping engine not available')
    #
    #     if structure_mapping_engine == 'itasser' or structure_mapping_engine == 'all':
    #         if 'homology_raw_dir' not in kwargs:
    #             raise TypeError('Path to homology models needed as homology_raw_dir')
    #
    #     if not structure_mapping_engine == 'none':
    #         if 'always_use_homology' not in kwargs:
    #             raise TypeError('Flag always_use_homology needed')
    #
    #     if structure_mapping_engine == 'uniprot' or structure_mapping_engine == 'all':
    #         print('Mapping UniProt IDs to the PDB...')
    #         self.map_uniprot_to_pdb()
    #
    #     if structure_mapping_engine == 'itasser' or structure_mapping_engine == 'all':
    #         print('Copying I-TASSER models...')
    #         self.get_itasser_models(homology_raw_dir=kwargs['homology_raw_dir'],
    #                                 custom_itasser_name_mapping=kwargs['custom_itasser_name_mapping'])
    #
    #     if not structure_mapping_engine == 'none':
    #         print('Setting representative structures...')
    #         self.set_representative_structure(always_use_homology=kwargs['always_use_homology'])


if __name__ == '__main__':
    pass
    # TODO: Run the GEM-PRO pipeline!