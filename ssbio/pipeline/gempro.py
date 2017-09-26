import logging
import os
import os.path as op
import shutil
from copy import copy

import pandas as pd
from Bio import SeqIO
from bioservices import KEGG
from bioservices import UniProt
from cobra.core import DictList
from six.moves.urllib.error import HTTPError
from slugify import Slugify

import ssbio.core.modelpro
import ssbio.databases.kegg
import ssbio.databases.pdb
import ssbio.databases.uniprot
import ssbio.protein.sequence.properties.residues
import ssbio.protein.sequence.properties.tmhmm
import ssbio.protein.sequence.utils.fasta
import ssbio.protein.structure.properties.msms
import ssbio.protein.structure.properties.quality
import ssbio.protein.structure.properties.residues
from ssbio import utils
from ssbio.core.genepro import GenePro
from ssbio.core.modelpro import ModelPro
from ssbio.core.object import Object
from ssbio.databases.kegg import KEGGProp
from ssbio.databases.uniprot import UniProtProp
from ssbio.protein.sequence.properties.scratch import SCRATCH

if utils.is_ipynb():
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm

custom_slugify = Slugify(safe_chars='-_.')
logging.getLogger("requests").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.ERROR)
log = logging.getLogger(__name__)
date = utils.Date()
bs_unip = UniProt()
bs_kegg = KEGG()


class GEMPRO(Object):

    """Generic class to represent all information for a GEM-PRO project.

    Initialize the GEM-PRO project with a genome-scale model, a list of genes, or a dict of genes and sequences.
    Specify the name of your project, along with the root directory where a folder with that name will be created.

    Main methods provided are:

    #. Mapping of sequence IDs
        * With KEGG mapper
        * With UniProt mapper
        * Allowing manual gene ID --> protein sequence entry
        * Allowing manual gene ID --> UniProt ID

    #. Consolidating sequence IDs and setting a representative sequence
        * Currently these are set based on available PDB IDs

    #. Mapping of representative sequence --> structures
        * With UniProt --> ranking of PDB structures
        * BLAST representative sequence --> PDB database

    #. Homology modeling
        * Mapping to existing models
        * Preparation for running I-TASSER
        * Parsing I-TASSER runs

    #. Running QC/QA on structures and setting a representative structure
        * Various cutoffs (mutations, insertions, deletions) can be set to filter structures

    Each step will also generate reports as Pandas DataFrames and also output various logging messages.

    Args:
        gem_name (str): The name of your GEM or just your project in general.
            This will be the name of the main folder that is created in root_dir.
        pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` - 
            choose a file type for files downloaded from the PDB
        root_dir (str): Path to where the folder named after ``gem_name`` will be created. If not provided,
            directories will not be created and output directories need to be specified for some steps
        genome_path (str): Simple reference link to the genome FASTA file (CDS)
        gem (Model): COBRApy Model object
        gem_file_path (str): Path to GEM file
        gem_file_type (str): GEM model type - ``sbml`` (or ``xml``), ``mat``, or ``json`` formats
        genes_list (list): List of gene IDs that you want to map
        genes_and_sequences (dict): Dictionary of gene IDs and their amino acid sequence strings
        description (str): Optional string to describe your project

    """

    def __init__(self, gem_name, pdb_file_type='cif', root_dir=None, genome_path=None,
                 gem=None,
                 gem_file_path=None, gem_file_type=None,
                 genes_list=None,
                 genes_and_sequences=None,
                 description=None,
                 custom_spont_id=None):
        Object.__init__(self, id=gem_name, description=description)

        self.custom_spont_id = custom_spont_id
        self.pdb_file_type = pdb_file_type
        self.genome_path = genome_path

        # Create directories
        self._root_dir = None
        if root_dir:
            self.root_dir = root_dir

        # Load a model object
        self.model = None
        if gem:
            self.load_cobra_model(gem)

        # Or, load a GEM file
        elif gem_file_path and gem_file_type:
            gem = ssbio.core.modelpro.model_loader(gem_file_path=gem_file_path,
                                                   gem_file_type=gem_file_type)
            self.load_cobra_model(gem)

        # Or, load a list of gene IDs
        elif genes_list:
            self.genes = genes_list

        # Or, load a dictionary of genes and their sequences
        elif genes_and_sequences:
            self.genes = list(genes_and_sequences.keys())

        # If neither a model or genes are input, you can still add IDs with method add_genes_by_id later
        else:
            self.genes = []
            log.warning('No model or genes input')

        if genes_and_sequences:
            self.manual_seq_mapping(genes_and_sequences)

        log.info('{}: number of genes'.format(len(self.genes)))

    @property
    def root_dir(self):
        """str: Directory where GEM-PRO project folder named after the attribute ``base_dir`` is located"""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.info('Changing root directory of GEM-PRO project "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)):
                raise IOError('GEM-PRO project "{}" does not exist in folder {}'.format(self.id, path))
        else:
            log.info('Creating GEM-PRO project directory in folder {}'.format(path))

        self._root_dir = path

        for d in [self.base_dir, self.model_dir, self.data_dir, self.genes_dir]:
            ssbio.utils.make_dir(d)

        log.info('{}: GEM-PRO project location'.format(self.base_dir))

        # Propagate changes to gene
        if hasattr(self, 'genes'):
            for g in self.genes:
                g.root_dir = self.genes_dir

    @property
    def base_dir(self):
        """str: GEM-PRO project folder"""
        if self.root_dir:
            return op.join(self.root_dir, self.id)
        else:
            return None

    @property
    def model_dir(self):
        """str: Directory where original GEMs and GEM-related files are stored"""
        if self.base_dir:
            return op.join(self.base_dir, 'model')
        else:
            return None

    @property
    def data_dir(self):
        """str: Directory where all data (dataframes and more) will be stored"""
        if self.base_dir:
            return op.join(self.base_dir, 'data')
        else:
            return None

    @property
    def genes_dir(self):
        """str: Directory where all gene specific information is stored"""
        if self.base_dir:
            return op.join(self.base_dir, 'genes')
        else:
            return None

    def load_cobra_model(self, model):
        """Load a COBRApy Model object into the GEM-PRO project.

        Args:
            model (Model): COBRApy ``Model`` object

        """
        self.model = ModelPro(model)
        for g in self.model.genes:
            if self.genes_dir:
                g.root_dir = self.genes_dir
            g.protein.pdb_file_type = self.pdb_file_type
        self.genes = self.model.genes

        log.info('{}: loaded model'.format(model.id))
        log.info('{}: number of reactions'.format(len(self.model.reactions)))
        log.info('{}: number of reactions linked to a gene'.format(ssbio.core.modelpro.true_num_reactions(self.model)))
        log.info('{}: number of genes (excluding spontaneous)'.format(ssbio.core.modelpro.true_num_genes(self.model,
                                                                                                         custom_spont_id=self.custom_spont_id)))
        log.info('{}: number of metabolites'.format(len(self.model.metabolites)))
        log.warning('IMPORTANT: All Gene objects have been transformed into GenePro objects, and will be for any new ones')

    @property
    def genes_with_structures(self):
        """DictList: All genes with any available protein structures"""
        return DictList(x for x in self.genes if x.protein.num_structures > 0)

    @property
    def genes_with_experimental_structures(self):
        """DictList: All genes that have at least one experimental structure"""
        return DictList(x for x in self.genes_with_structures if x.protein.num_structures_experimental > 0)

    @property
    def genes_with_homology_models(self):
        """DictList: All genes that have at least one homology model"""
        return DictList(x for x in self.genes_with_structures if x.protein.num_structures_homology > 0)

    @property
    def genes_with_a_representative_sequence(self):
        """DictList: All genes with a representative sequence"""
        tmp = DictList(x for x in self.genes if x.protein.representative_sequence)
        return DictList(y for y in tmp if y.protein.representative_sequence.seq)

    @property
    def genes_with_a_representative_structure(self):
        """DictList: All genes with a representative protein structure"""
        tmp = DictList(x for x in self.genes if x.protein.representative_structure)
        return DictList(y for y in tmp if y.protein.representative_structure.structure_file)

    @property
    def genes(self):
        return ssbio.core.modelpro.filter_out_spontaneous_genes(self._genes, custom_spont_id=self.custom_spont_id)

    @genes.setter
    def genes(self, genes_list):
        """Set the genes attribute to be a DictList of GenePro objects.

        A "protein" attribute will be added to each Gene.

        Args:
            genes_list: DictList of COBRApy Gene objects, or list of gene IDs

        """

        if not isinstance(genes_list, DictList):
            tmp_list = []
            for x in list(set(genes_list)):
                x = str(x)
                new_gene = GenePro(id=x, pdb_file_type=self.pdb_file_type, root_dir=self.genes_dir)
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
            new_gene = GenePro(id=x, pdb_file_type=self.pdb_file_type, root_dir=self.genes_dir)
            new_genes.append(new_gene)

        orig_num_genes = len(self.genes)

        # Add unique genes only
        if self.model:
            self.model.genes.union(new_genes)
        else:
            self.genes.union(new_genes)

        log.info('Added {}/{} genes to GEM-PRO project'.format(len(self.genes)-orig_num_genes, len(genes_list)))

    ####################################################################################################################
    ### SEQUENCE RELATED METHODS ###
    def kegg_mapping_and_metadata(self, kegg_organism_code, custom_gene_mapping=None, outdir=None,
                                  set_as_representative=False, force_rerun=False):
        """Map all genes in the model to KEGG IDs using the KEGG service.

        Steps:
            1. Download all metadata and sequence files in the sequences directory
            2. Creates a KEGGProp object in the protein.sequences attribute
            3. Returns a Pandas DataFrame of mapping results

        Args:
            kegg_organism_code (str): The three letter KEGG code of your organism
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model gene IDs.
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            set_as_representative (bool): If mapped KEGG IDs should be set as representative sequences
            force_rerun (bool): If you want to overwrite any existing mappings and files

        """

        # First map all of the organism's KEGG genes to UniProt
        kegg_to_uniprot = ssbio.databases.kegg.map_kegg_all_genes(organism_code=kegg_organism_code, target_db='uniprot')

        successfully_mapped_counter = 0

        for g in tqdm(self.genes):
            if custom_gene_mapping:
                kegg_g = custom_gene_mapping[g.id]
            else:
                kegg_g = g.id

            # Download both FASTA and KEGG metadata files
            kegg_prop = g.protein.load_kegg(kegg_id=kegg_g, kegg_organism_code=kegg_organism_code,
                                            download=True, outdir=outdir, set_as_representative=set_as_representative,
                                            force_rerun=force_rerun)

            # Update potentially old UniProt ID
            if kegg_g in kegg_to_uniprot.keys():
                kegg_prop.uniprot = kegg_to_uniprot[kegg_g]
                if g.protein.representative_sequence:
                    if g.protein.representative_sequence.kegg == kegg_prop.kegg:
                        g.protein.representative_sequence.uniprot = kegg_to_uniprot[kegg_g]

            # Keep track of missing mappings - missing is defined by no available sequence
            if kegg_prop.sequence_file:
                successfully_mapped_counter += 1

            log.debug('{}: loaded KEGG information for gene'.format(g.id))

        log.info('{}/{}: number of genes mapped to KEGG'.format(successfully_mapped_counter, len(self.genes)))
        log.info('Completed ID mapping --> KEGG. See the "df_kegg_metadata" attribute for a summary dataframe.')

    @property
    def df_kegg_metadata(self):
        kegg_pre_df = []
        df_cols = ['gene', 'kegg', 'refseq', 'uniprot', 'num_pdbs', 'pdbs', 'seq_len', 'sequence_file', 'metadata_file']

        for g in self.genes:
            kegg_mappings = g.protein.filter_sequences(KEGGProp)
            for kegg_prop in kegg_mappings:
                kegg_dict = kegg_prop.get_dict(df_format=True, only_keys=df_cols)
                kegg_dict['gene'] = g.id
                kegg_pre_df.append(kegg_dict)

        # Save a dataframe of the file mapping info
        df = pd.DataFrame.from_records(kegg_pre_df, columns=df_cols).set_index('gene')
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df)

    @property
    def missing_kegg_mapping(self):
        kegg_missing = []
        for g in self.genes:
            keggs = g.protein.filter_sequences(KEGGProp)
            no_sequence_file_available = True
            for k in keggs:
                if k.sequence_file:
                    no_sequence_file_available = False
                    break
            if no_sequence_file_available:
                kegg_missing.append(g.id)
        return list(set(kegg_missing))

    def uniprot_mapping_and_metadata(self, model_gene_source, custom_gene_mapping=None, outdir=None,
                                     set_as_representative=False, force_rerun=False):
        """Map all genes in the model to UniProt IDs using the UniProt mapping service.
        Also download all metadata and sequences.

        Args:
            model_gene_source (str): the database source of your model gene IDs.
                See: http://www.uniprot.org/help/programmatic_access.
                Common model gene sources are:
                    * Ensembl Genomes - ``ENSEMBLGENOME_ID`` (i.e. E. coli b-numbers)
                    * Entrez Gene (GeneID) - ``P_ENTREZGENEID``
                    * RefSeq Protein - ``P_REFSEQ_AC``
            custom_gene_mapping (dict): If your model genes differ from the gene IDs you want to map,
                custom_gene_mapping allows you to input a dictionary which maps model gene IDs to new ones.
                Dictionary keys must match model genes.
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            set_as_representative (bool): If mapped UniProt IDs should be set as representative sequences
            force_rerun (bool): If you want to overwrite any existing mappings and files

        """

        # Allow model gene --> custom ID mapping ({'TM_1012':'TM1012'})
        if custom_gene_mapping:
            genes_to_map = list(custom_gene_mapping.values())
        else:
            genes_to_map = [x.id for x in self.genes]

        # Map all IDs first to available UniProts
        genes_to_uniprots = bs_unip.mapping(fr=model_gene_source, to='ACC', query=genes_to_map)

        successfully_mapped_counter = 0
        for g in tqdm(self.genes):
            if custom_gene_mapping and g.id in custom_gene_mapping.keys():
                uniprot_gene = custom_gene_mapping[g.id]
            else:
                uniprot_gene = g.id

            if uniprot_gene not in list(genes_to_uniprots.keys()):
                log.debug('{}: unable to map to UniProt'.format(g.id))
            else:
                for mapped_uniprot in genes_to_uniprots[uniprot_gene]:
                    try:
                        uniprot_prop = g.protein.load_uniprot(uniprot_id=mapped_uniprot, download=True, outdir=outdir,
                                                              set_as_representative=set_as_representative,
                                                              force_rerun=force_rerun)
                    except HTTPError as e:
                        log.error('{}, {}: unable to complete web request'.format(g.id, mapped_uniprot))
                        print(e)
                        continue

                    if uniprot_prop.sequence_file or uniprot_prop.metadata_file:
                        successfully_mapped_counter += 1

        log.info('{}/{}: number of genes mapped to UniProt'.format(successfully_mapped_counter, len(self.genes)))
        log.info('Completed ID mapping --> UniProt. See the "df_uniprot_metadata" attribute for a summary dataframe.')

    def manual_uniprot_mapping(self, gene_to_uniprot_dict, outdir=None, set_as_representative=True):
        """Read a manual dictionary of model gene IDs --> UniProt IDs. By default sets them as representative.

        This allows for mapping of the missing genes, or overriding of automatic mappings.
        
        Input a dictionary of::
            
            {
                <gene_id1>: <uniprot_id1>,
                <gene_id2>: <uniprot_id2>,
            }

        Args:
            gene_to_uniprot_dict: Dictionary of mappings as shown above
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            set_as_representative (bool): If mapped UniProt IDs should be set as representative sequences

        """
        for g, u in tqdm(gene_to_uniprot_dict.items()):
            g = str(g)
            gene = self.genes.get_by_id(g)

            try:
                uniprot_prop = gene.protein.load_uniprot(uniprot_id=u,
                                                         outdir=outdir, download=True,
                                                         set_as_representative=set_as_representative)
            except HTTPError as e:
                log.error('{}, {}: unable to complete web request'.format(g, u))
                print(e)
                continue

        log.info('Completed manual ID mapping --> UniProt. See the "df_uniprot_metadata" attribute for a summary dataframe.')

    @property
    def df_uniprot_metadata(self):
        uniprot_pre_df = []
        df_cols = ['gene', 'uniprot', 'reviewed', 'gene_name', 'kegg', 'refseq', 'num_pdbs', 'pdbs', 'ec_number',
                   'pfam', 'seq_len', 'description', 'entry_date', 'entry_version', 'seq_date', 'seq_version',
                   'sequence_file', 'metadata_file']

        for g in self.genes:
            uniprot_mappings = g.protein.filter_sequences(UniProtProp)
            for uniprot_prop in uniprot_mappings:
                uniprot_dict = uniprot_prop.get_dict(df_format=True, only_keys=df_cols)
                uniprot_dict['gene'] = g.id
                uniprot_pre_df.append(uniprot_dict)

        df = pd.DataFrame.from_records(uniprot_pre_df, columns=df_cols).set_index('gene')
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df)

    @property
    def missing_uniprot_mapping(self):
        uniprot_missing = []
        for g in self.genes:
            ups = g.protein.filter_sequences(UniProtProp)

            no_sequence_file_available = True
            for u in ups:
                if u.sequence_file or u.metadata_file:
                    no_sequence_file_available = False
                    break
            if no_sequence_file_available:
                uniprot_missing.append(g.id)

        return list(set(uniprot_missing))

    # TODO: should also have a seq --> uniprot id function (has to be 100% match) (also needs organism)
    def manual_seq_mapping(self, gene_to_seq_dict, outdir=None, set_as_representative=True):
        """Read a manual input dictionary of model gene IDs --> protein sequences. By default sets them as representative.

        Args:
            gene_to_seq_dict (dict): Mapping of gene IDs to their protein sequence strings
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            set_as_representative (bool): If mapped sequences should be set as representative

        """
        if outdir:
            outdir_set = True
        else:
            outdir_set = False

        # Save the sequence information in individual FASTA files
        for g, s in gene_to_seq_dict.items():
            g = str(g)
            gene = self.genes.get_by_id(g)

            if not outdir_set:
                outdir = gene.protein.sequence_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')

            manual_info = gene.protein.load_manual_sequence(ident=g, seq=s, outdir=outdir,
                                                            write_fasta_file=True,
                                                            set_as_representative=set_as_representative)
            log.debug('{}: loaded manually defined sequence information'.format(g))

        log.info('Loaded in {} sequences'.format(len(gene_to_seq_dict)))

    def set_representative_sequence(self, force_rerun=False):
        """Combine information from KEGG, UniProt, and manual mappings. Saves a DataFrame of results.

            Manual mappings override all existing mappings. UniProt mappings override KEGG mappings except
            when KEGG mappings have PDBs associated with them and UniProt doesn't.

        """

        # TODO: rethink use of multiple database sources - may lead to inconsistency with genome sources

        sequence_missing = []
        successfully_mapped_counter = 0
        for g in tqdm(self.genes):
            repseq = g.protein.set_representative_sequence(force_rerun=force_rerun)

            if not repseq:
                sequence_missing.append(g.id)
            elif not repseq.sequence_file:
                sequence_missing.append(g.id)
            else:
                successfully_mapped_counter += 1

        log.info('{}/{}: number of genes with a representative sequence'.format(len(self.genes_with_a_representative_sequence),
                                                                                len(self.genes)))
        log.info('See the "df_representative_sequences" attribute for a summary dataframe.')

    @property
    def df_representative_sequences(self):
        seq_mapping_pre_df = []
        df_cols = ['gene', 'uniprot', 'kegg', 'num_pdbs', 'pdbs', 'seq_len', 'sequence_file', 'metadata_file']

        for g in self.genes_with_a_representative_sequence:
            gene_dict = g.protein.representative_sequence.get_dict(df_format=True, only_keys=df_cols)
            gene_dict['gene'] = g.id
            seq_mapping_pre_df.append(gene_dict)

        df = pd.DataFrame.from_records(seq_mapping_pre_df, columns=df_cols).set_index('gene')
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df)

    @property
    def missing_representative_sequence(self):
        return [x.id for x in self.genes if not self.genes_with_a_representative_sequence.has_id(x.id)]

    def write_representative_sequences_file(self, outname, outdir=None, set_ids_from_model=True):
        """Write all the model's sequences as a single FASTA file. By default, sets IDs to model gene IDs.

        Args:
            outname (str): Name of the output FASTA file without the extension
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            set_ids_from_model (bool): If the gene ID source should be the model gene IDs, not the original sequence ID

        """

        if not outdir:
            outdir = self.data_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        outfile = op.join(outdir, outname + '.faa')

        tmp = []
        for x in self.genes_with_a_representative_sequence:
            repseq = x.protein.representative_sequence
            copied_seq_record = copy(repseq)
            if set_ids_from_model:
                copied_seq_record.id = x.id
            tmp.append(copied_seq_record)

        SeqIO.write(tmp, outfile, "fasta")

        log.info('{}: wrote all representative sequences to file'.format(outfile))
        self.genome_path = outfile
        return self.genome_path

    def get_sequence_properties(self, representatives_only=True):
        """Run Biopython ProteinAnalysis and EMBOSS pepstats to summarize basic statistics of the protein sequences.
        Annotations are stored in the protein's sequence at ``.annotations``

        """
        for g in tqdm(self.genes_with_a_representative_sequence):
            g.protein.get_sequence_properties(representative_only=representatives_only)

    def get_scratch_predictions(self, path_to_scratch, results_dir, scratch_basename='scratch', num_cores=1,
                                exposed_buried_cutoff=25, custom_gene_mapping=None):
        """Run and parse ``SCRATCH`` results to predict secondary structure and solvent accessibility.
        Annotations are stored in the gene protein's representative sequence at:
            * ``.annotations``
            * ``.letter_annotations``

        Args:
            path_to_scratch (str): Path to SCRATCH executable
            results_dir (str): Path to SCRATCH results folder, which will have the files (scratch.ss, scratch.ss8,
                scratch.acc, scratch.acc20)
            scratch_basename (str): Basename of the SCRATCH results ('scratch' is default)
            exposed_buried_cutoff (int): Cutoff of exposed/buried for the acc20 predictions
            custom_gene_mapping (dict): Default parsing of SCRATCH output files is to look for the model gene IDs. If
                your output files contain IDs which differ from the model gene IDs, use this dictionary to map model
                gene IDs to result file IDs. Dictionary keys must match model genes.

        """
        if not self.genome_path:
            # Write all sequences as one file
            all_seqs = self.write_representative_sequences_file()

        # Runs SCRATCH or loads existing results in results_dir
        scratch = SCRATCH(project_name=scratch_basename, seq_file=self.genome_path)
        scratch.run_scratch(path_to_scratch=path_to_scratch, num_cores=num_cores, outdir=results_dir)

        counter = 0

        # Adding the scratch annotations to the representative_sequences letter_annotations
        for g in tqdm(self.genes_with_a_representative_sequence):
            if custom_gene_mapping:
                g_id = custom_gene_mapping[g.id]
            else:
                g_id = g.id

            if g_id in scratch.sspro_summary():
                # Secondary structure
                g.protein.representative_sequence.annotations.update(scratch.sspro_summary()[g_id])
                g.protein.representative_sequence.annotations.update(scratch.sspro8_summary()[g_id])
                g.protein.representative_sequence.letter_annotations['SS-sspro'] = scratch.sspro_results()[g_id]
                g.protein.representative_sequence.letter_annotations['SS-sspro8'] = scratch.sspro8_results()[g_id]

                # Solvent accessibility
                g.protein.representative_sequence.annotations.update(scratch.accpro_summary()[g_id])
                g.protein.representative_sequence.annotations.update(scratch.accpro20_summary(exposed_buried_cutoff)[g_id])
                g.protein.representative_sequence.letter_annotations['RSA-accpro'] = scratch.accpro_results()[g_id]
                g.protein.representative_sequence.letter_annotations['RSA-accpro20'] = scratch.accpro20_results()[g_id]

                counter += 1
            else:
                log.error('{}: missing SCRATCH results'.format(g.id))

        log.info('{}/{}: number of genes with SCRATCH predictions loaded'.format(counter, len(self.genes)))

    def get_tmhmm_predictions(self, tmhmm_results, custom_gene_mapping=None):
        """Parse TMHMM results and store in the representative sequences.
        
        This is a basic function to parse pre-run TMHMM results. Run TMHMM from the 
        web service (http://www.cbs.dtu.dk/services/TMHMM/) by doing the following:
            1. Write all representative sequences in the GEM-PRO using the function ``write_representative_sequences_file``
            2. Upload the file to ``http://www.cbs.dtu.dk/services/TMHMM/`` and choose "Extensive, no graphics" as the output
            3. Copy and paste the results (ignoring the top header and above "HELP with output formats") into a file and save it
            4. Run this function on that file

        Args:
            tmhmm_results (str): Path to TMHMM results (long format)
            custom_gene_mapping (dict): Default parsing of TMHMM output is to look for the model gene IDs. If
                your output file contains IDs which differ from the model gene IDs, use this dictionary to map model
                gene IDs to result file IDs. Dictionary keys must match model genes.

        """
        # TODO: refactor to Protein class
        tmhmm_dict = ssbio.protein.sequence.properties.tmhmm.parse_tmhmm_long(tmhmm_results)

        counter = 0
        for g in tqdm(self.genes_with_a_representative_sequence):
            if custom_gene_mapping:
                g_id = custom_gene_mapping[g.id]
            else:
                g_id = g.id

            if g_id in tmhmm_dict:
                log.debug('{}: loading TMHMM results'.format(g.id))
                if not tmhmm_dict[g_id]:
                    log.error("{}: missing TMHMM results".format(g.id))
                g.protein.representative_sequence.annotations['num_tm_helix-tmhmm'] = tmhmm_dict[g_id]['num_tm_helices']
                g.protein.representative_sequence.letter_annotations['TM-tmhmm'] = tmhmm_dict[g_id]['sequence']
                counter += 1
            else:
                log.error("{}: missing TMHMM results".format(g.id))

        log.info('{}/{}: number of genes with TMHMM predictions loaded'.format(counter, len(self.genes)))

    ### END SEQUENCE RELATED METHODS ###
    ####################################################################################################################

    ####################################################################################################################
    ### STRUCTURE RELATED METHODS ###
    def map_uniprot_to_pdb(self, seq_ident_cutoff=0.0, outdir=None, force_rerun=False):
        """Map UniProt IDs to a ranked list of PDB structures available.

        Uses the "Best structures" API availble from https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
        The list of PDB structures mapping to a UniProt accession sorted by coverage of the protein and,
        if the same, resolution.

        Args:
            seq_ident_cutoff (float): Cutoff results based on percent coverage (in decimal form)
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            force_rerun (bool): Obtain best structures mapping ignoring previously downloaded results

        """

        # First get all UniProt IDs and check if they have PDBs
        all_representative_uniprots = []
        for g in self.genes_with_a_representative_sequence:
            uniprot_id = g.protein.representative_sequence.uniprot
            if uniprot_id:
                # TODO: add warning or something for isoform ids?
                if '-' in uniprot_id:
                    uniprot_id = uniprot_id.split('-')[0]
                all_representative_uniprots.append(uniprot_id)
        log.info('Mapping UniProt IDs --> PDB IDs...')
        uniprots_to_pdbs = bs_unip.mapping(fr='ACC', to='PDB_ID', query=all_representative_uniprots)

        counter = 0
        # Now run the best_structures API for all genes
        for g in tqdm(self.genes_with_a_representative_sequence):
            uniprot_id = g.protein.representative_sequence.uniprot
            if uniprot_id:
                if '-' in uniprot_id:
                    uniprot_id = uniprot_id.split('-')[0]
                if uniprot_id in uniprots_to_pdbs:
                    best_structures = g.protein.map_uniprot_to_pdb(seq_ident_cutoff=seq_ident_cutoff, outdir=outdir, force_rerun=force_rerun)
                    if best_structures:
                        counter += 1
                        log.debug('{}: {} PDBs mapped'.format(g.id, len(best_structures)))
                else:
                    log.debug('{}, {}: no PDBs available'.format(g.id, uniprot_id))

        log.info('{}/{}: number of genes with at least one experimental structure'.format(len(self.genes_with_experimental_structures),
                                                                                          len(self.genes)))
        log.info('Completed UniProt --> best PDB mapping. See the "df_pdb_ranking" attribute for a summary dataframe.')

    @property
    def df_pdb_ranking(self):
        df = pd.DataFrame()
        for g in self.genes_with_experimental_structures:
            protein_df = g.protein.df_pdb_ranking.copy().reset_index()
            if not protein_df.empty:
                protein_df['gene'] = g.id
                df = df.append(protein_df)
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df.set_index('gene'))

    def blast_seqs_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, all_genes=False, display_link=False,
                          outdir=None, force_rerun=False):
        """BLAST each gene sequence to the PDB. Saves raw BLAST results (XML files).

        Also creates a summary dataframe accessible by the attribute "df_pdb_blast".

        Args:
            seq_ident_cutoff (float, optional): Cutoff results based on percent coverage (in decimal form)
            evalue (float, optional): Cutoff for the E-value - filters for significant hits. 0.001 is liberal, 0.0001 is stringent (default).
            all_genes (bool): If all genes should be BLASTed, or only those without any structures currently mapped
            display_link (bool, optional): Set to True if links to the HTML results should be displayed
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            force_rerun (bool, optional): If existing BLAST results should not be used, set to True. Default is False

        """
        counter = 0

        for g in tqdm(self.genes_with_a_representative_sequence):
            # If all_genes=False, BLAST only genes without a uniprot -> pdb mapping
            if g.protein.num_structures_experimental > 0 and not all_genes and not force_rerun:
                log.debug('{}: skipping BLAST, {} experimental structures already mapped '
                          'and all_genes flag is False'.format(g.id,
                                                               g.protein.num_structures_experimental))
                continue

            # BLAST the sequence to the PDB
            new_pdbs = g.protein.blast_representative_sequence_to_pdb(seq_ident_cutoff=seq_ident_cutoff,
                                                                   evalue=evalue,
                                                                   display_link=display_link,
                                                                   outdir=outdir,
                                                                   force_rerun=force_rerun)

            if new_pdbs:
                counter += 1
                log.debug('{}: {} PDBs BLASTed'.format(g.id, len(new_pdbs)))
            else:
                log.debug('{}: no BLAST results'.format(g.id))

        log.info('Completed sequence --> PDB BLAST. See the "df_pdb_blast" attribute for a summary dataframe.')
        log.info('{}: number of genes with additional structures added from BLAST'.format(counter))

    @property
    def df_pdb_blast(self):
        df = pd.DataFrame()
        for g in self.genes_with_experimental_structures:
            protein_df = g.protein.df_pdb_blast.copy().reset_index()
            if not protein_df.empty:
                protein_df['gene'] = g.id
                df = df.append(protein_df)
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df.set_index('gene'))

    def get_manual_homology_models(self, input_dict, outdir=None, clean=True, force_rerun=False):
        """Copy homology models to the GEM-PRO project.

        Requires an input of a dictionary formatted like so::
        
            {
                model_gene: {
                                homology_model_id1: {
                                                        'model_file': '/path/to/homology/model.pdb',
                                                        'file_type': 'pdb'
                                                        'additional_info': info_value
                                                    },
                                homology_model_id2: {   
                                                        'model_file': '/path/to/homology/model.pdb'
                                                        'file_type': 'pdb'
                                                    }
                            }
            }

        Args:
            input_dict (dict): Dictionary of dictionaries of gene names to homology model IDs and other information
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            clean (bool): If homology files should be cleaned and saved as a new PDB file
            force_rerun (bool): If homology files should be copied again even if they exist in the GEM-PRO directory

        """
        if outdir:
            outdir_set = True
        else:
            outdir_set = False

        counter = 0
        for g in tqdm(self.genes):
            if g.id not in input_dict:
                continue

            if not outdir_set:
                outdir = g.protein.structure_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')

            for hid, hdict in input_dict[g.id].items():
                if 'model_file' not in hdict or 'file_type' not in hdict:
                    raise KeyError('"model_file" and "file_type" must be keys in the manual input dictionary.')

                new_homology = g.protein.load_pdb(pdb_id=hid, pdb_file=hdict['model_file'],
                                                  file_type=hdict['file_type'], is_experimental=False)

                if clean:
                    new_homology.load_structure_path(new_homology.clean_structure(outdir=outdir, force_rerun=force_rerun),
                                                     hdict['file_type'])
                else:
                    if ssbio.utils.force_rerun(force_rerun, op.basename(hdict['model_file'])):
                        # Just copy the file to the structure directory and store the file name
                        shutil.copy2(hdict['model_file'], outdir)
                        new_homology.load_structure_path(op.join(outdir, hdict['model_file']), hdict['file_type'])

                # TODO: need to better handle other info in the provided dictionary, if any
                new_homology.update(hdict)

                log.debug('{}: updated homology model information and copied model file.'.format(g.id))
            counter += 1

        log.info('Updated homology model information for {} genes.'.format(counter))

    def get_itasser_models(self, homology_raw_dir, custom_itasser_name_mapping=None, outdir=None, force_rerun=False):
        """Copy generated I-TASSER models from a directory to the GEM-PRO directory.

        Args:
            homology_raw_dir (str): Root directory of I-TASSER folders.
            custom_itasser_name_mapping (dict): Use this if your I-TASSER folder names differ from your model gene names.
                Input a dict of {model_gene: ITASSER_folder}.
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            force_rerun (bool): If homology files should be copied again even if they exist in the GEM-PRO directory

        """
        counter = 0

        for g in tqdm(self.genes):
            if custom_itasser_name_mapping and g.id in custom_itasser_name_mapping:
                hom_id = custom_itasser_name_mapping[g.id]
                if not op.exists(op.join(homology_raw_dir, hom_id)):
                    hom_id = g.id
            else:
                hom_id = g.id

            # The name of the actual pdb file will be $GENEID_model1.pdb
            new_itasser_name = hom_id + '_model1'
            orig_itasser_dir = op.join(homology_raw_dir, hom_id)

            try:
                itasser_prop = g.protein.load_itasser_folder(ident=hom_id, itasser_folder=orig_itasser_dir,
                                                             organize=True, outdir=outdir, organize_name=new_itasser_name,
                                                             create_dfs=True, force_rerun=force_rerun)
            except OSError:
                log.debug('{}: homology model folder unavailable'.format(g.id))
                continue

            if itasser_prop.model_file:
                counter += 1
            else:
                log.debug('{}: homology model file unavailable, perhaps modelling did not finish'.format(g.id))

        log.info('Completed copying of {} I-TASSER models to GEM-PRO directory. See the "df_homology_models" attribute for a summary dataframe.'.format(counter))

    @property
    def df_homology_models(self):
        df = pd.DataFrame()
        for g in self.genes_with_homology_models:
            protein_df = g.protein.df_homology_models.copy().reset_index()
            if not protein_df.empty:
                protein_df['gene'] = g.id
                df = df.append(protein_df)
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df.set_index('gene'))

    def set_representative_structure(self, seq_outdir=None, struct_outdir=None, pdb_file_type=None,
                                     engine='needle', always_use_homology=False, rez_cutoff=0.0,
                                     seq_ident_cutoff=0.5, allow_missing_on_termini=0.2,
                                     allow_mutants=True, allow_deletions=False,
                                     allow_insertions=False, allow_unresolved=True,
                                     force_rerun=False):
        """Set the representative protein structure for a gene.
        Each gene can have a combination of the following, which will be analyzed to set a representative structure.
            * Homology model(s)
            * Ranked PDBs
            * BLASTed PDBs

        If the ``always_use_homology`` flag is true, homology models are always set as representative when they exist.
        If there are multiple homology models, we rank by the percent sequence coverage.

        Args:
            seq_outdir (str): Path to output directory of sequence alignment files, must be set if GEM-PRO directories
                were not created initially
            struct_outdir (str): Path to output directory of structure files, must be set if GEM-PRO directories
                were not created initially
            pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``,
                ``mmtf.gz`` - PDB structure file type that should be downloaded
            engine (str): ``needle`` or ``biopython`` - which pairwise sequence alignment engine should be used
                needle is the standard EMBOSS tool to run pairwise alignments
                biopython is Biopython's implementation of needle. Results can differ!
            always_use_homology (bool): If homology models should always be set as the representative structure
            rez_cutoff (float): Resolution cutoff, in Angstroms (only if experimental structure)
            seq_ident_cutoff (float): Percent sequence identity cutoff, in decimal form
            allow_missing_on_termini (float): Percentage of the total length of the reference sequence which will be ignored
                when checking for modifications. Example: if 0.1, and reference sequence is 100 AA, then only residues
                5 to 95 will be checked for modifications.
            allow_mutants (bool): If mutations should be allowed or checked for
            allow_deletions (bool): If deletions should be allowed or checked for
            allow_insertions (bool): If insertions should be allowed or checked for
            allow_unresolved (bool): If unresolved residues should be allowed or checked for
            force_rerun (bool): If sequence to structure alignment should be rerun

        """
        for g in tqdm(self.genes):
            repstruct = g.protein.set_representative_structure(seq_outdir=seq_outdir,
                                                               struct_outdir=struct_outdir,
                                                               pdb_file_type=pdb_file_type,
                                                               engine=engine,
                                                               rez_cutoff=rez_cutoff,
                                                               seq_ident_cutoff=seq_ident_cutoff,
                                                               always_use_homology=always_use_homology,
                                                               allow_missing_on_termini=allow_missing_on_termini,
                                                               allow_mutants=allow_mutants,
                                                               allow_deletions=allow_deletions,
                                                               allow_insertions=allow_insertions,
                                                               allow_unresolved=allow_unresolved,
                                                               force_rerun=force_rerun)

        log.info('{}/{}: number of genes with a representative structure'.format(len(self.genes_with_a_representative_structure),
                                                                                 len(self.genes)))
        log.info('See the "df_representative_structures" attribute for a summary dataframe.')

    @property
    def df_representative_structures(self):
        rep_struct_pre_df = []
        df_cols = ['gene', 'id', 'is_experimental', 'file_type', 'structure_file']

        for g in self.genes_with_a_representative_structure:
            repdict = g.protein.representative_structure.get_dict(df_format=True, only_keys=df_cols)
            repdict['gene'] = g.id
            rep_struct_pre_df.append(repdict)

        df = pd.DataFrame.from_records(rep_struct_pre_df, columns=df_cols).set_index('gene')
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df)

    @property
    def missing_representative_structure(self):
        return [x.id for x in self.genes if not self.genes_with_a_representative_structure.has_id(x.id)]

    def prep_itasser_models(self, itasser_installation, itlib_folder, runtype, create_in_dir=None,
                            execute_from_dir=None, all_genes=False, print_exec=False, **kwargs):
        """Prepare to run I-TASSER homology modeling for genes without structures, or all genes.

        Args:
            itasser_installation (str): Path to I-TASSER folder, i.e. ``~/software/I-TASSER4.4``
            itlib_folder (str): Path to ITLIB folder, i.e. ``~/software/ITLIB``
            runtype: How you will be running I-TASSER - local, slurm, or torque
            create_in_dir (str): Local directory where folders will be created
            execute_from_dir (str): Optional path to execution directory - use this if you are copying the homology
                models to another location such as a supercomputer for running
            all_genes (bool): If all genes should be prepped, or only those without any mapped structures
            print_exec (bool): If the execution statement should be printed to run modelling
            **kwargs:

        Todo:
            * kwargs - extra options for SLURM or Torque execution

        """
        # TODO: kwargs for slurm/torque options

        if not create_in_dir:
            if not self.data_dir:
                raise ValueError('Output directory must be specified')
            self.homology_models_dir = op.join(self.data_dir, 'homology_models')
        else:
            self.homology_models_dir = create_in_dir

        ssbio.utils.make_dir(self.homology_models_dir)

        if not execute_from_dir:
            execute_from_dir = self.homology_models_dir

        counter = 0
        for g in self.genes_with_a_representative_sequence:
            repstruct = g.protein.representative_structure
            if repstruct and not all_genes:
                log.debug('{}: representative structure set, skipping homology modeling'.format(g.id))
                continue

            g.protein.prep_itasser_modeling(itasser_installation=itasser_installation,
                                            itlib_folder=itlib_folder, runtype=runtype,
                                            create_in_dir=self.homology_models_dir,
                                            execute_from_dir=execute_from_dir,
                                            print_exec=print_exec, **kwargs)
            counter += 1

        log.info('Prepared I-TASSER modeling folders for {} genes in folder {}'.format(counter,
                                                                                       self.homology_models_dir))

    def pdb_downloader_and_metadata(self, outdir=None, pdb_file_type=None, force_rerun=False):
        """Download all structures which have been mapped to our genes.

        Args:
            outdir (str): Path to output directory of downloaded files, must be set if GEM-PRO directories
                were not created initially
            pdb_file_type (str): pdb, pdb.gz, mmcif, cif, cif.gz, xml.gz, mmtf, mmtf.gz - PDB structure file type that
                should be downloaded
            force_rerun (bool): If you want to overwrite any existing mappings and files

        """

        if not pdb_file_type:
            pdb_file_type = self.pdb_file_type

        counter = 0
        for g in tqdm(self.genes):
            pdbs = g.protein.pdb_downloader_and_metadata(outdir=outdir, pdb_file_type=pdb_file_type, force_rerun=force_rerun)

            if pdbs:
                counter += len(pdbs)

        log.info('Updated PDB metadata dataframe. See the "df_pdb_metadata" attribute for a summary dataframe.')
        log.info('Saved {} structures total'.format(counter))

    @property
    def df_pdb_metadata(self):
        """DataFrame: Parse PDB file headers and create a metadata table"""
        df = pd.DataFrame()
        for g in self.genes_with_experimental_structures:
            # Get per protein DataFrame
            protein_df = g.protein.df_pdb_metadata.copy().reset_index()
            protein_df['gene'] = g.id
            df = df.append(protein_df)
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df.set_index('gene'))

    @property
    def df_proteins(self):
        """DataFrame: Get a summary dataframe of all proteins in the model"""
        pre_df = []
        df_cols = ['gene', 'id', 'sequences', 'num_sequences', 'representative_sequence', 'num_structures',
                   'experimental_structures', 'num_experimental_structures',
                   'homology_models', 'num_homology_models',
                   'representative_structure', 'representative_chain', 'representative_chain_seq_coverage',
                   'num_sequence_alignments', 'num_structure_alignments']

        for g in self.genes:
            # Get per protein DataFrame
            protein_dict = g.protein.protein_statistics
            protein_dict['gene'] = g.id
            pre_df.append(protein_dict)

        df = pd.DataFrame.from_records(pre_df, columns=df_cols).set_index('gene')
        if df.empty:
            log.warning('Empty dataframe')
            return df
        else:
            return ssbio.utils.clean_df(df)

    def get_dssp_annotations(self, representatives_only=True, force_rerun=False):
        """Run DSSP on all protein structures and store calculations.
        Annotations are stored in each protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-dssp']``
        
        """
        for g in tqdm(self.genes):
            g.protein.get_dssp_annotations(representative_only=representatives_only, force_rerun=force_rerun)

    def get_msms_annotations(self, representatives_only=True, force_rerun=False):
        """Run MSMS on all protein structures and store calculations.
        Annotations are stored in each protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-msms']``

        """
        for g in tqdm(self.genes):
            g.protein.get_msms_annotations(representative_only=representatives_only, force_rerun=force_rerun)

    def get_freesasa_annotations(self, include_hetatms=False, representatives_only=True, force_rerun=False):
        """Run freesasa on all protein structures and store calculations.
        Annotations are stored in each protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-freesasa']``

        """
        for g in tqdm(self.genes):
            g.protein.get_freesasa_annotations(include_hetatms=include_hetatms,
                                               representative_only=representatives_only,
                                               force_rerun=force_rerun)


    def get_disulfide_bridges(self, representatives_only=True):
        """Run Biopython's disulfide bridge finder and store found bridges.
        Annotations are stored in each protein structure's chain sequence at:
        ``seq_record.annotations['SSBOND-biopython']``

        """
        for g in tqdm(self.genes):
            g.protein.find_disulfide_bridges(representative_only=representatives_only)


    ### END STRUCTURE RELATED METHODS ###
    ####################################################################################################################

    def __json_encode__(self):
        to_return = {}
        # Don't save properties, methods in the JSON
        for x in [a for a in dir(self) if not a.startswith('__') and not isinstance(getattr(type(self), a, None), property) and not callable(getattr(self,a))]:
            if self.model and x == 'genes':
                continue
            to_return.update({x: getattr(self, x)})
        return to_return

    def __json_decode__(self, **attrs):
        for k, v in attrs.items():
            setattr(self, k, v)
        if not self.model:
            self.genes = DictList(self.genes)
        else:
            self.genes = self.model.genes


if __name__ == '__main__':
    pass
    # TODO: Run the GEM-PRO pipeline!
