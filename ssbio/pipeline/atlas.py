import os
import os.path as op
import cobra.flux_analysis
import cobra.manipulation
import numpy as np
import pandas as pd
import ssbio.sequence.utils.blast
from Bio import SeqIO
import ssbio.cobra.utils
import ssbio.databases.ncbi
import ssbio.databases.patric
import ssbio.sequence.utils.fasta
from ssbio import utils
import ssbio.sequence.utils.alignment
date = utils.Date()
from ssbio.pipeline.gempro import GEMPRO
import sys
import logging
import copy
import shutil
from cobra.core import DictList
from cobra.core import Gene
from ssbio.core.genepro import GenePro
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from IPython.display import clear_output
    have_ipython = True
    from tqdm import tqdm_notebook as tqdm
except ImportError:
    have_ipython = False
    from tqdm import tqdm


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

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


class ATLAS():
    """Class to represent an ATLAS workflow to carry out multi-strain comparisons

    Main steps are:
    1. Strain-specific model construction based on orthologous genes & systems modeling
    2. Phylogenetic analysis to pick out important genes
    3. GEM-PRO of the "base strain"
    4. Structure property calculation & integrated structural systems analysis

    Each step may generate a report and also request additional files if something is missing
    """

    def __init__(self, base_gempro, base_genome_path=None):
        """Prepare a GEM-PRO model for ATLAS analysis

        Args:
            base_gempro (GEMPRO): Completed GEM-PRO model
            base_genome_path (str): Simple reference link to the genome FASTA file (CDS)
        """

        # Load the GEM-PRO (could be a model, could just be a list of genes)
        self.base_strain_gempro = base_gempro

        # Check if there is a genome file associated with this model - if not, write all sequences
        if not base_genome_path and self.base_strain_gempro.genome_path:
            self.base_strain_gempro.genome_path = self.base_strain_gempro.write_representative_sequences_file(outname='BASE_CDS')
        else:
            self.base_strain_gempro.genome_path = base_genome_path
            # TODO: must check if base_genome gene IDs can be matched to the base model

        ### Prepare ATLAS directories
        list_of_dirs = []

        # Everything will be stored where the GEM-PRO folder is
        base_gempro_dir = self.base_strain_gempro.base_dir

        # atlas_dir - directory where all ATLAS analysis will be carried out
        self.atlas_dir = op.join(base_gempro_dir, 'atlas')
        list_of_dirs.append(self.atlas_dir)

        # atlas_model_dir - directory where new strain specific models will be stored
        self.atlas_model_dir = op.join(self.atlas_dir, 'models')
        list_of_dirs.append(self.atlas_model_dir)

        # atlas_data_dir - directory where all data will be stored
        self.atlas_data_dir = op.join(self.atlas_dir, 'data')
        list_of_dirs.append(self.atlas_data_dir)

        # atlas_sequence_dir - sequence related files are stored here
        self.atlas_sequence_dir = op.join(self.atlas_dir, 'sequences')
        self.atlas_seq_genomes_dir = op.join(self.atlas_sequence_dir, 'by_organism')
        self.atlas_seq_genes_dir = op.join(self.atlas_sequence_dir, 'by_gene')
        list_of_dirs.extend([self.atlas_sequence_dir, self.atlas_seq_genomes_dir, self.atlas_seq_genes_dir])

        # Make the dirs
        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('Created directory: {}'.format(directory))
            else:
                log.debug('Directory already exists: {}'.format(directory))

        # Other initializations
        self.atlas_strains = DictList()

    def copy_base_model(self, new_id):
        """Copy the base strain model or the genes list (excluding GEM-PRO information) into a new model with a specified ID.

        Appends the model to the atlas_strains attribute.

        Args:
            new_id (str): New ID to be assigned to the copied model

        """
        copied_model = self.base_strain_gempro.model.copy()
        copied_model.id = new_id
        self.atlas_strains.append(copied_model)
        log.debug('{}: new model ID copied from base model'.format(new_id))

    def load_patric_genomes(self, ids_to_genome_file):
        """Load a dictionary of PATRIC IDs which point to the path of their respective genome files.

        Creates initial copies of the base strain model for each strain ID and stores them in self.atlas_strains

        Args:
            ids_to_genome_file (dict): Keys are PATRIC IDs (which will become your strain IDs) and
                values are absolute paths to the FASTA file containing CDS regions
        """
        pass
        # TODO: code to load IDs and link them to the genome files

    def download_genome_cds_patric(self, ids, force_rerun=False):
        """Download genome files from PATRIC give a list of IDs"""
        pass
        # TODO: run the below, and then load_patric_genomes
        # strains_to_fasta_file = {}
        # ids = ssbio.utils.force_list(ids)
        #
        # for patric_id in tqdm(ids):
        #     f = ssbio.databases.patric.download_genome_sequence(patric_id=patric_id, seqtype='protein',
        #                                                         outdir=self.atlas_seq_genomes_dir,
        #                                                         force_rerun=force_rerun)
        #     if f:
        #         strains_to_fasta_file[patric_id] = f
        #         log.debug('{}: Downloaded sequence'.format(patric_id))
        #     else:
        #         log.warning('{}: Unable to download sequence'.format(patric_id))
        #
        # if len(strains_to_fasta_file) > 0:
        #     log.info('Downloaded sequences of coding genes.')
        #
        # return strains_to_fasta_file

    def get_orthology_matrix(self):
        """Run run_makeblastdb, run_bidirectional_blast, and calculate_bbh for protein sequences.
        """
        # TODO: probably best to split this function up to allow for force_rerunning of different parts
        # Get the path to the reference genome
        r_file = self.base_strain_gempro_genome_path
        r_folder, r_name, r_ext = utils.split_folder_and_path(r_file)

        bbh_files = {}

        for strain_model in tqdm(self.atlas_strains):
            g_file = strain_model.get_genome_file_path(self.atlas_seq_genomes_dir)
            # TODO: replace usages of g_name with strain_model.id instead for consistency
            g_folder, g_name, g_ext = utils.split_folder_and_path(g_file)

            # Run bidirectional BLAST
            log.debug('{} vs {}: Running bidirectional BLAST'.format(r_name, g_name))
            r_vs_g, g_vs_r = ssbio.sequence.utils.blast.run_bidirectional_blast(reference=r_file, other_genome=g_file,
                                                                                dbtype='prot',
                                                                                outdir=self.atlas_seq_genomes_dir)

            # Using the BLAST files, find the BBH
            log.debug('{} vs {}: Finding BBHs'.format(r_name, g_name))
            bbh = ssbio.sequence.utils.blast.calculate_bbh(blast_results_1=r_vs_g, blast_results_2=g_vs_r,
                                                           r_name=r_name, g_name=g_name,
                                                           outdir=self.atlas_seq_genomes_dir)
            bbh_files[strain_model.id] = bbh

        # Make the orthologous genes matrix
        log.debug('Creating orthology matrix')
        ortho_matrix = ssbio.sequence.utils.blast.create_orthology_matrix(r_name=r_name,
                                                                          genome_to_bbh_files=bbh_files,
                                                                          outname='{}_{}_orthology.csv'.format(r_name, 'prot'),
                                                                          outdir=self.atlas_data_dir)

        log.info('Saved orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))
        self.df_orthology_matrix = pd.read_csv(ortho_matrix, index_col=0)

        # Remove extraneous strains from our analysis
        to_remove = []
        for i, strain_model in enumerate(self.atlas_strains):
            strain_id = strain_model.id

            if strain_id not in self.df_orthology_matrix.columns:
                log.warning('No orthologous genes found for organism {}. Will remove this strain from analysis.'.format(strain_id))
                to_remove.append(i)
                continue

            not_in_strain = self.df_orthology_matrix[pd.isnull(self.df_orthology_matrix[strain_id])][strain_id].index.tolist()
            if len(not_in_strain) == 0:
                log.info('{}: Strain has no differences from the base. Will remove this strain from analysis.')
                to_remove.append(i)
                continue
        if to_remove:
            # Remove strains with no differences
            for x in to_remove:
                self.atlas_strains.pop(x)
            log.info('Removed {} strains from analysis'.format(len(to_remove)))

        # Also need to check for potential gene IDs that are in the model and not in the orthology matrix
        # This is probably because: the CDS FASTA file for the base strain did not contain the correct ID
        # for the gene and consequently was not included in the orthology matrix
        # We should add these genes to the orthology matrix and set them as not present in the other strains
        base_strain_gene_ids = [x.id for x in self.base_strain_gempro.model.genes]
        in_base_strain_but_not_orthology = [x for x in base_strain_gene_ids if x not in self.df_orthology_matrix.index.tolist()]
        for notinorth in in_base_strain_but_not_orthology:
            # Add empty row
            self.df_orthology_matrix.ix[notinorth] = np.nan

        # Filter the matrix for genes within our base model only
        self.df_orthology_matrix_filtered = self.df_orthology_matrix[self.df_orthology_matrix.index.map(lambda x: x in base_strain_gene_ids)]
        log.info('Created base model specific "df_orthology_matrix_filtered" attribute.'.format(ortho_matrix))

        log.info('{} strains to be analyzed'.format(len(self.atlas_strains)))

    def pare_down_model(self, model_to_be_modified, genes_to_remove):
        """Remove genes from a model. Directly modifies the model.

        Args:
            list_of_strain_ids:

        """
        model_to_be_modified._trimmed = False
        model_to_be_modified._trimmed_genes = []
        model_to_be_modified._trimmed_reactions = {}

        # Filter out genes in genes_to_remove which do not show up in the model
        model_genes = [x.id for x in model_to_be_modified.genes]
        genes_to_remove = list(set(genes_to_remove).intersection(set(model_genes)))

        if len(genes_to_remove) == 0:
            log.debug('No genes marked for removal')
        else:
            log.debug('{} genes marked for removal'.format(len(genes_to_remove)))

            # Delete genes!
            cobra.manipulation.delete_model_genes(model_to_be_modified, genes_to_remove)

            if model_to_be_modified._trimmed:
                log.info('{}: deleted {} reactions, {} genes'.format(model_to_be_modified.id,
                                                                     len(model_to_be_modified._trimmed_reactions),
                                                                     len(model_to_be_modified._trimmed_genes)))

    def build_strain_specific_models(self):
        """Using the orthologous genes matrix, write strain specific models.

        """
        if len(self.df_orthology_matrix_filtered) == 0:
            raise RuntimeError('Empty orthology matrix')

        # For each genome, create the strain specific model
        for strain_model in tqdm(self.strain_models):
            strain_id = strain_model.id

            # Get a list of genes which do not have orthology in the strain
            not_in_strain = self.df_orthology_matrix_filtered[pd.isnull(self.df_orthology_matrix_filtered[strain_id])][strain_id].index.tolist()

            self.pare_down_model(model_to_be_modified=strain_model, genes_to_remove=not_in_strain)

            # Save the strain specific file
            outfile = op.join(self.atlas_model_dir, '{}.json'.format(strain_id))
            cobra.io.save_json_model(strain_model, outfile)
            log.debug('{}: saved model at {}'.format(strain_id, outfile))

        log.info('Created {} new strain-specific models'.format(len(self.strain_models)))

    def align_orthologous_genes_pairwise(self):
        """For each gene in the base strain, run a pairwise alignment for all orthologous gene sequences to it.

        """
        info_pre_df = []

        for base_gene in tqdm(self.base_strain_gempro.genes):
            # Get base strain gene fasta file path
            base_gene_id = base_gene.id
            base_gene_seq_file = base_gene.annotation['sequence']['representative']['sequence_file']
            base_gene_seq_len = base_gene.annotation['sequence']['representative']['seq_len']

            if not base_gene_seq_file:
                log.warning('{}: No representative sequence set in base strain'.format(base_gene_id))
                continue

            base_gene_seq_path = op.join(self.base_strain_gempro.sequence_dir, base_gene_id, base_gene_seq_file)

            gene_dir = op.join(self.atlas_seq_genes_dir, base_gene_id)

            info_dict = {'gene': base_gene_id}
            mutation_count = 0
            num_strains_with_gene = 0

            # Get gene file in all strains if it shows up as functional
            for strain_model in self.strain_models:
                strain_id = strain_model.id
                strain_gene = strain_model.genes.get_by_id(base_gene_id)

                if not strain_gene.functional:
                    log.debug('{}: gene not present in strain {}, not aligning'.format(base_gene_id, strain_id))
                    continue

                strain_gene_seq_path = strain_model.get_gene_sequence_path(base_gene_id, gene_dir)

                alignment_file = ssbio.sequence.utils.alignment.run_needle_alignment_on_files(id_a=base_gene_id,
                                                                                              id_b=strain_id,
                                                                                              faa_a=base_gene_seq_path,
                                                                                              faa_b=strain_gene_seq_path,
                                                                                              outdir=gene_dir)

                # Save in strain gene annotation
                strain_gene.annotation['sequence']['alignment_file'] = op.basename(alignment_file)

                # Some stuff you can count
                alignment_df = ssbio.sequence.utils.alignment.get_alignment_df_from_file(alignment_file)

                # TODO: also count: number of unique mutations (have to consider position, amino acid change)
                # indels? (average length of insertion, deletion)
                # TODO: can i also "save" the alignment somehow? not the dataframe?
                # TODO: keep track of strain with most mutations, least mutations
                # TODO: keep track of strains that conserve the length of the protein, others that extend or truncate it
                # need statistical test for that too (how long is "extended"/"truncated"?)
                # TODO: number of strains with at least 1 mutations
                # TODO: number of strains with <5% mutated, 5-10%, etc

                num_mutations = len(alignment_df[alignment_df.type == 'mutation'])
                mutation_count += num_mutations
                num_strains_with_gene += 1

            if num_strains_with_gene == 0:
                continue

            info_dict['num_strains_with_gene'] = num_strains_with_gene
            info_dict['mutation_count'] = mutation_count
            info_dict['mutation_percent'] = mutation_count / float(num_strains_with_gene*base_gene_seq_len)
            info_pre_df.append(info_dict)

        # Save a dataframe of the file mapping info
        cols = ['gene', 'num_strains_with_gene', 'mutation_count', 'mutation_percent']
        self.df_alignment_stats = pd.DataFrame.from_records(info_pre_df, columns=cols)
        log.info('Created alignment statistics dataframe. See the "df_alignment_stats" attribute.')

    def per_gene_total_num_mutations(self, gene_id):
        """Simply return the total number of mutations for a gene, summed over all strains.

        Returns:
            int: Total number of point mutations that show up in a gene throughout all strains

        """
        total_num_mutations = 0

        gene_dir = op.join(self.atlas_seq_genes_dir, gene_id)

        if not op.exists(gene_dir):
            raise ValueError('{}: Gene does not exist'.format(gene_id))

        for strain_model in tqdm(self.strain_models):
            strain_id = strain_model.id
            strain_gene = strain_model.genes.get_by_id(gene_id)

            if not strain_gene.functional:
                log.debug('{}: gene not present in strain {}, not counting'.format(gene_id, strain_id))
                continue

            alignment_file = op.join(gene_dir, strain_gene.annotation['sequence']['alignment_file'])
            alignment_df = ssbio.sequence.utils.alignment.get_alignment_df_from_file(alignment_file)

            num_mutations = len(alignment_df[alignment_df.type == 'mutation'])

            total_num_mutations += num_mutations

        return total_num_mutations