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

try:
    from IPython.display import clear_output
    have_ipython = True
    from tqdm import tqdm_notebook as tqdm
except ImportError:
    have_ipython = False
    from tqdm import tqdm

from ssbio import utils
import ssbio.sequence.utils.alignment

date = utils.Date()
from ssbio.pipeline.gempro import GEMPRO
from cobra.core import Model
import sys
import logging
import copy
import shutil

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from cobra.core import DictList
from cobra.core import Gene
# Casting Gene objects into GenePro objects
# This replaces any new instance of Genes with GenePros. Even when you load a model later using COBRApy methods
# Use with caution!
# See http://stackoverflow.com/questions/3464061/cast-base-class-to-derived-class-python-or-more-pythonic-way-of-extending-class
from ssbio.core.genepro import GenePro
def __new__(cls, *args, **kwargs):
    if cls == Gene:
        return object.__new__(GenePro)
    return object.__new__(cls)
Gene.__new__ = staticmethod(__new__)

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)

class AtlasModel(Model):
    """Class to represent a model within the ATLAS analysis. Just adding some fields to the annotation attributes.
    """

    def __init__(self, model, genome_file_name):
        Model.__init__(self, model)
        self.annotation['genome_file'] = op.basename(genome_file_name)

        for g in model.genes:
            g.annotation['sequence'] = {'seq_len'       : 0,
                                        'sequence_file'      : None,
                                        'alignment_file': None,
                                        'properties'    : {}
                                        }
            g.annotation['structure'] = {'base_structure': g.annotation['structure']['representative']}

    def get_genome_file_path(self, genome_file_dir):
        return op.join(genome_file_dir, self.annotation['genome_file'])

    def get_gene_sequence_path(self, gene_id, gene_dir):
        g = self.genes.get_by_id(gene_id)
        if not g.annotation['sequence']['sequence_file']:
            raise IOError('{}: Sequence does not exist'.format(gene_id))
        g_seq_path = op.join(gene_dir, g.annotation['sequence']['sequence_file'])
        return g_seq_path


class ATLAS():
    """Class to represent an ATLAS workflow to carry out multi-strain comparisons

    Main steps are:
    1. Strain-specific model construction based on orthologous genes & systems modeling
    2. Phylogenetic analysis to pick out important genes
    3. GEM-PRO of the "base strain"
    4. Structure property calculation & integrated structural systems analysis

    Each step may generate a report and also request additional files if something is missing
    """

    def __init__(self,
                 base_gempro_name,
                 base_dir,
                 base_gempro_file,
                 base_genome_id,
                 list_of_ncbi_ids=None,
                 email_for_ncbi='',
                 list_of_patric_ids=None,
                 dict_of_genomes=None,
                 seq_type='protein'):
        """Prepare for ATLAS analysis.

        """

        # Load the GEM-PRO
        self.base_strain_gempro = GEMPRO(gem_name=base_gempro_name, root_dir=base_dir,
                                         gem_file_path=base_gempro_file, gem_file_type='json')

        # Prepare ATLAS directories
        list_of_dirs = []
        self.base_gempro_dir = self.base_strain_gempro.base_dir
        # atlas_dir - directory where all ATLAS analysis will be carried out
        self.atlas_dir = op.join(self.base_gempro_dir, 'atlas')
        # model_dir - directory where base strain GEM and new models will be stored
        self.atlas_model_files = op.join(self.atlas_dir, 'models')
        # atlas_data_dir - directory where all data will be stored
        self.atlas_data_dir = op.join(self.atlas_dir, 'data')
        # figure_dir - directory where all figure_dir will be stored
        self.atlas_figure_dir = op.join(self.atlas_dir, 'figures')
        # sequence_dir - sequence related files are stored here
        self.atlas_sequence_dir = op.join(self.atlas_dir, 'sequences')

        list_of_dirs.extend([self.atlas_dir, self.atlas_data_dir, self.atlas_figure_dir,
                             self.atlas_model_files, self.atlas_sequence_dir])

        # Check to see what analysis we'll be doing
        self.seq_type = seq_type
        if seq_type == 'dna':
            self.seq_atlas_dir = op.join(self.atlas_sequence_dir, 'dna')
            self.seq_atlas_org_dir = op.join(self.seq_atlas_dir, 'by_organism')
            self.seq_atlas_gene_dir = op.join(self.seq_atlas_dir, 'by_gene')
            self.fasta_extension = 'fna'
            self.blast_seq_type = 'nucl'
        elif seq_type == 'protein':
            self.seq_atlas_dir = op.join(self.atlas_sequence_dir, 'protein')
            self.seq_atlas_org_dir = op.join(self.seq_atlas_dir, 'by_organism')
            self.seq_atlas_gene_dir = op.join(self.seq_atlas_dir, 'by_gene')
            self.fasta_extension = 'faa'
            self.blast_seq_type = 'prot'
        else:
            raise ValueError('seqtype must be "dna" or "protein"')

        # Make more dirs based on this analysis
        list_of_dirs.extend([self.seq_atlas_dir, self.seq_atlas_org_dir, self.seq_atlas_gene_dir])

        # Make the dirs
        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('Created directory: {}'.format(directory))
            else:
                log.debug('Directory already exists: {}'.format(directory))

        strains_to_fasta_file = {}
        # Loading the strains
        if list_of_ncbi_ids:
            if not email_for_ncbi:
                raise ValueError('Input an email if downloading from NCBI')
            log.info('Downloading genomes from NCBI...')
            strains_to_fasta_file.update(self.download_genome_cds_ncbi(list_of_patric_ids, email_for_ncbi))
        if list_of_patric_ids:
            log.info('Downloading genomes from PATRIC...')
            strains_to_fasta_file.update(self.download_genome_cds_patric(list_of_patric_ids))
        if dict_of_genomes:
            log.info('Loading existing genomes to ATLAS sequence directory...')
            strains_to_fasta_file.update(dict_of_genomes)
            for v in dict_of_genomes.values():
                if not op.exists(op.join(self.seq_atlas_org_dir, op.basename(v))):
                    shutil.copy(v, self.seq_atlas_org_dir)

        # Create initial strain models that are just copies of the base strain with a resetted annotation
        self.strain_ids = []
        self.strain_models = DictList([])
        log.info('Creating initial strain specific models based on the base model...')
        for strain_id, strain_fasta in tqdm(strains_to_fasta_file.items()):
            self.strain_ids.append(strain_id)

            if strain_id == base_genome_id:
                log.debug('Not making strain specific model for the base model')
                continue

            base_strain_model_copy = copy.deepcopy(self.base_strain_gempro.model)
            base_strain_model_copy.id = strain_id
            self.strain_models.append(AtlasModel(model=base_strain_model_copy,
                                                 genome_file_name=strain_fasta))

        # Set the base, reference genome
        self._reference_genome = base_genome_id
        self.base_strain_gempro.model.annotation['genome_file'] = op.basename(strains_to_fasta_file[base_genome_id])
        self.base_strain_gempro_genome_file = op.basename(strains_to_fasta_file[base_genome_id])
        self.base_strain_gempro_genome_path = op.join(self.seq_atlas_org_dir, self.base_strain_gempro_genome_file)

    def download_genome_cds_patric(self, ids, force_rerun=False):
        """Download genome files from PATRIC

        """
        strains_to_fasta_file = {}

        for x in tqdm(ids):
            f = ssbio.databases.patric.download_genome_sequence(patric_id=x, seqtype=self.seq_type,
                                                                outdir=self.seq_atlas_org_dir,
                                                                force_rerun=force_rerun)
            if f:
                strains_to_fasta_file[x] = f
                log.debug('{}: Downloaded sequence'.format(x))
            else:
                log.warning('{}: Unable to download sequence'.format(x))

        if len(strains_to_fasta_file) > 0:
            log.info('Downloaded sequences of coding genes.')

        return strains_to_fasta_file

    def download_genome_cds_ncbi(self, ids, email, force_rerun=False):
        """Loads a list of NCBI GIs or RefSeq complete genome IDs and downloads the CDS FASTA files

        Args:
            ids (list): List of strain_genome_ids
            email (str): your email

        """
        strains_to_fasta_file = {}

        for x in tqdm(ids):
            f = ssbio.databases.ncbi.download_genome_sequence(genome_accession_or_id=x, seqtype=self.seq_type,
                                                              email=email,
                                                              outdir=self.seq_atlas_org_dir, force_rerun=force_rerun)
            strains_to_fasta_file[x] = f
            log.debug('Downloaded sequence for {}'.format(x))

        if len(strains_to_fasta_file) > 0:
            log.info('Downloaded sequences of coding genes.')

        return strains_to_fasta_file

    @property
    def reference_genome(self):
        return self._reference_genome

    @reference_genome.setter
    def reference_genome(self, reference_genome_id):
        if reference_genome_id in self.strain_ids:
            self._reference_genome = reference_genome_id
        else:
            raise ValueError('Reference genome ID not in list of genomes.')
        log.info('Set reference genome to ID {}'.format(reference_genome_id))

    def get_orthology_matrix(self):
        """Run run_makeblastdb, run_bidirectional_blast, and calculate_bbh for both DNA and protein sequences.

        """
        # TODO: probably best to split this function up to allow for force_rerunning of different parts
        # Get the path to the reference genome
        r_file = self.base_strain_gempro_genome_path
        r_folder, r_name, r_ext = utils.split_folder_and_path(r_file)

        bbh_files = {}

        for strain_model in tqdm(self.strain_models):
            g_file = strain_model.get_genome_file_path(self.seq_atlas_org_dir)
            # TODO: replace usages of g_name with strain_model.id instead for consistency
            g_folder, g_name, g_ext = utils.split_folder_and_path(g_file)

            # Run bidirectional BLAST
            log.debug('{} vs {}: Running bidirectional BLAST'.format(r_name, g_name))
            r_vs_g, g_vs_r = ssbio.sequence.utils.blast.run_bidirectional_blast(reference=r_file, other_genome=g_file,
                                                                                dbtype=self.blast_seq_type,
                                                                                outdir=self.seq_atlas_org_dir)

            # Using the BLAST files, find the BBH
            log.debug('{} vs {}: Finding BBHs'.format(r_name, g_name))
            bbh = ssbio.sequence.utils.blast.calculate_bbh(blast_results_1=r_vs_g, blast_results_2=g_vs_r,
                                                           r_name=r_name, g_name=g_name,
                                                           outdir=self.seq_atlas_org_dir)
            bbh_files[strain_model.id] = bbh

        # Make the orthologous genes matrix
        log.debug('Creating orthology matrix')
        ortho_matrix = ssbio.sequence.utils.blast.create_orthology_matrix(r_name=r_name,
                                                                          genome_to_bbh_files=bbh_files,
                                                                          outname='{}_{}_orthology.csv'.format(r_name, self.blast_seq_type),
                                                                          outdir=self.atlas_data_dir)

        log.info('Saved orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))
        self.df_orthology_matrix = pd.read_csv(ortho_matrix, index_col=0)

        # Remove extraneous strains from our analysis
        to_remove = []
        for i, strain_model in enumerate(self.strain_models):
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
                self.strain_models.pop(x)
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

        log.info('{} strains to be analyzed'.format(len(self.strain_models)))

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
            outfile = op.join(self.atlas_model_files, '{}.json'.format(strain_id))
            cobra.io.save_json_model(strain_model, outfile)
            log.debug('{}: saved model at {}'.format(strain_id, outfile))

        log.info('Created {} new strain-specific models'.format(len(self.strain_models)))

    def write_orthologous_gene_sequences(self):
        """For each organism, write their orthologous gene files named like <STRAIN>_<BASESTRAIN_GENE>.faa

        """
        if len(self.df_orthology_matrix_filtered) == 0:
            raise RuntimeError('Empty orthology matrix')

        for strain_model in tqdm(self.strain_models):
            strain_id = strain_model.id

            # Load the strain genome file
            strain_fasta = strain_model.get_genome_file_path(self.seq_atlas_org_dir)
            strain_sequences = SeqIO.index(strain_fasta, 'fasta')
            log.debug('Loaded {}'.format(strain_fasta))

            # Get the list of orthologous genes, ignoring genes outside our context of the base model
            base_to_strain = self.df_orthology_matrix_filtered[pd.notnull(self.df_orthology_matrix_filtered[strain_id])][strain_id].to_dict()
            for base_g_id, strain_g_id in base_to_strain.items():

                # Make the gene directory
                gene_dir = op.join(self.seq_atlas_gene_dir, base_g_id)
                if not op.exists(gene_dir):
                    os.mkdir(gene_dir)

                # Get and write the strain sequence
                strain_seq_record = strain_sequences[strain_g_id]
                outfile = '{}_{}.{}'.format(base_g_id, strain_id, self.fasta_extension)

                # Save the filename in the strain model's gene annotation
                strain_model.genes.get_by_id(base_g_id).annotation['sequence']['sequence_file'] = outfile
                strain_model.genes.get_by_id(base_g_id).annotation['sequence']['seq_len'] = len(strain_seq_record.seq)

                if not op.exists(op.join(gene_dir, outfile)):
                    with open(op.join(gene_dir, outfile), 'w') as f:
                        SeqIO.write(strain_seq_record, f, "fasta")

            log.debug('Wrote all sequences for strain {}'.format(strain_id))

        log.info('Wrote all individual orthologous gene sequences for all strains.')

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

            gene_dir = op.join(self.seq_atlas_gene_dir, base_gene_id)

            info_dict = {}
            info_dict['gene'] = base_gene_id
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

    def align_orthologous_genes_multiple(self):
        """For each gene in the base strain, run a multiple alignment to all orthologous strain genes

        Returns:

        """
        pass

    def per_gene_total_num_mutations(self, gene_id):
        """Simply return the total number of mutations for a gene, summed over all strains.

        Returns:
            int: Total number of point mutations that show up in a gene throughout all strains

        """
        total_num_mutations = 0

        gene_dir = op.join(self.seq_atlas_gene_dir, gene_id)

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


    def predict_auxotrophies(self):
        """Find what nutrients must be added to allow them to grow in minimal media

        """
        pass
        # for strain_id, strain_model in self.genome_id_to_strain_model.items():
        #     model = ssbio.cobra.utils.model_loader(gem_file_path=strain_model, gem_file_type='json')
        #
        #     # Optimize the strain model for growth
        #     model.optimize()
        #
        #     if model.solution.f < 0.00001:
        #         log.info('Gap filling {}'.format(strain_id))
        #         # Run GrowMatch with only exchange reactions
        #         solution = cobra.flux_analysis.growMatch(model, dm_rxns=False, ex_rxns=True, iterations=1)
        #
        #         # Open the exchange reactions determined using GrowMatch
        #         for rxn in solution[0]:
        #             rxn = rxn.id.replace('_reverse', '')
        #             model.reactions[model.reactions.index(rxn)].lower_bound = -1

