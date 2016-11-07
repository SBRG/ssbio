import os
import os.path as op
import pandas as pd
import ssbio.cobra.utils
import ssbio.sequence.blast
import ssbio.sequence.fasta
import ssbio.databases.ncbi
import ssbio.databases.patric

import cobra.manipulation
import cobra.flux_analysis

from Bio import SeqIO
from Bio import Entrez

try:
    from IPython.display import clear_output
    have_ipython = True
    from tqdm import tqdm_notebook as tqdm
except ImportError:
    have_ipython = False
    from tqdm import tqdm

from ssbio import utils
import ssbio.sequence.alignment
from collections import defaultdict
date = utils.Date()
from ssbio.pipeline.gempro import GEMPRO
from cobra.core import DictList
from cobra.core import Model
import sys
import logging
import copy
import shutil

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
                                        'seq_file'      : None,
                                        'pairwise_alignment_file': None,
                                        'properties'    : {}
                                        }
            g.annotation['structure'] = {'base_structure': g.annotation['structure']['representative']}

    def get_genome_file_path(self, genome_file_dir):
        return op.join(genome_file_dir, self.annotation['genome_file'])


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
        log.info('Creating strain specific models based on the base model...')
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

        # Create the orthology matrix
        self.df_orthology_matrix = pd.DataFrame()

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
        # Get the path to the reference genome
        r_file = self.base_strain_gempro_genome_path
        r_folder, r_name, r_ext = utils.split_folder_and_path(r_file)

        bbh_files = {}

        for strain_model in self.strain_models:
            g_file = strain_model.get_genome_file_path(self.seq_atlas_org_dir)
            g_folder, g_name, g_ext = utils.split_folder_and_path(g_file)

            # Run bidirectional BLAST
            b1, b2 = ssbio.sequence.blast.run_bidirectional_blast(reference=r_file, other_genome=g_file,
                                                                  dbtype=self.blast_seq_type, outdir=self.seq_atlas_org_dir)

            # Using the BLAST files, find the BBH
            bbh = ssbio.sequence.blast.calculate_bbh(blast_results_1=b1, blast_results_2=b2,
                                                     r_name=r_name, g_name=g_name,
                                                     outdir=self.seq_atlas_org_dir)
            bbh_files[g_name] = bbh

        # Make the orthologous genes matrix
        ortho_matrix = ssbio.sequence.blast.create_orthology_matrix(r_name=r_name,
                                                                    genome_to_bbh_files=bbh_files,
                                                                    outname='{}_{}_orthology.csv'.format(r_name, self.blast_seq_type),
                                                                    outdir=self.atlas_data_dir)

        log.info('Saved orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))

        # TODO: check what the column names are!!
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

        # Remove strains with no differences
        for x in to_remove:
            self.strain_models.pop(x)

        log.info('Removed {} strains from analysis'.format(len(to_remove)))
        log.info('{} strains remain in analysis'.format(len(self.strain_models)))

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
        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        # For each genome, create the strain specific model
        for strain_model in tqdm(self.strain_models):
            strain_id = strain_model.id

            # Get a list of genes which do not have orthology in the strain
            not_in_strain = self.df_orthology_matrix[pd.isnull(self.df_orthology_matrix[strain_id])][strain_id].index.tolist()

            self.pare_down_model(model_to_be_modified=strain_model, genes_to_remove=not_in_strain)

            # Save the strain specific file
            outfile = op.join(self.atlas_model_files, '{}.json'.format(strain_id))
            cobra.io.save_json_model(strain_model, outfile)
            log.debug('{}: saved model at {}'.format(strain_id, outfile))

        log.info('Created {} new strain-specific models'.format(len(self.strain_models)))

    def write_orthologous_gene_sequences(self):
        """For each organism, write their orthologous gene files named like <STRAIN>_<BASESTRAIN_GENE>.faa

        """
        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        for strain_model in tqdm(self.strain_models):
            strain_id = strain_model.id

            # Load the strain genome file
            strain_fasta = strain_model.get_genome_file_path(self.seq_atlas_org_dir)
            strain_sequences = SeqIO.index(strain_fasta, 'fasta')
            log.debug('Loaded {}'.format(strain_fasta))

            # Get the list of orthologous genes
            base_to_strain = self.df_orthology_matrix[pd.notnull(self.df_orthology_matrix[strain_id])][strain_id].to_dict()
            for base_g_id, strain_g_id in base_to_strain.items():

                # Make the gene directory
                gene_dir = op.join(self.seq_atlas_gene_dir, base_g_id)
                if not op.exists(gene_dir):
                    os.mkdir(gene_dir)

                # Get and write the strain sequence
                strain_seq_record = strain_sequences[strain_g_id]
                outfile = '{}_{}.{}'.format(base_g_id, strain_id, self.fasta_extension)

                # Save the filename in the strain model's gene annotation
                strain_model.genes.get_by_id(strain_g_id).annotation['seq_file'] = outfile

                if not op.exists(op.join(gene_dir, outfile)):
                    with open(op.join(gene_dir, outfile), 'w') as f:
                        SeqIO.write(strain_seq_record, f, "fasta")

            log.debug('Wrote all sequences for strain {}'.format(strain_id))

        log.info('Wrote all individual orthologous gene sequences for all strains.')

    def align_orthologous_genes_pairwise(self):
        """For each gene in the base strain, run a pairwise alignment for all orthologous gene sequences to it.

        """
        pass
        # for base_gene in self.base_strain_gempro.genes:
        #     # Get base strain gene fasta file path
        #     base_gene_id = base_gene.id
        #     base_gene_seq_file = base_gene.annotation['sequence']['representative']['seq_file']
        #     base_gene_seq_path = op.join(self.base_strain_gempro.sequence_dir, base_gene_id, base_gene_seq_file)
        #
        #     # Get gene file in all strains if it shows up as functional
        #     for strain_model in self.strain_models:
        #         strain_id = strain_model.id
        #         strain_gene = strain_model.get_by_id(base_gene_id)
        #
        #         if not strain_gene.functional:
        #             log.debug('{}: gene not functional in strain {}, not aligning'.format(base_gene_id, strain_id))
        #             continue
        #
        #         strain_gene_seq_file = strain_gene.annotation['sequence']['seq_file']
        #         ssbio.sequence.alignment.run_needle_alignment_on_files()

    def align_orthologous_genes_multiple(self):
        """For each gene in the base strain, run a

        Returns:

        """
        pass



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

