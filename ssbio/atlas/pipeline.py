import os
import os.path as op
import pandas as pd
import ssbio.cobra.utils
import ssbio.sequence.blast
import ssbio.sequence.fasta
import ssbio.databases.ncbi

import cobra.manipulation
import cobra.flux_analysis

from Bio import SeqIO
from Bio import Entrez
from tqdm import tqdm

from ssbio import utils

date = utils.Date()

import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
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

    def __init__(self, base_strain_name, root_dir, seq_type='prot', base_gem_file_path=None, base_gem_file_type=None,
                 reference_genome=None):
        """Prepare for ATLAS analysis.

        Args:
            base_strain_name (str): Name of your base strain, and folder which will be created
            root_dir (str): Path to folder in which base_strain_name folder will be created
            seq_type (str): Whether to run analysis on "nucl" (DNA) or "prot" (amino acid) sequences
        """

        self._reference_genome = reference_genome
        self.seq_type = seq_type

        # Load the base strain model
        if base_gem_file_path and base_gem_file_type:
            self.model = ssbio.cobra.utils.model_loader(base_gem_file_path, base_gem_file_type)
            log.info('Loaded COBRA model from {}'.format(base_gem_file_path))
            # Log information on the number of things
            log.info('Number of reactions: {}'.format(len(self.model.reactions)))
            log.info(
                    'Number of reactions linked to a gene: {}'.format(ssbio.cobra.utils.true_num_reactions(self.model)))
            log.info(
                    'Number of genes (excluding spontaneous): {}'.format(ssbio.cobra.utils.true_num_genes(self.model)))
            log.info('Number of metabolites: {}'.format(len(self.model.metabolites)))

        list_of_dirs = []
        self.root_dir = root_dir
        self.base_dir = op.join(self.root_dir, base_strain_name)
        # notebooks_dir - directory where ipython notebooks will be stored for manual analyses
        self.notebooks_dir = op.join(self.base_dir, 'notebooks')
        # atlas_dir - directory where all ATLAS analysis will be carried out
        self.atlas_dir = op.join(self.base_dir, 'atlas')
        # model_files - directory where base strain GEM and new models will be stored
        self.model_files = op.join(self.atlas_dir, 'model_files')
        # data_dir - directory where all data_dir will be stored
        self.data_dir = op.join(self.atlas_dir, 'data')
        # figures - directory where all figures will be stored
        self.figures_dir = op.join(self.atlas_dir, 'figures')
        # seq_dir - sequence related files are stored here
        self.seq_dir = op.join(self.atlas_dir, 'sequence_files')

        list_of_dirs.extend(
                [self.base_dir, self.notebooks_dir, self.atlas_dir, self.data_dir, self.figures_dir, self.model_files,
                 self.seq_dir])

        if seq_type == 'nucl':
            self.seq_atlas_dir = op.join(self.seq_dir, 'dna')
            self.seq_atlas_org_dir = op.join(self.seq_atlas_dir, 'by_organism')
            self.seq_atlas_gene_dir = op.join(self.seq_atlas_dir, 'by_gene')
            self.fasta_extension = 'fna'

        if seq_type == 'prot':
            self.seq_atlas_dir = op.join(self.seq_dir, 'protein')
            self.seq_atlas_org_dir = op.join(self.seq_atlas_dir, 'by_organism')
            self.seq_atlas_gene_dir = op.join(self.seq_atlas_dir, 'by_gene')
            self.fasta_extension = 'faa'

        list_of_dirs.extend([self.seq_atlas_dir, self.seq_atlas_org_dir, self.seq_atlas_gene_dir])

        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('Created directory: {}'.format(directory))
            else:
                log.debug('Directory already exists: {}'.format(directory))

        self.orthology_matrix = pd.DataFrame()
        self.genome_id_to_fasta_file = {}
        self.genome_id_to_strain_model = {}

    def download_genome_cds(self, ids, email, force_rerun=False):
        """Loads a list of NCBI GIs or RefSeq complete genome IDs and downloads the CDS FASTA files

        Args:
            ids (list): List of strain_genome_ids
            email (str): your email

        """

        for x in tqdm(ids):
            f = ssbio.databases.ncbi.download_genome_sequence(genome_accession_or_id=x, seqtype=self.seq_type,
                                                              email=email,
                                                              outdir=self.seq_atlas_org_dir, force_rerun=force_rerun)
            self.genome_id_to_fasta_file[x] = f
            log.debug('Downloaded sequence for {}'.format(x))

        if len(self.genome_id_to_fasta_file) > 0:
            log.info('Downloaded sequences of coding genes.')

    @property
    def reference_genome(self):
        return self._reference_genome

    @reference_genome.setter
    def reference_genome(self, reference_genome_id):
        if reference_genome_id in self.genome_id_to_fasta_file.keys():
            self._reference_genome = reference_genome_id
        else:
            raise ValueError('Reference genome ID not in list of genomes.')
        log.info('Set reference genome to ID {}'.format(reference_genome_id))

    def find_bbh(self):
        """Run run_makeblastdb, run_bidirectional_blast, and calculate_bbh for both DNA and protein sequences.

        """
        genomes = list(self.genome_id_to_fasta_file.values())
        r_file = self.genome_id_to_fasta_file[self.reference_genome]

        bbh_files = {}

        for g_file in tqdm(genomes):
            g_folder, g_name, g_ext = utils.split_folder_and_path(g_file)

            # Run bidirectional BLAST
            b1, b2 = ssbio.sequence.blast.run_bidirectional_blast(reference=r_file, other_genome=g_file,
                                                                  dbtype=self.seq_type, outdir=self.seq_atlas_org_dir)

            # Using the BLAST files, find the BBH
            bbh = ssbio.sequence.blast.calculate_bbh(blast_results_1=b1, blast_results_2=b2,
                                                     r_name=self.reference_genome, g_name=g_name,
                                                     outdir=self.seq_atlas_org_dir)
            bbh_files[g_name] = bbh

        # Make the orthologous genes matrix
        ortho_matrix = ssbio.sequence.blast.create_orthology_matrix(r_name=self.reference_genome,
                                                                    genome_to_bbh_files=bbh_files,
                                                                    outname='{}_{}_orthology.csv'.format(
                                                                        self.reference_genome, self.seq_type),
                                                                    outdir=self.data_dir)

        log.info('Saved orthology matrix at {}'.format(ortho_matrix))

        self.orthology_matrix = pd.read_csv(ortho_matrix, index_col=0)

    def write_orthologous_gene_sequences(self):
        """For each organism, write their orthologous gene files named like <STRAIN>_<BASESTRAIN_GENE>.faa

        """
        if len(self.orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        genomes = self.genome_id_to_fasta_file

        for strain_id, fasta_file in tqdm(genomes.items()):
            if strain_id not in self.orthology_matrix.columns:
                log.warning('No orthologous genes found for organism {}'.format(strain_id))
                continue

            # Load the strain FASTA file
            strain_sequences = SeqIO.index(fasta_file, 'fasta')
            log.debug('Loaded {}'.format(fasta_file))

            # Get the list of orthologous genes
            base_to_strain = self.orthology_matrix[pd.notnull(self.orthology_matrix[strain_id])][strain_id].to_dict()
            for base_g_id, strain_g_id in base_to_strain.items():

                # Make the gene directory
                gene_dir = op.join(self.seq_atlas_gene_dir, base_g_id)
                if not op.exists(gene_dir):
                    os.mkdir(gene_dir)

                # Get and write the strain sequence
                strain_seq_record = strain_sequences[strain_g_id]
                outfile = '{}_{}.{}'.format(base_g_id, strain_id, self.fasta_extension)
                if not op.exists(op.join(gene_dir, outfile)):
                    with open(op.join(gene_dir, outfile), 'w') as f:
                        SeqIO.write(strain_seq_record, f, "fasta")

            log.debug('Wrote all sequences for strain {}'.format(strain_id))

        log.info('Wrote all individual gene sequences for all strains.')

    def build_strain_specific_models(self):
        """Using the orthologous genes matrix, write strain specific models.

        """
        if len(self.orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        if not hasattr(self, 'model'):
            raise RuntimeError('No GEM loaded')

        # For each genome, create the strain specific model
        for strain_id in self.genome_id_to_fasta_file.keys():
            # Get a list of genes which do not have orthology in the strain
            not_in_strain = self.orthology_matrix[pd.isnull(self.orthology_matrix[strain_id])][strain_id].index.tolist()

            # Make a copy of the base strain
            my_new_strain_model = self.model.copy()
            my_new_strain_model._trimmed = False
            my_new_strain_model._trimmed_genes = []
            my_new_strain_model._trimmed_reactions = {}

            # Filter out genes which do not show up in the base strain model
            model_genes = [x.id for x in self.model.genes]
            genes_to_remove = list(set(not_in_strain).intersection(set(model_genes)))

            if len(genes_to_remove) == 0:
                log.debug('No genes marked for removal from base strain')
            else:
                log.debug('{} genes marked for removal from base strain'.format(len(genes_to_remove)))

                # Change the model's name to correspond with the strain
                my_new_strain_model.id = strain_id

                # TODO: allow input of model name as well (use a name:id mapping dict or something)
                # my_new_strain_model.name = strain_name

                # Delete genes!
                cobra.manipulation.delete_model_genes(my_new_strain_model, genes_to_remove)

                if my_new_strain_model._trimmed == True:
                    log.info('{}: deleted {} reactions, {} genes'.format(strain_id,
                                                                         len(my_new_strain_model._trimmed_reactions),
                                                                         len(my_new_strain_model._trimmed_genes)))

            # Save the strain specific file
            outfile = op.join(self.model_files, '{}.json'.format(strain_id))
            cobra.io.save_json_model(my_new_strain_model, outfile)
            self.genome_id_to_strain_model[strain_id] = outfile
            log.debug('{}: saved model at {}'.format(strain_id, outfile))

            del my_new_strain_model

    def predict_auxotrophies(self):
        """Find what nutrients must be added to allow them to grow in minimal media

        """
        for strain_id, strain_model in self.genome_id_to_strain_model.items():
            model = ssbio.cobra.utils.model_loader(gem_file_path=strain_model, gem_file_type='json')

            # Optimize the strain model for growth
            model.optimize()

            if model.solution.f < 0.00001:
                log.info('Gap filling {}'.format(strain_id))
                # Run GrowMatch with only exchange reactions
                solution = cobra.flux_analysis.growMatch(model, dm_rxns=False, ex_rxns=True, iterations=1)

                # Open the exchange reactions determined using GrowMatch
                for rxn in solution[0]:
                    rxn = rxn.id.replace('_reverse', '')
                    model.reactions[model.reactions.index(rxn)].lower_bound = -1

