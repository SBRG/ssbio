import os
import os.path as op
import cobra.flux_analysis
import cobra.manipulation
import pandas as pd
import ssbio.sequence.utils.blast
import ssbio.cobra.utils
import ssbio.databases.ncbi
import ssbio.databases.patric
import ssbio.sequence.utils.fasta
from ssbio.core.object import Object
from ssbio import utils
import ssbio.sequence.utils.alignment
from ssbio.pipeline.gempro import GEMPRO
from Bio import SeqIO
import sys
import logging
from slugify import Slugify
from cobra.core import DictList
from cobra.core import Gene
from ssbio.core.genepro import GenePro
import numpy as np
import ssbio.sequence.properties.residues
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

date = utils.Date()
custom_slugify = Slugify(safe_chars='-_.')
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


class ATLAS(Object):
    """Class to represent an ATLAS workflow to carry out multi-strain comparisons

    Main steps are:
    1. Strain-specific model construction based on orthologous genes & systems modeling
    2. Phylogenetic analysis to pick out important genes
    3. GEM-PRO of the "base strain"
    4. Structure property calculation & integrated structural systems analysis

    Each step may generate a report and also request additional files if something is missing
    """

    def __init__(self, atlas_name, base_gempro, base_genome_path=None, description=None):
        """Prepare a GEM-PRO model for ATLAS analysis

        Args:
            atlas_name (str): Name of your ATLAS project
            base_gempro (GEMPRO): Completed GEM-PRO model
            base_genome_path (str): Simple reference link to the genome FASTA file (CDS)
            description (str): Optional string to describe your project
        """
        Object.__init__(self, id=atlas_name, description=description)
        self.atlas_strains = DictList()

        # Load the GEM-PRO (could be a model, could just be a list of genes)
        self.base_strain_gempro = base_gempro

        # Check if there is a genome file associated with this model - if not, write all sequences and use that
        if not base_genome_path and not self.base_strain_gempro.genome_path:
            self.base_strain_gempro.genome_path = self.base_strain_gempro.write_representative_sequences_file(outname='BASE_CDS')
        else:
            self.base_strain_gempro.genome_path = base_genome_path
            # TODO: must check if base_genome gene IDs can be matched to the base model

        # Prepare ATLAS directories
        self.atlas_dir = None
        self.atlas_model_dir = None
        self.atlas_data_dir = None
        self.atlas_sequence_dir = None
        self.atlas_seq_genomes_dir = None
        self.atlas_seq_genes_dir = None
        self._create_dirs(self.base_strain_gempro.base_dir)

        # Other attributes
        self.df_orthology_matrix = pd.DataFrame()
        # Mark if the orthology matrix has gene IDs (thus we need to retrieve seqs from the genome file) or if
        # it is in the orthology matrix itself
        self._orthology_matrix_has_sequences = False

    def _create_dirs(self, root_dir):
        """Internal method to create ATLAS directories.

        Args:
            root_dir (str): Root directory to create a single new ATLAS project folder in

        """
        list_of_dirs = []

        if not root_dir:
            root_dir = os.getcwd()

        if not op.exists(root_dir):
            raise ValueError('{}: folder does not exist'.format(root_dir))

        # atlas_dir - directory where all ATLAS analysis will be carried out
        self.atlas_dir = op.join(root_dir, 'atlas')
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

        log.info('{}: ATLAS project files location'.format(self.atlas_dir))

    def load_manual_orthology_matrix(self, df, clean_names=True,
                                     remove_strains_with_no_orthology=True,
                                     remove_strains_with_no_differences=False,
                                     remove_genes_not_in_base_model=True):
        """Load a manually curated orthology matrix to use in ATLAS. Genes = rows, strains = columns.

        Args:
            df (DataFrame): Pandas DataFrame with genes as the rows and strains as the columns
            clean_names (bool): Remove unwanted characters from gene names and strain IDs
            remove_strains_with_no_orthology (bool): Remove strains which have no orthologous genes found
            remove_strains_with_no_differences (bool): Remove strains which have all the same genes as the base model.
                Default is False because since orthology is found using a PID cutoff, all genes may be present but
                differences may be on the sequence level.
            remove_genes_not_in_base_model (bool): Remove genes from the orthology matrix which are not present in our
                base model. This happens if we use a genome file for our model that has other genes in it.

        """
        self._orthology_matrix_has_sequences = True

        if clean_names:
            new_rows = [custom_slugify(x) for x in df.index]
            new_cols = [custom_slugify(y) for y in df.columns]
            df.index = new_rows
            df.columns = new_cols

        self.df_orthology_matrix = df

        # Make the copies of the base model
        for strain_id in tqdm(self.df_orthology_matrix.columns):
            self._copy_base_gempro(new_id=strain_id)

        # Filter the strains and orthology matrix
        self._filter_with_orthology_matrix(remove_strains_with_no_orthology=remove_strains_with_no_orthology,
                                           remove_strains_with_no_differences=remove_strains_with_no_differences,
                                           remove_genes_not_in_base_model=remove_genes_not_in_base_model)

    def _copy_base_gempro(self, new_id):
        """Copy the base strain GEM-PRO into a new GEM-PRO with a specified ID.

        Appends the model to the atlas_strains attribute.

        Args:
            new_id (str): New ID to be assigned to the copied model

        Returns:
            GEMPRO: copied GEM-PRO to represent the new strain

        """
        # Temporarily disable logging messages
        logging.disable(logging.WARNING)
        if self.base_strain_gempro.model:
            # If there is a SBML model associated with the GEMPRO, copy that model
            copied_model = GEMPRO(gem_name=new_id, create_dirs=False, gem=self.base_strain_gempro.model.copy())
            copied_model.model.id = new_id

            # Reset the GenePro attributes
            for x in copied_model.genes:
                x.reset_protein()
        else:
            # Otherwise, just copy the list of genes over and rename the IDs
            strain_genes = [x.id for x in self.base_strain_gempro.genes]
            copied_model = GEMPRO(gem_name=new_id, create_dirs=False, genes_list=strain_genes)
        # Re-enable logging
        logging.disable(logging.NOTSET)

        self.atlas_strains.append(copied_model)
        log.debug('{}: new model copied from base model'.format(new_id))

        return self.atlas_strains.get_by_id(new_id)

    def load_strain_ids_and_genomes(self, ids_to_genome_file):
        """Load a dictionary of PATRIC IDs which point to the path of their respective genome files.

        Creates initial copies of the base strain model for each strain ID and stores them in self.atlas_strains

        Args:
            ids_to_genome_file (dict): Keys are PATRIC IDs (which will become your strain IDs) and
                values are absolute paths to the FASTA file containing CDS regions

        """
        for strain_id, genome_path in tqdm(ids_to_genome_file.items()):
            strain_gp = self._copy_base_gempro(new_id=strain_id)
            strain_gp.genome_path = genome_path

    def download_genome_cds_patric(self, ids, force_rerun=False):
        """Download genome files from PATRIC give a list of IDs

        Args:
            ids (str, list): PATRIC ID or list of PATRIC IDs
            force_rerun: If genome files should be downloaded again even if they exist

        """
        strains_to_fasta_file = {}
        ids = ssbio.utils.force_list(ids)

        log.info('Downloading sequences from PATRIC...')
        for patric_id in tqdm(ids):
            f = ssbio.databases.patric.download_coding_sequences(patric_id=patric_id, seqtype='protein',
                                                                 outdir=self.atlas_seq_genomes_dir,
                                                                 force_rerun=force_rerun)
            if f:
                strains_to_fasta_file[patric_id] = f
                log.debug('{}: downloaded sequence'.format(patric_id))
            else:
                log.warning('{}: unable to download sequence'.format(patric_id))

        log.info('Copying base strain to create initial strain specific models...')
        self.load_strain_ids_and_genomes(strains_to_fasta_file)
        log.info('Created {} new strain GEM-PROs, accessible at "atlas_strains"'.format(len(strains_to_fasta_file)))

    def get_orthology_matrix(self, pid_cutoff=None, bitscore_cutoff=None, evalue_cutoff=None, filter_condition='OR',
                             remove_strains_with_no_orthology=True,
                             remove_strains_with_no_differences=False,
                             remove_genes_not_in_base_model=True):
        """Create the orthology matrix by finding best bidirectional BLAST hits. Genes = rows, strains = columns

        Runs run_makeblastdb, run_bidirectional_blast, and calculate_bbh for protein sequences.

        Args:
            pid_cutoff (float): Minimum percent identity between BLAST hits to filter for in the range [0, 100]
            bitscore_cutoff (float): Minimum bitscore allowed between BLAST hits
            evalue_cutoff (float): Maximum E-value allowed between BLAST hits
            filter_condition (str): 'OR' or 'AND', how to combine cutoff filters. 'OR' gives more results since it
            is less stringent, as you will be filtering for hits with (>80% PID or >30 bitscore or <0.0001 evalue).
            remove_strains_with_no_orthology (bool): Remove strains which have no orthologous genes found
            remove_strains_with_no_differences (bool): Remove strains which have all the same genes as the base model.
                Default is False because since orthology is found using a PID cutoff, all genes may be present but
                differences may be on the sequence level.
            remove_genes_not_in_base_model (bool): Remove genes from the orthology matrix which are not present in our
                base model. This happens if we use a genome file for our model that has other genes in it.

        Returns:
            DataFrame: Orthology matrix calculated from best bidirectional BLAST hits.

        """
        # TODO: document and test other cutoffs

        # Get the path to the reference genome
        r_file = self.base_strain_gempro.genome_path

        bbh_files = {}

        for strain_model in tqdm(self.atlas_strains):
            g_file = strain_model.genome_path

            # Run bidirectional BLAST
            log.debug('{} vs {}: Running bidirectional BLAST'.format(self.base_strain_gempro.id, strain_model.id))
            r_vs_g, g_vs_r = ssbio.sequence.utils.blast.run_bidirectional_blast(reference=r_file, other_genome=g_file,
                                                                                dbtype='prot',
                                                                                outdir=self.atlas_seq_genomes_dir)

            # Using the BLAST files, find the BBH
            log.debug('{} vs {}: Finding BBHs'.format(self.base_strain_gempro.id, strain_model.id))
            bbh = ssbio.sequence.utils.blast.calculate_bbh(blast_results_1=r_vs_g, blast_results_2=g_vs_r,
                                                           outdir=self.atlas_seq_genomes_dir)
            bbh_files[strain_model.id] = bbh

        # Make the orthologous genes matrix
        log.debug('Creating orthology matrix')
        ortho_matrix = ssbio.sequence.utils.blast.create_orthology_matrix(r_name=self.base_strain_gempro.id,
                                                                          genome_to_bbh_files=bbh_files,
                                                                          pid_cutoff=pid_cutoff,
                                                                          bitscore_cutoff=bitscore_cutoff,
                                                                          evalue_cutoff=evalue_cutoff,
                                                                          filter_condition=filter_condition,
                                                                          outname='{}_{}_orthology.csv'.format(self.base_strain_gempro.id, 'prot'),
                                                                          outdir=self.atlas_data_dir)

        log.info('Saved orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))
        self.df_orthology_matrix = pd.read_csv(ortho_matrix, index_col=0)

        # Filter the matrix to genes only in our analysis, and also check for strains with no differences or no orthologous genes
        self._filter_with_orthology_matrix(remove_strains_with_no_orthology=remove_strains_with_no_orthology,
                                           remove_strains_with_no_differences=remove_strains_with_no_differences,
                                           remove_genes_not_in_base_model=remove_genes_not_in_base_model)

        return self.df_orthology_matrix

    def _filter_with_orthology_matrix(self,
                                      remove_strains_with_no_orthology=True,
                                      remove_strains_with_no_differences=False,
                                      remove_genes_not_in_base_model=True):
        """Filters the orthology matrix by removing genes not in our base model, and also
            removes strains from the analysis which have: 0 orthologous genes or no difference from the base strain.

        Args:
            remove_strains_with_no_orthology (bool): Remove strains which have no orthologous genes found
            remove_strains_with_no_differences (bool): Remove strains which have all the same genes as the base model.
                Default is False because since orthology is found using a PID cutoff, all genes may be present but
                differences may be on the sequence level.
            remove_genes_not_in_base_model (bool): Remove genes from the orthology matrix which are not present in our
                base model. This happens if we use a genome file for our model that has other genes in it.

        """

        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        initial_num_strains = len(self.atlas_strains)

        # Adding names to the row and column of the orthology matrix
        self.df_orthology_matrix = self.df_orthology_matrix.rename_axis('gene').rename_axis("strain", axis="columns")

        # Gene filtering (of the orthology matrix)
        if remove_genes_not_in_base_model:
            # Check for gene IDs that are in the model and not in the orthology matrix
            # This is probably because: the CDS FASTA file for the base strain did not contain the correct ID
            # for the gene and consequently was not included in the orthology matrix
            # Save these and report them
            base_strain_gene_ids = [x.id for x in self.base_strain_gempro.genes]
            self.missing_in_orthology_matrix = [x for x in base_strain_gene_ids if x not in self.df_orthology_matrix.index.tolist()]
            self.missing_in_base_strain = [y for y in self.df_orthology_matrix.index.tolist() if y not in base_strain_gene_ids]

            # Filter the matrix for genes within our base model only
            self.df_orthology_matrix = self.df_orthology_matrix[self.df_orthology_matrix.index.map(lambda x: x in base_strain_gene_ids)]
            log.info('Filtered orthology matrix for genes present in base model')
            log.warning('{} genes are in your base model but not your orthology matrix, see the attribute "missing_in_orthology_matrix"'.format(len(self.missing_in_orthology_matrix)))
            log.warning('{} genes are in the orthology matrix but not your base model, see the attribute "missing_in_base_strain"'.format(len(self.missing_in_base_strain)))

        # Strain filtering
        for strain_model in self.atlas_strains.copy():
            if remove_strains_with_no_orthology:
                if strain_model.id not in self.df_orthology_matrix.columns:
                    self.atlas_strains.remove(strain_model)
                    log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_model.id))
                    continue
                elif self.df_orthology_matrix[strain_model.id].isnull().all():
                    self.atlas_strains.remove(strain_model)
                    log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_model.id))
                    continue

            if remove_strains_with_no_differences:
                not_in_strain = self.df_orthology_matrix[pd.isnull(self.df_orthology_matrix[strain_model.id])][strain_model.id].index.tolist()
                if len(not_in_strain) == 0:
                    self.atlas_strains.remove(strain_model)
                    log.info('{}: strain has no differences from the base, removed from analysis.')
                    continue

        log.info('{} strains to be analyzed, {} strains removed'.format(len(self.atlas_strains), initial_num_strains-len(self.atlas_strains)))

    def _pare_down_model(self, gempro, genes_to_remove):
        """Mark genes as non-functional in a GEM-PRO. If there is a COBRApy model associated with it, the
            COBRApy method delete_model_genes is utilized to delete genes.

        Args:
            gempro (GEMPRO): GEMPRO object
            genes_to_remove (list): List of gene IDs to remove from the model

        """
        # Filter out genes in genes_to_remove which do not show up in the model
        model_genes = [x.id for x in gempro.genes]
        genes_to_remove = list(set(genes_to_remove).intersection(set(model_genes)))
        genes_to_remove.extend(self.missing_in_orthology_matrix)

        if len(genes_to_remove) == 0:
            log.info('{}: no genes marked for removal'.format(gempro.id))
            return
        else:
            log.debug('{}: {} genes marked for removal'.format(gempro.id, len(genes_to_remove)))

        # If a COBRApy model exists, utilize the delete_model_genes method
        if gempro.model:
            gempro.model._trimmed = False
            gempro.model._trimmed_genes = []
            gempro.model._trimmed_reactions = {}

            # Delete genes!
            cobra.manipulation.delete_model_genes(gempro.model, genes_to_remove)

            if gempro.model._trimmed:
                log.info('{}: deleted {} reactions, {} genes'.format(gempro.model.id,
                                                                     len(gempro.model._trimmed_reactions),
                                                                     len(gempro.model._trimmed_genes)))
        # Otherwise, just mark the genes as non-functional
        else:
            for g in genes_to_remove:
                my_gene = gempro.genes.get_by_id(g)
                my_gene.functional = False
            log.info('{}: deleted {} genes'.format(gempro.id, len(genes_to_remove)))

    def build_strain_specific_models(self):
        """Using the orthologous genes matrix, modify the strain specific models based on if orthologous genes exist.

        Also stores the sequences directly in the GEM-PRO protein sequences attribute for the strains.
        """

        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix')

        # For each genome, create the strain specific model
        for strain_model in tqdm(self.atlas_strains):
            log.debug('{}: building strain specific model'.format(strain_model.id))

            # Get a list of genes which do not have orthology in the strain
            not_in_strain = self.df_orthology_matrix[pd.isnull(self.df_orthology_matrix[strain_model.id])][strain_model.id].index.tolist()

            # Mark genes non-functional
            self._pare_down_model(gempro=strain_model, genes_to_remove=not_in_strain)

            # Load sequences into the base and strain models
            self._load_strain_sequences(strain_model.id)

        log.info('Created {} new strain-specific models and loaded in sequences'.format(len(self.atlas_strains)))

    def _load_strain_sequences(self, strain_id):
        """Load strain sequences from the orthology matrix into the base model for comparisons, and into the
            strain-specific model itself.

        """
        strain_model = self.atlas_strains.get_by_id(strain_id)

        if self._orthology_matrix_has_sequences:  # Load directly from the orthology matrix if it contains sequences
            strain_sequences = self.df_orthology_matrix[strain_id].to_dict()
        else:  # Otherwise load from the genome file if the orthology matrix contains gene IDs
            # Load the genome FASTA file
            log.debug('{}: loading strain genome CDS file'.format(strain_model.genome_path))
            strain_sequences = SeqIO.index(strain_model.genome_path, 'fasta')

        for strain_gene in strain_model.genes:
            if strain_gene.functional:
                if self._orthology_matrix_has_sequences:
                    strain_gene_key = strain_gene.id
                else:
                    # Pull the gene ID of the strain from the orthology matrix
                    strain_gene_key = self.df_orthology_matrix.loc[strain_gene.id, strain_model.id]
                    log.debug('{}: original gene ID in strain fasta file'.format(strain_gene_key))

                # Load into the base strain for comparisons
                base_gene = self.base_strain_gempro.genes.get_by_id(strain_gene.id)
                new_id = '{}_{}'.format(strain_gene.id, strain_id)
                if base_gene.protein.sequences.has_id(new_id):
                    log.debug('{}: sequence already loaded into base model'.format(new_id))
                    continue
                base_gene.protein.load_manual_sequence(seq=strain_sequences[strain_gene_key], ident=new_id,
                                                       set_as_representative=False)
                log.debug('{}: loaded sequence into base model'.format(new_id))

                # Load into the strain GEM-PRO
                strain_gene.protein.load_manual_sequence(seq=strain_sequences[strain_gene_key], ident=new_id,
                                                         set_as_representative=True)
                log.debug('{}: loaded sequence into strain model'.format(new_id))

    def align_orthologous_genes_pairwise(self, gapopen=10, gapextend=0.5):
        """For each gene in the base strain, run a pairwise alignment for all orthologous gene sequences to it."""
        for base_gene in tqdm(self.base_strain_gempro.genes):
            if len(base_gene.protein.sequences) > 1:
                alignment_dir = op.join(self.atlas_seq_genes_dir, base_gene.id)
                if not op.exists(alignment_dir):
                    os.mkdir(alignment_dir)
                base_gene.protein.pairwise_align_sequences_to_representative(gapopen=gapopen, gapextend=gapextend,
                                                                             outdir=alignment_dir, parse=True)

    def align_orthologous_genes_multiple(self):
        """For each gene in the base strain, run a multiple alignment to all orthologous strain genes"""
        pass

    def get_atlas_summary_df(self):
        """Create a single data frame which summarizes all genes per row.

        Returns:
            DataFrame: Pandas DataFrame of the results

        """
        all_info = []
        for g in self.base_strain_gempro.genes_with_a_representative_sequence:
            info = {}
            info['Gene_ID'] = g.id

            # Protein object
            p = g.protein
            info['Protein_sequences'] = len(p.sequences)
            info['Protein_structures'] = len(p.structures)

            # SeqProp
            rseq = p.representative_sequence
            if rseq.seq_record:
                info['RepSeq_sequence_length'] = rseq.seq_len
                info['RepSeq_num_sequence_alignments'] = len(rseq.sequence_alignments)
                info['RepSeq_num_structure_alignments'] = len(rseq.structure_alignments)

                # SeqRecord annotations (properties calculated that summarize the whole sequence)
                rseq_sr = rseq.seq_record
                for annotation_name, annotation in rseq_sr.annotations.items():
                    info['RepSeq_' + annotation_name] = annotation

                # SeqRecord alignment annotations
                all_num_mutations = []
                all_num_deletions = []
                all_len_deletions = []
                all_num_insertions = []
                all_len_insertions = []
                all_percent_identity = []
                all_percent_similarity = []
                for aln in rseq.sequence_alignments:
                    # Gather the strain speicific stuff
                    if '{}_'.format(p.id) not in aln.annotations['b_seq']:
                        continue
                    info[aln.annotations['b_seq'].split('{}_'.format(p.id))[1]] = aln.annotations['percent_identity']

                    # Gather the percent identities/similarities
                    all_percent_identity.append(aln.annotations['percent_identity'])
                    all_percent_similarity.append(aln.annotations['percent_similarity'])

                    # Gather the number of residues that are mutated (filter for different mutations of same residue)
                    num_mutations = len(list(set([x[1] for x in aln.annotations['mutations']])))
                    all_num_mutations.append(num_mutations)

                    # Gather the number of deletions as well as the length of the deletion
                    if not aln.annotations['deletions']:
                        num_deletions = 0
                        len_deletions = [0]
                    else:
                        num_deletions = len(aln.annotations['deletions'])
                        len_deletions = [x[1] for x in aln.annotations['deletions']]
                    all_num_deletions.append(num_deletions)
                    # Get the total length of the deletion for this one strain
                    avg_len_deletions = np.sum(len_deletions)
                    all_len_deletions.append(avg_len_deletions)

                    # Gather the number of insertions as well as the length of the insertion
                    if not aln.annotations['insertions']:
                        num_insertions = 0
                        len_insertions = [0]
                    else:
                        num_insertions = len(aln.annotations['insertions'])
                        len_insertions = [x[1] for x in aln.annotations['insertions']]
                    all_num_insertions.append(num_insertions)
                    # Get the total length of insertion for this one strain
                    avg_len_insertions = np.sum(len_insertions)
                    all_len_insertions.append(avg_len_insertions)

                info['ATLAS_mean_num_mutations'] = np.mean(all_num_mutations)
                info['ATLAS_mean_num_deletions'] = np.mean(all_num_deletions)
                info['ATLAS_mean_len_deletions'] = np.mean(all_len_deletions)
                info['ATLAS_mean_num_insertions'] = np.mean(all_num_insertions)
                info['ATLAS_mean_len_insertions'] = np.mean(all_len_insertions)
                info['ATLAS_mean_percent_identity'] = np.mean(all_percent_identity)
                info['ATLAS_mean_percent_similarity'] = np.mean(all_percent_similarity)

                # Other mutation analysis
                single, fingerprint = rseq.sequence_mutation_summary()

                # Mutations that show up in more than 10% of strains
                singles = []
                for k, v in single.items():
                    k = [str(x) for x in k]
                    if len(v) / len(rseq.sequence_alignments) >= 0.01:
                        singles.append(''.join(k))  # len(v) is the number of strains which have this mutation
                info['ATLAS_popular_mutations'] = ';'.join(singles)

                # Mutation groups that show up in more than 10% of strains
                allfingerprints = []
                for k, v in fingerprint.items():
                    if len(v) / len(rseq.sequence_alignments) >= 0.01:
                        fingerprints = []
                        for m in k:
                            y = [str(x) for x in m]
                            fingerprints.append(''.join(y))
                        allfingerprints.append('-'.join(fingerprints))
                info['ATLAS_popular_mutation_groups'] = ';'.join(allfingerprints)

            # StructProp
            rstruct = p.representative_structure
            if rstruct:
                if rstruct.structure_file:
                    info['RepStruct_ID'] = rstruct.id
                    info['RepStruct_is_experimental'] = rstruct.is_experimental
                    info['RepStruct_description'] = rstruct.description
                    info['RepStruct_repseq_coverage'] = rstruct.reference_seq_top_coverage

                    # ChainProp
                    rchain = rstruct.representative_chain
                    info['RepChain_ID'] = rchain.id

                    # ChainProp SeqRecord annotations
                    rchain_sr = rchain.seq_record
                    for annotation_name, annotation in rchain_sr.annotations.items():
                        info['RepChain_' + annotation_name] = annotation

            all_info.append(info)

        cols = ['Gene_ID', 'Protein_sequences', 'Protein_structures',
                'RepSeq_num_sequence_alignments', 'RepSeq_num_structure_alignments', 'RepSeq_num_tm_helix-tmhmm',
                'RepSeq_sequence_length',
                'RepStruct_ID', 'RepStruct_description', 'RepStruct_is_experimental', 'RepStruct_repseq_coverage',
                'RepChain_ID',
                'ATLAS_mean_percent_identity', 'ATLAS_mean_percent_similarity', 'ATLAS_mean_num_mutations',
                'ATLAS_popular_mutations', 'ATLAS_popular_mutation_groups', 'ATLAS_mean_num_deletions',
                'ATLAS_mean_num_insertions', 'ATLAS_mean_len_deletions', 'ATLAS_mean_len_insertions',
                'RepSeq_aromaticity', 'RepSeq_instability_index', 'RepSeq_isoelectric_point', 'RepSeq_molecular_weight',
                'RepSeq_monoisotopic', 'RepSeq_percent_acidic', 'RepSeq_percent_aliphatic', 'RepSeq_percent_aromatic',
                'RepSeq_percent_B-sspro8', 'RepSeq_percent_basic', 'RepSeq_percent_buried-accpro',
                'RepSeq_percent_buried-accpro20', 'RepSeq_percent_C-sspro', 'RepSeq_percent_C-sspro8',
                'RepSeq_percent_charged', 'RepSeq_percent_E-sspro', 'RepSeq_percent_E-sspro8',
                'RepSeq_percent_exposed-accpro', 'RepSeq_percent_exposed-accpro20', 'RepSeq_percent_G-sspro8',
                'RepSeq_percent_H-sspro', 'RepSeq_percent_H-sspro8', 'RepSeq_percent_helix_naive',
                'RepSeq_percent_I-sspro8', 'RepSeq_percent_non-polar', 'RepSeq_percent_polar',
                'RepSeq_percent_S-sspro8', 'RepSeq_percent_small', 'RepSeq_percent_strand_naive',
                'RepSeq_percent_T-sspro8', 'RepSeq_percent_tiny', 'RepSeq_percent_turn_naive',
                'RepChain_percent_B-dssp', 'RepChain_percent_C-dssp', 'RepChain_percent_E-dssp',
                'RepChain_percent_G-dssp', 'RepChain_percent_H-dssp', 'RepChain_percent_I-dssp',
                'RepChain_percent_S-dssp', 'RepChain_percent_T-dssp', 'RepChain_SSBOND-biopython']
        cols.extend([x.id for x in self.atlas_strains])
        df_atlas_summary = pd.DataFrame(all_info, columns=cols)
        # Drop columns that don't have anything in them
        df_atlas_summary.dropna(axis=1, how='all', inplace=True)

        return df_atlas_summary

    def get_atlas_per_gene_mutation_df(self, gene_id):
        """Create a single data frame which summarizes a gene and its mutations.

        Args:
            gene_id (str): Gene ID in the base model

        Returns:
            DataFrame: Pandas DataFrame of the results

        """
        # TODO: also count: number of unique mutations (have to consider position, amino acid change)
        # TODO: keep track of strain with most mutations, least mutations
        # TODO: keep track of strains that conserve the length of the protein, others that extend or truncate it
        # need statistical test for that too (how long is "extended"/"truncated"?)
        # TODO: number of strains with at least 1 mutations
        # TODO: number of strains with <5% mutated, 5-10%, etc

        g = self.base_strain_gempro.genes.get_by_id(gene_id)

        single, fingerprint = g.protein.representative_sequence.sequence_mutation_summary()

        structure_type_suffix = 'NA'
        appender = []

        for k, strains in single.items():

            # Mutations in the strain
            to_append = {}
            orig_res = k[0]
            resnum = int(k[1])
            mutated_res = k[2]
            num_strains = len(strains)
            strain_ids = [str(x.split(g.id + '_')[1]) for x in strains]
            to_append['base_residue'] = orig_res
            to_append['base_resnum'] = resnum
            to_append['strain_residue'] = mutated_res
            to_append['num_strains'] = num_strains
            to_append['strains_mutated'] = ';'.join(strain_ids)
            to_append['at_disulfide_bridge'] = False

            # Residue properties
            origres_props = ssbio.sequence.properties.residues.residue_biochemical_definition(orig_res)
            mutres_props = ssbio.sequence.properties.residues.residue_biochemical_definition(mutated_res)
            to_append['base_residue_prop'] = origres_props
            to_append['strain_residue_prop'] = mutres_props

            # Grantham score - score a mutation based on biochemical properties
            grantham_s, grantham_txt = ssbio.sequence.properties.residues.grantham_score(orig_res, mutated_res)
            to_append['grantham_score'] = grantham_s
            to_append['grantham_annotation'] = grantham_txt

            # Structure properties - predicted from sequence
            repseq_letter_annotations = g.protein.representative_sequence.seq_record.letter_annotations
            for k in repseq_letter_annotations.keys():
                if k != 'repchain_resnums':
                    to_append['{}_PRED'.format(k)] = repseq_letter_annotations[k][resnum - 1]

            if g.protein.representative_structure:
                if g.protein.representative_structure.is_experimental:
                    structure_type_suffix = 'EXP'
                else:
                    structure_type_suffix = 'HOM'

                # Structure properties - calculated
                mapped = g.protein.representative_structure._map_repseq_resnums_to_repchain_index(resnum)
                if resnum in mapped:
                    repchain_index = mapped[resnum]

                    if not g.protein.representative_structure.representative_chain.seq_record:
                        appender.append(to_append)
                        continue

                    repchain_letter_annotations = g.protein.representative_structure.representative_chain.seq_record.letter_annotations
                    for k in repchain_letter_annotations.keys():
                        to_append['{}_{}'.format(k, structure_type_suffix)] = repchain_letter_annotations[k][repchain_index]

                    # At disulfide bond?
                    repchain_annotations = g.protein.representative_structure.representative_chain.seq_record.annotations
                    if 'SSBOND-biopython' in repchain_annotations:
                        structure_resnum = g.protein.representative_structure.map_repseq_resnums_to_structure_resnums(
                            resnum)
                        if resnum in structure_resnum:
                            ssbonds = repchain_annotations['SSBOND-biopython']
                            ssbonds_res = []
                            for x in ssbonds:
                                ssbonds_res.append(x[0])
                                ssbonds_res.append(x[1])

                            if structure_resnum in ssbonds_res:
                                to_append['at_disulfide_bridge'] = True

            appender.append(to_append)

        cols = ['base_residue', 'base_resnum', 'strain_residue', 'num_strains',
                'base_residue_prop', 'strain_residue_prop', 'grantham_score', 'grantham_annotation',
                'at_disulfide_bridge',
                'SS-sspro_PRED', 'SS-sspro8_PRED', 'SS-dssp_{}'.format(structure_type_suffix), 'TM-tmhmm_PRED',
                'RSA-accpro_PRED', 'RSA-accpro20_PRED', 'RSA-dssp_{}'.format(structure_type_suffix), 'ASA-dssp_{}'.format(structure_type_suffix),
                'CA_DEPTH-msms_{}'.format(structure_type_suffix), 'RES_DEPTH-msms_{}'.format(structure_type_suffix),
                'PHI-dssp_{}'.format(structure_type_suffix), 'PSI-dssp_{}'.format(structure_type_suffix), 'structure_resnums_{}'.format(structure_type_suffix),
                'strains_mutated']

        df_gene_summary = pd.DataFrame.from_records(appender, columns=cols)

        # Drop columns that don't have anything in them
        df_gene_summary.dropna(axis=1, how='all', inplace=True)

        df_gene_summary.sort_values(by='base_resnum', inplace=True)
        df_gene_summary = df_gene_summary.set_index('base_resnum')
        return df_gene_summary

    def download_mutation_images(self, outdir):
        # TODO: dunno if this works
        import ipywidgets
        import math

        views = []
        for g in self.base_strain_gempro.genes:
            if g.protein.representative_structure:
                view = g.protein.view_all_mutations(grouped=False, structure_opacity=0.5,
                                                    opacity_range=(0.6, 1), scale_range=(.5, 5))
                view._remote_call("setSize", target='Widget', args=['300px', '300px'])
                view.download_image(filename='{}_mutations.png'.format(g.id))
                views.append(view)

        hboxes = [ipywidgets.HBox(views[i * 3:i * 3 + 3])
                  for i in range(int(math.ceil(len(views) / 3.0)))]
        vbox = ipywidgets.VBox(hboxes)
        return vbox