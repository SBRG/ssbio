import logging
import os
import os.path as op
import sys

import cobra.flux_analysis
import cobra.manipulation
import numpy as np
import pandas as pd
from Bio import SeqIO
from cobra.core import DictList
from slugify import Slugify

import ssbio.core.modelpro
import ssbio.databases.ncbi
import ssbio.databases.patric
import ssbio.protein.sequence.properties.residues
import ssbio.protein.sequence.utils.alignment
import ssbio.protein.sequence.utils.blast
import ssbio.protein.sequence.utils.fasta
from ssbio import utils
from ssbio.core.object import Object
from ssbio.pipeline.gempro import GEMPRO

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


class ATLAS2(Object):

    """REDO ATLAS workflow"""

    def __init__(self, atlas_name, spark_context, root_dir, reference_gempro, reference_genome_path=None, description=None):
        """Prepare a GEM-PRO model for ATLAS analysis

        Args:
            atlas_name (str): Name of your ATLAS project
            spark_context (SparkContext): Spark Context for parallelization
            root_dir (str): Path to where the folder named after ``atlas_name`` will be created.
            reference_gempro (GEMPRO): GEM-PRO model to use as the reference genome
            reference_genome_path (str): Path to reference genome FASTA file
            description (str): Optional string to describe your project

        """
        Object.__init__(self, id=atlas_name, description=description)

        # Create directories
        self._root_dir = None
        self.root_dir = root_dir

        self.sc = spark_context

        self.strains = DictList()
        self.df_orthology_matrix = pd.DataFrame()
        # Mark if the orthology matrix has gene IDs (thus we need to retrieve seqs from the genome file) or if
        # it is in the orthology matrix itself
        self._orthology_matrix_has_sequences = False

        # Load the GEM-PRO (could be a model, could just be a list of genes)
        # Check if there is a genome file associated with this model - if not, write all sequences and use that
        self.reference_gempro = reference_gempro
        if not reference_genome_path and not self.reference_gempro.genome_path:
            self.reference_gempro.genome_path = self.reference_gempro.write_representative_sequences_file(outname=self.reference_gempro.id)
        else:
            self.reference_gempro.genome_path = reference_genome_path
            # TODO: must also check if reference_genome_path gene IDs can be matched to the reference_gempro

        self._empty_reference_gempro = None

    def set_empty_reference_gempro(self):
        self._empty_reference_gempro = None
        if self.reference_gempro.model:
            # If there is a SBML model associated with the GEMPRO, copy that model
            logging.disable(logging.WARNING)
            self._empty_reference_gempro = GEMPRO(gem_name='Copied reference GEM-PRO',
                                                  gem=self.reference_gempro.model.copy())
            logging.disable(logging.NOTSET)
            # Reset the GenePro attributes
            for x in self._empty_reference_gempro.genes:
                x.reset_protein()
        else:
            # Otherwise, just copy the list of genes over and rename the IDs
            strain_genes = [x.id for x in self.reference_gempro.genes]
            if len(strain_genes) == 0:
                raise ValueError('GEM-PRO has no genes, unable to run multi-strain analysis')
            self._empty_reference_gempro = GEMPRO(gem_name='Copied reference GEM-PRO', genes_list=strain_genes)

    @property
    def root_dir(self):
        """str: Directory where ATLAS project folder named after the attribute ``base_dir`` is located"""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.info('Changing root directory of project "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)):
                raise IOError('Project "{}" does not exist in folder {}'.format(self.id, path))
        else:
            log.info('Creating project directory in folder {}'.format(path))

        self._root_dir = path

        for d in [self.base_dir, self.model_dir, self.data_dir,
                  self.sequences_dir, self.sequences_by_gene_dir, self.sequences_by_organism_dir]:
            ssbio.utils.make_dir(d)

        log.info('{}: project location'.format(self.base_dir))

    @property
    def base_dir(self):
        """str: ATLAS project folder"""
        if self.root_dir:
            return op.join(self.root_dir, self.id)
        else:
            return None

    @property
    def model_dir(self):
        """str: Directory where strain-specific GEMs are stored"""
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
    def sequences_dir(self):
        """str: Base directory for genome protein sequences and alignments"""
        if self.base_dir:
            return op.join(self.base_dir, 'sequences')
        else:
            return None

    @property
    def sequences_by_gene_dir(self):
        """str: Directory where all gene specific information and pairwise alignments are stored"""
        if self.sequences_dir:
            return op.join(self.sequences_dir, 'by_gene')
        else:
            return None

    @property
    def sequences_by_organism_dir(self):
        """str: Directory where all strain specific genome and BLAST files are stored"""
        if self.sequences_dir:
            return op.join(self.sequences_dir, 'by_organism')
        else:
            return None

    @property
    def strains_rdd(self):
        if len(self.strains) == 0:
            raise ValueError('No strains to create RDD from!')
        return self.sc.parallelize(self.strains)

    def load_strains(self, strains_to_fasta_files):
        """Load in a dictionary of strain IDs to their protein FASTA files"""
        logging.disable(logging.WARNING)
        for strain_id, strain_fasta_file in strains_to_fasta_files.items():
            strain_gp = GEMPRO(gem_name=strain_id)
            strain_gp.genome_path = strain_fasta_file
            self.strains.append(strain_gp)
        logging.disable(logging.NOTSET)

    def get_orthology_matrix(self, outfile, outdir=None,
                             pid_cutoff=None, bitscore_cutoff=None, evalue_cutoff=None, filter_condition='OR',
                             force_rerun=False):
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

        ################################################################################################################
        # BIDIRECTIONAL BLAST
        def run_bidirectional_blast(strain_gempro,
                                    r_file=self.reference_gempro.genome_path,
                                    outdir=self.sequences_by_organism_dir):
            import ssbio.protein.sequence.utils.blast
            g_file = strain_gempro.genome_path

            # Run bidirectional BLAST
            r_vs_g, g_vs_r = ssbio.protein.sequence.utils.blast.run_bidirectional_blast(reference=r_file,
                                                                                        other_genome=g_file,
                                                                                        dbtype='prot',
                                                                                        outdir=outdir)

            # Using the BLAST files, find the BBH
            bbh = ssbio.protein.sequence.utils.blast.calculate_bbh(blast_results_1=r_vs_g, blast_results_2=g_vs_r,
                                                                   outdir=outdir)
            return strain_gempro.id, bbh

        log.info('Running bidirectional BLAST and finding best bidirectional hits (BBH)...')
        result = self.strains_rdd.map(run_bidirectional_blast).collect()
        bbh_files = dict(result)
        ################################################################################################################

        ################################################################################################################
        # ORTHOLOGY MATRIX
        if not outdir:
            outdir = self.data_dir
        log.info('Creating orthology matrix from BBHs...')
        ortho_matrix = ssbio.protein.sequence.utils.blast.create_orthology_matrix(r_name=self.reference_gempro.id,
                                                                                  genome_to_bbh_files=bbh_files,
                                                                                  pid_cutoff=pid_cutoff,
                                                                                  bitscore_cutoff=bitscore_cutoff,
                                                                                  evalue_cutoff=evalue_cutoff,
                                                                                  filter_condition=filter_condition,
                                                                                  outname=outfile,
                                                                                  outdir=outdir,
                                                                                  force_rerun=force_rerun)

        log.info('Saved orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))
        self.df_orthology_matrix = pd.read_csv(ortho_matrix, index_col=0)
        self.df_orthology_matrix = self.df_orthology_matrix.rename_axis('gene').rename_axis("strain", axis="columns")
        ################################################################################################################

    def filter_genes_and_strains(self, remove_genes_not_in_reference_model=True,
                                 remove_strains_with_no_orthology=True, remove_strains_with_no_differences=False,
                                 custom_keep_strains=None, custom_keep_genes=None):
        """Filters the analysis by keeping a subset of strains or genes based on certain criteria.

        Args:
            remove_genes_not_in_reference_model (bool): Remove genes from reference model not in orthology matrix
            remove_strains_with_no_orthology (bool): Remove strains which have no orthologous genes found
            remove_strains_with_no_differences (bool): Remove strains which have all the same genes as the base model.
                Default is False because since orthology is found using a PID cutoff, all genes may be present but
                differences may be on the sequence level.
            custom_keep_genes (list): List of gene IDs to keep in analysis
            custom_keep_strains (list): List of strain IDs to keep in analysis

        """

        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix, please calculate first!')

        reference_strain_gene_ids = [x.id for x in self.reference_gempro.genes]
        initial_num_genes = len(reference_strain_gene_ids)

        strain_ids = [x.id for x in self.strains]
        initial_num_strains = len(strain_ids)

        # Gene filtering
        if custom_keep_genes:
            to_remove = [x for x in reference_strain_gene_ids if x not in custom_keep_genes]
            if self.reference_gempro.model:
                cobra.manipulation.delete_model_genes(self.reference_gempro.model, to_remove)
            else:
                for g_id in to_remove:
                    self.reference_gempro.genes.remove(g_id)
        if remove_genes_not_in_reference_model:
            to_remove = [x for x in reference_strain_gene_ids if x not in self.df_orthology_matrix.index.tolist()]
            if self.reference_gempro.model:
                cobra.manipulation.delete_model_genes(self.reference_gempro.model, to_remove)
            else:
                for g_id in to_remove:
                    self.reference_gempro.genes.remove(g_id)

        # Strain filtering
        if custom_keep_strains or remove_strains_with_no_orthology or remove_strains_with_no_differences:
            for strain_id in strain_ids:
                if custom_keep_strains:
                    if strain_id not in custom_keep_strains:
                        self.strains.remove(strain_id)
                        continue

                if remove_strains_with_no_orthology:
                    if strain_id not in self.df_orthology_matrix.columns:
                        self.strains.remove(strain_id)
                        log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_id))
                        continue
                    elif self.df_orthology_matrix[strain_id].isnull().all():
                        self.strains.remove(strain_id)
                        log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_id))
                        continue

                if remove_strains_with_no_differences:
                    not_in_strain = self.df_orthology_matrix[pd.isnull(self.df_orthology_matrix[strain_id])][strain_id].index.tolist()
                    if len(not_in_strain) == 0:
                        self.strains.remove(strain_id)
                        log.info('{}: strain has no differences from the base, removed from analysis.')
                        continue

        log.info('{} genes to be analyzed, originally {}'.format(len(self.reference_gempro.functional_genes), initial_num_genes))
        log.info('{} strains to be analyzed, originally {}'.format(len(self.strains), initial_num_strains))

    def build_strain_specific_models(self, save_models=False):
        """Using the orthologous genes matrix, create and modify the strain specific models based on if orthologous
        genes exist.
        """

        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix, please calculate first!')

        # Make sure to set a copy of the original model to a empty reference one
        self.set_empty_reference_gempro()

        # Create an emptied copy of the reference GEM-PRO
        def build_strain_specific_model(strain_gempro, empty_reference_gempro=self._empty_reference_gempro,
                                        orth_matrix=self.df_orthology_matrix,
                                        save_model=save_models, model_dir=self.model_dir):
            # For each genome, load the metabolic model or genes from the reference GEM-PRO
            if empty_reference_gempro.model:
                logging.disable(logging.WARNING)
                strain_gempro.load_cobra_model(empty_reference_gempro.model)
                logging.disable(logging.NOTSET)
                strain_gempro.model._trimmed = False
                strain_gempro.model._trimmed_genes = []
                strain_gempro.model._trimmed_reactions = {}
            elif empty_reference_gempro.genes:
                strain_gempro.genes = empty_reference_gempro.genes.copy()

            # Get a list of genes which do not have orthology in the strain
            genes_to_remove = orth_matrix[pd.isnull(orth_matrix[strain_gempro.id])][strain_gempro.id].index.tolist()

            # Mark genes non-functional
            strain_genes = [x.id for x in strain_gempro.genes]
            genes_to_remove = list(set(genes_to_remove).intersection(set(strain_genes)))

            if len(genes_to_remove) > 0:
                # If a COBRApy model exists, utilize the delete_model_genes method
                if strain_gempro.model:
                    cobra.manipulation.delete_model_genes(strain_gempro.model, genes_to_remove)
                # Otherwise, just mark the genes as non-functional
                else:
                    for g in genes_to_remove:
                        strain_gempro.genes.get_by_id(g).functional = False

            if save_model:
                cobra.io.save_json_model(model=strain_gempro.model,
                                         filename=op.join(model_dir, '{}.json'.format(strain_gempro.id)))

            return strain_gempro

        self.strains = DictList(self.strains_rdd.map(build_strain_specific_model).collect())
        # for s in self.strains:
        #     build_strain_specific_model(s)

    def load_strains_sequences(self, save_models=False):
        """Load all strains sequences from the orthology matrix into the base model for comparisons, and into the
        strain-specific model itself.

        """

        def load_sequences_to_strain(strain_gempro, orth_matrix=self.df_orthology_matrix,
                                     save_model=save_models, model_dir=self.model_dir):
            strain_sequences = SeqIO.index(strain_gempro.genome_path, 'fasta')
            # strain_sequences = SeqIO.to_dict(SeqIO.parse(strain_gempro.genome_path, 'fasta'))
            for strain_gene in strain_gempro.functional_genes:
                # Pull the gene ID of the strain from the orthology matrix
                strain_gene_key = orth_matrix.at[strain_gene.id, strain_gempro.id]
                # Load into the strain GEM-PRO
                new_id = '{}_{}'.format(strain_gene.id, strain_gempro.id)
                strain_gene.protein.load_manual_sequence(seq=strain_sequences[strain_gene_key], ident=new_id,
                                                         set_as_representative=True)
            if save_model:
                strain_gempro.save_json(outfile=op.join(model_dir, '{}_gp.json'.format(strain_gempro.id)))
                strain_gempro.save_pickle(outfile=op.join(model_dir, '{}_gp.pckl'.format(strain_gempro.id)))
            return strain_gempro

        self.strains = DictList(self.strains_rdd.map(load_sequences_to_strain).collect())
        # for s in self.strains:
        #     load_sequences_to_strain(s)

        def load_sequences_to_reference_gene(g, strain_gempros=self.strains,
                                             orth_matrix=self.df_orthology_matrix):
            for strain_gempro in strain_gempros:
                strain_sequences = SeqIO.index(strain_gempro.genome_path, 'fasta')
                strain_gene = strain_gempro.genes.get_by_id(g.id)
                if strain_gene.functional:
                    # Pull the gene ID of the strain from the orthology matrix
                    strain_gene_key = orth_matrix.at[strain_gene.id, strain_gempro.id]
                    new_id = '{}_{}'.format(strain_gene.id, strain_gempro.id)
                    if g.protein.sequences.has_id(new_id):
                        continue
                    g.protein.load_manual_sequence(seq=strain_sequences[strain_gene_key], ident=new_id,
                                                   set_as_representative=False)
            return g

        # for g in self.reference_gempro.functional_genes:
        #     load_sequences_to_reference_gene(g)
        genes_rdd = self.sc.parallelize(self.reference_gempro.functional_genes)
        result = genes_rdd.map(load_sequences_to_reference_gene).collect()
        for modified_g in result:
            original_gene = self.reference_gempro.genes.get_by_id(modified_g.id)
            original_gene.copy_modified_gene(modified_g)

    def align_orthologous_genes_pairwise(self, gapopen=10, gapextend=0.5):
        """For each gene in the base strain, run a pairwise alignment for all orthologous gene sequences to it."""
        def run_alignments(g, outdir=self.sequences_by_gene_dir):
            if len(g.protein.sequences) > 1:
                alignment_dir = op.join(outdir, g.id)
                if not op.exists(alignment_dir):
                    os.mkdir(alignment_dir)
                g.protein.pairwise_align_sequences_to_representative(gapopen=gapopen, gapextend=gapextend,
                                                                     outdir=alignment_dir, parse=True)
            return g

        genes_rdd = self.sc.parallelize(self.reference_gempro.functional_genes)
        result = genes_rdd.map(run_alignments).collect()
        for modified_g in result:
            original_gene = self.reference_gempro.genes.get_by_id(modified_g.id)
            original_gene.copy_modified_gene(modified_g)

def get_atlas_summary_df(self):
    """Create a single data frame which summarizes all genes per row.

    Returns:
        DataFrame: Pandas DataFrame of the results

    """
    all_info = []
    for g in self.reference_gempro.genes_with_a_representative_sequence:
        info = {}
        info['Gene_ID'] = g.id
        info['Gene_name'] = g.name

        # Protein object
        p = g.protein
        info['Protein_sequences'] = len(p.sequences)
        info['Protein_structures'] = len(p.structures)

        # SeqProp
        rseq = p.representative_sequence
        info['RepSeq_ID'] = rseq.id
        info['RepSeq_sequence_length'] = rseq.seq_len
        info['RepSeq_num_sequence_alignments'] = len([x for x in p.sequence_alignments if x.annotations['ssbio_type'] == 'seqalign'])
        info['RepSeq_num_structure_alignments'] = len([x for x in p.sequence_alignments if x.annotations['ssbio_type'] == 'structalign'])

        # SeqRecord annotations (properties calculated that summarize the whole sequence)
        for annotation_name, annotation in rseq.annotations.items():
            info['RepSeq_' + annotation_name] = annotation

        # SeqRecord alignment annotations
        all_num_mutations = []
        all_num_deletions = []
        all_len_deletions = []
        all_num_insertions = []
        all_len_insertions = []
        all_percent_identity = []
        all_percent_similarity = []
        for aln in p.sequence_alignments:
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
        single, fingerprint = p.sequence_mutation_summary()

        # Mutations that show up in more than 10% of strains
        singles = []
        for k, v in single.items():
            k = [str(x) for x in k]
            if len(v) / len(p.sequence_alignments) >= 0.01:
                singles.append(''.join(k))  # len(v) is the number of strains which have this mutation
        info['ATLAS_popular_mutations'] = ';'.join(singles)

        # Mutation groups that show up in more than 10% of strains
        allfingerprints = []
        for k, v in fingerprint.items():
            if len(v) / len(p.sequence_alignments) >= 0.01:
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
                info['RepStruct_repseq_coverage'] = p.representative_chain_seq_coverage

                # ChainProp
                rchain = p.representative_chain
                info['RepChain_ID'] = rchain

                # ChainProp SeqRecord annotations
                rchain_sr = rstruct.chains.get_by_id(rchain).seq_record
                for annotation_name, annotation in rchain_sr.annotations.items():
                    info['RepChain_' + annotation_name] = annotation

        all_info.append(info)

    cols = ['Gene_ID', 'Gene_name', 'Protein_sequences', 'Protein_structures',
            'RepSeq_ID', 'RepSeq_sequence_length',
            'RepSeq_num_sequence_alignments', 'RepSeq_num_structure_alignments',
            'RepStruct_ID', 'RepChain_ID', 'RepStruct_description',
            'RepStruct_is_experimental', 'RepStruct_repseq_coverage',
            'ATLAS_mean_percent_identity', 'ATLAS_mean_percent_similarity', 'ATLAS_mean_num_mutations',
            'ATLAS_popular_mutations', 'ATLAS_popular_mutation_groups', 'ATLAS_mean_num_deletions',
            'ATLAS_mean_num_insertions', 'ATLAS_mean_len_deletions', 'ATLAS_mean_len_insertions',
            'RepSeq_aromaticity', 'RepSeq_instability_index', 'RepSeq_isoelectric_point', 'RepSeq_molecular_weight',
            'RepSeq_monoisotopic', 'RepSeq_num_tm_helix-tmhmm', 'RepSeq_percent_acidic', 'RepSeq_percent_aliphatic',
            'RepSeq_percent_aromatic', 'RepSeq_percent_B-sspro8', 'RepSeq_percent_basic',
            'RepSeq_percent_buried-accpro', 'RepSeq_percent_buried-accpro20', 'RepSeq_percent_C-sspro',
            'RepSeq_percent_C-sspro8', 'RepSeq_percent_charged', 'RepSeq_percent_E-sspro',
            'RepSeq_percent_E-sspro8', 'RepSeq_percent_exposed-accpro', 'RepSeq_percent_exposed-accpro20',
            'RepSeq_percent_G-sspro8', 'RepSeq_percent_H-sspro', 'RepSeq_percent_H-sspro8',
            'RepSeq_percent_helix_naive', 'RepSeq_percent_I-sspro8', 'RepSeq_percent_non-polar',
            'RepSeq_percent_polar', 'RepSeq_percent_S-sspro8', 'RepSeq_percent_small',
            'RepSeq_percent_strand_naive', 'RepSeq_percent_T-sspro8', 'RepSeq_percent_tiny',
            'RepSeq_percent_turn_naive', 'RepChain_percent_B-dssp', 'RepChain_percent_C-dssp',
            'RepChain_percent_E-dssp', 'RepChain_percent_G-dssp', 'RepChain_percent_H-dssp',
            'RepChain_percent_I-dssp', 'RepChain_percent_S-dssp', 'RepChain_percent_T-dssp',
            'RepChain_SSBOND-biopython']
    cols.extend([x.id for x in self.strains])

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

    g = self.reference_gempro.genes.get_by_id(gene_id)

    single, fingerprint = g.protein.sequence_mutation_summary(alignment_type='seqalign')

    structure_type_suffix = 'NA'
    appender = []

    for k, strains in single.items():
        # Mutations in the strain
        to_append = {}
        orig_res = k[0]
        resnum = int(k[1])
        mutated_res = k[2]
        num_strains_mutated = len(strains)
        strain_ids = [str(x.split(g.id + '_')[1]) for x in strains]
        to_append['ref_residue'] = orig_res
        to_append['ref_resnum'] = resnum
        to_append['strain_residue'] = mutated_res
        to_append['num_strains_mutated'] = num_strains_mutated
        to_append['strains_mutated'] = ';'.join(strain_ids)
        to_append['at_disulfide_bridge'] = False

        # Residue properties
        origres_props = ssbio.protein.sequence.properties.residues.residue_biochemical_definition(orig_res)
        mutres_props = ssbio.protein.sequence.properties.residues.residue_biochemical_definition(mutated_res)
        to_append['ref_residue_prop'] = origres_props
        to_append['strain_residue_prop'] = mutres_props

        # Grantham score - score a mutation based on biochemical properties
        grantham_s, grantham_txt = ssbio.protein.sequence.properties.residues.grantham_score(orig_res, mutated_res)
        to_append['grantham_score'] = grantham_s
        to_append['grantham_annotation'] = grantham_txt

        # Get all per residue annotations - predicted from sequence and calculated from structure
        to_append.update(g.protein.get_residue_annotations(seq_resnum=resnum, use_representatives=True))

        # Check structure type
        if g.protein.representative_structure:
            if g.protein.representative_structure.is_experimental:
                to_append['structure_type'] = 'EXP'
            else:
                to_append['structure_type'] = 'HOM'

            # At disulfide bond?
            repchain = g.protein.representative_chain
            repchain_annotations = g.protein.representative_structure.chains.get_by_id(repchain).seq_record.annotations
            if 'SSBOND-biopython' in repchain_annotations:
                structure_resnum = g.protein.map_seqprop_resnums_to_structprop_resnums(resnums=resnum,
                                                                                       use_representatives=True)
                if resnum in structure_resnum:
                    ssbonds = repchain_annotations['SSBOND-biopython']
                    ssbonds_res = []
                    for x in ssbonds:
                        ssbonds_res.append(x[0])
                        ssbonds_res.append(x[1])

                    if structure_resnum in ssbonds_res:
                        to_append['at_disulfide_bridge'] = True

        appender.append(to_append)

    if not appender:
        return pd.DataFrame()

    cols = ['ref_residue', 'ref_resnum', 'strain_residue', 'num_strains_mutated', 'strains_mutated',
            'ref_residue_prop', 'strain_residue_prop', 'grantham_score', 'grantham_annotation',
            'at_disulfide_bridge',
            'seq_SS-sspro', 'seq_SS-sspro8', 'seq_RSA-accpro', 'seq_RSA-accpro20', 'seq_TM-tmhmm',
            'struct_SS-dssp', 'struct_RSA-dssp', 'struct_ASA-dssp',
            'struct_CA_DEPTH-msms', 'struct_RES_DEPTH-msms',
            'struct_PHI-dssp', 'struct_PSI-dssp',
            'struct_resnum', 'struct_residue'
            'strains_mutated']

    df_gene_summary = pd.DataFrame.from_records(appender, columns=cols)

    # Drop columns that don't have anything in them
    df_gene_summary.dropna(axis=1, how='all', inplace=True)

    df_gene_summary.sort_values(by='ref_resnum', inplace=True)
    df_gene_summary = df_gene_summary.set_index('ref_resnum')
    return df_gene_summary

def download_mutation_images(self):
    # TODO: dunno if this works
    import ipywidgets
    import math

    views = []
    for g in self.reference_gempro.genes:
        if g.protein.representative_structure:
            view = g.protein.view_all_mutations(alignment_type='seqalign', grouped=False, structure_opacity=0.5,
                                                opacity_range=(0.6, 1), scale_range=(.5, 5))
            view._remote_call("setSize", target='Widget', args=['300px', '300px'])
            view.download_image(filename='{}_{}_mutations.png'.format(g.id, g.name))
            views.append(view)

    hboxes = [ipywidgets.HBox(views[i * 3:i * 3 + 3])
              for i in range(int(math.ceil(len(views) / 3.0)))]
    vbox = ipywidgets.VBox(hboxes)
    return vbox