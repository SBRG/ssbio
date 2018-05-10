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
import json
import ssbio.core.modelpro
import ssbio.databases.ncbi
import ssbio.databases.patric
import ssbio.protein.sequence.properties.residues
import ssbio.protein.sequence.utils.alignment
import ssbio.protein.sequence.utils.blast
import ssbio.protein.sequence.utils.fasta
from joblib import Parallel, delayed
import ssbio.utils
import ssbio.io
from copy import deepcopy
from ssbio.core.object import Object
from ssbio.pipeline.gempro import GEMPRO
from collections import defaultdict
from more_itertools import locate

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

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


class ATLAS2(Object):

    """REDO ATLAS workflow"""

    def __init__(self, atlas_name, root_dir, reference_gempro, reference_genome_path,
                 strains_to_fasta_files, description=None):
        """Prepare a GEM-PRO model for ATLAS analysis

        Args:
            atlas_name (str): Name of your ATLAS project
            root_dir (str): Path to where the folder named after ``atlas_name`` will be created.
            reference_gempro (GEMPRO): GEM-PRO model to use as the reference genome
            reference_genome_path (str): Path to reference genome FASTA file
            strains_to_fasta_files (dict): Strain IDs to their genome FASTA files
            description (str): Optional string to describe your project

        """
        Object.__init__(self, id=atlas_name, description=description)
        self._root_dir = None
        self.root_dir = root_dir

        # Load the GEM-PRO (could be a model, could just be a list of genes)
        self.reference_gempro = reference_gempro
        self.reference_gempro.genome_path = reference_genome_path

        self.gene_protein_pickles = deepcopy(self.reference_gempro.gene_protein_pickles)
        if not self.gene_protein_pickles:
            self.reference_gempro.save_protein_pickles()
            self.gene_protein_pickles = deepcopy(self.reference_gempro.gene_protein_pickles)

        self.strain_ids = []
        """list: Strain IDs to analyze"""
        self.strain_infodict = defaultdict(dict)
        """dict: Strain genome paths and functional gene information dictionary """

        for strain_id, strain_genome_path in strains_to_fasta_files.items():
            self.strain_ids.append(strain_id)
            self.strain_infodict[strain_id]['genome_path'] = strain_genome_path

        self.df_orthology_matrix = pd.DataFrame()
        """DataFrame: Pandas Dataframe representation of the orthology matrix, containing strain FASTA sequence IDs"""

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
            log.info('Project directory in folder {}'.format(path))

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

    def get_orthology_matrix(self, outfile, sc, outdir=None,
                             pid_cutoff=None, bitscore_cutoff=None, evalue_cutoff=None,
                             force_rerun=False):
        """Create the orthology matrix by finding best bidirectional BLAST hits. Genes = rows, strains = columns

        Runs run_makeblastdb, run_bidirectional_blast, and calculate_bbh for protein sequences.

        Args:
            outfile (str): Filename with extension of the orthology matrix (ie. df_orthology.csv)
            outdir (str): Path to output of orthology matrix, default is ATLAS data_dir
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
        # Output directory for orthology matrix
        if not outdir:
            outdir = self.data_dir

        ortho_matrix = op.join(outdir, outfile)
        if ssbio.utils.force_rerun(flag=force_rerun, outfile=ortho_matrix):
            if not sc:
                raise ValueError('Please initialize SparkContext')
            ################################################################################################################
            # BIDIRECTIONAL BLAST
            def run_bidirectional_blast(strain_id, strain_genome_path,
                                        r_file=self.reference_gempro.genome_path,
                                        outdir2=self.sequences_by_organism_dir):
                import ssbio.protein.sequence.utils.blast

                # Run bidirectional BLAST
                r_vs_g, g_vs_r = ssbio.protein.sequence.utils.blast.run_bidirectional_blast(reference=r_file,
                                                                                            other_genome=strain_genome_path,
                                                                                            dbtype='prot',
                                                                                            outdir=outdir2)

                # Using the BLAST files, find the BBH
                bbh = ssbio.protein.sequence.utils.blast.calculate_bbh(blast_results_1=r_vs_g, blast_results_2=g_vs_r,
                                                                       outdir=outdir2)
                return strain_id, bbh

            log.info('Running bidirectional BLAST and finding best bidirectional hits (BBH)...')
            for_strains_rdd = [(k, v['genome_path']) for k,v in self.strain_infodict.items()]
            from random import shuffle
            shuffle(for_strains_rdd)
            strains_rdd = sc.parallelize(for_strains_rdd)
            result = strains_rdd.map(lambda x: run_bidirectional_blast(strain_id=x[0], strain_genome_path=x[1])).collect()
            bbh_files = dict(result)
            ################################################################################################################

            ################################################################################################################
            # ORTHOLOGY MATRIX
            log.info('Creating orthology matrix from BBHs...')
            ortho_matrix = ssbio.protein.sequence.utils.blast.create_orthology_matrix(r_name=self.reference_gempro.id,
                                                                                      genome_to_bbh_files=bbh_files,
                                                                                      pid_cutoff=pid_cutoff,
                                                                                      bitscore_cutoff=bitscore_cutoff,
                                                                                      evalue_cutoff=evalue_cutoff,
                                                                                      outname=outfile,
                                                                                      outdir=outdir,
                                                                                      force_rerun=force_rerun)

        log.info('Orthology matrix at {}. See the "df_orthology_matrix" attribute.'.format(ortho_matrix))
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

        initial_num_strains = len(self.strain_ids)

        # Gene filtering
        to_remove_genes = []
        if custom_keep_genes:
            to_remove_genes.extend([x for x in reference_strain_gene_ids if x not in custom_keep_genes])
        if remove_genes_not_in_reference_model:
            to_remove_genes.extend([x for x in reference_strain_gene_ids if x not in self.df_orthology_matrix.index.tolist()])

        to_remove_genes = list(set(to_remove_genes))
        if self.reference_gempro.model:
            cobra.manipulation.delete_model_genes(self.reference_gempro.model, to_remove_genes)
        else:
            for g_id in to_remove_genes:
                self.reference_gempro.genes.get_by_id(g_id).functional = False

        # Create new orthology matrix with only our genes of interest
        new_gene_subset = [x.id for x in self.reference_gempro.functional_genes]
        tmp_new_orthology_matrix = self.df_orthology_matrix[self.df_orthology_matrix.index.isin(new_gene_subset)]

        # Strain filtering
        if custom_keep_strains or remove_strains_with_no_orthology or remove_strains_with_no_differences:
            for strain_id in self.strain_ids:
                if custom_keep_strains:
                    if strain_id not in custom_keep_strains:
                        self.strain_ids.remove(strain_id)
                        continue

                if remove_strains_with_no_orthology:
                    if strain_id not in tmp_new_orthology_matrix.columns:
                        self.strain_ids.remove(strain_id)
                        log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_id))
                        continue
                    elif tmp_new_orthology_matrix[strain_id].isnull().all():
                        self.strain_ids.remove(strain_id)
                        log.info('{}: no orthologous genes found for this strain, removed from analysis.'.format(strain_id))
                        continue

                if remove_strains_with_no_differences:
                    not_in_strain = tmp_new_orthology_matrix[pd.isnull(tmp_new_orthology_matrix[strain_id])][strain_id].index.tolist()
                    if len(not_in_strain) == 0:
                        self.strain_ids.remove(strain_id)
                        log.info('{}: strain has no differences from the base, removed from analysis.')
                        continue

        log.info('{} genes to be analyzed, originally {}'.format(len(self.reference_gempro.functional_genes), initial_num_genes))
        log.info('{} strains to be analyzed, originally {}'.format(len(self.strain_ids), initial_num_strains))

    def _write_strain_functional_genes(self, strain_id, ref_functional_genes, orth_matrix, force_rerun=False):
        """Create strain functional genes json file"""
        func_genes_path = op.join(self.model_dir, '{}_funcgenes.json'.format(strain_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=func_genes_path):
            gene_to_func = {k:True for k in ref_functional_genes}
            # Get a list of genes which do not have orthology in the strain
            genes_to_remove = orth_matrix[pd.isnull(orth_matrix[strain_id])][strain_id].index.tolist()

            # Mark genes non-functional
            genes_to_remove = list(set(genes_to_remove).intersection(set(ref_functional_genes)))

            if len(genes_to_remove) > 0:
                for g in genes_to_remove:
                    gene_to_func[g] = False

            with open(func_genes_path, 'w') as f:
                json.dump(gene_to_func, f)
        else:
            with open(func_genes_path, 'r') as f:
                gene_to_func = json.load(f)

        return strain_id, gene_to_func

    def write_strain_functional_genes(self, force_rerun=False):
        """Wrapper function for _write_strain_functional_genes"""
        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix, please calculate first!')
        ref_functional_genes = [g.id for g in self.reference_gempro.functional_genes]
        log.info('Building strain specific models...')
        result = []
        for s in tqdm(self.strain_ids):
            result.append(self._write_strain_functional_genes(s, ref_functional_genes, self.df_orthology_matrix, force_rerun=force_rerun))

        for strain_id, functional_genes in result:
            self.strain_infodict[strain_id]['functional_genes'] = functional_genes

    def _build_strain_specific_model(self, strain_id, ref_functional_genes, orth_matrix, force_rerun=False):
        """Create strain GEMPRO, set functional genes"""
        gp_noseqs_path = op.join(self.model_dir, '{}_gp.pckl'.format(strain_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=gp_noseqs_path):
            logging.disable(logging.WARNING)

            strain_gp = GEMPRO(gem_name=strain_id)

            # if self.reference_gempro.model:
            #     strain_gp.load_cobra_model(deepcopy(self.reference_gempro.model))
            #     # Reset the GenePro attributes
            #     for x in strain_gp.genes:
            #         x.reset_protein()
            # else:
                # Otherwise, just copy the list of genes over and rename the IDs
            strain_genes = [x for x in ref_functional_genes]
            strain_gp.add_gene_ids(strain_genes)

            logging.disable(logging.NOTSET)

            # Get a list of genes which do not have orthology in the strain
            genes_to_remove = orth_matrix[pd.isnull(orth_matrix[strain_id])][strain_id].index.tolist()

            # Mark genes non-functional
            strain_genes = [x.id for x in strain_gp.genes]
            genes_to_remove = list(set(genes_to_remove).intersection(set(strain_genes)))

            if len(genes_to_remove) > 0:
                # If a COBRApy model exists, utilize the delete_model_genes method
                # if strain_gp.model:
                #     cobra.manipulation.delete_model_genes(strain_gp.model, genes_to_remove)
                # # Otherwise, just mark the genes as non-functional
                # else:
                for g in genes_to_remove:
                    strain_gp.genes.get_by_id(g).functional = False

            strain_gp.save_pickle(outfile=gp_noseqs_path)

        return strain_id, gp_noseqs_path

    def build_strain_specific_models(self, joblib=False, cores=1, force_rerun=False):
        """Wrapper function for _build_strain_specific_model"""
        if len(self.df_orthology_matrix) == 0:
            raise RuntimeError('Empty orthology matrix, please calculate first!')
        ref_functional_genes = [g.id for g in self.reference_gempro.functional_genes]
        log.info('Building strain specific models...')
        if joblib:
            result = DictList(Parallel(n_jobs=cores)(delayed(self._build_strain_specific_model)(s, ref_functional_genes, self.df_orthology_matrix, force_rerun=force_rerun) for s in self.strain_ids))
        # if sc:
        #     strains_rdd = sc.parallelize(self.strain_ids)
        #     result = strains_rdd.map(self._build_strain_specific_model).collect()
        else:
            result = []
            for s in tqdm(self.strain_ids):
                result.append(self._build_strain_specific_model(s, ref_functional_genes, self.df_orthology_matrix, force_rerun=force_rerun))

        for strain_id, gp_noseqs_path in result:
            self.strain_infodict[strain_id]['gp_noseqs_path'] = gp_noseqs_path

    def _load_sequences_to_strain(self, strain_id, force_rerun=False):
        """Load strain GEMPRO with functional genes defined, load sequences to it, save as new GEMPRO"""
        gp_seqs_path = op.join(self.model_dir, '{}_gp_withseqs.pckl'.format(strain_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=gp_seqs_path):
            gp_noseqs = ssbio.io.load_pickle(self.strain_infodict[strain_id]['gp_noseqs_path'])
            strain_sequences = SeqIO.index(self.strain_infodict[strain_id]['genome_path'], 'fasta')
            for strain_gene in gp_noseqs.functional_genes:
                # Pull the gene ID of the strain from the orthology matrix
                strain_gene_key = self.df_orthology_matrix.at[strain_gene.id, strain_id]
                # Load into the strain GEM-PRO
                new_id = '{}_{}'.format(strain_gene.id, strain_id)
                if strain_gene.protein.sequences.has_id(new_id):
                    continue
                strain_gene.protein.load_manual_sequence(seq=strain_sequences[strain_gene_key], ident=new_id,
                                                         set_as_representative=True)
            gp_noseqs.save_pickle(outfile=gp_seqs_path)

        return strain_id, gp_seqs_path

    def load_sequences_to_strains(self, joblib=False, cores=1, force_rerun=False):
        """Wrapper function for _load_sequences_to_strain"""
        log.info('Loading sequences to strain GEM-PROs...')
        if joblib:
            result = DictList(Parallel(n_jobs=cores)(delayed(self._load_sequences_to_strain)(s, force_rerun) for s in self.strain_ids))
        else:
            result = []
            for s in tqdm(self.strain_ids):
                result.append(self._load_sequences_to_strain(s, force_rerun))

        for strain_id, gp_seqs_path in result:
            self.strain_infodict[strain_id]['gp_seqs_path'] = gp_seqs_path

    def _load_sequences_to_reference_gene(self, g_id, force_rerun=False):
        """Load orthologous strain sequences to reference Protein object, save as new pickle"""
        protein_seqs_pickle_path = op.join(self.sequences_by_gene_dir, '{}_protein_withseqs.pckl'.format(g_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_pickle_path):
            protein_pickle_path = self.gene_protein_pickles[g_id]
            protein_pickle = ssbio.io.load_pickle(protein_pickle_path)

            for strain, info in self.strain_infodict.items():
                strain_sequences = SeqIO.index(info['genome_path'], 'fasta')
                strain_gene_functional = info['functional_genes'][g_id]
                if strain_gene_functional:
                    # Pull the gene ID of the strain from the orthology matrix
                    strain_gene_key = self.df_orthology_matrix.at[g_id, strain]
                    new_id = '{}_{}'.format(g_id, strain)
                    if protein_pickle.sequences.has_id(new_id):
                        continue
                    protein_pickle.load_manual_sequence(seq=strain_sequences[strain_gene_key],
                                                        ident=new_id,
                                                        set_as_representative=False)
            protein_pickle.save_pickle(outfile=protein_seqs_pickle_path)

        return g_id, protein_seqs_pickle_path

    def load_sequences_to_reference(self, sc=None, force_rerun=False):
        """Wrapper for _load_sequences_to_reference_gene"""
        log.info('Loading sequences to reference GEM-PRO...')
        from random import shuffle
        g_ids = [g.id for g in self.reference_gempro.functional_genes]
        shuffle(g_ids)

        def _load_sequences_to_reference_gene_sc(g_id, outdir=self.sequences_by_gene_dir,
                                                 g_to_pickle=self.gene_protein_pickles,
                                                 strain_infodict=self.strain_infodict,
                                                 orth_matrix=self.df_orthology_matrix, force_rerun=force_rerun):
            """Load orthologous strain sequences to reference Protein object, save as new pickle"""
            import ssbio.utils
            import ssbio.io
            from Bio import SeqIO
            import os.path as op

            protein_seqs_pickle_path = op.join(outdir, '{}_protein_withseqs.pckl'.format(g_id))

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_pickle_path):
                protein_pickle_path = g_to_pickle[g_id]
                protein_pickle = ssbio.io.load_pickle(protein_pickle_path)

                for strain, info in strain_infodict.items():
                    strain_sequences = SeqIO.index(info['genome_path'], 'fasta')
                    strain_gene_functional = info['functional_genes'][g_id]
                    if strain_gene_functional:
                        # Pull the gene ID of the strain from the orthology matrix
                        strain_gene_key = orth_matrix.at[g_id, strain]
                        new_id = '{}_{}'.format(g_id, strain)
                        if protein_pickle.sequences.has_id(new_id):
                            continue
                        protein_pickle.load_manual_sequence(seq=strain_sequences[strain_gene_key],
                                                            ident=new_id,
                                                            set_as_representative=False)
                protein_pickle.save_pickle(outfile=protein_seqs_pickle_path)

            return g_id, protein_seqs_pickle_path

        if sc:
            genes_rdd = sc.parallelize(g_ids)
            result = genes_rdd.map(_load_sequences_to_reference_gene_sc).collect()
        else:
            result = []
            for g in tqdm(g_ids):
                result.append(self._load_sequences_to_reference_gene(g, force_rerun))

        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        updated = []
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle
            updated.append(g_id)
        not_updated = set(list(self.gene_protein_pickles.keys())).difference(updated)
        log.info('No change to {} genes, removing from gene_protein_pickles'.format(len(not_updated)))
        log.debug(not_updated)
        for rem in not_updated:
            del self.gene_protein_pickles[rem]

    def store_disorder(self, sc=None, force_rerun=False):
        """Wrapper for _store_disorder"""
        log.info('Loading sequences to reference GEM-PRO...')
        from random import shuffle
        g_ids = [g.id for g in self.reference_gempro.functional_genes]
        shuffle(g_ids)

        def _store_disorder_sc(g_id, outdir=self.sequences_by_gene_dir,
                               g_to_pickle=self.gene_protein_pickles, force_rerun=force_rerun):
            """Load orthologous strain sequences to reference Protein object, save as new pickle"""
            import ssbio.utils
            import ssbio.io
            import os.path as op

            protein_seqs_pickle_path = op.join(outdir, '{}_protein_withseqs_dis.pckl'.format(g_id))

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_pickle_path):
                protein_pickle_path = g_to_pickle[g_id]
                protein_pickle = ssbio.io.load_pickle(protein_pickle_path)
                protein_pickle.get_all_disorder_predictions(representative_only=False)
                protein_pickle.save_pickle(outfile=protein_seqs_pickle_path)

            return g_id, protein_seqs_pickle_path

        if sc:
            genes_rdd = sc.parallelize(g_ids)
            result = genes_rdd.map(_store_disorder_sc).collect()
        else:
            result = []
            for g in tqdm(g_ids):
                result.append(self._load_sequences_to_reference_gene(g, force_rerun))

        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle

    def _align_orthologous_gene_pairwise(self, g_id, gapopen=10, gapextend=0.5, engine='needle', parse=True, force_rerun=False):
        """Align orthologous strain sequences to representative Protein sequence, save as new pickle"""
        protein_seqs_aln_pickle_path = op.join(self.sequences_by_gene_dir, '{}_protein_withseqs_dis_aln.pckl'.format(g_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_aln_pickle_path):
            protein_seqs_pickle_path = self.gene_protein_pickles[g_id]
            protein_pickle = ssbio.io.load_pickle(protein_seqs_pickle_path)

            if not protein_pickle.representative_sequence:
                log.error('{}: no representative sequence to align to'.format(g_id))
                return

            if len(protein_pickle.sequences) < 1:
                log.error('{}: no other sequences to align to'.format(g_id))
                return

            alignment_dir = op.join(self.sequences_by_gene_dir, g_id)
            ssbio.utils.make_dir(alignment_dir)
            protein_pickle.pairwise_align_sequences_to_representative(gapopen=gapopen, gapextend=gapextend,
                                                                      engine=engine, outdir=alignment_dir,
                                                                      parse=parse, force_rerun=force_rerun)
            protein_pickle.save_pickle(outfile=protein_seqs_aln_pickle_path)

        return g_id, protein_seqs_aln_pickle_path

    def align_orthologous_genes_pairwise(self, sc=None, joblib=False, cores=1, gapopen=10, gapextend=0.5,
                                         engine='needle', parse=True, force_rerun=False):
        """Wrapper for _align_orthologous_gene_pairwise"""
        log.info('Aligning sequences to reference GEM-PRO...')
        from random import shuffle
        g_ids = [g.id for g in self.reference_gempro.functional_genes]
        shuffle(g_ids)

        def _align_orthologous_gene_pairwise_sc(g_id, g_to_pickle=self.gene_protein_pickles,
                                                gapopen=gapopen, gapextend=gapextend, engine=engine, parse=parse,
                                                outdir=self.sequences_by_gene_dir,
                                                force_rerun=force_rerun):
            """Align orthologous strain sequences to representative Protein sequence, save as new pickle"""
            import ssbio.utils
            import ssbio.io
            import os.path as op

            protein_seqs_aln_pickle_path = op.join(outdir, '{}_protein_withseqs_dis_aln.pckl'.format(g_id))

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_aln_pickle_path):
                protein_seqs_pickle_path = g_to_pickle[g_id]
                protein_pickle = ssbio.io.load_pickle(protein_seqs_pickle_path)

                if not protein_pickle.representative_sequence:
                    return

                if len(protein_pickle.sequences) < 1:
                    return

                alignment_dir = op.join(outdir, g_id)
                ssbio.utils.make_dir(alignment_dir)
                protein_pickle.pairwise_align_sequences_to_representative(gapopen=gapopen, gapextend=gapextend,
                                                                          engine=engine, outdir=alignment_dir,
                                                                          parse=parse, force_rerun=force_rerun)
                protein_pickle.save_pickle(outfile=protein_seqs_aln_pickle_path)

            return g_id, protein_seqs_aln_pickle_path

        if sc:
            genes_rdd = sc.parallelize(g_ids)
            result_raw = genes_rdd.map(_align_orthologous_gene_pairwise_sc).collect()
        else:
            result_raw = []
            for g in tqdm(g_ids):
                result_raw.append(self._align_orthologous_gene_pairwise(g, gapopen=gapopen, gapextend=gapextend,
                                                                    engine=engine, parse=parse,
                                                                    force_rerun=force_rerun))

        result = [x for x in result_raw if x is not None]
        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle

    def loadseqstoref_alignorth(self, sc, start=None, end=None, force_rerun=False):
        from random import shuffle
        log.info('Loading sequences to reference GEM-PRO and aligning sequences...')

        g_ids = [g.id for g in self.reference_gempro.functional_genes]

        if start and end:
            g_ids = g_ids[start:end]

        shuffle(g_ids)

        def _do_all_shit(g_id, outdir=self.sequences_by_gene_dir, g_to_pickle=self.gene_protein_pickles,
                         strain_infodict=self.strain_infodict,
                         orth_matrix=self.df_orthology_matrix, force_rerun=force_rerun):
            import ssbio.utils
            import ssbio.io
            from Bio import SeqIO
            import os.path as op

            protein_seqs_aln_pickle_path = op.join(outdir, '{}_protein_withseqs_aln.pckl'.format(g_id))

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_seqs_aln_pickle_path):
                # protein_seqs_pickle_path = op.join(outdir, '{}_protein_withseqs.pckl'.format(g_id))
                protein_pickle_path = g_to_pickle[g_id]
                protein_pickle = ssbio.io.load_pickle(protein_pickle_path)
                for strain, info in strain_infodict.items():
                    strain_sequences = SeqIO.index(info['genome_path'], 'fasta')
                    strain_gene_functional = info['functional_genes'][g_id]
                    if strain_gene_functional:
                        # Pull the gene ID of the strain from the orthology matrix
                        strain_gene_key = orth_matrix.at[g_id, strain]
                        new_id = '{}_{}'.format(g_id, strain)
                        if protein_pickle.sequences.has_id(new_id):
                            continue
                        protein_pickle.load_manual_sequence(seq=strain_sequences[strain_gene_key],
                                                            ident=new_id,
                                                            set_as_representative=False)
                # protein_pickle.save_pickle(outfile=protein_seqs_pickle_path)

                # protein_seqs_dis_pickle_path = op.join(outdir, '{}_protein_withseqs_dis.pckl'.format(g_id))
                # protein_pickle.get_all_disorder_predictions(representative_only=False)
                # protein_pickle.save_pickle(outfile=protein_seqs_dis_pickle_path)

                protein_seqs_aln_pickle_path = op.join(outdir, '{}_protein_withseqs_aln.pckl'.format(g_id))
                if not protein_pickle.representative_sequence:
                    return
                if len(protein_pickle.sequences) < 1:
                    return
                alignment_dir = op.join(outdir, g_id)
                ssbio.utils.make_dir(alignment_dir)
                protein_pickle.pairwise_align_sequences_to_representative(gapopen=10, gapextend=0.5,
                                                                          engine='needle', outdir=alignment_dir,
                                                                          parse=True, force_rerun=force_rerun)
                protein_pickle.save_pickle(outfile=protein_seqs_aln_pickle_path)

            return g_id, protein_seqs_aln_pickle_path

        if sc:
            genes_rdd = sc.parallelize(g_ids)
            result_raw = genes_rdd.map(_do_all_shit).collect()
        else:
            result_raw = []
            for g in tqdm(g_ids):
                result_raw.append(_do_all_shit(g))

        result = [x for x in result_raw if x is not None]
        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        updated = []
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle
            updated.append(g_id)
        not_updated = set(list(self.gene_protein_pickles.keys())).difference(updated)
        log.info('No change to {} genes, removing from gene_protein_pickles'.format(len(not_updated)))
        log.debug(not_updated)
        for rem in not_updated:
            del self.gene_protein_pickles[rem]

    def load_protein_pickles(self):
        log.info('Loading new Protein objects back to reference GEM-PRO...')
        for g_id, protein in tqdm(self.gene_protein_pickles.items()):
            g = self.reference_gempro.genes.get_by_id(g_id)
            g.protein = ssbio.io.load_pickle(protein)


def calculate_residue_counts_perstrain(protein_pickle_path, outdir, pdbflex_keys_file, wt_pid_cutoff=None, force_rerun=False):
    """Writes out a feather file for a PROTEIN counting amino acid occurences for ALL STRAINS along with SUBSEQUENCES"""
    from collections import defaultdict
    from ssbio.protein.sequence.seqprop import SeqProp
    from ssbio.protein.sequence.properties.residues import _aa_property_dict_one
    log = logging.getLogger(__name__)

    protein_id = op.splitext(op.basename(protein_pickle_path))[0].split('_')[0]

    protein_df_outfile = op.join(outdir, '{}_protein_strain_properties.fthr'.format(protein_id))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=protein_df_outfile):

        protein = ssbio.io.load_pickle(protein_pickle_path)

        # First calculate disorder cuz i forgot to
        protein.get_all_disorder_predictions(representative_only=True)

        # Then get all subsequences
        all_protein_subseqs = protein.get_all_subsequences(pdbflex_keys_file=pdbflex_keys_file)
        if not all_protein_subseqs:
            log.error('{}: cannot run subsequence calculator'.format(protein.id))
            return

        # Each strain gets a dictionary
        strain_to_infodict = defaultdict(dict)

        for seqprop_to_analyze in protein.sequences:
            if seqprop_to_analyze.id == protein.representative_sequence.id:
                strain_id = 'K12'
            elif type(seqprop_to_analyze) == SeqProp and seqprop_to_analyze.id != protein.id:  # This is to filter out other KEGGProps or UniProtProps
                strain_id = seqprop_to_analyze.id.split('_', 1)[1]  # This split should work for all strains
            else:
                continue

            ## Additional filtering for genes marked as orthologous but actually have large deletions or something
            ## TODO: experiment with other cutoffs?
            if wt_pid_cutoff:
                aln = protein.sequence_alignments.get_by_id('{0}_{0}_{1}'.format(protein.id, seqprop_to_analyze.id))
                if aln.annotations['percent_identity'] < wt_pid_cutoff:
                    continue

            ###### Calculate "all" properties ######
            seqprop_to_analyze.get_biopython_pepstats()

            # [ALL] aa_count
            if 'amino_acids_percent-biop' not in seqprop_to_analyze.annotations:  # May not run if weird amino acids in the sequence
                log.warning('Protein {}, sequence {}: skipping, unable to run Biopython ProteinAnalysis'.format(protein.id,
                                                                                                            seqprop_to_analyze.id))
                continue
            strain_to_infodict[strain_id].update({'aa_count_{}'.format(k): v for k, v in seqprop_to_analyze.annotations['amino_acids_content-biop'].items()})

            # [ALL] aa_count_total
            strain_to_infodict[strain_id]['aa_count_total'] = seqprop_to_analyze.seq_len

            ###### Calculate subsequence properties ######
            for prop, propdict in all_protein_subseqs.items():
                strain_to_infodict[strain_id].update(protein.get_subseq_props(property_dict=propdict, property_name=prop,
                                                                              seqprop=seqprop_to_analyze))

        protein_df = pd.DataFrame(strain_to_infodict)
        protein_df.reset_index().to_feather(protein_df_outfile)

    return protein_pickle_path, protein_df_outfile


import ssbio.io
import pandas as pd
from collections import OrderedDict


def make_mutation_characterization_perstrain_df(protein_pickle_path, pdbflex_keys_file, outdir, wt_pid_cutoff=None, force_rerun=False):
    from ssbio.protein.sequence.properties.residues import characterize_residue_mutation

    protein = ssbio.io.load_pickle(protein_pickle_path)

    # First calculate disorder cuz i forgot to
    protein.get_all_disorder_predictions(representative_only=True)

    # Then get all subsequences
    all_protein_subseqs = protein.get_all_subsequences(pdbflex_keys_file=pdbflex_keys_file)
    if not all_protein_subseqs:
        log.error('{}: cannot run subsequence calculator'.format(protein.id))
        return

    strains_and_dfs = []
    for sa in protein.sequence_alignments:
        if sa.annotations['ssbio_type'] != 'seqalign':
            continue

        try:
            # strain_id = sa.annotations['b_seq'].split('_')[1]
            strain_id = sa.annotations['b_seq'].split('_', 1)[1]
        except IndexError:
            continue

        ## Additional filtering for genes marked as orthologous but actually have large deletions or something
        ## TODO: experiment with other cutoffs?
        if wt_pid_cutoff:
            if sa.annotations['percent_identity'] < wt_pid_cutoff:
                continue

        mutations = sa.annotations['mutations']
        insertions = sa.annotations['insertions']
        deletions = sa.annotations['deletions']
        # XTODO: simple characterization of insertions and deletions locations

        protein_dir = ssbio.utils.make_dir(op.join(outdir, protein.id + '_mutchar'))
        outpath = op.join(protein_dir, strain_id + '_mutchar.fthr')

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
            strain_pre_df = []
            if mutations:
                for change in mutations:
                    infodict = OrderedDict()
                    wt = change[0]
                    resnum = change[1]
                    mut = change[2]
                    infodict['wt'] = wt
                    infodict['resnum'] = resnum
                    infodict['mut'] = mut

                    props = characterize_residue_mutation(wt, mut, use_extended_def=True)

                    ### Grantham
                    infodict['prot_mutdiff_grantham_2D'] = ssbio.protein.sequence.properties.residues.grantham_score(wt, mut)[0]

                    ### General characterization of mutation
                    infodict['general_polar_to_nonpolar'] = False
                    if 'Polar' in props[wt] and 'Non-polar' in props[mut]:
                        infodict['general_polar_to_nonpolar'] = True
                    infodict['general_nonpolar_to_polar'] = False
                    if 'Non-polar' in props[wt] and 'Polar' in props[mut]:
                        infodict['general_nonpolar_to_polar'] = True
                    infodict['general_chrg_to_notchrg'] = False
                    if 'Charged' in props[wt] and 'Charged' not in props[mut]:
                        infodict['general_chrg_to_notchrg'] = True
                    infodict['general_notchrg_to_chrg'] = False
                    if 'Charged' not in props[wt] and 'Charged' in props[mut]:
                        infodict['general_notchrg_to_chrg'] = True
                    infodict['general_notsmall_to_small'] = False
                    if ('Small' not in props[wt] or 'Tiny' not in props[wt]) and ('Small' in props[mut] and 'Tiny' in props[mut]):
                        infodict['general_notsmall_to_small'] = True
                    infodict['general_notbulk_to_bulk'] = False
                    if 'Bulky' not in props[wt] and 'Bulky' not in props[mut]:
                        infodict['general_notbulk_to_bulk'] = True
                    infodict['general_dis_to_notdis'] = False
                    if 'Disorder promoting' in props[wt] and 'Disorder promoting' not in props[mut]:
                        infodict['general_dis_to_notdis'] = True
                    infodict['general_notord_to_ord'] = False
                    if 'Order promoting' not in props[wt] and 'Order promoting' in props[mut]:
                        infodict['general_notord_to_ord'] = True
                    infodict['general_dis_to_ord'] = False
                    if 'Disorder promoting' in props[wt] and 'Order promoting' in props[mut]:
                        infodict['general_dis_to_ord'] = True
                    infodict['general_anyres_to_gprs'] = False
                    if 'Pathogen enriched' not in props[wt] and 'Pathogen enriched' in props[mut]:
                        infodict['general_anyres_to_gprs'] = True
                    # Removed since this is only relevant in TM domains, calculated below
                    # infodict['general_tm_stabilizing'] = False
                    # if 'TM to Thr stabilizing' in props[wt] and 'TM stabilizing' in props[mut]:
                    #     infodict['general_tm_stabilizing'] = True
                    infodict['general_notcarb_to_carb'] = False
                    if 'Carbonylation susceptible' not in props[wt] and 'Carbonylation susceptible' in props[mut]:
                        infodict['general_notcarb_to_carb'] = True
                    infodict['general_carb_to_notcarb'] = False
                    if 'Carbonylation susceptible' in props[wt] and 'Carbonylation susceptible' not in props[mut]:
                        infodict['general_notcarb_to_carb'] = True

                    for aa in ['C', 'M', 'Y']:
                        infodict['general_{0}_to_not{0}'.format(aa)] = False
                        if wt == aa and mut != aa:
                            infodict['general_{0}_to_not{0}'.format(aa)] = True
                        infodict['general_not{0}_to_{0}'.format(aa)] = False
                        if wt != aa and mut == aa:
                            infodict['general_not{0}_to_{0}'.format(aa)] = True

                    ### Characterization within disordered region
                    for region in ['disorder_2D', 'disorder_3D', 'ss_disorder_2D', 'ss_disorder_3D']:
                        infodict['mut_dis_to_ord_{}'.format(region)] = None
                        infodict['mut_notord_to_ord_{}'.format(region)] = None
                        infodict['mut_dis_to_notdis_{}'.format(region)] = None

                        # If there exists no region, do not characterize this mutation as something related to it
                        if region in all_protein_subseqs:
                            # If there exists a region, check to see if the residue number exists in it
                            if resnum in all_protein_subseqs[region]['subseq_resnums']:
                                # mut_dis_to_notdis
                                prefix = 'mut_dis_to_notdis'
                                key = '{}_{}'.format(prefix, region)
                                if key not in infodict: raise KeyError('Forgot to initialize key')
                                infodict[key] = False
                                if 'Disorder promoting' in props[wt] and 'Disorder promoting' not in props[mut]:
                                    infodict[key] = True

                                # mut_notord_to_ord
                                prefix = 'mut_notord_to_ord'
                                key = '{}_{}'.format(prefix, region)
                                if key not in infodict: raise KeyError('Forgot to initialize key')
                                infodict[key] = False
                                if 'Order promoting' not in props[wt] and 'Order promoting' in props[mut]:
                                    infodict[key] = True

                                # mut_dis_to_ord
                                prefix = 'mut_dis_to_ord'
                                key = '{}_{}'.format(prefix, region)
                                if key not in infodict: raise KeyError('Forgot to initialize key')
                                infodict[key] = False
                                if 'Disorder promoting' in props[wt] and 'Order promoting' in props[mut]:
                                    infodict[key] = True

                    ### Characterization within surface exposed residues
                    for region in ['acc_2D', 'acc_3D', 'surface_3D']:
                        infodict['mut_carb_to_notcarb_{}'.format(region)] = None
                        infodict['mut_cys_to_notcys_{}'.format(region)] = None
                        infodict['mut_notmet_to_met_{}'.format(region)] = None
                        infodict['mut_chrg_to_notchrg_{}'.format(region)] = None
                        infodict['mut_negchrg_to_notnegchrg_{}'.format(region)] = None
                        infodict['mut_poschrg_to_notposchrg_{}'.format(region)] = None

                        # If there exists no region, do not characterize this mutation as something related to it
                        if region in all_protein_subseqs:
                            # If there exists a region, check to see if the residue number exists in it
                            if resnum in all_protein_subseqs[region]['subseq_resnums']:
                                # mut_carb_to_notcarb
                                prefix = 'mut_carb_to_notcarb'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Carbonylation susceptible' in props[wt] and 'Carbonylation susceptible' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                for aa in ['C', 'M', 'Y']:
                                    prefix = 'mut_{0}_to_not{0}'.format(aa)
                                    key = '{}_{}'.format(prefix, region)
                                    infodict[key] = False
                                    if wt == aa and mut != aa:
                                        infodict[key] = True

                                    prefix = 'mut_not{0}_to_{0}'.format(aa)
                                    key = '{}_{}'.format(prefix, region)
                                    infodict[key] = False
                                    if wt != aa and mut == aa:
                                        infodict[key] = True

                                # mut_chrg_to_notchrg
                                prefix = 'mut_chrg_to_notchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Charged' in props[wt] and 'Charged' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_notposchrg_to_poschrg
                                prefix = 'mut_poschrg_to_notposchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Basic' in props[wt] and 'Basic' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_notnegchrg_to_negchrg
                                prefix = 'mut_negchrg_to_notnegchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Acidic' in props[wt] and 'Acidic' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                    ### Characterization within metal binding regions
                    for region in ['metal_2_5D', 'metal_3D', 'csa_2_5D', 'sites_2_5D']:
                        infodict['mut_notbulk_to_bulk_{}'.format(region)] = None
                        infodict['mut_carb_to_notcarb_{}'.format(region)] = None
                        infodict['mut_cys_to_notcys_{}'.format(region)] = None
                        infodict['mut_notmet_to_met_{}'.format(region)] = None
                        infodict['mut_chrg_to_notchrg_{}'.format(region)] = None
                        infodict['mut_notposchrg_to_poschrg_{}'.format(region)] = None
                        infodict['mut_notnegchrg_to_negchrg_{}'.format(region)] = None
                        infodict['mut_poschrg_to_notposchrg_{}'.format(region)] = None
                        infodict['mut_negchrg_to_notnegchrg_{}'.format(region)] = None

                        # If there exists no region, do not characterize this mutation as something related to it
                        if region in all_protein_subseqs:
                            # If there exists a region, check to see if the residue number exists in it
                            if resnum in all_protein_subseqs[region]['subseq_resnums']:
                                # mut_notbulk_to_bulk
                                prefix = 'mut_notbulk_to_bulk'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Bulky' not in props[wt] and 'Bulky' in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_carb_to_notcarb
                                prefix = 'mut_carb_to_notcarb'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Carbonylation susceptible' in props[wt] and 'Carbonylation susceptible' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                for aa in ['C', 'M']:
                                    prefix = 'mut_{0}_to_not{0}'.format(aa)
                                    key = '{}_{}'.format(prefix, region)
                                    infodict[key] = False
                                    if wt == aa and mut != aa:
                                        infodict[key] = True

                                    prefix = 'mut_not{0}_to_{0}'.format(aa)
                                    key = '{}_{}'.format(prefix, region)
                                    infodict[key] = False
                                    if wt != aa and mut == aa:
                                        infodict[key] = True

                                # mut_chrg_to_notchrg
                                prefix = 'mut_chrg_to_notchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Charged' in props[wt] and 'Charged' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_notposchrg_to_poschrg
                                prefix = 'mut_notposchrg_to_poschrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Basic' not in props[wt] and 'Basic' in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_notnegchrg_to_negchrg
                                prefix = 'mut_notnegchrg_to_negchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Acidic' not in props[wt] and 'Acidic' in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_poschrg_to_notposchrg
                                prefix = 'mut_poschrg_to_notposchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Basic' in props[wt] and 'Basic' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_negchrg_to_notnegchrg
                                prefix = 'mut_negchrg_to_notnegchrg'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'Acidic' in props[wt] and 'Acidic' not in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                    ### Characterization within TM domain
                    for region in ['tm_2D', 'tm_3D']:
                        infodict['mut_tmunstab_to_tmstab_{}'.format(region)] = None
                        infodict['mut_notmet_to_met_{}'.format(region)] = None

                        # If there exists no region, do not characterize this mutation as something related to it
                        if region in all_protein_subseqs:
                            # If there exists a region, check to see if the residue number exists in it
                            if resnum in all_protein_subseqs[region]['subseq_resnums']:
                                # mut_tmunstab_to_tmstab
                                prefix = 'mut_tmunstab_to_tmstab'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if 'TM to Thr stabilizing' in props[wt] and 'TM stabilizing' in props[mut]:
                                    infodict['{}_{}'.format(prefix, region)] = True

                                # mut_notmet_to_met
                                prefix = 'mut_notmet_to_met'
                                key = '{}_{}'.format(prefix, region)
                                infodict[key] = False
                                if wt != 'M' and mut == 'M':
                                    infodict['{}_{}'.format(prefix, region)] = True


                    ### Singular defined domains
                    #### DNA binding domain
                    dna_region = 'dna_2_5D'
                    infodict['mut_notdis_to_dis_{}'.format(dna_region)] = None
                    infodict['mut_ord_to_notord_{}'.format(dna_region)] = None
                    infodict['mut_ord_to_dis_{}'.format(dna_region)] = None
                    infodict['mut_at_{}'.format(dna_region)] = False
                    if dna_region in all_protein_subseqs:
                        if resnum in all_protein_subseqs[dna_region]['subseq_resnums']:
                            infodict['mut_at_{}'.format(dna_region)] = True
                            if 'Disorder promoting' not in props[wt] and 'Disorder promoting' in props[mut]:
                                infodict['mut_notdis_to_dis_{}'.format(dna_region)] = True
                            else:
                                infodict['mut_notdis_to_dis_{}'.format(dna_region)] = False
                            if 'Order promoting' in props[wt] and 'Order promoting' not in props[mut]:
                                infodict['mut_ord_to_notord_{}'.format(dna_region)] = True
                            else:
                                infodict['mut_ord_to_notord_{}'.format(dna_region)] = False
                            if 'Order promoting' in props[wt] and 'Disorder promoting' in props[mut]:
                                infodict['mut_ord_to_dis_{}'.format(dna_region)] = True
                            else:
                                infodict['mut_ord_to_dis_{}'.format(dna_region)] = False

                    # Now contained in the metal one above
                    # #### Catalytic site
                    # csa_region = 'csa_2_5D'
                    # infodict['mut_at_{}'.format(csa_region)] = False
                    # if csa_region in all_protein_subseqs:
                    #     if resnum in all_protein_subseqs[csa_region]['subseq_resnums']:
                    #         infodict['mut_at_{}'.format(csa_region)] = True
                    #
                    # #### UniProt site
                    # uni_region = 'sites_2_5D'
                    # infodict['mut_at_{}'.format(uni_region)] = False
                    # if uni_region in all_protein_subseqs:
                    #     if resnum in all_protein_subseqs[uni_region]['subseq_resnums']:
                    #         infodict['mut_at_{}'.format(uni_region)] = True

                    ## XTODO: need to add these to all_protein_subseqs and also search within vicinity of them
                    ## XTODO: do not use CarSPred, use CarbonylDB
                    #### Experimentally determined carbonylatable residues
                    # carbonyldb_region = 'carbonyldb_exp'
                    # infodict['mut_at_{}'.format(carbonyldb_region)] = False
                    # if 'experimental_carb_ox-carbonyldb_exp' in protein.representative_sequence.annotations:
                    #     if resnum in protein.representative_sequence.annotations['experimental_cys_ox-carbonyldb_exp']:
                    #         infodict['mut_at_{}'.format(carbonyldb_region)] = True

                    #### CarSPred residue
                    cars_region = 'predcarb_2D'
                    infodict['mut_at_{}'.format(cars_region)] = False
                    if 'predicted_carbonylated-carspred' in protein.representative_sequence.letter_annotations:
                        subfeat_indices = list(locate(protein.representative_sequence.letter_annotations['predicted_carbonylated-carspred'],
                                                      lambda x: x != (False, 0.0)))
                        subfeat_resnums = [x+1 for x in subfeat_indices]
                        if resnum in subfeat_resnums:
                            infodict['mut_at_{}'.format(cars_region)] = True

                    #### Experimentally determined oxidizable cysteines
                    redoxdb_region = 'redoxdb_exp'
                    infodict['mut_at_{}'.format(redoxdb_region)] = False
                    if 'experimental_cys_ox-redoxdb' in protein.representative_sequence.annotations:
                        if resnum in protein.representative_sequence.annotations['experimental_cys_ox-redoxdb']:
                            infodict['mut_at_{}'.format(redoxdb_region)] = True

                    #### Predicted disulfide bridge
                    disbrdg_region = 'disbrdg_3D'
                    infodict['mut_at_{}'.format(disbrdg_region)] = False
                    if protein.representative_structure:
                        if 'SSBOND-biopython' in protein.get_representative_chain.seq_record.annotations:
                            struct_ssbonds = []
                            for x in protein.get_representative_chain.seq_record.annotations['SSBOND-biopython']:
                                for y in x:
                                    struct_ssbonds.append(y[1])
                            mapped = protein.map_structprop_resnums_to_seqprop_resnums(struct_ssbonds, use_representatives=True)
                            if resnum in list(mapped.values()):
                                infodict['mut_at_{}'.format(disbrdg_region)] = True

                    #### Any available FoldX info
                    foldx_suffix = 'foldx_3D'
                    infodict['prop_stability_{}'.format(foldx_suffix)] = None
                    if 'ecoref_foldx_predictions' in protein.representative_sequence.annotations:
                        for k,v in protein.representative_sequence.annotations['ecoref_foldx_predictions'].items():
                            foldx_calculated = ssbio.utils.split_mutation_string(k)
                            if wt == foldx_calculated[0] and mut == foldx_calculated[2] and str(resnum) == foldx_calculated[1]:
                                infodict['stability_{}'.format(foldx_suffix)] = v

                    strain_pre_df.append(infodict)

            strain_df = pd.DataFrame(strain_pre_df)
            strain_df.index = [protein.id] * len(strain_df)

            # If you want to save them use the code below but remember to require an outdir as well as write the force rerun code
            strain_df.reset_index().to_feather(outpath)

        else:
            strain_df = pd.read_feather(outpath)
            strain_df = strain_df.set_index(strain_df.columns[0])

        strains_and_dfs.append((strain_id, strain_df))

    return strains_and_dfs



mut_color_mapping_human = {'mut_dis_to_ord_disorder_2D': 'deep aqua',
'mut_notord_to_ord_disorder_2D': 'blue blue',
'mut_dis_to_notdis_disorder_2D': 'deep aqua',
'mut_dis_to_ord_disorder_3D': 'deep aqua',
'mut_notord_to_ord_disorder_3D': 'blue blue',
'mut_dis_to_notdis_disorder_3D': 'deep aqua',
'mut_dis_to_ord_ss_disorder_2D': 'deep aqua',
'mut_notord_to_ord_ss_disorder_2D': 'blue blue',
'mut_dis_to_notdis_ss_disorder_2D': 'deep aqua',
'mut_dis_to_ord_ss_disorder_3D': 'deep aqua',
'mut_notord_to_ord_ss_disorder_3D': 'blue blue',
'mut_dis_to_notdis_ss_disorder_3D': 'deep aqua',
'mut_carb_to_notcarb_acc_2D': 'lime green',
'mut_C_to_notC_acc_2D': 'violet',
'mut_Y_to_notY_acc_2D': 'black',
'mut_notM_to_M_acc_2D': 'orange',
'mut_chrg_to_notchrg_acc_2D': 'hot pink',
'mut_carb_to_notcarb_acc_3D': 'lime green',
'mut_C_to_notC_acc_3D': 'violet',
'mut_Y_to_notY_acc_3D': 'black',
'mut_notM_to_M_acc_3D': 'orange',
'mut_chrg_to_notchrg_acc_3D': 'hot pink',
'mut_carb_to_notcarb_surface_3D': 'lime green',
'mut_C_to_notC_surface_3D': 'violet',
'mut_Y_to_notY_surface_3D': 'black',
'mut_notM_to_M_surface_3D': 'orange',
'mut_chrg_to_notchrg_surface_3D': 'hot pink',
'mut_notbulk_to_bulk_metal_2_5D': 'pale yellow',
'mut_carb_to_notcarb_metal_2_5D': 'lime green',
'mut_C_to_notC_metal_2_5D': 'violet',
'mut_notM_to_M_metal_2_5D': 'orange',
'mut_chrg_to_notchrg_metal_2_5D': 'hot pink',
'mut_notposchrg_to_poschrg_metal_2_5D': 'light pink',
'mut_notnegchrg_to_negchrg_metal_2_5D': 'baby blue',
'mut_negchrg_to_notnegchrg_metal_2_5D': 'light pink',
'mut_negchrg_to_notnegchrg_metal_3D': 'light pink',
'mut_poschrg_to_notposchrg_metal_2_5D': 'baby blue',
'mut_poschrg_to_notposchrg_metal_3D': 'baby blue',
'mut_notbulk_to_bulk_metal_3D': 'pale yellow',
'mut_carb_to_notcarb_metal_3D': 'lime green',
'mut_C_to_notC_metal_3D': 'violet',
'mut_notM_to_M_metal_3D': 'orange',
'mut_chrg_to_notchrg_metal_3D': 'hot pink',
'mut_notposchrg_to_poschrg_metal_3D': 'light pink',
'mut_notnegchrg_to_negchrg_metal_3D': 'baby blue',
'mut_tmunstab_to_tmstab_tm_2D': 'bright yellow',
'mut_notM_to_M_tm_2D': 'orange',
'mut_tmunstab_to_tmstab_tm_3D': 'bright yellow',
'mut_notM_to_M_tm_3D': 'orange',
'mut_notdis_to_dis_dna_2_5D': 'bright teal',
'mut_ord_to_notord_dna_2_5D': 'bright teal',
'mut_ord_to_dis_dna_2_5D': 'bright teal',
'mut_at_dna_2_5D': 'black',
'mut_at_csa_2_5D': 'black',
'mut_at_sites_2_5D': 'black',
'mut_at_predcarb_2D': 'black',
'mut_at_redoxdb_exp': 'black',
'mut_at_disbrdg_3D': 'black',
                           'general_notcarb_to_carb': 'red',
                           'general_C_to_notC'  : 'violet',
'general_Y_to_notY'  : 'black',
                           'general_notM_to_M'  : 'orange',
                           'general_notC_to_C':'black',
                           'general_dis_to_notdis': 'white',
                           'general_notord_to_ord': 'white'}

import seaborn as sns
mut_color_mapping = {k: sns.xkcd_palette([v]).as_hex()[0] for k,v in mut_color_mapping_human.items()}