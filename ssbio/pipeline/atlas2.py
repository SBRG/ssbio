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
                             pid_cutoff=None, bitscore_cutoff=None, evalue_cutoff=None, filter_condition='OR',
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

    def load_sequences_to_reference(self, sc=None, joblib=False, cores=1, force_rerun=False):
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
        # elif joblib:
        #     result = Parallel(n_jobs=cores)(delayed(self._load_sequences_to_reference_gene)(g, force_rerun) for g in g_ids)
        else:
            result = []
            for g in tqdm(g_ids):
                result.append(self._load_sequences_to_reference_gene(g, force_rerun))

        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle

    def _align_orthologous_gene_pairwise(self, g_id, gapopen=10, gapextend=0.5, engine='needle', parse=True, force_rerun=False):
        """Align orthologous strain sequences to representative Protein sequence, save as new pickle"""
        protein_seqs_aln_pickle_path = op.join(self.sequences_by_gene_dir, '{}_protein_withseqs_aln.pckl'.format(g_id))

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

            protein_seqs_aln_pickle_path = op.join(outdir, '{}_protein_withseqs_aln.pckl'.format(g_id))

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
            result = [x for x in result_raw if x is not None]
        # elif joblib:
        #     result = Parallel(n_jobs=cores)(delayed(self._align_orthologous_gene_pairwise)(g,
        #                                                                                    gapopen=gapopen,
        #                                                                                    gapextend=gapextend,
        #                                                                                    engine=engine,
        #                                                                                    parse=parse,
        #                                                                                    force_rerun=force_rerun) for g in g_ids)
        else:
            result = []
            for g in tqdm(g_ids):
                result.append(self._align_orthologous_gene_pairwise(g, gapopen=gapopen, gapextend=gapextend,
                                                                    engine=engine, parse=parse,
                                                                    force_rerun=force_rerun))

        log.info('Storing paths to new Protein objects in self.gene_protein_pickles...')
        for g_id, protein_pickle in result:
            self.gene_protein_pickles[g_id] = protein_pickle

    def load_protein_pickles(self):
        log.info('Loading new Protein objects back to reference GEM-PRO...')
        for g_id, protein in tqdm(self.gene_protein_pickles.items()):
            g = self.reference_gempro.genes.get_by_id(g_id)
            g.protein = ssbio.io.load_pickle(protein)

    def get_atlas_summary_df(self):
        """Create a single data frame which summarizes all genes per row.

        Returns:
            DataFrame: Pandas DataFrame of the results

        """
        all_info = []
        for g in self.reference_gempro.functional_genes:
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
        cols.extend([x.id for x in self.strain_ids])

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