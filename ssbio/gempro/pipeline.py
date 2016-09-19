import argparse
import os
import os.path as op
import shutil
import warnings

from cobra.io import load_matlab_model
from cobra.io.sbml import create_cobra_model_from_sbml_file

from ssbio import utils
from ssbio.cobra.utils import is_spontaneous
from ssbio.databases.identifiers import bioservices_uniprot_mapper
date = utils.Date()

from ssbio.databases.kegg import get_kegg_aa_seq
from tqdm import tqdm
import requests

from ssbio.databases.pdb import top_pdb_blast_hit


class GEMPRO():
    """Generic class to represent all information of a GEM-PRO for a GEM

    Main steps are:
    1. Preparation of folders
    2. Automatic ID mapping
    3. QC/QA
    4. Homology modeling
    5. Representative structures

    Each step may generate a report and also request additional files if mappings are incomplete
    """

    def __init__(self, gem_file, gem_name, root_dir):
        self.root_dir = root_dir
        self.gem_file_temp = gem_file
        self.model_dir = op.join(root_dir, gem_name)
        # model_files - directory where original gems and gem-related files are stored
        self.model_files = op.join(self.model_dir, 'model_files')

    def prep_folders(self):
        """Take in a sequence string and prepares the folder for it to run ITASSER
        """
        # data_frames - directory where all data frames will be stored (all stages)
        self.data_frames = op.join(self.model_dir, 'data_frames')

        # notebooks - directory where ipython notebooks will be stored for manual analyses
        self.notebooks = op.join(self.model_dir, 'notebooks')

        # # figures - directory where all figures will be stored (all stages)
        self.figures = op.join(self.model_dir, 'figures')

        # missing - directory where missing reports go
        self.missing = op.join(self.model_dir, 'missing')

        # struct_files - directory where structure related files will be downloaded/are located
        self.struct_files = op.join(self.model_dir, 'structure_files')
        self.struct_exp_files = op.join(self.struct_files, 'experimental')
        self.struct_homology_files = op.join(self.struct_files, 'homology_models')

        # seq_files - sequence related files are stored here
        self.seq_files = op.join(self.model_dir, 'sequence_files')

        self.seq_main_files = op.join(self.seq_files, 'amino_acid_sequences')
        self.seq_main_metadata = op.join(self.seq_files, 'amino_acid_metadata')

        self.seq_align_files = op.join(self.seq_files, 'alignment')

        self.seq_pdb_files = op.join(self.seq_files, 'pdb_sequences')
        self.seq_homology_files = op.join(self.seq_files, 'homology_model_sequences')

        for directory in [self.model_dir, self.data_frames, self.notebooks, self.figures, self.missing,
                          self.model_files, self.struct_files, self.struct_exp_files,
                          self.struct_homology_files, self.seq_files, self.seq_main_files,
                          self.seq_main_metadata, self.seq_align_files, self.seq_pdb_files,
                          self.seq_homology_files]:
            if not op.exists(directory):
                os.makedirs(directory)
                print('[PREP] Created directory: {}'.format(directory))
            else:
                print('[PREP] Directory already exists: {}'.format(directory))

        if not op.exists(op.join(self.model_files, op.basename(self.gem_file_temp))):
            shutil.copy(self.gem_file_temp, self.model_files)
        self.gem_file = op.join(self.model_files, op.basename(self.gem_file_temp))

    def load_model(self):
        """Load the GEM using COBRApy. Accept SBML or MAT files as input

        Returns: None - loads the model in the GEMPRO object instance

        """
        extension = op.splitext(self.gem_file)[1]
        # TODO: auto parsing should not be done -- explicit is better than implicit
        # TODO: improved parsing of input file?
        if extension.lower() == '.xml':
            self.model = create_cobra_model_from_sbml_file(self.gem_file, print_time=False)
            print('[PREP] Loaded model: {}'.format(self.gem_file))
        elif extension.lower() == '.mat':
            self.model = load_matlab_model(self.gem_file)
            print('[PREP] Loaded model: {}'.format(self.gem_file))

        # obtain list of all gene ids, excluding spontaneous things
        # TODO: account for isoforms
        self.genes = [x.id for x in self.model.genes]
        self.genes_filtered = [y for y in self.genes if not is_spontaneous(y)]

    def kegg_download_all_sequences(self, kegg_organism):
        """Download all sequences from KEGG for all genes in the model

        Args:
            kegg_organism:

        Returns:

        """
        self.gene_to_protein_fasta = {}
        self.genes_to_kegg_missing = []
        for g in tqdm(self.genes_filtered):
            outfile = op.join(self.seq_main_files, '{}-{}.faa'.format(kegg_organism, g))
            if op.exists(outfile):
                self.gene_to_protein_fasta[g] = outfile
                continue
            else:
                try:
                    outfile = get_kegg_aa_seq(kegg_organism, g, out_dir=self.seq_main_files)
                    self.gene_to_protein_fasta[g] = outfile
                except requests.exceptions.HTTPError:
                    self.genes_to_kegg_missing.append(g)

        if self.genes_to_kegg_missing:
            err_file = op.join(self.missing, '{}_kegg_seq_missing.err'.format(date.short_date))
            warnings.warn(
                '[INFO] Some genes were not able to be mapped to a Kegg entry. Please see {} for a list of these genes'.format(
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.genes_to_kegg_missing))

    def blast_pdb_hits(self):
        """Create a dataframe of best PDB hits to each gene's sequence

        Returns: Path to results dataframe

        """
        for g,kg in self.gene_to_protein_fasta:
            fasta_file = op.join(self.seq_main_files, kg)
            top_pdb_blast_hit()

    def set_gene_source(self, db, uniprot_mapping_from=None, isoform_presence=False):

        self.uniprot_mapping_from = uniprot_mapping_from
        self.isoform_presence = isoform_presence
        if self.isoform_presence == True:
            warnings.warn('Isoform automatic mapping is in beta!')

    def gene_id_mapping(self):
        """Take all genes present in the model and map to UniProt and PDB IDs

        Currenty uses the UniProt mapping service through the bioservices package.

        Returns: Path to mapping dataframe

        """
        # allow for input of manual mapping files
        # format should be csv of
        # MODEL_NAME, MAPPED_SOURCE1[, MAPPED_SOURCE2, ...]
        # MODEL_GENE, MAPPED_GENE1[, MAPPED_GENEID2, ...]



        # mapping gene IDs to uniprot IDs
        self.genes_to_uniprots = bioservices_uniprot_mapper(self.uniprot_mapping_from, 'ACC', genes_filtered)
        self.genes_to_uniprots_missing = list(set(genes_filtered).difference(self.genes_to_uniprots.keys()))

        # report genes unable to be mapped
        if self.genes_to_uniprots_missing:
            err_file = op.join(self.missing, '{}_uniprot_mapping_missing.err'.format(date.short_date))
            warnings.warn(
                '[INFO] Some genes were not able to be mapped to a UniProt entry. Please see {} for a list of these genes'.format(
                    err_file))
            with open(err_file, mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(self.genes_to_uniprots_missing))

            # TODO: accept manual input of a file with the correct mappings?

    def run_pipeline(self):
        """Run the entire GEM-PRO pipel ine.

        Options include:
        ...

        Returns:

        """
        self.prep_folders()
        self.load_model()
        # map gene IDs
        # check for missing maps, request manual mapping file?
        pass


if __name__ == '__main__':
    # run the GEM-PRO pipeline!

    # parse arguments
    p = argparse.ArgumentParser(description='Runs the GEM-PRO pipeline')
    p.add_argument('gemfile', help='Path to the GEM file')
    p.add_argument('gemname', help='Name you would like to use to refer to this GEM')
    p.add_argument('rootdir', help='Directory where GEM-PRO files should be stored')

    args = p.parse_args()

    my_gem_pro = GEMPRO(args.gemfile, args.gemname, args.rootdir)
    my_gem_pro.run_pipeline()
