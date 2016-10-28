import os.path as op
import unittest
from ssbio.atlas.pipeline import ATLAS


class TestATLAS(unittest.TestCase):
    """Unit tests for the GEMPRO pipeline
    """

    @classmethod
    def setUpClass(self):
        # self.tempdir = tempfile.TemporaryDirectory()
        # self.ROOT_DIR = self.tempdir.name
        self.ROOT_DIR = 'C:\\Users\\nathan\\Dropbox (UCSD SBRG)\\projects\\'

        # Prepare your ATLAS analysis
        self.atlas = ATLAS(base_strain_name='iJO1366', root_dir=self.ROOT_DIR,
                           seq_type='prot',
                           base_gem_file_path=op.join(self.ROOT_DIR, 'iJO1366', 'atlas\\model_dir\\iJO1366.xml'),
                           base_gem_file_type='sbml',
                           reference_genome='U00096')

        # Input the list of strains you would like to compare
        with open(op.join(self.ROOT_DIR, 'iJO1366\\atlas\\data_dir\\161016-55strains.in')) as f:
            genomes = f.read().splitlines()
        self.atlas.download_genome_cds(genomes, email='nmih@ucsd.edu')

        # Using renamed sequence files for the base strain
        self.atlas.genome_id_to_fasta_file['U00096'] = op.join(self.ROOT_DIR, 'iJO1366\\atlas\\sequence_files\\protein\\by_organism\\U00096_renamed.faa')

        self.atlas.find_bbh()

    def test_build_strain_specific_models(self):
        self.atlas.build_strain_specific_models()

    @classmethod
    def tearDownClass(self):
        # self.tempdir.cleanup()
        pass
