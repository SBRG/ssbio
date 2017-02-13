import os.path as op
import unittest
from ssbio.pipeline.gempro import GEMPRO
from ssbio.pipeline.atlas import ATLAS


class TestATLAS(unittest.TestCase):
    """Unit tests for the ATLAS pipeline"""

    @classmethod
    def setUpClass(self):
        GEM_NAME = 'ecoli_test'
        ROOT_DIR = op.join('test_files', 'out')
        gem_file = op.join('test_files', 'Ec_core_flux1.xml')
        my_gempro = GEMPRO(GEM_NAME, ROOT_DIR, gem_file_path=gem_file, gem_file_type='sbml')

        # Prepare your ATLAS analysis
        self.atlas = ATLAS(atlas_name='atlas_test', base_gempro=my_gempro)

    def test_copy_base_gempro(self):
        self.atlas.copy_base_gempro('copied')
        # atlas_strains should now have a GEMPRO with the id 'copied'
        self.assertTrue(self.atlas.atlas_strains.has_id('copied'))
        # The GEMPRO should also append 'copied' to the gene IDs
        self.atlas.atlas_strains.get_by_id('copied').genes[0].id.endswith('_copied')
        self.atlas.atlas_strains.get_by_id('copied').genes[0].protein.id.endswith('_copied')

