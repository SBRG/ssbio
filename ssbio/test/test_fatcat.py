import unittest
import os.path as op
import ssbio.protein.structure.properties.fatcat as fatcat


class TestFATCAT(unittest.TestCase):
    """Unit tests for FATCAT"""

    def test_run_fatcat_all_by_all(self):
        OUT_DIR = op.join('test_files', 'out')
        FATCAT_SH = '/home/nathan/software/fatcat/runFATCAT.sh'

        structs = [op.join('test_files', 'structures', '12as-A_clean.pdb'),
                   op.join('test_files', 'structures', '1af6-A_clean.pdb'),
                   op.join('test_files', 'structures', '1a9x-A_clean.pdb')]

        tm_scores = fatcat.run_fatcat_all_by_all(structs, fatcat_sh=FATCAT_SH, outdir=OUT_DIR)
