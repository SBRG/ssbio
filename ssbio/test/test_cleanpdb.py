import os.path as op
import unittest
import tempfile

from ssbio.tools.cleanpdb import CleanPDB as clean


class TestCleanPDB(unittest.TestCase):
    """Unit tests for CleanPDB
    """

    def test_clean_pdb1(self):
        infile = 'structures/1kf6.pdb'
        out_dir = 'structures'
        outfile_basename = '1kf6_clean.pdb'

        my_cleaner = clean(infile)
        default_cleaned_pdb = my_cleaner.clean_pdb(out_dir=tempfile.gettempdir())

        default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

        # test if the filename is correct
        self.assertEqual(default_cleaned_pdb_basename, outfile_basename)

        # test if the file contents are equal
        self.assertEqual(open(op.join(out_dir, '1kf6_clean_tester.pdb'),'r').read(),
                         open(default_cleaned_pdb,'r').read())

        # test that the file does not equal the original file
        self.assertNotEqual(open(op.join(out_dir, '1kf6.pdb'),'r').read(),
                            open(default_cleaned_pdb,'r').read())

    def test_clean_pdb2(self):
        infile = 'structures/PHOE_ECOLI_model1.pdb'
        out_dir = 'structures'
        outfile_basename = 'PHOE_ECOLI_model1_clean.pdb'

        my_cleaner = clean(infile)
        default_cleaned_pdb = my_cleaner.clean_pdb(out_dir=tempfile.gettempdir())

        default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

        # test if the filename is correct
        self.assertEqual(default_cleaned_pdb_basename, outfile_basename)

        # test if the file contents are equal
        self.assertEqual(open(op.join(out_dir, 'PHOE_ECOLI_model1_clean_tester.pdb'), 'r').read(),
                         open(default_cleaned_pdb, 'r').read())

        # test that the file does not equal the original file
        self.assertNotEqual(open(op.join(out_dir, 'PHOE_ECOLI_model1.pdb'), 'r').read(),
                            open(default_cleaned_pdb, 'r').read())

    def test_clean_pdb3(self):
        infile = 'structures/E04142.pdb'
        out_dir = 'structures'
        outfile_basename = 'E04142_clean.pdb'

        my_cleaner = clean(infile)
        default_cleaned_pdb = my_cleaner.clean_pdb(out_dir=tempfile.gettempdir())

        default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

        # test if the filename is correct
        self.assertEqual(default_cleaned_pdb_basename, outfile_basename)

        # test if the file contents are equal
        self.assertEqual(open(op.join(out_dir, 'E04142_clean_tester.pdb'), 'r').read(),
                         open(default_cleaned_pdb, 'r').read())

        # test that the file does not equal the original file
        self.assertNotEqual(open(op.join(out_dir, 'E04142.pdb'), 'r').read(),
                            open(default_cleaned_pdb, 'r').read())

    def test_clean_pdb4(self):
        infile = 'structures/1cbn.pdb'
        out_dir = 'structures'
        outfile_basename = '1cbn_clean.pdb'

        my_cleaner = clean(infile)
        default_cleaned_pdb = my_cleaner.clean_pdb(out_dir=tempfile.gettempdir())

        default_cleaned_pdb_basename = op.basename(default_cleaned_pdb)

        # test if the filename is correct
        self.assertEqual(default_cleaned_pdb_basename, outfile_basename)

        # test if the file contents are equal
        self.assertEqual(open(op.join(out_dir, '1cbn_clean_tester.pdb'), 'r').read(),
                         open(default_cleaned_pdb, 'r').read())

        # test that the file does not equal the original file
        self.assertNotEqual(open(op.join(out_dir, '1cbn.pdb'), 'r').read(),
                            open(default_cleaned_pdb, 'r').read())