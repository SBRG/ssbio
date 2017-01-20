import os.path as op
import tempfile
import unittest

from ssbio.structure.utils.cleanpdb import CleanPDB
from ssbio.structure.utils.mutatepdb import MutatePDB
from ssbio.structure.utils.pdbioext import PDBIOExt


class TestMutatePDB(unittest.TestCase):
    """Unit tests for MutatePDB
    """

    def test_mutate_pdb(self):
        files = [('1kf6.pdb', '1kf6_mutated_tester.pdb')]

        working_dir = 'test_structures'
        out_suffix = '_mutated'
        muts = [('A',0,'Y')]

        for infile, outfile in files:
            outfile_new = '{}{}.pdb'.format(op.splitext(infile)[0], out_suffix)
            infile_path = op.join(working_dir, infile)

            my_pdb = PDBIOExt(infile_path, file_type='pdb')
            custom_mutate = MutatePDB(muts)
            default_mutated_pdb = my_pdb.write_pdb(custom_selection=custom_mutate, out_suffix=out_suffix,
                                                   out_dir=tempfile.gettempdir())
            default_mutated_pdb_basename = op.basename(default_mutated_pdb)

            # test if the filename is correct
            self.assertEqual(default_mutated_pdb_basename, outfile_new)

            # test if the file contents are equal
            self.assertEqual(open(default_mutated_pdb, 'r').read(),
                             open(op.join(working_dir, outfile), 'r').read())

            # test that the file does not equal the original file
            self.assertNotEqual(open(default_mutated_pdb, 'r').read(),
                                open(infile_path, 'r').read())

    def test_clean_and_mutate_pdb(self):
        files = [('1kf6.pdb', '1kf6_clean_mutated_tester.pdb')]

        working_dir = 'test_structures'
        out_suffix = '_clean_mutated'
        muts = [('A', 0, 'Y')]

        for infile, outfile in files:
            outfile_new = '{}{}.pdb'.format(op.splitext(infile)[0], out_suffix)
            infile_path = op.join(working_dir, infile)

            my_pdb = PDBIOExt(infile_path, file_type='pdb')
            my_cleaner = CleanPDB(keep_chains=[m[0] for m in muts])
            my_clean_pdb = my_pdb.write_pdb(out_suffix='_clean', out_dir=tempfile.gettempdir(),
                                            custom_selection=my_cleaner, force_rerun=True)

            my_pdb = PDBIOExt(my_clean_pdb, file_type='pdb')
            custom_mutate = MutatePDB(muts)
            default_mutated_pdb = my_pdb.write_pdb(custom_selection=custom_mutate, out_suffix='_mutated',
                                                   out_dir=tempfile.gettempdir(), force_rerun=True)
            default_mutated_pdb_basename = op.basename(default_mutated_pdb)

            # test if the filename is correct
            self.assertEqual(default_mutated_pdb_basename, outfile_new)

            # test if the file contents are equal
            self.assertEqual(open(default_mutated_pdb, 'r').read(),
                             open(op.join(working_dir, outfile), 'r').read())

            # test that the file does not equal the original file
            self.assertNotEqual(open(default_mutated_pdb, 'r').read(),
                                open(infile_path, 'r').read())