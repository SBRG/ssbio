import unittest

import ssbio.utils

class TestUtils(unittest.TestCase):
    """Unit tests for utils
    """

    def test_split_folder_and_path(self):
        test_path = '/this/is/a/hypothetical/path/to_a_file.myfile'
        result = ('/this/is/a/hypothetical/path', 'to_a_file', '.myfile')

        self.assertEqual(result, ssbio.utils.split_folder_and_path(test_path))

    def test_outfile_name_maker(self):
        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta')
        correct_out = 'P00001.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta', outext='.mao')
        correct_out = 'P00001.mao'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta', outext='.new', outfile='P00001_aligned')
        correct_out = 'P00001_aligned.new'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta', outfile='P00001_aligned')
        correct_out = 'P00001_aligned.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta', outfile='P00001_aligned', outdir='/my/dir/')
        correct_out = '/my/dir/P00001_aligned.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_name_maker(infile='P00001.fasta', outfile='P00001_aligned', outdir='/my/dir/', outext='.test')
        correct_out = '/my/dir/P00001_aligned.test'
        self.assertEqual(test_out, correct_out)