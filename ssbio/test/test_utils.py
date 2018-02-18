import pytest
import unittest
import os.path as op
import ssbio.utils


class Dummy(object):
    def __init__(self):
        self.outdir_set = 'an_output_directory'
        self.outdir_notset = None


@pytest.fixture(scope='class')
def dummy_object():
    return Dummy()


def test_double_check_attribute():
    dummy = dummy_object()
    ssbio.utils.double_check_attribute(object=dummy, setter=None, backup_attribute='outdir_set')
    ssbio.utils.double_check_attribute(object=dummy, setter='newoutdir', backup_attribute='outdir_set')
    with pytest.raises(ValueError):
        ssbio.utils.double_check_attribute(object=dummy, setter=None, backup_attribute='outdir_notset')


class TestUtils(unittest.TestCase):
    """Unit tests for utils
    """

    def test_split_folder_and_path(self):
        test_path = op.join('this','is','a','hypothetical','path','to_a_file.myfile')
        result = (op.join('this','is','a','hypothetical','path'), 'to_a_file', '.myfile')

        self.assertEqual(result, ssbio.utils.split_folder_and_path(test_path))

    def test_force_rerun(self):
        # TODO: need to use op.join for the files, for tests to work properly on windows
        self.assertTrue(ssbio.utils.force_rerun(flag=True, outfile=op.join('not','existing','file.txt')))

        self.assertTrue(ssbio.utils.force_rerun(flag=False, outfile=op.join('not','existing','file.txt')))

        self.assertTrue(ssbio.utils.force_rerun(flag=True, outfile=op.join('test_files','Ec_core_flux1.xml')))

        self.assertFalse(ssbio.utils.force_rerun(flag=False, outfile=op.join('test_files','Ec_core_flux1.xml')))

    def test_outfile_name_maker(self):
        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta')
        correct_out = 'P00001.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', append_to_name='_new')
        correct_out = 'P00001_new.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outext='.mao')
        correct_out = 'P00001.mao'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outext='.mao', append_to_name='_new')
        correct_out = 'P00001_new.mao'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outext='.new', outname='P00001_aligned')
        correct_out = 'P00001_aligned.new'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outname='P00001_aligned')
        correct_out = 'P00001_aligned.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outname='P00001_aligned', append_to_name='_new')
        correct_out = 'P00001_aligned_new.out'
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outname='P00001_aligned', outdir=op.join('my','dir'))
        correct_out = op.join('my','dir','P00001_aligned.out')
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname='P00001.fasta', outname='P00001_aligned', outdir=op.join('my','dir'), outext='.test')
        correct_out = op.join('my','dir','P00001_aligned.test')
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname=op.join('test','other','dir','P00001.fasta'), append_to_name='_new')
        correct_out = op.join('test','other','dir','P00001_new.out')
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname=op.join('test','other','dir','P00001.fasta'), outname='P00001_aligned')
        correct_out = op.join('test','other','dir','P00001_aligned.out')
        self.assertEqual(test_out, correct_out)

        test_out = ssbio.utils.outfile_maker(inname=op.join('test','other','dir','P00001.fasta'), outname='P00001_aligned', outdir=op.join('my','dir'))
        correct_out = op.join('my','dir','P00001_aligned.out')
        self.assertEqual(test_out, correct_out)