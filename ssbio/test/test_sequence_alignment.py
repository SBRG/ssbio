import os.path as op
import unittest
import tempfile
import ssbio.sequence.alignment

class TestSequenceAlignment(unittest.TestCase):
    """Unit tests for ssbio.sequence.alignment
    """

    def test_run_needle_alignment_on_files(self):
        """Test if needle alignment runs correctly on 2 input files

        """
        working_dir = 'test_sequences'
        infile_a = op.join(working_dir, 'P9WGE7.fasta')
        infile_b = op.join(working_dir, '1gn3_A.faa')
        test_outfile = op.join(working_dir, 'P9WGE7_1gn3_A_align.txt')

        raw_txt = ssbio.sequence.alignment.run_needle_alignment_on_files(id_a='P9WGE7', faa_a=infile_a,
                                                                         id_b='1gn3_A', faa_b=infile_b,
                                                                         outdir=tempfile.gettempdir(),
                                                                         outfile='test_alignment.txt')

        outfile = op.join(tempfile.gettempdir(), 'test_alignment.txt')

        # test if the file content is correct
        # ignore the first 12 lines because it is dated
        with open(outfile) as f:
            test_lines = f.readlines()[12:]
        with open(test_outfile) as f:
            correct_lines = f.readlines()[12:]

        self.assertEqual(test_lines, correct_lines)


    def test_run_needle_alignment_on_str(self):
        """Test if needle alignment runs correctly on 2 input strings

        """
        str1 = 'MAEYTLPDLDWDYGALEPHISGQINELHHSKHHATYVKGANDAVAKLEEA'
        str2 = 'XAEYTLPDLDWDYGALEPHISGQINELYHSKHHANDVKGANDAVAKLEEA'

        result = ['########################################',
                 '',
                 '#=======================================',
                 '#',
                 '# Aligned_sequences: 2',
                 '# 1: asis',
                 '# 2: asis',
                 '# Matrix: EBLOSUM62',
                 '# Gap_penalty: 10.0',
                 '# Extend_penalty: 0.5',
                 '#',
                 '# Length: 50',
                 '# Identity:      46/50 (92.0%)',
                 '# Similarity:    47/50 (94.0%)',
                 '# Gaps:           0/50 ( 0.0%)',
                 '# Score: 245.0',
                 '# ',
                 '#',
                 '#=======================================',
                 '',
                 'asis               1 MAEYTLPDLDWDYGALEPHISGQINELHHSKHHATYVKGANDAVAKLEEA     50',
                 '                     .||||||||||||||||||||||||||:||||||..||||||||||||||',
                 'asis               1 XAEYTLPDLDWDYGALEPHISGQINELYHSKHHANDVKGANDAVAKLEEA     50',
                 '',
                 '',
                 '#---------------------------------------',
                 '#---------------------------------------']

        raw_txt = ssbio.sequence.alignment.run_needle_alignment_on_str(id_a='a', seq_a=str1,
                                                                       id_b='b', seq_b=str2).splitlines()[12:]

        # test if the result is the same
        # ignore the first 12 lines because it is dated
        self.assertEqual(raw_txt, result)

    def test_run_needle_alignment_on_str_plus_outfile(self):
        """Test if needle alignment runs correctly on 2 input strings

        """
        str1 = 'MAEYTLPDLDWDYGALEPHISGQINELHHSKHHATYVKGANDAVAKLEEA'
        str2 = 'XAEYTLPDLDWDYGALEPHISGQINELYHSKHHANDVKGANDAVAKLEEA'

        result = ['',
                  '#=======================================',
                  '#',
                  '# Aligned_sequences: 2',
                  '# 1: asis',
                  '# 2: asis',
                  '# Matrix: EBLOSUM62',
                  '# Gap_penalty: 10.0',
                  '# Extend_penalty: 0.5',
                  '#',
                  '# Length: 50',
                  '# Identity:      46/50 (92.0%)',
                  '# Similarity:    47/50 (94.0%)',
                  '# Gaps:           0/50 ( 0.0%)',
                  '# Score: 245.0',
                  '# ',
                  '#',
                  '#=======================================',
                  '',
                  'asis               1 MAEYTLPDLDWDYGALEPHISGQINELHHSKHHATYVKGANDAVAKLEEA     50',
                  '                     .||||||||||||||||||||||||||:||||||..||||||||||||||',
                  'asis               1 XAEYTLPDLDWDYGALEPHISGQINELYHSKHHANDVKGANDAVAKLEEA     50',
                  '',
                  '',
                  '#---------------------------------------',
                  '#---------------------------------------']

        outfile = op.join(tempfile.gettempdir(), 'test_alignment2.txt')

        raw_txt = ssbio.sequence.alignment.run_needle_alignment_on_str(id_a='a', seq_a=str1,
                                                                       id_b='b', seq_b=str2, outfile=outfile).splitlines()[12:]

        # test if the result is the same
        # ignore the first 12 lines because it is dated
        self.assertEqual(raw_txt, result)

        self.assertTrue(op.exists(outfile))