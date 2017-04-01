import os.path as op
import unittest

from Bio.Align import MultipleSeqAlignment

import ssbio.protein.sequence.utils.alignment


class TestAlignment(unittest.TestCase):
    """Unit tests for Alignment"""

    def test_pairwise_sequence_alignment_str(self):
        base_id = 'a_test'
        base = 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVAVVASASA'
        muts_id = 'b_test'
        muts = 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVARRVAASSSAS'

        alignment = ssbio.protein.sequence.utils.alignment.pairwise_sequence_alignment(a_seq_id=base_id, a_seq=base,
                                                                                       b_seq_id=muts_id, b_seq=muts,
                                                                                       engine='biopython')
        self.assertTrue(isinstance(alignment, MultipleSeqAlignment))
        self.assertEqual(alignment[0].id, base_id)
        self.assertEqual(alignment[0].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVAV--VA-SASA-')
        self.assertEqual(alignment[1].id, muts_id)
        self.assertEqual(alignment[1].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVA-RRVAASSSAS')
        self.assertEqual(alignment.annotations, {'end': 51, 'percent_identity': 88.23529411764706, 'score': 224.0, 'start': 0})

        alignment2 = ssbio.protein.sequence.utils.alignment.pairwise_sequence_alignment(a_seq_id=base_id, a_seq=base,
                                                                                        b_seq_id=muts_id, b_seq=muts,
                                                                                        engine='needle', force_rerun=True)
        self.assertTrue(isinstance(alignment2, MultipleSeqAlignment))
        self.assertEqual(alignment2[0].id, base_id)
        self.assertEqual(alignment2[0].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVA--VVASASA-')
        self.assertEqual(alignment2[1].id, muts_id)
        self.assertEqual(alignment2[1].seq, 'MSPRVGVTLSGRYRLQRLIATGGMGQVWEAVDNRLGRRVARRVAASSSAS')
        self.assertEqual(alignment2.annotations, {'score': 213.5, 'percent_gaps': 6.0, 'percent_similarity': 92.0, 'percent_identity': 90.0})

    def test_run_needle_alignment_on_files(self):
        """Test if needle alignment runs correctly on 2 input files

        """
        outdir = op.join('test_files', 'out')
        working_dir = op.join('test_files', 'sequences')
        infile_a = op.join(working_dir, 'P9WGE7.fasta')
        infile_b = op.join(working_dir, '1gn3_A.faa')
        test_outfile = op.join(working_dir, 'P9WGE7_1gn3_A_align.txt')

        outfile = ssbio.protein.sequence.utils.alignment.run_needle_alignment_on_files(id_a='P9WGE7', faa_a=infile_a,
                                                                                       id_b='1gn3_A', faa_b=infile_b,
                                                                                       outdir=outdir,
                                                                                       outfile='test_alignment.txt',
                                                                                       force_rerun=True)

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
        outdir = op.join('test_files', 'out')
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

        outfile = op.join(outdir, 'test_alignment2.txt')

        outpath = ssbio.protein.sequence.utils.alignment.run_needle_alignment(seq_a=str1, seq_b=str2, outfile=outfile)

        # test if the result is the same
        # ignore the first 12 lines because it is dated
        with open(outpath, 'r') as f:
            raw_txt = f.read().splitlines()[12:]

        self.assertEqual(raw_txt, result)

        self.assertTrue(op.exists(outfile))