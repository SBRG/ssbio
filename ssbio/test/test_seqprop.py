import unittest
import os.path as op
from ssbio.sequence import SeqProp

class TestSeqProp(unittest.TestCase):
    """Unit tests for SeqProp"""

    @classmethod
    def setUpClass(self):
        self.sp = SeqProp(ident='P0ABP8', sequence_path='test_files/sequences/P0ABP8.fasta',
                          metadata_path='test_files/sequences/P0ABP8.txt')

    def test_load_seq_file(self):
        self.assertEqual(self.sp.sequence_file, 'P0ABP8.fasta')
        self.assertEqual(self.sp.seq_len, 239)

    def test_load_metadata_file(self):
        self.assertEqual(self.sp.metadata_file, 'P0ABP8.txt')

    def test_get_seq_str(self):
        seq = 'MATPHINAEMGDFADVVLMPGDPLRAKYIAETFLEDAREVNNVRGMLGFTGTYKGRKISVMGHGMGIPSCSIYTKELITDFGVKKIIRVGSCGAVLPHVKLRDVVIGMGACTDSKVNRIRFKDHDFAAIADFDMVRNAVDAAKALGIDARVGNLFSADLFYSPDGEMFDVMEKYGILGVEMEAAGIYGVAAEFGAKALTICTVSDHIRTHEQTTAAERQTTFNDMIKIALESVLLGDKE'
        self.assertTrue(self.sp.seq_str)
        self.assertEqual(self.sp.seq_str, seq)

    def test_equal_to(self):
        newsp = SeqProp(ident='P0ABP8', sequence_path='test_files/sequences/P0ABP8.fasta',
                        metadata_path='test_files/sequences/P0ABP8.txt')
        self.assertTrue(self.sp.equal_to(newsp))

    def test_equal_to_fasta(self):
        self.assertTrue(self.sp.equal_to_fasta('test_files/sequences/P0ABP8.fasta'))

    def test_num_pdbs(self):
        self.assertEqual(self.sp.num_pdbs, 0)
        self.sp.pdbs = ['1abc', '2abc']
        self.assertEqual(self.sp.num_pdbs, 2)

    def test_write_fasta_file(self):
        self.sp.write_fasta_file(outname='P0ABP8_new', outdir='test_files/out/', force_rerun=True)
        self.assertTrue(op.exists('test_files/out/P0ABP8_new.faa'))
        self.assertEqual(self.sp.sequence_file, 'P0ABP8_new.faa')