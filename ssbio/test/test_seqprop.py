import unittest

from ssbio.sequence import SeqProp

class TestSeqProp(unittest.TestCase):
    """Unit tests for SeqProp"""

    @classmethod
    def setUpClass(self):
        self.sp = SeqProp(ident='P0ABP8',
                          sequence_file='test_files/sequences/P0ABP8.fasta',
                          metadata_file='test_files/sequences/P0ABP8.txt')

    def test_load_seq_file(self):
        self.assertEqual(self.sp.sequence_file, 'P0ABP8.fasta')
        self.assertEqual(self.sp.sequence_len, 239)

    def test_load_metadata_file(self):
        self.assertEqual(self.sp.metadata_file, 'P0ABP8.txt')

    def test_get_seq_str(self):
        seq = 'MATPHINAEMGDFADVVLMPGDPLRAKYIAETFLEDAREVNNVRGMLGFTGTYKGRKISVMGHGMGIPSCSIYTKELITDFGVKKIIRVGSCGAVLPHVKLRDVVIGMGACTDSKVNRIRFKDHDFAAIADFDMVRNAVDAAKALGIDARVGNLFSADLFYSPDGEMFDVMEKYGILGVEMEAAGIYGVAAEFGAKALTICTVSDHIRTHEQTTAAERQTTFNDMIKIALESVLLGDKE'
        my_seq_str = self.sp.get_seq_str()
        self.assertTrue(self.sp.seq_str)
        self.assertEqual(self.sp.seq_str, seq)

    def test_get_seq_record(self):
        my_seq_record = self.sp.get_seq_record()
        self.assertTrue(self.sp.seq_record)

    def test_equal_to(self):
        newsp = SeqProp(ident='P0ABP8',
                          sequence_file='test_files/sequences/P0ABP8.fasta',
                          metadata_file='test_files/sequences/P0ABP8.txt')
        self.assertTrue(self.sp.equal_to(newsp))

    def test_equal_to_fasta(self):
        self.assertTrue(self.sp.equal_to_fasta('test_files/sequences/P0ABP8.fasta'))

    def test_num_pdbs(self):
        self.assertEqual(self.sp.num_pdbs(), 0)
        self.sp.pdbs = ['1abc', '2abc']
        self.assertEqual(self.sp.num_pdbs(), 2)