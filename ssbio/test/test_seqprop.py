import pytest
import os.path as op
from Bio import SeqIO
from BCBio import GFF

from ssbio.protein.sequence.seqprop import SeqProp

@pytest.fixture(scope='module')
def sequence_id():
    return 'P0ABP8'

@pytest.fixture(scope='module')
def sequence_file():
    return 'P0ABP8.fasta'

@pytest.fixture(scope='module')
def metadata_file():
    return 'P0ABP8.txt'

@pytest.fixture(scope='module')
def feature_file():
    return 'P0ABP8.gff'

@pytest.fixture(scope='module')
def sequence_path(test_files_sequences, sequence_file):
    return op.join(test_files_sequences, sequence_file)

@pytest.fixture(scope='module')
def metadata_path(test_files_sequences, metadata_file):
    return op.join(test_files_sequences, metadata_file)

@pytest.fixture(scope='module')
def feature_path(test_files_sequences, feature_file):
    return op.join(test_files_sequences, feature_file)

@pytest.fixture(scope='class')
def seqprop_with_id(sequence_id):
    return SeqProp(ident=sequence_id)

@pytest.fixture(scope='class')
def seqprop_with_id_seq(sequence_id):
    return SeqProp(ident=sequence_id,
                   seq='MTESTERSEQUENCE')

@pytest.fixture(scope='class')
def seqprop_with_id_seqpath(sequence_id, sequence_path):
    return SeqProp(ident=sequence_id,
                   sequence_path=sequence_path)

@pytest.fixture(scope='class')
def seqprop_with_id_seqpath_metapath(sequence_id, sequence_path, metadata_path):
    return SeqProp(ident=sequence_id,
                   sequence_path=sequence_path,
                   metadata_path=metadata_path)


class TestSeqPropJustID():
    def test_init(self, seqprop_with_id):
        assert seqprop_with_id.seq_record == None

        assert seqprop_with_id.annotations == {}
        assert seqprop_with_id.letter_annotations == {}
        assert seqprop_with_id.features == []

        assert seqprop_with_id.sequence_file == None
        with pytest.raises(IOError):
            seqprop_with_id.sequence_dir
        with pytest.raises(IOError):
            seqprop_with_id.sequence_path

        assert seqprop_with_id.metadata_file == None
        with pytest.raises(IOError):
            seqprop_with_id.metadata_dir
        with pytest.raises(IOError):
            seqprop_with_id.metadata_path

    def test_load_sequence_path(self, seqprop_with_id, sequence_path):
        # Test file paths
        seqprop_with_id.load_sequence_path(sequence_path)
        assert seqprop_with_id.sequence_file == op.basename(sequence_path)
        assert seqprop_with_id.sequence_path == sequence_path

        # Test SeqRecord
        orig_sr = SeqIO.read(sequence_path, 'fasta')
        assert seqprop_with_id.description == orig_sr.description
        assert seqprop_with_id.annotations == orig_sr.annotations
        assert seqprop_with_id.letter_annotations == orig_sr.letter_annotations
        assert seqprop_with_id.features == orig_sr.features
        assert seqprop_with_id.seq_str == str(orig_sr.seq)

    def test_load_metadata_path(self, seqprop_with_id, metadata_path):
        # Test file paths
        seqprop_with_id.load_metadata_path(metadata_path)
        assert seqprop_with_id.metadata_file == op.basename(metadata_path)
        assert seqprop_with_id.metadata_path == metadata_path

    def test_load_feature_path(self, seqprop_with_id, feature_path):
        # Test file paths
        seqprop_with_id.load_feature_path(feature_path)
        assert seqprop_with_id.feature_file == op.basename(feature_path)
        assert seqprop_with_id.feature_path == feature_path

        # Test features
        assert seqprop_with_id.features != []
        with open(feature_path) as handle:
            feats = list(GFF.parse(handle))
        assert len(seqprop_with_id.features) == len(feats[0].features)

    def test_change_root_dir(self, tmpdir, seqprop_with_id):
        # Test changing the sequence and metadata dirs
        t = tmpdir.strpath
        seqprop_with_id.sequence_dir = t
        seqprop_with_id.metadata_dir = t
        seqprop_with_id.feature_dir = t

        assert seqprop_with_id.sequence_dir == t
        with pytest.raises(OSError):  # Sequence path should throw an error since the file was never moved there
            seqprop_with_id.sequence_path

        assert seqprop_with_id.metadata_dir == t
        with pytest.raises(OSError):  # Metadata path should throw an error since the file was never moved there
            seqprop_with_id.metadata_path

        assert seqprop_with_id.feature_dir == t
        with pytest.raises(OSError):  # Metadata path should throw an error since the file was never moved there
            seqprop_with_id.feature_path



class TestSeqPropWithSeqFile():
    pass

# class TestSeqProp(unittest.TestCase):
#     """Unit tests for SeqProp"""
#
#     @classmethod
#     def setUpClass(self):
#         self.sp = SeqProp(ident='P0ABP8', sequence_path='test_files/sequences/P0ABP8.fasta',
#                           metadata_path='test_files/sequences/P0ABP8.txt')
#
#     def test_load_seq_file(self):
#         self.assertEqual(self.sp.sequence_file, 'P0ABP8.fasta')
#         self.assertEqual(self.sp.seq_len, 239)
#
#     def test_load_metadata_file(self):
#         self.assertEqual(self.sp.metadata_file, 'P0ABP8.txt')
#
#     def test_get_seq_str(self):
#         seq = 'MATPHINAEMGDFADVVLMPGDPLRAKYIAETFLEDAREVNNVRGMLGFTGTYKGRKISVMGHGMGIPSCSIYTKELITDFGVKKIIRVGSCGAVLPHVKLRDVVIGMGACTDSKVNRIRFKDHDFAAIADFDMVRNAVDAAKALGIDARVGNLFSADLFYSPDGEMFDVMEKYGILGVEMEAAGIYGVAAEFGAKALTICTVSDHIRTHEQTTAAERQTTFNDMIKIALESVLLGDKE'
#         self.assertTrue(self.sp.seq_str)
#         self.assertEqual(self.sp.seq_str, seq)
#
#     def test_equal_to(self):
#         newsp = SeqProp(ident='P0ABP8', sequence_path='test_files/sequences/P0ABP8.fasta',
#                         metadata_path='test_files/sequences/P0ABP8.txt')
#         self.assertTrue(self.sp.equal_to(newsp))
#
#     def test_equal_to_fasta(self):
#         self.assertTrue(self.sp.equal_to_fasta('test_files/sequences/P0ABP8.fasta'))
#
#     def test_num_pdbs(self):
#         self.assertEqual(self.sp.num_pdbs, 0)
#         self.sp.pdbs = ['1abc', '2abc']
#         self.assertEqual(self.sp.num_pdbs, 2)
#
#     def test_write_fasta_file(self):
#         self.sp.write_fasta_file(outfile='test_files/out/P0ABP8_new.faa', force_rerun=True)
#         self.assertTrue(op.exists('test_files/out/P0ABP8_new.faa'))
#         self.assertEqual(self.sp.sequence_file, 'P0ABP8_new.faa')