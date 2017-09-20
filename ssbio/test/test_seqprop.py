import pytest
import os.path as op
from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from ssbio.protein.sequence.seqprop import SeqProp

@pytest.fixture(scope='class')
def seq_str_example():
    """Dummy sequence string to load"""
    return "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"

@pytest.fixture(scope='class')
def seq_record_example():
    """Dummy SeqRecord to load"""
    return SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
                     IUPAC.protein),
                     id="YP_025292.1", name="HokC",
                     description="toxic membrane protein, small",
                     annotations={'hello':'world'})

@pytest.fixture(scope='class')
def seq_record_loaded_from_file_example(sequence_path):
    """Original SeqRecord loaded from sequence file"""
    return SeqIO.read(sequence_path, "fasta")

@pytest.fixture(scope='class')
def features_loaded_from_file_example(feature_path):
    """Original list of features"""
    with open(feature_path) as handle:
        feats = list(GFF.parse(handle))
        return feats[0].features

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
    """Path to the FASTA file"""
    return op.join(test_files_sequences, sequence_file)

@pytest.fixture(scope='module')
def metadata_path(test_files_sequences, metadata_file):
    """Path to the metadata file"""
    return op.join(test_files_sequences, metadata_file)

@pytest.fixture(scope='module')
def feature_path(test_files_sequences, feature_file):
    """Path to the GFF file"""
    return op.join(test_files_sequences, feature_file)

@pytest.fixture(scope='class')
def seqprop_with_i(sequence_id):
    """SeqProp with ID"""
    return SeqProp(id=sequence_id, seq=None)

@pytest.fixture(scope='class')
def seqprop_with_i_s_m_f(sequence_id, sequence_path, metadata_path, feature_path):
    """SeqProp with ID + sequence file + metadata file + feature file"""
    return SeqProp(id=sequence_id,
                   seq=None,
                   sequence_path=sequence_path,
                   metadata_path=metadata_path,
                   feature_path=feature_path)


class TestSeqPropWithId():
    """Class to test a bare SeqProp object with just an ID"""

    def test_init(self, seqprop_with_i, sequence_id):
        """Test initializing with just an ID"""
        assert seqprop_with_i.id == sequence_id

        # If just an ID initialized, everything should be empty
        assert seqprop_with_i.seq == None
        assert seqprop_with_i.name == '<unknown name>'
        assert seqprop_with_i.description == '<unknown description>'
        assert len(seqprop_with_i.annotations) == 0
        assert len(seqprop_with_i.letter_annotations) == 0
        assert len(seqprop_with_i.features) == 0

        # Files should not exist and raise errors if accessed
        assert seqprop_with_i.sequence_file == None
        with pytest.raises(IOError):
            seqprop_with_i.sequence_dir
        with pytest.raises(IOError):
            seqprop_with_i.sequence_path
        assert seqprop_with_i.metadata_file == None
        with pytest.raises(IOError):
            seqprop_with_i.metadata_dir
        with pytest.raises(IOError):
            seqprop_with_i.metadata_path
        assert seqprop_with_i.feature_file == None
        with pytest.raises(IOError):
            seqprop_with_i.feature_dir
        with pytest.raises(IOError):
            seqprop_with_i.feature_path

    def test_set_seq_with_str(self, seqprop_with_i, seq_str_example):
        """Test setting the seq attribute with a sequence string"""
        seqprop_with_i.seq = seq_str_example
        assert type(seqprop_with_i.seq) == Seq
        assert str(seqprop_with_i.seq) == seq_str_example

    def test_set_seq_with_seqrecord(self, seqprop_with_i, seq_record_example):
        """Test setting the seq attribute with a SeqRecord"""
        seqprop_with_i.seq = seq_record_example
        assert type(seqprop_with_i.seq) == Seq
        assert seqprop_with_i.seq == seq_record_example.seq
        assert seqprop_with_i.name == seq_record_example.name
        assert seqprop_with_i.description == seq_record_example.description
        assert seqprop_with_i.annotations == seq_record_example.annotations

    def test_get_emboss_pepstats_failure(self, seqprop_with_i):
        """Test that EMBOSS pepstats does not run when no file has been written"""
        with pytest.raises(IOError):
            seqprop_with_i.get_emboss_pepstats()

    def test_write_fasta_file(self, seqprop_with_i, tmpdir, test_files_outputs, seq_record_example):
        """Test that everything behaves properly when writing the SeqProp to a FASTA file"""
        # Add dummy annotations to the SeqProp - check to see if they stay in the SeqProp even after Seq is written
        seqprop_with_i.letter_annotations.update({'test_la_key': 'X' * len(seqprop_with_i.seq)})
        seqprop_with_i.features.append(SeqFeature(FeatureLocation(1, 3)))

        # Write the Seq to a FASTA file
        outpath = tmpdir.join('test_seqprop_with_i_write_fasta_file.fasta').strpath
        seqprop_with_i.write_fasta_file(outfile=outpath, force_rerun=True)

        # Test that the file was written
        assert op.exists(outpath)
        assert op.getsize(outpath) > 0

        # Test that file paths are correct
        assert seqprop_with_i.sequence_path == outpath
        assert seqprop_with_i.sequence_file == 'test_seqprop_with_i_write_fasta_file.fasta'
        assert seqprop_with_i.sequence_dir == tmpdir

        # Once a file is written, the annotations should not be lost, even though the sequence now
            # loads from the written file as a Seq
        assert seqprop_with_i.description == seq_record_example.description
        assert seqprop_with_i.annotations == seq_record_example.annotations
        assert seqprop_with_i.letter_annotations == {'test_la_key': 'X' * len(seq_record_example.seq)}
        assert len(seqprop_with_i.features) == 1

        # Test that sequence cannot be changed
        with pytest.raises(ValueError):
            seqprop_with_i.seq = 'THISWILLNOTBETHESEQ'
        assert seqprop_with_i.seq == seq_record_example.seq

    def test_get_residue_annotations(self, seqprop_with_i):
        """Test retrieval of residue letter_annotations"""
        stuff = seqprop_with_i.get_residue_annotations(start_resnum=1, end_resnum=34)
        assert stuff == {'test_la_key': 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'}

    def test_get_biopython_pepstats(self, seqprop_with_i):
        """Test storing Biopython pepstats and consistency of results"""
        seqprop_with_i.get_biopython_pepstats()

        results = {'instability_index': 27.172727272727272, 'aromaticity': 0.022727272727272728,
                   'percent_turn_naive': 0.022727272727272728,  'percent_strand_naive': 0.2954545454545454,
                   'monoisotopic': False, 'isoelectric_point': 8.84234619140625, 'molecular_weight': 4820.8507,
                   'percent_helix_naive': 0.38636363636363635}
        for k, v in results.items():
            assert seqprop_with_i.annotations[k] == pytest.approx(v)

    def test_get_emboss_pepstats_success(self, seqprop_with_i):
        """Test that EMBOSS pepstats does run when a file has been written"""
        seqprop_with_i.get_emboss_pepstats()
        assert 'percent_charged' in seqprop_with_i.annotations

    def test_set_metadata_path(self, seqprop_with_i, metadata_path, metadata_file, test_files_sequences):
        """Test setting the metadata file"""
        seqprop_with_i.metadata_path = metadata_path

        # Test that file paths are correct
        assert seqprop_with_i.metadata_path == metadata_path
        assert seqprop_with_i.metadata_file == metadata_file
        assert seqprop_with_i.metadata_dir == test_files_sequences

    def test_set_sequence_path(self, seqprop_with_i, seq_record_loaded_from_file_example,
                               sequence_path, sequence_file, test_files_sequences):
        """Test setting the seq attribute with a sequence file, and that seq is now loaded from file"""
        seqprop_with_i.sequence_path = sequence_path

        # Test that file paths are correct
        assert seqprop_with_i.sequence_path == sequence_path
        assert seqprop_with_i.sequence_file == sequence_file
        assert seqprop_with_i.sequence_dir == test_files_sequences

        # Test that the loaded sequence is the same as the original sequence
        assert seqprop_with_i.seq == seq_record_loaded_from_file_example.seq

        # Test that sequence cannot be changed
        with pytest.raises(ValueError):
            seqprop_with_i.seq = 'THISWILLNOTBETHESEQ'

    def test_set_features(self, seqprop_with_i, features_loaded_from_file_example):
        """Test setting the features attribute in memory"""
        seqprop_with_i.features = features_loaded_from_file_example[:5]
        assert seqprop_with_i.features == features_loaded_from_file_example[:5]
        assert seqprop_with_i.feature_file == None

    def test_write_gff_file(self, seqprop_with_i, tmpdir):
        """Test writing the features, and that features are now loaded from a file"""
        outpath = tmpdir.join('test_seqprop_with_i_write_gff_file.gff').strpath
        seqprop_with_i.write_gff_file(outfile=outpath, force_rerun=True)

        # Test that the file was written
        assert op.exists(outpath)
        assert op.getsize(outpath) > 0

        # Test that file paths are correct
        assert seqprop_with_i.feature_path == outpath
        assert seqprop_with_i.feature_file == 'test_seqprop_with_i_write_gff_file.gff'
        assert seqprop_with_i.feature_dir == tmpdir

        # Test that features cannot be changed
        with pytest.raises(ValueError):
            seqprop_with_i.features = ['NOFEATURES']

    def test_set_feature_path(self, seqprop_with_i, features_loaded_from_file_example,
                              feature_path, feature_file, test_files_sequences):
        """Test loading a feature file, and that old features are overwritten"""
        # Test that the existing feature set is not the same as the new one to be loaded
        assert len(seqprop_with_i.features) != len(features_loaded_from_file_example)

        seqprop_with_i.feature_path = feature_path

        # Test that file paths are correct
        assert seqprop_with_i.feature_path == feature_path
        assert seqprop_with_i.feature_file == feature_file
        assert seqprop_with_i.feature_dir == test_files_sequences

        # Test that features cannot be changed
        with pytest.raises(ValueError):
            seqprop_with_i.features = ['NOFEATURES']

        # Test that number of features stored is same
        assert len(seqprop_with_i.features) == len(features_loaded_from_file_example)


class TestSeqPropWithIdAndFiles():
    """Class to test a SeqProp object assigned files"""

    def test_init(self, seqprop_with_i_s_m_f, sequence_id,
                  sequence_path, metadata_path, feature_path, test_files_sequences,
                  sequence_file, metadata_file, feature_file,
                  seq_record_loaded_from_file_example, features_loaded_from_file_example):
        """Test initializing with assigned files"""
        assert seqprop_with_i_s_m_f.id == sequence_id
        assert seqprop_with_i_s_m_f.seq == seq_record_loaded_from_file_example.seq
        assert seqprop_with_i_s_m_f.name == seq_record_loaded_from_file_example.name
        assert seqprop_with_i_s_m_f.description == seq_record_loaded_from_file_example.description
        assert seqprop_with_i_s_m_f.annotations == seq_record_loaded_from_file_example.annotations
        assert seqprop_with_i_s_m_f.letter_annotations == seq_record_loaded_from_file_example.letter_annotations
        assert len(seqprop_with_i_s_m_f.features) == len(features_loaded_from_file_example)

        # Files should not exist and raise errors if accessed
        assert seqprop_with_i_s_m_f.sequence_file == sequence_file
        assert seqprop_with_i_s_m_f.sequence_dir == test_files_sequences
        assert seqprop_with_i_s_m_f.sequence_path == sequence_path
        assert seqprop_with_i_s_m_f.metadata_file == metadata_file
        assert seqprop_with_i_s_m_f.metadata_dir == test_files_sequences
        assert seqprop_with_i_s_m_f.metadata_path == metadata_path
        assert seqprop_with_i_s_m_f.feature_file == feature_file
        assert seqprop_with_i_s_m_f.feature_dir == test_files_sequences
        assert seqprop_with_i_s_m_f.feature_path == feature_path

    def test_change_root_dir(self, seqprop_with_i_s_m_f, tmpdir):
        # Test changing the sequence, metadata and feature dirs
        t = tmpdir.strpath
        seqprop_with_i_s_m_f.sequence_dir = t
        seqprop_with_i_s_m_f.metadata_dir = t
        seqprop_with_i_s_m_f.feature_dir = t

        assert seqprop_with_i_s_m_f.sequence_dir == t
        with pytest.raises(OSError):  # Sequence path should throw an error since the file was never moved there
            seqprop_with_i_s_m_f.sequence_path

        assert seqprop_with_i_s_m_f.metadata_dir == t
        with pytest.raises(OSError):  # Metadata path should throw an error since the file was never moved there
            seqprop_with_i_s_m_f.metadata_path

        assert seqprop_with_i_s_m_f.feature_dir == t
        with pytest.raises(OSError):  # Feature path should throw an error since the file was never moved there
            seqprop_with_i_s_m_f.feature_path