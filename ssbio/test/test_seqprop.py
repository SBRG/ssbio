import pytest
import os.path as op
from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

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
def seqprop_with_id_seq(sequence_id, seq_record_example):
    return SeqProp(ident=sequence_id,
                   seq=seq_record_example)

@pytest.fixture(scope='class')
def seqprop_with_id_seqpath(sequence_path):
    return SeqProp(ident='TESTID',
                   sequence_path=sequence_path)

@pytest.fixture(scope='class')
def seqprop_with_id_seqpath_metapath(sequence_id, sequence_path, metadata_path):
    return SeqProp(ident=sequence_id,
                   sequence_path=sequence_path,
                   metadata_path=metadata_path)

@pytest.fixture(scope='class')
def seqprop_with_id_seqpath_metapath_featpath(sequence_id, sequence_path, metadata_path, feature_path):
    return SeqProp(ident=sequence_id,
                   sequence_path=sequence_path,
                   metadata_path=metadata_path,
                   feature_path=feature_path)

@pytest.fixture(scope='module')
def seq_record_example():
    return SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
                     IUPAC.protein),
                     id="YP_025292.1", name="HokC",
                     description="toxic membrane protein, small",
                     annotations={'hello':'world'})

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
        assert seqprop_with_id.seq_len == len(orig_sr.seq)
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

    def test_change_root_dir(self, seqprop_with_id, tmpdir):
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
        with pytest.raises(OSError):  # Feature path should throw an error since the file was never moved there
            seqprop_with_id.feature_path


class TestSeqPropWithSeq():
    def test_init(self, seqprop_with_id_seq, seq_record_example):
        assert seqprop_with_id_seq.seq_record.seq == seq_record_example.seq

        # Identifiers will be different here..
        assert seqprop_with_id_seq.id != seqprop_with_id_seq.seq_record.id

        # If a sequence is initialized, the F/A/LA/D should be copied from the sequence
        assert seqprop_with_id_seq.description == seqprop_with_id_seq.seq_record.description
        assert seqprop_with_id_seq.annotations == seqprop_with_id_seq.seq_record.annotations
        assert seqprop_with_id_seq.letter_annotations == seqprop_with_id_seq.seq_record.letter_annotations
        assert seqprop_with_id_seq.features == seqprop_with_id_seq.seq_record.features

        assert seqprop_with_id_seq.sequence_file == None
        with pytest.raises(IOError):
            seqprop_with_id_seq.sequence_dir
        with pytest.raises(IOError):
            seqprop_with_id_seq.sequence_path

        assert seqprop_with_id_seq.metadata_file == None
        with pytest.raises(IOError):
            seqprop_with_id_seq.metadata_dir
        with pytest.raises(IOError):
            seqprop_with_id_seq.metadata_path

    def test_add_annotations(self, seqprop_with_id_seq):
        # Add dummy annotations to the SeqProp
        seqprop_with_id_seq.annotations.update({'test_key': 'test_annotation'})
        seqprop_with_id_seq.letter_annotations.update({'test_la1': 'X'*44})

        # Test that annotations and letter_annotations are copied over to the SeqRecord
        assert seqprop_with_id_seq.seq_record.annotations == seqprop_with_id_seq.annotations
        assert seqprop_with_id_seq.seq_record.letter_annotations == seqprop_with_id_seq.letter_annotations

        # Test that SeqProp letter_annotations is a RestrictedDict
        with pytest.raises(TypeError):
            seqprop_with_id_seq.letter_annotations.update({'test_la2': 'MTESTER'})

    def test_write_fasta_file(self, seqprop_with_id_seq, tmpdir, test_files_outputs):
        outpath = tmpdir.join('test_write_fasta_file.fasta').strpath
        seqprop_with_id_seq.write_fasta_file(outfile=outpath, force_rerun=True)

        # Test that the file was written
        assert op.exists(outpath)

        # Test that the file contents are correct
        with open(op.join(test_files_outputs, 'test_write_fasta_file.fasta')) as orig_file:
            with open(outpath) as f:
                assert f.read() == orig_file.read()

        # Once a file is written, the annotations should not be lost, even though the sequence now
            # loads from the written file as a SeqRecord
        assert seqprop_with_id_seq.description == seqprop_with_id_seq.seq_record.description
        assert seqprop_with_id_seq.annotations == seqprop_with_id_seq.seq_record.annotations
        assert seqprop_with_id_seq.letter_annotations == seqprop_with_id_seq.seq_record.letter_annotations
        assert seqprop_with_id_seq.features == seqprop_with_id_seq.seq_record.features

    def test_add_annotations2(self, seqprop_with_id_seq):
        # Another test_add_annotations, now to see if things change when the sequence is written
        # Add dummy annotations to the SeqProp
        seqprop_with_id_seq.annotations.update({'test_key': 'test_annotation2'})
        seqprop_with_id_seq.letter_annotations.update({'test_la1': '2'*44})

        # Test that annotations and letter_annotations are copied over to the SeqRecord
        assert seqprop_with_id_seq.seq_record.annotations == seqprop_with_id_seq.annotations
        assert seqprop_with_id_seq.seq_record.letter_annotations == seqprop_with_id_seq.letter_annotations

        # Test that SeqProp letter_annotations is a RestrictedDict
        with pytest.raises(TypeError):
            seqprop_with_id_seq.letter_annotations.update({'test_la2': 'MTESTER'})

    def test_get_biopython_pepstats(self, seqprop_with_id_seq):
        seqprop_with_id_seq.get_biopython_pepstats()

        results = {'instability_index': 27.172727272727272, 'aromaticity': 0.022727272727272728,
                   'percent_turn_naive': 0.022727272727272728,  'percent_strand_naive': 0.2954545454545454,
                   'monoisotopic': False, 'isoelectric_point': 8.84234619140625, 'molecular_weight': 4820.8507,
                   'percent_helix_naive': 0.38636363636363635}
        for k, v in results.items():
            assert seqprop_with_id_seq.annotations[k] == v
            assert seqprop_with_id_seq.seq_record.annotations[k] == v


class TestSeqPropWithSeqFile():
    pass