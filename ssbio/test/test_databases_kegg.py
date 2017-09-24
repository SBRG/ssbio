import pytest
import os.path as op
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from BCBio import GFF
from ssbio.databases.kegg import KEGGProp

@pytest.fixture(scope='class')
def seq_record_loaded_from_file_example(fasta_path):
    """Original SeqRecord loaded from sequence file"""
    return SeqIO.read(fasta_path, "fasta")

@pytest.fixture(scope='module')
def kegg_id():
    return 'mtu:Rv0417'

@pytest.fixture(scope='module')
def fasta_file():
    return 'mtu-Rv0417.faa'

@pytest.fixture(scope='module')
def txt_file():
    return 'mtu-Rv0417.kegg'

@pytest.fixture(scope='module')
def fasta_path(test_files_sequences, fasta_file):
    return op.join(test_files_sequences, fasta_file)

@pytest.fixture(scope='module')
def txt_path(test_files_sequences, txt_file):
    return op.join(test_files_sequences, txt_file)

@pytest.fixture(scope='class')
def keggprop_with_i(kegg_id):
    return KEGGProp(id=kegg_id,
                    seq=None)

@pytest.fixture(scope='class')
def keggprop_with_i_s_m_f(kegg_id, fasta_path, txt_path):
    return KEGGProp(id=kegg_id,
                    seq=None,
                    fasta_path=fasta_path,
                    txt_path=txt_path)


class TestKEGGPropWithId():

    """Class to test a bare KEGGProp object with just an ID"""

    def test_init(self, keggprop_with_i, kegg_id):
        """Test initializing with just an ID"""
        assert keggprop_with_i.id == kegg_id

        # If just an ID initialized, everything should be empty
        assert keggprop_with_i.seq == None
        assert keggprop_with_i.name == '<unknown name>'
        assert keggprop_with_i.description == '<unknown description>'
        assert len(keggprop_with_i.annotations) == 0
        assert len(keggprop_with_i.letter_annotations) == 0
        assert len(keggprop_with_i.features) == 0

        # Files should not exist and raise errors if accessed
        assert keggprop_with_i.sequence_file == None
        with pytest.raises(IOError):
            keggprop_with_i.sequence_dir
        with pytest.raises(IOError):
            keggprop_with_i.sequence_path
        assert keggprop_with_i.metadata_file == None
        with pytest.raises(IOError):
            keggprop_with_i.metadata_dir
        with pytest.raises(IOError):
            keggprop_with_i.metadata_path
        assert keggprop_with_i.feature_file == None
        with pytest.raises(IOError):
            keggprop_with_i.feature_dir
        with pytest.raises(IOError):
            keggprop_with_i.feature_path

    def test_set_sequence_path(self, keggprop_with_i, fasta_path, fasta_file, test_files_sequences):
        """Test setting the seq attribute with a sequence file"""
        keggprop_with_i.sequence_path = fasta_path

        # Test that file paths are correct
        assert keggprop_with_i.sequence_path == fasta_path
        assert keggprop_with_i.sequence_file == fasta_file
        assert keggprop_with_i.sequence_dir == test_files_sequences

    def test_set_feature_path(self, keggprop_with_i, features_loaded_from_file_example,
                              gff_path, gff_file, test_files_sequences):
        """Test loading a feature file, and that old features are overwritten"""
        # Test that the existing feature set is not the same as the new one to be loaded
        assert len(keggprop_with_i.features) != len(features_loaded_from_file_example)

        keggprop_with_i.feature_path = gff_path

        # Test that file paths are correct
        assert keggprop_with_i.feature_path == gff_path
        assert keggprop_with_i.feature_file == gff_file
        assert keggprop_with_i.feature_dir == test_files_sequences

        # Test that features cannot be changed
        with pytest.raises(ValueError):
            keggprop_with_i.features = ['NOFEATURES']

        # Test that number of features stored is same
        assert len(keggprop_with_i.features) == len(features_loaded_from_file_example)

    def test_set_metadata_path(self, keggprop_with_i, txt_path, txt_file, test_files_sequences,
                               txt_record_loaded_from_file_example):
        keggprop_with_i.metadata_path = txt_path

        # Unset sequence and feature paths
        keggprop_with_i.sequence_path = None
        keggprop_with_i.feature_path = None

        # Test that file paths are correct
        assert keggprop_with_i.metadata_path == txt_path
        assert keggprop_with_i.metadata_file == txt_file
        assert keggprop_with_i.metadata_dir == test_files_sequences

        # Test loaded information
        assert keggprop_with_i.description == txt_record_loaded_from_file_example.description
        assert keggprop_with_i.bigg == None
        for k in ['ecj:JW4347', 'eco:b4384']:
            assert k in keggprop_with_i.kegg
        for r in ['NP_418801.1', 'WP_000224877.1']:
            assert r in keggprop_with_i.refseq
        assert keggprop_with_i.kegg == 'mtu:Rv0417'
        assert keggprop_with_i.gene_name == 'deoD'
        for p in ['1A69', '1ECP', '1K9S', '1OTX', '1OTY', '1OU4', '1OUM', '1OV6', '1OVG',
                    '3ONV', '3OOE', '3OOH', '3OPV', '3UT6', '4TS3', '4TS9', '4TTA', '4TTI',
                    '4TTJ', '5I3C', '5IU6']:
            assert p in keggprop_with_i.pdbs
        for g in ['GO:0004731', 'GO:0005829', 'GO:0006152', 'GO:0006974', 'GO:0016020', 'GO:0019686', 'GO:0042802']:
            assert g in keggprop_with_i.go
        assert keggprop_with_i.pfam == ['PF01048']
        assert keggprop_with_i.ec_number == None  ## TODO: parse
        assert keggprop_with_i.reviewed == False  ## TODO: parse
        for u in ['Q2M5T3', 'P09743']:
            assert u in keggprop_with_i.alt_keggs
        assert keggprop_with_i.taxonomy == 'Escherichia coli (strain K12)'
        assert keggprop_with_i.seq_version == 2
        assert keggprop_with_i.seq_date == '2007-01-23'
        assert keggprop_with_i.entry_version == 106
        assert keggprop_with_i.entry_date == '2017-08-30'

        # Test that features are loaded directly from this metadata file
        assert len(keggprop_with_i.features) == len(txt_record_loaded_from_file_example.features)


class TestKEGGPropWithIdAndFiles():

    """Class to test a bare KEGGProp object with just an ID"""

    def test_init(self, keggprop_with_i_s_m_f, kegg_id,
                  fasta_path, txt_path, gff_path, test_files_sequences,
                  fasta_file, txt_file, gff_file,
                  seq_record_loaded_from_file_example,
                  features_loaded_from_file_example,
                  txt_record_loaded_from_file_example):
        """Test initializing with assigned files"""
        assert keggprop_with_i_s_m_f.id == kegg_id
        assert keggprop_with_i_s_m_f.seq == seq_record_loaded_from_file_example.seq
        assert keggprop_with_i_s_m_f.name == seq_record_loaded_from_file_example.name
        assert keggprop_with_i_s_m_f.description == txt_record_loaded_from_file_example.description
        assert keggprop_with_i_s_m_f.annotations == {}  # No annotations will be loaded from files
        assert keggprop_with_i_s_m_f.letter_annotations == txt_record_loaded_from_file_example.letter_annotations
        assert len(keggprop_with_i_s_m_f.features) == len(features_loaded_from_file_example)

        # Files should exist
        assert keggprop_with_i_s_m_f.sequence_file == fasta_file
        assert keggprop_with_i_s_m_f.sequence_dir == test_files_sequences
        assert keggprop_with_i_s_m_f.sequence_path == fasta_path
        assert keggprop_with_i_s_m_f.metadata_file == txt_file
        assert keggprop_with_i_s_m_f.metadata_dir == test_files_sequences
        assert keggprop_with_i_s_m_f.metadata_path == txt_path
        assert keggprop_with_i_s_m_f.feature_file == gff_file
        assert keggprop_with_i_s_m_f.feature_dir == test_files_sequences
        assert keggprop_with_i_s_m_f.feature_path == gff_path