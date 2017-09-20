import pytest
import os.path as op
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from BCBio import GFF
from ssbio.databases.uniprot import UniProtProp

@pytest.fixture(scope='class')
def seq_record_loaded_from_file_example(fasta_path):
    """Original SeqRecord loaded from sequence file"""
    return SeqIO.read(fasta_path, "fasta")

@pytest.fixture(scope='class')
def xml_record_loaded_from_file_example(xml_path):
    """Original SeqRecord loaded from sequence file"""
    return SeqIO.read(xml_path, "uniprot-xml")

@pytest.fixture(scope='class')
def features_loaded_from_file_example(gff_path):
    """Original list of features"""
    with open(gff_path) as handle:
        feats = list(GFF.parse(handle))
        return feats[0].features

@pytest.fixture(scope='module')
def uniprot_id():
    return 'P0ABP8'

@pytest.fixture(scope='module')
def fasta_file():
    return 'P0ABP8.fasta'

@pytest.fixture(scope='module')
def xml_file():
    return 'P0ABP8.xml'

@pytest.fixture(scope='module')
def gff_file():
    return 'P0ABP8.gff'

@pytest.fixture(scope='module')
def fasta_path(test_files_sequences, fasta_file):
    return op.join(test_files_sequences, fasta_file)

@pytest.fixture(scope='module')
def xml_path(test_files_sequences, xml_file):
    return op.join(test_files_sequences, xml_file)

@pytest.fixture(scope='module')
def gff_path(test_files_sequences, gff_file):
    return op.join(test_files_sequences, gff_file)

@pytest.fixture(scope='class')
def uniprotprop_with_i(uniprot_id):
    return UniProtProp(id=uniprot_id,
                       seq=None)

@pytest.fixture(scope='class')
def uniprotprop_with_i_s_m_f(uniprot_id, fasta_path, xml_path, gff_path):
    return UniProtProp(id=uniprot_id,
                       seq=None,
                       fasta_path=fasta_path,
                       xml_path=xml_path,
                       gff_path=gff_path)


class TestUniProtPropWithID():
    """Class to test a bare UniProtProp object with just an ID"""

    def test_init(self, uniprotprop_with_i, uniprot_id):
        """Test initializing with just an ID"""
        assert uniprotprop_with_i.id == uniprot_id

        # If just an ID initialized, everything should be empty
        assert uniprotprop_with_i.seq == None
        assert uniprotprop_with_i.name == '<unknown name>'
        assert uniprotprop_with_i.description == '<unknown description>'
        assert len(uniprotprop_with_i.annotations) == 0
        assert len(uniprotprop_with_i.letter_annotations) == 0
        assert len(uniprotprop_with_i.features) == 0

        # Files should not exist and raise errors if accessed
        assert uniprotprop_with_i.sequence_file == None
        with pytest.raises(IOError):
            uniprotprop_with_i.sequence_dir
        with pytest.raises(IOError):
            uniprotprop_with_i.sequence_path
        assert uniprotprop_with_i.metadata_file == None
        with pytest.raises(IOError):
            uniprotprop_with_i.metadata_dir
        with pytest.raises(IOError):
            uniprotprop_with_i.metadata_path
        assert uniprotprop_with_i.feature_file == None
        with pytest.raises(IOError):
            uniprotprop_with_i.feature_dir
        with pytest.raises(IOError):
            uniprotprop_with_i.feature_path

    def test_set_sequence_path(self, uniprotprop_with_i, fasta_path, fasta_file, test_files_sequences):
        """Test setting the seq attribute with a sequence file"""
        uniprotprop_with_i.sequence_path = fasta_path

        # Test that file paths are correct
        assert uniprotprop_with_i.sequence_path == fasta_path
        assert uniprotprop_with_i.sequence_file == fasta_file
        assert uniprotprop_with_i.sequence_dir == test_files_sequences

    def test_set_feature_path(self, uniprotprop_with_i, features_loaded_from_file_example,
                              gff_path, gff_file, test_files_sequences):
        """Test loading a feature file, and that old features are overwritten"""
        # Test that the existing feature set is not the same as the new one to be loaded
        assert len(uniprotprop_with_i.features) != len(features_loaded_from_file_example)

        uniprotprop_with_i.feature_path = gff_path

        # Test that file paths are correct
        assert uniprotprop_with_i.feature_path == gff_path
        assert uniprotprop_with_i.feature_file == gff_file
        assert uniprotprop_with_i.feature_dir == test_files_sequences

        # Test that features cannot be changed
        with pytest.raises(ValueError):
            uniprotprop_with_i.features = ['NOFEATURES']

        # Test that number of features stored is same
        assert len(uniprotprop_with_i.features) == len(features_loaded_from_file_example)

    def test_set_metadata_path(self, uniprotprop_with_i, xml_path, xml_file, test_files_sequences,
                               xml_record_loaded_from_file_example):
        uniprotprop_with_i.metadata_path = xml_path

        # Unset sequence and feature paths
        uniprotprop_with_i.sequence_path = None
        uniprotprop_with_i.feature_path = None

        # Test that file paths are correct
        assert uniprotprop_with_i.metadata_path == xml_path
        assert uniprotprop_with_i.metadata_file == xml_file
        assert uniprotprop_with_i.metadata_dir == test_files_sequences

        # Test loaded information
        assert uniprotprop_with_i.description == xml_record_loaded_from_file_example.description
        assert uniprotprop_with_i.bigg == None
        for k in ['ecj:JW4347', 'eco:b4384']:
            assert k in uniprotprop_with_i.kegg
        for r in ['NP_418801.1', 'WP_000224877.1']:
            assert r in uniprotprop_with_i.refseq
        assert uniprotprop_with_i.uniprot == 'P0ABP8'
        assert uniprotprop_with_i.gene_name == 'deoD'
        for p in ['1A69', '1ECP', '1K9S', '1OTX', '1OTY', '1OU4', '1OUM', '1OV6', '1OVG',
                    '3ONV', '3OOE', '3OOH', '3OPV', '3UT6', '4TS3', '4TS9', '4TTA', '4TTI',
                    '4TTJ', '5I3C', '5IU6']:
            assert p in uniprotprop_with_i.pdbs
        for g in ['GO:0004731', 'GO:0005829', 'GO:0006152', 'GO:0006974', 'GO:0016020', 'GO:0019686', 'GO:0042802']:
            assert g in uniprotprop_with_i.go
        assert uniprotprop_with_i.pfam == ['PF01048']
        assert uniprotprop_with_i.ec_number == None  ## TODO: parse
        assert uniprotprop_with_i.reviewed == False  ## TODO: parse
        for u in ['Q2M5T3', 'P09743']:
            assert u in uniprotprop_with_i.alt_uniprots
        assert uniprotprop_with_i.taxonomy == 'Escherichia coli (strain K12)'
        assert uniprotprop_with_i.seq_version == 2
        assert uniprotprop_with_i.seq_date == '2007-01-23'
        assert uniprotprop_with_i.entry_version == 106
        assert uniprotprop_with_i.entry_date == '2017-08-30'

        # Test that features are loaded directly from this metadata file
        assert len(uniprotprop_with_i.features) == len(xml_record_loaded_from_file_example.features)