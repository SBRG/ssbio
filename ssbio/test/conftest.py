import os.path as op
import pytest
import tempfile


@pytest.fixture(scope='module')
def test_files_tempdir(test_files):
    # /tmp/
    return tempfile.gettempdir()


@pytest.fixture(scope='module')
def test_files():
    return 'test_files'


@pytest.fixture(scope='module')
def test_files_sequences(test_files):
    """Sequence files"""
    # ssbio/test/test_files/sequences
    return op.join(test_files, 'sequences')


@pytest.fixture(scope='module')
def test_files_structures(test_files):
    """Structure files"""
    # ssbio/test/test_files/structures
    return op.join(test_files, 'structures')


@pytest.fixture(scope='module')
def test_files_models(test_files):
    """SBML model files"""
    # ssbio/test/test_files/models
    return op.join(test_files, 'models')


@pytest.fixture(scope='module')
def test_files_outputs(test_files):
    """ssbio local output directory - likely to use tempdir in most cases"""
    # ssbio/test/test_files/out
    return op.join(test_files, 'out')


@pytest.fixture(scope='module')
def test_gem_small_json(test_files_models):
    """Mini JSON GEM for testing (E. coli)"""
    # ssbio/test/test_files/models/mini.json
    return op.join(test_files_models, 'mini.json')


@pytest.fixture(scope='module')
def test_gem_large_json(test_files_models):
    """Large JSON GEM for testing (MTB)"""
    # ssbio/test/test_files/models/iNJ661.json
    return op.join(test_files_models, 'iNJ661.json')


@pytest.fixture(scope='module')
def list_of_genes_ecoli():
    """List of genes for testing (E. coli)"""
    return ['b0720','b1241','b1761','b0003']


@pytest.fixture(scope='module')
def dict_of_genes_seqs_ecoli():
    """Dictionary of genes and sequence strings for testing (E. coli)"""
    return {'b0870': 'MIDLRSDTVTRPSRAMLEAMMAAPVGDDVYGDDPTVNALQDYAAELSGKEAAIFLPTGTQANLVALLSHCERGEEYIVGQAAHNYLFEAGGAAVLGSIQPQPIDAAADGTLPLDKVAMKIKPDDIHFARTKLLSLENTHNGKVLPREYLKEAWEFTRERNLALHVDGARIFNAVVAYGCELKEITQYCDSFTICLSKGLGTPVGSLLVGNRDYIKRAIRWRKMTGGGMRQSGILAAAGIYALKNNVARLQEDHDNAAWMAEQLREAGADVMRQDTNMLFVRVGEENAAALGEYMKARNVLINASPIVRLVTHLDVSREQLAEVAAHWRAFLAR',
            'b3041': 'MNQTLLSSFGTPFERVENALAALREGRGVMVLDDEDRENEGDMIFPAETMTVEQMALTIRHGSGIVCLCITEDRRKQLDLPMMVENNTSAYGTGFTVTIEAAEGVTTGVSAADRITTVRAAIADGAKPSDLNRPGHVFPLRAQAGGVLTRGGHTEATIDLMTLAGFKPAGVLCELTNDDGTMARAPECIEFANKHNMALVTIEDLVAYRQAHERKAS'}


@pytest.fixture(scope='module')
def test_fasta_file_multseqs(test_files_sequences):
    # XTODO: create test file
    return op.join(test_files_sequences, 'multiple_sequences.fasta')


@pytest.fixture(scope='module')
def pdb_ids_working():
    return ['3bwm','1b4x','2XRN','4atp']


@pytest.fixture(scope='module')
def pdb_ids_obsolete():
    return ['4xdb','5aow','5vqc','5len']


@pytest.fixture(scope='module')
def pdb_ids_false():
    return ['soda','meow','1984','pycharm']