import os.path as op
import pytest
import tempfile


@pytest.fixture(scope='module')
def test_files():
    return 'test_files'


@pytest.fixture(scope='module')
def test_files_sequences(test_files):
    return op.join(test_files, 'sequences')


@pytest.fixture(scope='module')
def test_files_structures(test_files):
    return op.join(test_files, 'structures')


@pytest.fixture(scope='module')
def test_files_models(test_files):
    return op.join(test_files, 'models')


@pytest.fixture(scope='module')
def test_files_outputs(test_files):
    return op.join(test_files, 'out')


@pytest.fixture(scope='module')
def test_files_tempdir(test_files):
    return tempfile.gettempdir()


@pytest.fixture(scope='module')
def test_gem_small_sbml(test_files_models):
    return op.join(test_files_models, 'Ec_core_flux1.xml')


@pytest.fixture(scope='module')
def test_gem_large_json(test_files_models):
    return op.join(test_files_models, 'iNJ661.json')


@pytest.fixture(scope='module')
def list_of_genes_ecoli():
    return ['b0720','b1241','b1761','b0003']


@pytest.fixture(scope='module')
def dict_of_genes_seqs_ecoli():
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