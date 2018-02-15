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
def test_files_outputs(test_files):
    return op.join(test_files, 'out')

@pytest.fixture(scope='module')
def test_files_tempdir(test_files):
    return tempfile.gettempdir()

@pytest.fixture(scope='module')
def pdb_ids_working():
    return ['3bwm','1b4x','2XRN','4atp']

@pytest.fixture(scope='module')
def pdb_ids_obsolete():
    return ['4xdb','5aow','5vqc','5len']

@pytest.fixture(scope='module')
def pdb_ids_false():
    return ['soda','meow','1984','pycharm']