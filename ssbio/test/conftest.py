import os.path as op
import pytest

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