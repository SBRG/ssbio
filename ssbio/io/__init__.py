import logging
import pickle

from json_tricks import dump, load

log = logging.getLogger(__name__)


def save_json(obj, outfile, allow_nan=True, compression=False):
    """Save an ssbio object as a JSON file using json_tricks"""
    if compression:
        with open(outfile, 'wb') as f:
            dump(obj, f, allow_nan=allow_nan, compression=compression)
    else:
        with open(outfile, 'w') as f:
            dump(obj, f, allow_nan=allow_nan, compression=compression)
    log.info('Saved {} (id: {}) to {}'.format(type(obj), obj.id, outfile))


def load_json(file, new_root_dir=None, decompression=False):
    """Load a JSON file using json_tricks"""
    if decompression:
        with open(file, 'rb') as f:
            my_object = load(f, decompression=decompression)
    else:
        with open(file, 'r') as f:
            my_object = load(f, decompression=decompression)
    if new_root_dir:
        my_object.root_dir = new_root_dir

    return my_object


def save_pickle(obj, outfile, protocol=2):
    """Save the object as a pickle file

    Args:
        outfile (str): Filename
        protocol (int): Pickle protocol to use. Default is 2 to remain compatible with Python 2

    Returns:
        str: Path to pickle file

    """
    with open(outfile, 'wb') as f:
        pickle.dump(obj, f, protocol=protocol)

    return outfile


def load_pickle(file, encoding=None):
    """Load a pickle file.

    Args:
        file (str): Path to pickle file

    Returns:
        object: Loaded object from pickle file

    """
    # TODO: test set encoding='latin1' for 2/3 incompatibility
    if encoding:
        with open(file, 'rb') as f:
            return pickle.load(f, encoding=encoding)

    with open(file, 'rb') as f:
       return pickle.load(f)