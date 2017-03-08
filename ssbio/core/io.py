from json_tricks.np import dump, load
import logging
log = logging.getLogger(__name__)


def save_json(obj, outfile, allow_nan=True, compression=True):
    """Save an ssbio object as a JSON file using json_tricks"""
    if compression:
        with open(outfile, 'wb') as f:
            dump(obj, f, allow_nan=allow_nan, compression=compression)
    else:
        with open(outfile, 'w') as f:
            dump(obj, f, allow_nan=allow_nan, compression=compression)
    log.info('Saved {} (id: {}) to {}'.format(type(obj), obj.id, outfile))


def load_json(file, new_root_dir=None, decompression=True):
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