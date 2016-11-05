import datetime
import os.path as op
import glob
import os
import pandas as pd
from contextlib import contextmanager
import sys, os

import os
import sys
import time
import collections
from collections import defaultdict
from contextlib import contextmanager
from collections import OrderedDict, Callable
import distutils.spawn

try:
    from IPython.display import clear_output

    have_ipython = True
except ImportError:
    have_ipython = False


@contextmanager
def suppress_stdout():
    """Suppress stout messages.
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def split_folder_and_path(filepath):
    """Split a file path into its folder, filename, and extension

    Args:
        path (str): path to a file

    Returns:
        Tuple of (folder, filename (without extension), extension)

    """
    dirname = op.dirname(filepath)
    filename = op.basename(filepath)
    splitext = op.splitext(filename)
    filename_without_extension = splitext[0]
    extension = splitext[1]

    return dirname, filename_without_extension, extension


def outfile_name_maker(infile, outext='', outfile='', outdir=''):
    """Create a default name for an output file based on the infile name, unless a output name is specified.

    Args:
        infile: Path to input file.
        outfile: Optional specified name of output file.
        outdir: Optional path to output directory

    Returns:
        str: Path to final output destination

    Examples:

        >>> outfile_name_maker(infile='P00001.fasta')
        'P00001.out'

        >>> outfile_name_maker(infile='P00001.fasta', outext='.mao')
        'P00001.mao'

        >>> outfile_name_maker(infile='P00001.fasta', outext='.new', outfile='P00001_aligned')
        'P00001_aligned.new'

        >>> outfile_name_maker(infile='P00001.fasta', outfile='P00001_aligned')
        'P00001_aligned.out'

        >>> outfile_name_maker(infile='P00001.fasta', outfile='P00001_aligned', outdir='/my/dir/')
        '/my/dir/P00001_aligned'

    """
    # If extension not provided, default is "out"
    if not outext:
        outext = '.out'

    # If output filename not provided, default is to take name of infile
    if not outfile:
        orig_dir, outfile, orig_ext = split_folder_and_path(infile)

    # Join the output filename and output extension
    outfile = '{}{}'.format(outfile, outext)
    # Join the output directory and output filename
    outfile = op.join(outdir, outfile)

    return outfile


def force_rerun(flag, outfile):
    """Check if we should force rerunning of a command if an output file exists.

    Args:
        flag (bool): Flag to force rerun
        outfile (str): Path to output file which may already exist

    Returns:
        bool: If we should force rerunning of a command

    Examples:
        >>> force_rerun(flag=False, outfile='/not/existing/file.txt')
        True

        >>> force_rerun(flag=True, outfile='/existing/file.txt')
        True

        >>> force_rerun(flag=False, outfile='/existing/file.txt')
        False

    """
    # If flag is True, always run
    if flag:
        return True
    # If flag is False but file doesn't exist, also run
    elif not flag and not op.exists(outfile):
        return True
    # Otherwise, do not run
    else:
        return False


def program_exists(prog_name):
    """Check if a program is available as a command line executable on a system.

    Args:
        prog_name: Name of the program.

    Returns:
        bool: True if the program is available.

    """
    if distutils.spawn.find_executable(prog_name):
        return True
    else:
        return False


def dict_head(d, disp=5):
    """Return the head of a dictionary.

    Returns the first 5 key/value pairs in a dictionary

    Args:
        d: dictionary to get head
        disp: number of elements to display

    Returns:
        Dict

    """
    return {k: d[k] for k in list(d.keys())[:disp]}


def rank_dated_files(pattern, dir):
    """Search a directory for files that match a pattern. Return an ordered list of these files.

    Args:
        pattern: glob pattern to search for

    Returns:
        Rank-ordered list, usually by date (shortdate format, e.g. 161010).
        The most recent file will be in position 0.

    """
    files = glob.glob(op.join(dir, pattern))
    return sorted(files, reverse=True)


def input_parser(args):
    """Parse command line inputs
    """
    pass


def find(lst, a, case_sensitive=True):
    """Return indices of a list which have elements that match an object or list of objects

    Args:
        lst: list of values
        a: object(s) to check equality

    Returns:
        list: list of indicies of lst which equal a

    """
    a = force_list(a)

    if not case_sensitive:
        lst = [x.lower() for x in lst]
        a = [y.lower() for y in a]

    return [i for i, x in enumerate(lst) if x in a]


def not_find(lst, a, case_sensitive=True):
    """Return indices of a list which have elements that DO NOT match an object or list of objects

    Args:
        lst: list of values
        a: object(s) to check inequality

    Returns:
        list: list of indicies of lst which do not equal a

    """
    a = force_list(a)

    if not case_sensitive:
        lst = [x.lower() for x in lst]
        a = [y.lower() for y in a]

    return [i for i, x in enumerate(lst) if x not in a]


def filter_list(lst, takeout, case_sensitive=True):
    """Return a modified list removing items specified.

    Args:
        lst: Original list of values
        takeout: Object or objects to remove from lst

    Returns:
        list: Filtered list of values

    """
    takeout = force_list(takeout)

    if not case_sensitive:
        lst = [x.lower() for x in lst]
        takeout = [y.lower() for y in takeout]

    return [x for x in lst if x not in takeout]


def filter_list_by_indices(lst, indices):
    """Return a modified list containing only the indices indicated.

    Args:
        lst: Original list of values
        indices: List of indices to keep from the original list

    Returns:
        list: Filtered list of values

    """
    return [x for i, x in enumerate(lst) if i in indices]


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def force_string(val=None):
    """Force a string representation of an object

    Args:
        val (str): object to parse into a string

    Returns:
        String

    """
    if val is None:
        return ''
    return val if isinstance(val, str) else ';'.join(val)


def force_list(val=None):
    """Force a list representation of an object

    Args:
        val: object to parse into a list

    Returns:

    """
    if val is None:
        return []
    if isinstance(val, pd.Series):
        return val.tolist()
    return val if isinstance(val, list) else [val]


def force_lower_list(val=None):
    """Force a lowercase list representation of strings

    Args:
        val: string or strings to parse into a list

    Returns:
        list: with lowercase values

    """
    return [x.lower() for x in force_list(val)]


def force_upper_list(val=None):
    """Force a UPPERCASE list representation of strings

    Args:
        val: string or strings to parse into a list

    Returns:
        list: with UPPERCASE values

    """
    return [x.upper() for x in force_list(val)]


def split_list(a_list):
    half = len(a_list) // 2
    return a_list[:half], a_list[half:]


def input_list_parser(instring, filetype=''):
    """Always return a list of files with varying input

    1. /path/to/folder -> list of files in folder (full paths)
    2. /path/to/file -> list of files (singular list)
    3. file1,file2 -> list of files
    4.

    Args:
        instring:

    Returns:

    """
    if filetype:
        # TODO
        searchstring = '*'
    else:
        searchstring = '*'

    if op.isdir(instring):
        os.chdir(instring)
        return glob.glob('*')

    if op.isfile(instring):
        return force_list(instring)

    else:
        return instring.split(',')


def flatlist_dropdup(list_of_lists):
    return list(set([str(item) for sublist in list_of_lists for item in sublist]))


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def combinations(iterable, r):
    """Calculate combinations
    combinations('ABCD', 2) --> AB AC AD BC BD CD
    combinations(range(4), 3) --> 012 013 023 123

    Args:
        iterable:
        r:

    Returns:

    """

    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield list(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield list(pool[i] for i in indices)


class Date():
    def __init__(self):
        self.short_date = self.date_prefix()

    def date_prefix(self):
        today = datetime.date.today()
        return today.strftime('%y%m%d')


class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
                not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))


def percentage_to_float(x):
    return float(x.strip('%')) / 100
