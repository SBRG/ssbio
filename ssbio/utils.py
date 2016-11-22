import os
import sys
import datetime
import os.path as op
import glob
import pandas as pd
import contextlib
import logging
import distutils.spawn
import subprocess
import shutil
from collections import OrderedDict
from collections import Callable
try:
    from IPython.display import clear_output
    have_ipython = True
except ImportError:
    have_ipython = False

log = logging.getLogger(__name__)


class Date():
    """Various methods to return formatted dates for today.
    """
    def __init__(self):
        self.short_date = self.short_date_prefix()

    def short_date_prefix(self):
        today = datetime.date.today()
        return today.strftime('%y%m%d')


def todays_short_date():
    today = Date()
    return today.short_date


class DefaultOrderedDict(OrderedDict):
    """Class to combine defaultdict and OrderedDict.
    Source: http://stackoverflow.com/a/6190500/562769

    """
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


@contextlib.contextmanager
def suppress_stdout():
    """Suppress the stdout of any function.

    Usage:
    with ssbio.utils.suppress_stdout():
        my_function_here()
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
        path (str): Path to a file

    Returns:
        tuple: of (folder, filename (without extension), extension)

    """
    dirname = op.dirname(filepath)
    filename = op.basename(filepath)
    splitext = op.splitext(filename)
    filename_without_extension = splitext[0]
    extension = splitext[1]

    return dirname, filename_without_extension, extension


def outfile_name_maker(inname, outext='.out', outfile='', outdir='', append_to_name=''):
    """Create a default name for an output file based on the inname name, unless a output name is specified.

    Args:
        inname: Path to input file
        outext: Optional specified extension for output file (with the "."). Default is ".out".
        outfile: Optional specified name of output file.
        outdir: Optional path to output directory.

    Returns:
        str: Path to final output destination.

    Examples:

        >>> outfile_name_maker(inname='P00001.fasta')
        'P00001.out'

        >>> outfile_name_maker(inname='P00001.fasta', append_to_name='_new')
        'P00001_new.out'

        >>> outfile_name_maker(inname='P00001.fasta', outext='.mao')
        'P00001.mao'

        >>> outfile_name_maker(inname='P00001.fasta', outext='.mao', append_to_name='_new')
        'P00001_new.mao'

        >>> outfile_name_maker(inname='P00001.fasta', outext='.new', outfile='P00001_aligned')
        'P00001_aligned.new'

        >>> outfile_name_maker(inname='P00001.fasta', outfile='P00001_aligned')
        'P00001_aligned.out'

        >>> outfile_name_maker(inname='P00001.fasta', outfile='P00001_aligned', append_to_name='_new')
        'P00001_aligned_new.out'

        >>> outfile_name_maker(inname='P00001.fasta', outfile='P00001_aligned', outdir='/my/dir/')
        '/my/dir/P00001_aligned.out'

    """

    # If output filename not provided, default is to take name of inname
    if not outfile:
        orig_dir, outfile, orig_ext = split_folder_and_path(inname)

    # Append additional stuff to the filename if specified
    if append_to_name:
        outfile = outfile + append_to_name

    # Join the output filename and output extension
    final_outfile = '{}{}'.format(outfile, outext)

    # Join the output directory and output filename
    final_outfile = op.join(outdir, final_outfile)

    return final_outfile


def force_rerun(flag, outfile):
    """Check if we should force rerunning of a command if an output file exists.

    Args:
        flag (bool): Flag to force rerun.
        outfile (str): Path to output file which may already exist.

    Returns:
        bool: If we should force rerunning of a command

    Examples:
        >>> force_rerun(flag=True, outfile='/not/existing/file.txt')
        True

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
    # If flag is False but filesize of output is 0, also run
    elif not flag and not os.stat(outfile).st_size != 0:
        return True
    # Otherwise, do not run
    else:
        return False


def copy_file_to_new_location(infile, copy_to_folder, rename_to=None):
    """Copy a file to a folder (which must exist) and optionally rename it to something else.

    Args:
        infile: Path to file
        copy_to_folder: Path to folder to copy to
        rename_to: Optional name of new file

    Returns:
        str: Path to new copied file

    """
    orig_name = op.basename(infile)
    if rename_to:
        new_path = op.join(copy_to_folder, rename_to)
    else:
        new_path = op.join(copy_to_folder, orig_name)

    shutil.copy2(infile, new_path)

    if op.isfile(new_path):
        log.debug('{}: copied {}'.format(new_path, orig_name))
    else:
        log.debug('{}: error copying {}'.format(new_path, orig_name))

    return new_path


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


def command_runner(program, args, force_rerun_flag, outfile):
    """Run a program with command-line arguments.

    Args:
        program: Name of the program.
        args: Any program flags as they would be formatted in the command-line (ie. "-i test.in -o test.out").
        force_rerun_flag: If the program should be rerun even if the output file exists.
        outfile: Name out the output file which may have been generated.
            This does not specify what the outfile will be, that should be done in the args.

    Returns:
        bool: If the program ran successfully.

    """
    # Check if pepstats is installed
    if not program_exists(program):
        raise OSError('{}: program not installed'.format(program))

    # Check for force rerunning
    if force_rerun(flag=force_rerun_flag, outfile=outfile):
        cmd = '{} {}'.format(program, args)
        command = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = command.communicate()
        log.debug('{}: Ran program, output to {}'.format(program, outfile))
    else:
        log.debug('{}: Output already exists'.format(outfile))


def dict_head(d, N=5):
    """Return the head of a dictionary. It will be random!

    Default is to return the first 5 key/value pairs in a dictionary.

    Args:
        d: Dictionary to get head.
        N: Number of elements to display.

    Returns:
        dict: the first N items of the dictionary.

    """
    return {k: d[k] for k in list(d.keys())[:N]}


def rank_dated_files(pattern, dir, descending=True):
    """Search a directory for files that match a pattern. Return an ordered list of these files by filename.

    Args:
        pattern: The glob pattern to search for.
        dir: Path to directory where the files will be searched for.
        descending: Default True, will sort alphabetically by descending order.

    Returns:
        list: Rank-ordered list by filename.

    """
    files = glob.glob(op.join(dir, pattern))
    return sorted(files, reverse=descending)


def find(lst, a, case_sensitive=True):
    """Return indices of a list which have elements that match an object or list of objects

    Args:
        lst: list of values
        a: object(s) to check equality
        case_sensitive: if the search should be case sensitive

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
        case_sensitive: if the search should be case sensitive

    Returns:
        list: list of indices of lst which do not equal a

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
        case_sensitive: if the search should be case sensitive

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


def split_list_by_n(l, n):
    """Split a list into lists of size n.

    Args:
        l: List of stuff.
        n: Size of new lists.

    Returns:
        list: List of lists each of size n derived from l.

    """
    n = max(1, n)
    return list(l[i:i+n] for i in range(0, len(l), n))


def split_list_into_n_lists(l, n):
    """Split a list into n lists.

    Args:
        l: List of stuff.
        n: Number of new lists to generate.

    Returns:
        list: List of n lists.

    """
    return [l[i::n] for i in range(n)]


def input_list_parser(infile_list):
    """Always return a list of files with varying input.

    >>> input_list_parser(['/path/to/folder/'])
    ['/path/to/folder/file1.txt', '/path/to/folder/file2.txt', '/path/to/folder/file3.txt']

    >>> input_list_parser(['/path/to/file.txt'])
    ['/path/to/file.txt']

    >>> input_list_parser(['file1.txt'])
    ['file1.txt']

    Args:
        infile_list: List of arguments

    Returns:
        list: Standardized list of files

    """

    final_list_of_files = []

    for x in infile_list:

        # If the input is a folder
        if op.isdir(x):
            os.chdir(x)
            final_list_of_files.extend(glob.glob('*'))

        # If the input is a file
        if op.isfile(x):
            final_list_of_files.append(x)

    return final_list_of_files


def flatlist_dropdup(list_of_lists):
    """Make a single list out of a list of lists, and drop all duplicates.

    Args:
        list_of_lists: List of lists.

    Returns:
        list: List of single objects.

    """
    return list(set([str(item) for sublist in list_of_lists for item in sublist]))


def combinations(iterable, r):
    """Calculate combinations

    >>> list(combinations('ABCD',2))
    [['A', 'B'], ['A', 'C'], ['A', 'D'], ['B', 'C'], ['B', 'D'], ['C', 'D']]

    >>> list(combinations(range(4), 3))
    [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]

    Args:
        iterable: Any iterable object.
        r: Size of combination.

    Yields:
        list: Combination of size r.

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


def percentage_to_float(x):
    """Convert a string representation of a percentage to float.

    >>> percentage_to_float('55%')
    0.55

    Args:
        x: String representation of a percentage

    Returns:
        float: Percentage in decimal form

    """
    return float(x.strip('%')) / 100
