from __future__ import print_function
import os
import sys
import datetime
import os.path as op
import glob
import pandas as pd
import numpy as np
import contextlib
import logging
import distutils.spawn
import subprocess
import shlex
import requests
import json
import warnings
import gzip
from collections import OrderedDict
from collections import Callable
import operator

log = logging.getLogger(__name__)


def is_ipynb():
    """Return True if the module is running in IPython kernel,
    False if in IPython shell or other Python shell.

    Copied from: http://stackoverflow.com/a/37661854/1592810
    There are other methods there too

    >>> is_ipynb()
    False

    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole?
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal running IPython?
            return False
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


class Date():
    """Various methods to return formatted dates for today.
    """
    def __init__(self):
        self.today = datetime.date.today()
        self.short_date = self.today.strftime('%y%m%d')
        self.long_date = self.today.strftime('%Y-%m-%d')


def todays_short_date():
    today = Date()
    return today.short_date


def todays_long_date():
    today = Date()
    return today.long_date


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


def clean_df(df, fill_nan=True, drop_empty_columns=True):
    """Clean a pandas dataframe by:
        1. Filling empty values with Nan
        2. Dropping columns with all empty values

    Args:
        df: Pandas DataFrame
        fill_nan (bool): If any empty values (strings, None, etc) should be replaced with NaN
        drop_empty_columns (bool): If columns whose values are all empty should be dropped

    Returns:
        DataFrame: cleaned DataFrame

    """
    if fill_nan:
        df = df.fillna(value=np.nan)
    if drop_empty_columns:
        df = df.dropna(axis=1, how='all')
    return df.sort_index()


def clean_single_dict(indict, prepend_to_keys=None, remove_keys_containing=None):
    """Clean a dict with values that contain single item iterators to single items

    Args:
        indict (dict): Dictionary to be cleaned
        prepend_to_keys (str): String to prepend to all keys
        remove_keys_containing (str): Text to check for in keys to ignore

    Returns:
         dict: Cleaned dictionary

    Examples:
        >>> clean_single_dict(indict={'test1': [1], 'test2': ['H']})
        {'test1': 1, 'test2': 'H'}

        >>> clean_single_dict(indict={'test1': [1], 'test2': ['H']}, prepend_to_keys='struct_')
        {'struct_test1': 1, 'struct_test2': 'H'}

        >>> clean_single_dict(indict={'test1': [1], 'ignore': ['H']}, prepend_to_keys='struct_', remove_keys_containing='ignore')
        {'struct_test1': 1}

    """
    if not prepend_to_keys:
        prepend_to_keys = ''

    outdict = {}
    for k, v in indict.items():
        if remove_keys_containing:
            if remove_keys_containing in k:
                continue
        outdict[prepend_to_keys + k] = v[0]

    return outdict


def double_check_attribute(object, setter, backup_attribute, custom_error_text=None):
    """Check if a parameter to be used is None, if it is, then check the specified backup attribute and throw
    an error if it is also None.

    Args:
        object: The original object
        setter: Any input object
        backup_attribute (str): Attribute in <object> to be double checked
        custom_error_text (str): If a custom string for the error should be raised

    Raises:
         ValueError: If both setter and backup_attribute are None

    """
    if not setter:
        next_checker = getattr(object, backup_attribute)
        if not next_checker:
            if custom_error_text:
                raise ValueError(custom_error_text)
            else:
                raise ValueError('Attribute replacing "{}" must be specified'.format(backup_attribute))


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


def is_non_zero_file(fpath):
    """Check if a file exists, and that it contains some kind of contents

    Args:
        fpath (str): Path to file

    Returns:
        bool: If file exists and contains contents

    """
    return op.isfile(fpath) and os.path.getsize(fpath) > 0


def outfile_maker(inname, outext='.out', outname='', outdir='', append_to_name=''):
    """Create a default name for an output file based on the inname name, unless a output name is specified.

    Args:
        inname: Path to input file
        outext: Optional specified extension for output file (with the "."). Default is ".out".
        outfile: Optional specified name of output file.
        outdir: Optional path to output directory.

    Returns:
        str: Path to final output destination.

    Examples:

        >>> outfile_maker(inname='P00001.fasta')
        'P00001.out'

        >>> outfile_maker(inname='P00001')
        'P00001.out'

        >>> outfile_maker(inname='P00001.fasta', append_to_name='_new')
        'P00001_new.out'

        >>> outfile_maker(inname='P00001.fasta', outext='.mao')
        'P00001.mao'

        >>> outfile_maker(inname='P00001.fasta', outext='.mao', append_to_name='_new')
        'P00001_new.mao'

        >>> outfile_maker(inname='P00001.fasta', outext='.new', outname='P00001_aligned')
        'P00001_aligned.new'

        >>> outfile_maker(inname='P00001.fasta', outname='P00001_aligned')
        'P00001_aligned.out'

        >>> outfile_maker(inname='P00001.fasta', outname='P00001_aligned', append_to_name='_new')
        'P00001_aligned_new.out'

        >>> outfile_maker(inname='P00001.fasta', outname='P00001_aligned', outdir='/my/dir/')
        '/my/dir/P00001_aligned.out'

        >>> outfile_maker(inname='/test/other/dir/P00001.fasta', append_to_name='_new')
        '/test/other/dir/P00001_new.out'

        >>> outfile_maker(inname='/test/other/dir/P00001.fasta', outname='P00001_aligned')
        '/test/other/dir/P00001_aligned.out'

        >>> outfile_maker(inname='/test/other/dir/P00001.fasta', outname='P00001_aligned', outdir='/my/dir/')
        '/my/dir/P00001_aligned.out'

    """

    # TODO: CHECK IF OUTNAME IS A VALID FILE NAME!
    orig_dir, orig_name, orig_ext = split_folder_and_path(inname)

    # If output filename not provided, default is to take name of inname
    if not outname:
        outname = orig_name

    # Create new path in the same directory of old path if a new one isn't specified
    if not outdir:
        outdir = orig_dir

    # Append additional stuff to the filename if specified
    if append_to_name:
        outname += append_to_name

    # Join the output filename and output extension
    final_outfile = op.join(outdir, '{}{}'.format(outname, outext))

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

        >>> force_rerun(flag=True, outfile='./utils.py')
        True

        >>> force_rerun(flag=False, outfile='./utils.py')
        False

    """
    # If flag is True, always run
    if flag:
        return True
    # If flag is False but file doesn't exist, also run
    elif not flag and not op.exists(outfile):
        return True
    # If flag is False but filesize of output is 0, also run
    elif not flag and not is_non_zero_file(outfile):
        return True
    # Otherwise, do not run
    else:
        return False


def gunzip_file(infile, outfile=None, outdir=None, delete_original=False, force_rerun_flag=False):
    """Decompress a gzip file and optionally set output values.

    Args:
        infile: Path to .gz file
        outfile: Name of output file
        outdir: Path to output directory
        delete_original: If original .gz file should be deleted
        force_rerun_flag: If file should be decompressed if outfile already exists

    Returns:
        str: Path to decompressed file

    """
    if not outfile:
        outfile = infile.replace('.gz', '')

    if not outdir:
        outdir = ''
    else:
        outdir = op.dirname(infile)
    outfile = op.join(outdir, op.basename(outfile))

    if force_rerun(flag=force_rerun_flag, outfile=outfile):
        gz = gzip.open(infile, "rb")
        decoded = gz.read()

        with open(outfile, "wb") as new_file:
            new_file.write(decoded)

        gz.close()
        log.debug('{}: file unzipped'.format(outfile))
    else:
        log.debug('{}: file already unzipped'.format(outfile))

    if delete_original:
        os.remove(infile)

    return outfile


def request_file(link, outfile, force_rerun_flag=False):
    """Download a file given a URL if the outfile does not exist already.

    Args:
        link (str): Link to download file.
        outfile (str): Path to output file, will make a new file if it does not exist. Will not download if it does
            exist, unless force_rerun_flag is True.
        force_rerun_flag (bool): Flag to force re-downloading of the file if it exists already.

    Returns:
        str: Path to downloaded file.

    """
    if force_rerun(flag=force_rerun_flag, outfile=outfile):
        req = requests.get(link)
        if req.status_code == 200:
            with open(outfile, 'w') as f:
                f.write(req.text)
            log.debug('Loaded and saved {} to {}'.format(link, outfile))
        else:
            log.error('{}: request error {}'.format(link, req.status_code))
    return outfile


def request_json(link, outfile, force_rerun_flag, outdir=None):
    """Download a file in JSON format from a web request

    Args:
        link: Link to web request
        outfile: Name of output file
        outdir: Directory of output file
        force_rerun_flag: If true, redownload the file

    Returns:
        dict: contents of the JSON request

    """
    if not outdir:
        outdir = ''
    outfile = op.join(outdir, outfile)

    if force_rerun(flag=force_rerun_flag, outfile=outfile):
        text_raw = requests.get(link)
        my_dict = text_raw.json()
        with open(outfile, 'w') as f:
            json.dump(my_dict, f)

        log.debug('Loaded and saved {} to {}'.format(link, outfile))
    else:
        with open(outfile, 'r') as f:
            my_dict = json.load(f)
        log.debug('Loaded {}'.format(outfile))

    return my_dict


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


def command_runner(shell_command, force_rerun_flag, outfile_checker, cwd=None, silent=False):
    """Run a shell command with subprocess, with additional options to check if output file exists and printing stdout.

    Args:
        shell_command (str): Command as it would be formatted in the command-line (ie. "program -i test.in -o test.out").
        force_rerun_flag: If the program should be rerun even if the output file exists.
        outfile_checker (str): Name out the output file which may have been generated. This does not specify what the outfile
            will be, that should be done in the program's args or predetermined.
        cwd (str): Path to working directory where command will be executed.
        silent (bool): If program STDOUT should be printed to the current shell.

    Returns:
        bool: If the program ran successfully.

    """
    program_and_args = shlex.split(shell_command)

    # Check if program is installed
    if not program_exists(program_and_args[0]):
        raise OSError('{}: program not installed'.format(program_and_args[0]))

    # Format outfile if working in cwd
    if cwd:
        # TODO: should this be done, or should user explicitly define whole outfile path?
        outfile_checker = op.join(cwd, op.basename(outfile_checker))

    # Check for force rerunning
    if force_rerun(flag=force_rerun_flag, outfile=outfile_checker):
        if silent:
            command = subprocess.Popen(program_and_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            out, err = command.communicate()
            ret = command.returncode
        else:
            # Prints output
            for path in execute(cmd=program_and_args, cwd=cwd):
                print(path, end="")

        # TODO: check return code and log properly
        log.debug('{}: Ran program, output to {}'.format(program_and_args[0], outfile_checker))
    else:
        log.debug('{}: Output already exists'.format(outfile_checker))


def execute(cmd, cwd=None):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, universal_newlines=True)
    # TODO: should also print stderr or give an option for that
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def write_torque_script(command, outfile, walltime, queue, name, out, err, print_exec=True):

    with open(outfile, 'w') as script:
        script.write('#PBS -l walltime={}\n'.format(walltime))
        script.write('#PBS -q regular\n')
        script.write('#PBS -N {}\n'.format(name))
        script.write('#PBS -o {}.o\n'.format(out))
        script.write('#PBS -e {}.e\n'.format(err))
        script.write('cd ${PBS_O_WORKDIR}\n')
        script.write(command)

    os.chmod(outfile, 0o755)

    if print_exec:
        print('qsub {};'.format(outfile))

    return outfile


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
        val: object to parse into a string

    Returns:
        str: String representation

    """
    if val is None:
        return ''
    if isinstance(val, list):
        newval = [str(x) for x in val]
        return ';'.join(newval)
    if isinstance(val, str):
        return val
    else:
        return str(val)


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


def conv_to_float(indata, inf_str=''):
    """Try to convert an arbitrary string to a float. Specify what will be replaced with "Inf".
    
    Args:
        indata (str): String which contains a float
        inf_str (str): If string contains something other than a float, and you want to replace it with float("Inf"), 
            specify that string here.

    Returns:
        float: Converted string representation

    """
    if indata.strip() == inf_str:
        outdata = float('Inf')
    else:
        try:
            outdata = float(indata)
        except:
            raise ValueError('Unable to convert {} to float'.format(indata))

    return outdata


def make_dir(path):
    """Make a directory if it does not already exist

    Args:
        path: Path to new directory

    """
    if not op.exists(path):
        os.mkdir(path)
    return path


def remap( x, oMin, oMax, nMin, nMax ):
    """Map to a 0 to 1 scale
        http://stackoverflow.com/questions/929103/convert-a-number-range-to-another-range-maintaining-ratio

    """

    #range check
    if oMin == oMax:
        log.warning("Zero input range, unable to rescale")
        return x

    if nMin == nMax:
        log.warning("Zero output range, unable to rescale")
        return x

    #check reversed input range
    reverseInput = False
    oldMin = min( oMin, oMax )
    oldMax = max( oMin, oMax )
    if not oldMin == oMin:
        reverseInput = True

    #check reversed output range
    reverseOutput = False
    newMin = min( nMin, nMax )
    newMax = max( nMin, nMax )
    if not newMin == nMin :
        reverseOutput = True

    portion = (x-oldMin)*(newMax-newMin)/(oldMax-oldMin)
    if reverseInput:
        portion = (oldMax-x)*(newMax-newMin)/(oldMax-oldMin)

    result = portion + newMin
    if reverseOutput:
        result = newMax - portion

    return result


def scale_calculator(multiplier, elements, rescale=None):
    """Get a dictionary of scales for each element in elements.

    Examples:
        >>> scale_calculator(1, [2,7,8])
        {8: 1, 2: 1, 7: 1}

        >>> scale_calculator(1, [2,2,2,3,4,5,5,6,7,8])
        {2: 3, 3: 1, 4: 1, 5: 2, 6: 1, 7: 1, 8: 1}

        >>> scale_calculator(1, [2,2,2,3,4,5,5,6,7,8], rescale=(0.5,1))
        {2: 1.0, 3: 0.5, 4: 0.5, 5: 0.75, 6: 0.5, 7: 0.5, 8: 0.5}

        >>> scale_calculator(1, {2:3, 3:1, 4:1, 5:2, 6:1, 7:1, 8:1}, rescale=(0.5,1))
        {2: 1.0, 3: 0.5, 4: 0.5, 5: 0.75, 6: 0.5, 7: 0.5, 8: 0.5}

        >>> scale_calculator(1, [(2,2,2),(3,),(4,),(5,),(5,),(6,7,8)], rescale=(0.5,1))
        {(2, 2, 2): 0.5, (5,): 1.0, (3,): 0.5, (6, 7, 8): 0.5, (4,): 0.5}

    Args:
        mutiplier (int, float): Base float to be multiplied
        elements (list, dict): Dictionary which contains object:count
            or list of objects that may have repeats which will be counted
        rescale (tuple): Min and max values to rescale to

    Returns:
        dict: Scaled values of mutiplier for each element in elements

    """

    # TODO: think about what happens when:
    # TODO: 1. there is only one (or n) of each element, and rescale is set to seomthing. what is the original min/max to scale from?
    # TODO: 2. can we normalize the scale based on other counts? (ie. other gene mutation frequencies)

    if isinstance(elements, list):
        unique_elements = list(set(elements))
        scales = {}
        for x in unique_elements:
            count = elements.count(x)
            scales[x] = multiplier * count
    elif isinstance(elements, dict):
        scales = {}
        for k,count in elements.items():
            scales[k] = multiplier * int(count)
    else:
        raise ValueError('Input list of elements or dictionary of elements & counts')

    if not rescale:
        return scales
    else:
        new_scales = {}
        for k,v in scales.items():
            new_scales[k] = remap(v, min(scales.values()), max(scales.values()), rescale[0], rescale[1])
        return new_scales


def label_sequential_regions(inlist):
    """Input a list of labeled tuples and return a dictionary of sequentially labeled regions.

    Args:
        inlist (list): A list of tuples with the first number representing the index and the second the index label.

    Returns:
        dict: Dictionary of labeled regions.

    Examples:

        >>> label_sequential_regions([(1, 'O'), (2, 'O'), (3, 'O'), (4, 'M'), (5, 'M'), (6, 'I'), (7, 'M'), (8, 'O'), (9, 'O')])
        {'O1': [1, 2, 3], 'M1': [4, 5], 'I1': [6], 'M2': [7], 'O2': [8, 9]}

    """
    import more_itertools as mit

    df = pd.DataFrame(inlist).set_index(0)

    labeled = {}
    for label in df[1].unique():
        iterable = df[df[1] == label].index.tolist()
        labeled.update({'{}{}'.format(label, i + 1): items for i, items in
                        enumerate([list(group) for group in mit.consecutive_groups(iterable)])})

    return labeled


def check_condition(left, condition, right):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq}
    if condition not in ops:
        raise KeyError('{}: condition not supported'.format(condition))
    return ops[condition](left, right)
