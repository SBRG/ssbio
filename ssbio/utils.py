import datetime
import os.path as op
import glob
import os
import pandas as pd

def force_string(val=None):
    """Force a string representation of an object

    Args:
        val: object to parse into a string

    Returns:

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

def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

def input_list_parser(instring):
    """Always return a list of files with varying input

    1. /path/to/folder -> list of files in folder (full paths)
    2. /path/to/file -> list of files (singular list)
    3. file1,file2 -> list of files
    4.

    Args:
        instring:

    Returns:

    """
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
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

class Date():
    def __init__(self):
        self.short_date = self.date_prefix()

    def date_prefix(self):
        today = datetime.date.today()
        return today.strftime('%y%m%d')

if __name__ == '__main__':
    d = Date()
    print(d.short_date)