# CLEARING OUTPUT

import os
import sys
import time
import collections
from collections import defaultdict
from contextlib import contextmanager

try:
    from IPython.display import clear_output
    have_ipython = True
except ImportError:
    have_ipython = False

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


# MATPLOTLIB FIGURES
# import seaborn as sns
# sns.set_style("dark")
# sns.set_context("poster")

# from matplotlib import rcParams
# rcParams['figure.dpi'] = 300
# rcParams['lines.linewidth'] = 2
# rcParams['axes.facecolor'] = 'white'
# rcParams['patch.edgecolor'] = 'white'

# from matplotlib import font_manager as fm
# proptease = fm.FontProperties()


def flatlist_dropdup(list_of_lists):
    return list(set([str(item) for sublist in list_of_lists for item in sublist]))


def de_unicodeify(data):
    if isinstance(data, basestring):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(de_unicodeify, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(de_unicodeify, data))
    else:
        return data