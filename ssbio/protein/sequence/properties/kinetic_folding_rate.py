__author__ = 'Ke Chen'
__email__ = "kec003@ucsd.edu"

# TODO: replace urllib usage with six library
try:
    from urllib.request import urlopen
    from urllib.request import build_opener
    from urllib.request import HTTPCookieProcessor
    from urllib.parse import urlparse
    from urllib.parse import urlencode
except ImportError:
    from urlparse import urlparse
    from urllib2 import urlopen
    from urllib import urlencode
    from urllib2 import build_opener
    from urllib2 import HTTPCookieProcessor

import math
import scipy.constants
import ssbio.protein.sequence.utils

# R (molar gas constant) in calories
r_j = scipy.constants.R
r_cal = scipy.constants.R / scipy.constants.calorie


def get_foldrate(seq, secstruct):
    """Submit sequence and structural class to FOLD-RATE calculator (http://www.iitm.ac.in/bioinfo/fold-rate/)
    to calculate kinetic folding rate.

    Args:
        seq (str, Seq, SeqRecord): Amino acid sequence
        secstruct (str): Structural class: `all-alpha``, ``all-beta``, ``mixed``, or ``unknown``

    Returns:
        float: Kinetic folding rate k_f

    """

    seq = ssbio.protein.sequence.utils.cast_to_str(seq)

    url = 'http://www.iitm.ac.in/bioinfo/cgi-bin/fold-rate/foldrateCalculator.pl'

    values = {'sequence': seq, 'eqn': secstruct}
    data = urlencode(values)
    data = data.encode('ASCII')
    response = urlopen(url, data)

    result = str(response.read())
    ind = str.find(result, 'The folding rate,')
    result2 = result[ind:ind + 70]
    ind1 = str.find(result2, '=')
    ind2 = str.find(result2, '/sec')
    rate = result2[ind1 + 2:ind2]

    return rate


def get_foldrate_at_temp(ref_rate, new_temp, ref_temp=37.0):
    """Scale the predicted kinetic folding rate of a protein to temperature T, based on the relationship ln(k_f)∝1/T

    Args:
        ref_rate (float): Kinetic folding rate calculated from the function :func:`~ssbio.protein.sequence.properties.kinetic_folding_rate.get_foldrate`
        new_temp (float): Temperature in degrees C
        ref_temp (float): Reference temperature, default to 37 C

    Returns:
        float: Kinetic folding rate k_f at temperature T

    """

    # Not much data available on this slope value, however its effect on growth rate in a model is very small
    slope = 22000

    # Get folding rate for the reference temperature
    preFactor = float(ref_rate) + slope / (float(ref_temp) + 273.15)

    # Calculate folding rate at desired temperature
    rate = math.exp(preFactor - slope / (float(new_temp) + 273.15))

    return rate