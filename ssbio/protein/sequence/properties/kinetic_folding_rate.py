"""This module provides a function to predict the kinetic folding rate (k_f) given an amino acid sequence and its 
    structure classficiation (alpha/beta/mixed). The main method is a web wrapper for the FOLD-RATE calculator 
    available at http://www.iitm.ac.in/bioinfo/fold-rate/

This method is adapted from:

    Gromiha, M. M., Thangakani, A. M., & Selvaraj, S. (2006). 'FOLD-RATE: prediction of protein folding rates from 
        amino acid sequence', Nucleic acids research, 34/Web Server issue: W70–4. DOI: 10.1093/nar/gkl043

For an example of usage of this parameter in a genome-scale model:

    Chen, K., Gao, Y., Mih, N., O'Brien, E., Yang, L., Palsson, B.O. (2017). 
        'Thermo-sensitivity of growth is determined by chaperone-mediated proteome re-allocation.',
        Submitted to PNAS.

"""

__author__ = 'Ke Chen'
__email__ = "kec003@eng.ucsd.edu"

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


def get_folding_rate_for_seq(seq, opn, temp, refT=37.0):
    """Predict the kinetic folding rate of a protein at temperature T.

    Args:
        seq (str, Seq, SeqRecord): Amino acid sequence
        opn (str): Structural class: 'all-alpha', 'all-beta', 'mixed', or 'unknown'
        temp (float): Temperature in degrees C
        refT (float): Reference temperature, default to 37 C

    Returns:
        float: Kinetic folding rate kf (kegg-1)

    """

    # Not much data available on this slope value, however its effect on growth rate in a model is very small
    slope = 22000

    # Get folding rate for the reference temperature
    ref_rate = get_folding_rate_from_url(seq, opn)
    preFactor = float(ref_rate) + slope / (float(refT) + 273.15)

    # Calculate folding rate at desired temperature
    rate = math.exp(preFactor - slope / (float(temp) + 273.15))

    return rate


def get_folding_rate_from_url(seq, opn):
    """Submit sequence and structural class to FOLD-RATE calculator (http://www.iitm.ac.in/bioinfo/fold-rate/)
        to calculate kinetic folding rate.

    Args:
        seq (str, Seq, SeqRecord): Amino acid sequence
        opn (str): Structural class: 'all-alpha', 'all-beta', 'mixed', or 'unknown'

    Returns:
        float: Kinetic folding rate ln(kf) ln(kegg-1)

    """

    seq = ssbio.protein.sequence.utils.cast_to_str(seq)

    url = 'http://www.iitm.ac.in/bioinfo/cgi-bin/fold-rate/foldrateCalculator.pl'

    values = {'sequence': seq, 'eqn': opn}
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