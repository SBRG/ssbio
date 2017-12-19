"""
FOLD-RATE - Kinetic Folding Rate
==================================

Description
-----------

Home page: FOLD-RATE_

This module provides a function to predict the kinetic folding rate (k_f) given an amino acid sequence and its
structural classficiation (alpha/beta/mixed).

Instructions
------------

#. Obtain your protein's sequence
#. Determine the main secondary structure composition of the protein (``all-alpha``, ``all-beta``, ``mixed``, or ``unknown``)
#. Input the sequence and secondary structure composition into the function ``get_foldrate``

FAQs
----

* What is the main secondary structure composition of my protein?

    * ``all-alpha`` = dominated by α-helices; α > 40% and β < 5%
    * ``all-beta`` = dominated by β-strands; β > 40% and α < 5%
    * ``mixed`` = contain both α-helices and β-strands; α > 15% and β > 10%

* What is the kinetic folding rate?

    * Protein folding rate is a measure of slow/fast folding of proteins from the unfolded state to native
      three-dimensional structure.

* What units is it in?

    * Number of proteins folded per second

* How can I install FOLD-RATE?

    * FOLD-RATE is only available as a web server. *ssbio* provides a wrapper for the web server and allows you to
      submit protein sequences to it along with caching the output files.

* How do I cite FOLD-RATE?

    * Gromiha MM, Thangakani AM & Selvaraj S (2006) FOLD-RATE: prediction of protein folding rates from amino acid
      sequence. Nucleic Acids Res. 34: W70–4 Available at: http://dx.doi.org/10.1093/nar/gkl043

* How can this parameter be used on a genome-scale?

    * See: Chen K, Gao Y, Mih N, O'Brien EJ, Yang L & Palsson BO (2017) Thermosensitivity of growth is determined by
      chaperone-mediated proteome reallocation. Proceedings of the National Academy of Sciences 114: 11548–11553
      Available at: http://www.pnas.org/content/114/43/11548.abstract

API
---
.. API will show up below...

.. Links
.. _FOLD-RATE: http://www.iitm.ac.in/bioinfo/fold-rate/

"""

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


def get_folding_rate_for_seq(seq, secstruct, temp, refT=37.0):
    """Scale the predicted kinetic folding rate of a protein to temperature T, based on the relationship ln(k_f)∝1/T

    Args:
        seq (str, Seq, SeqRecord): Amino acid sequence
        secstruct (str): Structural class: ``all-alpha``, ``all-beta``, ``mixed``, or ``unknown``
        temp (float): Temperature in degrees C
        refT (float): Reference temperature, default to 37 C

    Returns:
        float: Kinetic folding rate k_f at temperature T.

    """

    # Not much data available on this slope value, however its effect on growth rate in a model is very small
    slope = 22000

    # Get folding rate for the reference temperature
    ref_rate = get_foldrate(seq, secstruct)
    preFactor = float(ref_rate) + slope / (float(refT) + 273.15)

    # Calculate folding rate at desired temperature
    rate = math.exp(preFactor - slope / (float(temp) + 273.15))

    return rate