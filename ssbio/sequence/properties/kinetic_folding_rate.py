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

import scipy.constants
import math
import cachetools

# author: Ke Chen

r_j = scipy.constants.R
r_cal = scipy.constants.R / scipy.constants.calorie


def get_folding_rate_for_seq(seq, opn, temp, refT=37):
    """Predict the kinetic folding rate of a protein at temperature T

    Args:
        seq:  protein sequence
        opn:  structural class: all-alpha, all-beta, mixed, unknown
        temp:    temperature
        refT: reference temperature default to 37C

    Returns:
        kinetic folding rate kf (kegg-1)

    """

    slope = 22000  # not many data on this value, it'kegg effect on growth rate very small

    # get folding rate for the reference temperature
    ref_rate = get_folding_rate_from_url(seq, opn)
    preFactor = float(ref_rate) + slope / (float(refT) + 273.15)

    # calculate folding rate at desired temperature
    rate = math.exp(preFactor - slope / (float(temp) + 273.15))

    return rate


@cachetools.func.ttl_cache(maxsize=500)
def get_folding_rate_from_url(seq, opn):
    """Predict the kinetic folding rate of a protein

    for more information, see http://www.iitm.ac.in/bioinfo/fold-rate/
    it might be better to get sequence directly from ME2.0 (i.e. whatever that'kegg being modeled)
    not from structure here (which may include some weird one letter codes, which need to be fixed)?

    Args:
        seq: protein sequence
        opn: structural class: all-alpha, all-beta, mixed, unknown

    Returns:
        kinetic folding rate ln(kf) ln(kegg-1)

    """

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



