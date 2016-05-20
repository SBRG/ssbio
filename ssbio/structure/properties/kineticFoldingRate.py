import math
import requests
import urllib2
import urllib

R = 8.3144621 # J/K.mol
R = 1.9872041 # cal/K.mol

def get_folding_rate_at_T(gene,T,refT=37):
#def get_folding_rate_at_T(seq,opn,T,refT=37):

    """Predict the kinetic folding rate of a protein at temperature T

    Args:
        seq:  protein sequence
        opn:  structural class: all-alpha, all-beta, mixed, unknown
        T:    temperature
        refT: reference temperature default to 37C

    Returns:
        kinetic folding rate kf (s-1)

    """

    slope = 22000 # not many data on this value, but it's effect on growth rate very small
    # seq = f(gene) # a function linked to ME2.0 to get protein sequence?
    # TODO: SEE BELOW
    # opn = get_ss_class(gene)
    # testing scripts
    #seq='MPRSLKKGPFIDLHLLKKVEKAVESGDKKPLRTWSRRSTIFPNMIGLTIAVHNGRQHVPVFVTDEMVGHKLGEFAPTRTYRGHAADKKAKKK'
    #opn='mixed'

    # get folding rate for the reference temperature
    ref_rate = get_folding_rate_from_url(seq,opn)
    preFactor = float(ref_rate) + slope/(float(refT)+273.15)

    # calculate folding rate at desired temperature
    rate = math.exp(preFactor-slope/(float(T)+273.15))

    return rate

## for more information, see http://www.iitm.ac.in/bioinfo/fold-rate/
## it might be better to get sequence directly from ME2.0 (i.e. whatever that's being modeled)
## not from structure here (which may include some wired one letter codes, which need to be fixed)?

def get_folding_rate_from_url(seq,opn='mixed'):

    """Predict the kinetic folding rate of a protein

    Args:
        seq: protein sequence
        opn: structural class: all-alpha, all-beta, mixed, unknown

    Returns:
        kinetic folding rate ln(kf) ln(s-1)

    """

    url = 'http://www.iitm.ac.in/bioinfo/cgi-bin/fold-rate/foldrateCalculator.pl'

    values={'sequence':seq,'eqn':opn}
    data = urllib.urlencode(values)
    response = urllib2.urlopen(url,data)
    result = response.read()
    ind = str.find(result,'The folding rate,')
    result2 = result[ind:ind+70]
    ind1 = str.find(result2,'=')
    ind2 = str.find(result2,'/sec')
    rate = result2[ind1+2:ind2]
    return rate


## the final form should be folding_rate_at_T(gene,T,ref=37)

def folding_rate_at_T(seq,opn,T,refT=37):

    """Predict the kinetic folding rate of a protein at temperature T

    Args:
        seq:  protein sequence
        opn:  structural class: all-alpha, all-beta, mixed, unknown
        T:    temperature
        refT: reference temperature default to 37C

    Returns:
        kinetic folding rate kf (s-1)

    """

    slope = 22000 # not many data on this value, it's effect on growth rate very small
    # seq = f(gene) # a function linked to ME2.0 to get protein sequence?
    # opn = get_ss_class(gene)

    # get folding rate for the reference temperature
    ref_rate = get_folding_rate_from_url(seq,opn)
    preFactor = float(ref_rate) + slope/(float(refT)+273.15)

    # calculate folding rate at desired temperature
    rate = math.exp(preFactor-slope/(float(T)+273.15))

    return rate

## The following two procedure gets the secondary structure class
## i have it still with the original format you gave me
## please change accordingly if you have reorganized file format

# get secondary structure across chain
def get_dssp_ss_content_multiplechains(prag,chain):
    a=prag.getData('resnum')
    b=prag.getData('name')
    c=prag.getData('secondary')
    d=prag.getData('chain')
    resid=[]
    SSstr=[]
    idx = [i for i, x in enumerate(d) if x == chain]
    for i in idx:
        if b[i]=='CA':
            resid.append(i)
            SSstr.append(c[i])

    N = float(len(SSstr))
    helix_alpha=SSstr.count('H')/N
    helix_3_10=SSstr.count('G')/N
    extended=SSstr.count('E')/N
    return helix_alpha,helix_3_10,extended


# define secondary structure class
def get_ss_class(gene):

    # how do you store file now? need to change here...
    dsspfile = gene+'model1.dssp'
    pdbfile = gene+'model1.pdb'

    prag = pr.parsePDB(pdbfile)
    pr.parseDSSP(dsspfile,prag)
    alpha,threeTen,beta = get_dssp_ss_content_multiplechains(prag,chain)

    if alpha == 0 and beta > 0:
        classification = 'all-beta'
    elif beta == 0 and alpha > 0:
        classification = 'all-alpha'
    elif beta == 0 and alpha == 0:
        classification = 'mixed'
    elif float(alpha)/beta >= 20:
        classification = 'all-alpha'
    else:
        classification = 'mixed'

    # if structure not found, classification = 'unknown'

    return classification
