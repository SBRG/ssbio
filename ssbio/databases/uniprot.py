from tqdm import tqdm
from bioservices.uniprot import UniProt

bsup = UniProt()
import pandas as pd
import re
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import warnings

def uniprot_valid_id(instring):
    # regex from: http://www.uniprot.org/help/accession_numbers
    valid_id = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    if valid_id.match(instring):
        return True
    else:
        return False


# status can be "unreviewed" or "reviewed"
def uniprot_reviewed_checker(uniprot_ids):

    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]

    invalid_ids = [i for i in uniprot_ids if not uniprot_valid_id(i)]
    uniprot_ids = [i for i in uniprot_ids if uniprot_valid_id(i)]

    if invalid_ids:
        warnings.warn("Invalid UniProt IDs {} will be ignored".format(invalid_ids))

    # splitting query up into managable sizes (200 IDs each)
    Nmax = 200
    N, rest = divmod(len(uniprot_ids), Nmax)

    uni_rev_dict = {}

    if rest>0:
        N+=1
    for i in range(0,N):
        i1 = i*Nmax
        i2 = (i+1)*Nmax
        if i2>len(uniprot_ids):
            i2 = len(uniprot_ids)

        query = uniprot_ids[i1:i2]

        query_string = ''
        for x in query:
            query_string += 'id:' + x + '+OR+'
        query_string = query_string.strip('+OR+')

        uni_rev_raw = StringIO(bsup.search(query_string, columns='id,reviewed', frmt='tab'))
        uni_rev_df = pd.read_table(uni_rev_raw, sep='\t', index_col=0)
        uni_rev_df = uni_rev_df.fillna(False)

        # no_metadata = uni_rev_df[pd.isnull(uni_rev_df.Status)].index.tolist()
        # if no_metadata:
        #     warnings.warn("Unable to retrieve metadata for {}.".format(no_metadata))
        uni_rev_df = uni_rev_df[pd.notnull(uni_rev_df.Status)]

        uni_rev_df = uni_rev_df.replace(to_replace="reviewed", value=True)
        uni_rev_df = uni_rev_df.replace(to_replace="unreviewed", value=False)
        uni_rev_dict_adder = uni_rev_df.to_dict()['Status']
        uni_rev_dict.update(uni_rev_dict_adder)

    return uni_rev_dict