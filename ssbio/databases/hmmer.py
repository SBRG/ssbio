import urllib.parse as urllib
import urllib.request as urllib2
from urllib.request import urlopen, Request
from urllib.error import URLError
import json

import os.path as op


def manual_get_pfam_annotations(seq, outpath, searchtype='phmmer', force_rerun=False):
    """Retrieve and download PFAM results from the HMMER search tool.

    Args:
        seq:
        outpath:
        searchtype:
        force_rerun:

    Returns:

    Todo:
        * Document and test!

    """
    if op.exists(outpath):
        with open(outpath, 'r') as f:
            json_results = json.loads(json.load(f))

    else:
        fseq = '>Seq\n' + seq
        if searchtype == 'phmmer':
            parameters = {'seqdb': 'pdb', 'seq': fseq}
        if searchtype == 'hmmscan':
            parameters = {'hmmdb': 'pfam', 'seq': fseq}
        enc_params = urllib.urlencode(parameters).encode('utf-8')
        request = urllib2.Request('http://www.ebi.ac.uk/Tools/hmmer/search/{}'.format(searchtype), enc_params)
        url = (urllib2.urlopen(request).geturl() + '?output=json')
        request = str(url)
        request_read = urlopen(request).read().decode("utf-8")

        with open(outpath, 'w') as f:
            json.dump(request_read, f)

        json_results = json.loads(request_read)

    return json_results['results']['hits']