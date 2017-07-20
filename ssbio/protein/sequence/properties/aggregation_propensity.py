"""This module provides a function to predict the aggregation propensity of proteins, specifically the number 
    of aggregation-prone segments on an unfolded protein sequence. The main method is a web wrapper for the 
    AMYLPRED2 web server available at http://aias.biol.uoa.gr/AMYLPRED2/. In order to obtain the best balance between 
    sensitivity and specificity, we follow the author’s guidelines to consider every 5 consecutive residues agreed 
    among at least 5 methods contributing 1 to the aggregation propensity

This method is adapted from:

    Tsolis, A. C., Papandreou, N. C., Iconomidou, V. A., & Hamodrakas, S. J. (2013). 'A consensus method for the 
        prediction of "aggregation-prone" peptides in globular proteins', PloS one, 8/1: e54175. 
        DOI: 10.1371/journal.pone.0054175

For an example of usage of this parameter in a genome-scale model:

    Chen, K., Gao, Y., Mih, N., O'Brien, E., Yang, L., Palsson, B.O. (2017). 
        'Thermo-sensitivity of growth is determined by chaperone-mediated proteome re-allocation.',
        Submitted to PNAS.

"""
from __future__ import print_function

__author__ = 'Ke Chen'
__email__ = "kec003@eng.ucsd.edu"

import glob
import time
import requests
import os.path as op
import ssbio.protein.sequence.utils

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
try:
    from http.cookiejar import CookieJar
except ImportError:
    from cookielib import CookieJar


class AMYLPRED:
    def __init__(self, email, password):
        """AMYLPRED2 requires registration at http://aias.biol.uoa.gr/AMYLPRED2/. Set your email and password 
            used to login here.
        
        Args:
            email (str): Account email 
            password (str): Account password
             
        """

        self.email = email
        self.password = password

    def get_aggregation_propensity(self, seq, cutoff_v=5, cutoff_n=5):
        """Predict aggregation propensity.

        Args:
            seq (str, Seq, SeqRecord): Amino acid sequence
            cutoff_v (int): The minimal number of methods that agree on a residue being a aggregation-prone residue
            cutoff_n (int): The minimal number of consecutive residues to be considered as a 'stretch' of 
                aggregation-prone region

        Returns:
            int: Aggregation propensity

        """

        seq = ssbio.protein.sequence.utils.cast_to_str(seq)

        output = self.consensus_aggregation(seq)
        agg_index, agg_conf = self.write_propensity(len(seq), output, cutoff_v=cutoff_v, cutoff_n=cutoff_n)

        return agg_index

    def consensus_aggregation(self, seq, outfile, run_amylmuts=False):
        """Run the AMYLPRED2 web server for an amino acid sequence. 
        
        AMYLMUTS is an optional method to run as it is the most time consuming and generates a slightly different
            result every submission.
        
        Args:
            seq (str): Amino acid sequence as a string
            outdir (str): Directory to where output files should be saved
            run_amylmuts (bool): If AMYLMUTS method should be run, default False

        Returns:

        """
        url = "http://aias.biol.uoa.gr/AMYLPRED2/login.php"
        cj = CookieJar()
        opener = build_opener(HTTPCookieProcessor(cj))
        formdata = {"email": self.email, "password": self.password}
        data_encoded = urlencode(formdata)
        data_encoded = data_encoded.encode('ASCII')
        response = opener.open(url, data_encoded)

        # AMYLMUTS is most time consuming and generates a slightly different result every submission
        methods = ['AGGRESCAN', 'NETCSSP', 'PAFIG', 'APD', 'AMYLPATTERN',
                   'SECSTR', 'BSC', 'WALTZ', 'CONFENERGY', 'TANGO']

        if run_amylmuts:
            methods.append('AMYLMUTS')

        output = {}
        timeCounts = 0

        for met in methods:
            # TODO: utilize saved file as cached resultr

            # First check if there is an existing results file
            existing_results = glob.glob('*_{}.txt'.format(met))
            if existing_results:
                results_file = existing_results[0]
            else:
                values = {'seq_data': seq, 'method': met}
                data = urlencode(values)
                data = data.encode('ASCII')
                url_input = "http://aias.biol.uoa.gr/cgi-bin/AMYLPRED2/amylpred2.pl"
                response = opener.open(url_input, data)
                result = str(response.read())
                ind = str.find(result, 'Job ID')
                result2 = result[ind:ind + 50]
                ind1 = str.find(result2, ':')
                ind2 = str.find(result2, '<BR>')
                job_id = result2[ind1 + 2:ind2]

                # Waiting for the calculation to complete
                url_result = 'http://aias.biol.uoa.gr/AMYLPRED2/tmp/' + job_id + '.txt'
                print(url_result)
                print("Waiting for %s results" % met, end='.')
                while True:
                    result = urlopen(url_result).read()
                    if not result:
                        time.sleep(1)
                        timeCounts += 1
                        print('.', end='')
                    else:
                        response = requests.get(url_result)
                        break
                results_file = "{}_{}.txt".format(url_result.split('/')[-1].strip('.txt'), met)
                with open(results_file, "wb") as handle:
                    for data in response.iter_content():
                        handle.write(data)

            print("")
            method, hits = self.get_method_hits(results_file, met)
            # if method.lower() == met.lower():
            output[met] = hits
            # elif method == 'Beta-strand contiguity' and met == 'BSC':
            # output[met]=hits
            # elif method == 'Hexapeptide Conf. Energy' and met == 'CONFENERGY':
        print("Time Spent: %d seconds" % timeCounts)
        return output

    def get_method_hits(self, results_file, met):
        result = str(open(results_file).read())
        ind_s = str.find(result, 'HITS')
        ind_e = str.find(result, '**NOTE')
        tmp = result[ind_s + 10:ind_e].strip(" ")
        hits_resid = []
        method = None
        if ":" in tmp:
            method = tmp.split(":")[0]
            hits = tmp.split(":")[1]
            if "-" in hits:
                for ele in hits.split(","):
                    ele = ele.replace('\\r\\n\\r\\n', '')
                    res_s = ele.split("-")[0]
                    res_e = ele.split("-")[1]
                    for i in range(int(res_s), int(res_e) + 1):
                        hits_resid.append(i)
        if method:
            return method, hits_resid
        else:
            return met, hits_resid

    def get_aggregation_index(self, output, cutoff_v, cutoff_n):
        # all index for all aggregation potential > cutoff_v
        IDX = [i for i in range(0, len(output)) if output[i] >= cutoff_v]
        # hold number of continueous aa with aggregation potential > cutoff_v
        # for each continueous stretch
        Agg_index = []
        IDX_for_summation = []
        if IDX:
            # start index for each stretch
            start_IDX = [IDX[0]]
            i = 0
            counter = 0
            while i < len(IDX) - 1:
                if IDX[i] + 1 == IDX[i + 1]:
                    counter += 1
                else:
                    Agg_index.append(counter + 1)
                    start_IDX.append(IDX[i + 1])
                    counter = 0
                i += 1
                if i == len(IDX) - 1: Agg_index.append(counter + 1)
            IDX2 = [i for i in range(0, len(Agg_index)) if Agg_index[i] >= cutoff_n]
            for i in IDX2:
                IDX_for_summation.extend(range(start_IDX[i], Agg_index[i] + start_IDX[i]))
            sum = 0
            for i in IDX_for_summation:
                sum += output[i]
            sum_all = 0
        else:
            IDX2 = []
            sum = 0
        if len(IDX_for_summation) == 0:
            return len(IDX2), 0
        else:
            return len(IDX2), sum / float(len(IDX_for_summation))

    def write_propensity(self, N, output, cutoff_v=5, cutoff_n=5):
        Hits_res_count = {}
        for met in output.keys():
            for res in output[met]:
                if res in Hits_res_count.keys():
                    Hits_res_count[res] += 1
                else:
                    Hits_res_count[res] = 1

        AGG_Propensity = []
        output2 = []
        for i in range(1, int(N) + 1):
            if i in Hits_res_count.keys():
                output2.append(Hits_res_count[i])
            else:
                output2.append(0.0)
        agg_index, agg_conf = self.get_aggregation_index(output2, cutoff_v, cutoff_n)
        # These are parameters that you can explore, "5,5" means
        # 5 consecutive residues agreed among at least 5 methods
        # is considered contributing 1 to the aggregation propensity
        return agg_index, agg_conf


# TODO: clean up execution script to run on FASTA file(s)
# if __name__ == '__main__':
#     from ssbio import utils
#     date = utils.Date()
#     # load inputs from command line
#
#     p = argparse.ArgumentParser(description='Run AMYLPRED2 on a FASTA file or a folder of FASTA files.')
#     p.add_argument('email', help='http://aias.biol.uoa.gr/AMYLPRED2/login.php Email')
#     p.add_argument('password', help='Password')
#     p.add_argument('infile', help='FASTA file or directory of FASTA files.')
#     args = p.parse_args()
#
#     curr_dir = os.getcwd()
#     # initialize the class with your email and password for the site
#     agg_predictions = agg.AMYLPRED(args.email, args.password)
#
#     prop_dir = 'properties'
#     if not op.exists(prop_dir):
#         os.mkdir(prop_dir)
#
#     agg_prop_dir = op.join(prop_dir,'aggregation_propensity')
#     if not op.exists(agg_prop_dir):
#         os.mkdir(agg_prop_dir)
#
#     # TODO: improve arg parsing for files/dirs
#     # TODO: this was all done in a rush - current dir and infile should be improved
#     if len(args.infile) == 1 and op.isdir(args.infile[0]):
#         os.chdir(args.infile[0])
#         files = glob.glob('*')
#     else:
#         files = args.infile
#
#     results = []
#
#     for file in tqdm(files):
#         if op.isdir(file):
#             continue
#
#         # load the sequence file, also the ID
#         seq_records = fasta.load_fasta_file(file)
#         seq_id = op.splitext(op.basename(file))[0]
#
#         seq_folder = op.join(agg_prop_dir, seq_id)
#         if not op.exists(seq_folder):
#             os.mkdir(seq_folder)
#
#         os.chdir(seq_folder)
#         # TODO: seems useless to return seqrecords to just convert them to strings
#         for seq_record in seq_records:
#             agg_index = agg_predictions.get_aggregation_propensity(str(seq_record.seq))
#             result = {'id':seq_id, 'agg_index':agg_index}
#             results.append(result)
#         os.chdir(curr_dir)
#
#     agg_df = pd.DataFrame(results)
#     agg_df.to_csv(op.join(prop_dir, '{}_aggprop_results.csv'.format(date.short_date)))
#     print('Saved results in properties/aggregation_propensity and summarized in aggprop_results.csv')