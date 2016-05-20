import requests
import urllib2
import urllib
from cookielib import CookieJar
import time

def get_aggregation_propensity(gene,cutoff_v=5,cutoff_n=5):
    """Predict aggregation propensity

    Args:
        gene:     b number
        cutoff_v: the minimal number of methods that agree on 
                  a residue being a aggregation-prone residue 
        cutoff_n: minimal number of consecutive residues to be
                  considered as a 'stretch' of aggregation-prone region

    Returns:
        aggregation propensity

    """
 
    # seq=f(gene) # from ME2.0 or PDB structure
    output = consensus_aggregation(seq,gene)
    agg_index,agg_conf = write_propensity(gene,len(seq),output,cutoff_v=5,cutoff_n=5) 
   
    return agg_index



def consensus_aggregation(seq,gene):
    url = "http://aias.biol.uoa.gr/AMYLPRED2/login.php"
    cj = CookieJar()
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
    formdata = { "email" : "kec003@ucsd.edu", "password": "KoKa123456"}
    data_encoded = urllib.urlencode(formdata)
    response = opener.open(url, data_encoded)

    # To do: usually AMYLMUTS is most time consuming, 
    #        and generate a slightly different result every submission
    #        consider remove this from list to save computational time?
    #        will need redo statistics
    Methods = ['AGGRESCAN','NETCSSP','PAFIG','AMYLMUTS','APD','AMYLPATTERN', \
               'SECSTR','BSC','WALTZ','CONFENERGY','TANGO']

    output = {}
    timeCounts = 0
    for met in Methods:
        values= {'seq_name': gene, 'seq_data': seq, 'method':met}
        data = urllib.urlencode(values)
        #url_input = "http://aias.biol.uoa.gr/AMYLPRED2/input.php"
        url_input = "http://aias.biol.uoa.gr/cgi-bin/AMYLPRED2/amylpred2.pl"
        response = opener.open(url_input,data)
        result = response.read()
        ind = str.find(result,'Job ID')
        result2 = result[ind:ind+50]
        ind1 = str.find(result2,':')
        ind2 = str.find(result2,'<BR>')
        job_id = result2[ind1+2:ind2]

        # waiting for the calculation to complete
        url_result = 'http://aias.biol.uoa.gr/AMYLPRED2/tmp/'+job_id+'.txt'
        print url_result
        print "Waiting for %s results" % met
        while True:
            result = urllib2.urlopen(url_result).read()
            if not result:
                time.sleep(1)
                timeCounts += 1
                print ".",
            else:
                break
        print ""
        method,hits=get_method_hits(url_result,met)
        #if method.lower() == met.lower():
        output[met]=hits
        #elif method == 'Beta-strand contiguity' and met == 'BSC':
            #output[met]=hits
        #elif method == 'Hexapeptide Conf. Energy' and met == 'CONFENERGY':
    print "Time Spent on gene %s: %d seconds" %(gene,timeCounts)
    return output

def get_method_hits(url_result,met):
    result = urllib2.urlopen(url_result).read()
    ind_s = str.find(result,'HITS')
    ind_e = str.find(result,'**NOTE')
    tmp = result[ind_s+10:ind_e].strip(" ")
    hits_resid = []
    method = None
    if ":" in tmp:
        method = tmp.split(":")[0]
        hits = tmp.split(":")[1]
        if "-" in hits:
            for ele in hits.split(","):
                res_s = ele.split("-")[0]
                res_e = ele.split("-")[1]
                for i in range(int(res_s), int(res_e)+1):
                    hits_resid.append(i)
    if method:
        return method,hits_resid
    else:
        return met,hits_resid

def get_aggregation_index(output,cutoff_v,cutoff_n):
    # all index for all aggregation potential > cutoff_v
    IDX=[i for i in range(0,len(output)) if output[i] >= cutoff_v]
    # hold number of continueous aa with aggregation potential > cutoff_v
    # for each continueous stretch
    Agg_index = []
    IDX_for_summation = []
    if IDX:
        # start index for each stretch
        start_IDX=[IDX[0]]
        i = 0
        counter = 0
        while i < len(IDX)-1:
            if IDX[i]+1 == IDX[i+1]:
                counter += 1
            else:
                Agg_index.append(counter+1)
                start_IDX.append(IDX[i+1])
                counter = 0
            i += 1
            if i == len(IDX)-1: Agg_index.append(counter+1)
        IDX2=[i for i in range(0,len(Agg_index)) if Agg_index[i] >= cutoff_n]    
        for i in IDX2:
            IDX_for_summation.extend(range(start_IDX[i],Agg_index[i]+start_IDX[i]))
        sum = 0
        for i in IDX_for_summation:
            sum += output[i]
        sum_all = 0
    else:
        IDX2=[]
        sum = 0
    if len(IDX_for_summation) == 0:
        return len(IDX2),0
    else:
        return len(IDX2),sum/float(len(IDX_for_summation))

def write_propensity(gene,N,output,cutoff_v=5,cutoff_n=5):
    Hits_res_count = {}
    for met in output.keys():
        for res in output[met]:
            if res in Hits_res_count.keys():
                Hits_res_count[res] += 1
            else:
                Hits_res_count[res] = 1

    AGG_Propensity = []
    output2 = []
    for i in range(1,int(N)+1):
        if i in Hits_res_count.keys():
            output2.append(Hits_res_count[i])
        else:
            output2.append(0.0)
    agg_index, agg_conf = get_aggregation_index(output2,cutoff_v,cutoff_n) # these are parameters that you can explore, "5,5" means
                                                                           # 5 consecutive residues agreed among at least 5 methods
                                                                           # is considered contributing 1 to the aggregation propensity
    return agg_index,agg_conf
