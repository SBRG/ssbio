import pandas as pd
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import bioservices
bsup = bioservices.uniprot.UniProt()

from dateutil.parser import parse as dateparse
from tqdm import tqdm
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

# TODO: integrate these two functions into one
def uniprot_reviewed_checker_batch(uniprot_ids):
    """A UniProt ID can be be "unreviewed" or "reviewed". This checks the
    status and returns a dictionary indicating this status

    ids: a string or list of UniProt IDs

    """

    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]

    if not isinstance(uniprot_ids, list):
        return None

    uni_rev_dict = {}
    for uniprot_id in uniprot_ids:
        if not uniprot_valid_id(uniprot_id):
            warnings.warn(
                'Invalid UniProt ID found. Removing from results.')
            # uni_rev_dict[uniprot_id] = None
            uniprot_ids.remove(uniprot_id)

    # splitting query up into managable sizes (200 IDs each)
    Nmax = 200
    N, rest = divmod(len(uniprot_ids), Nmax)

    if rest > 0:
        N += 1
    for i in range(0, N):
        i1 = i * Nmax
        i2 = (i + 1) * Nmax
        if i2 > len(uniprot_ids):
            i2 = len(uniprot_ids)

        query = uniprot_ids[i1:i2]

        query_string = ''
        for x in query:
            query_string += 'id:' + x + '+OR+'
        query_string = query_string.strip('+OR+')

        uni_rev_raw = StringIO(bsup.search(
            query_string, columns='id,reviewed', frmt='tab'))
        uni_rev_df = pd.read_table(uni_rev_raw, sep='\t', index_col=0)
        uni_rev_dict_adder = uni_rev_df.to_dict()['Status']
        uni_rev_bool = {k: v.lower() == 'reviewed' for k,
                        v in uni_rev_dict_adder.items()}
        uni_rev_dict.update(uni_rev_bool)

    return uni_rev_dict





# TODO: 160501 - copied from ssbio_new - please ORGANIZE!
def uniprot_ec(uniprot_id):
    r = requests.post('http://www.uniprot.org/uniprot/?query=%s&columns=ec&format=tab' % uniprot_id)

    ec = r.content.splitlines()[1]

    if len(ec) == 0:
        ec = None

    return ec

def uniprot_sites(uniprot_id):
    r = requests.post('http://www.uniprot.org/uniprot/%s.gff' % uniprot_id)
    gff = sio.StringIO(r.content)

    gff_df = pd.read_table(gff, sep='\t',skiprows=2,header=None)
    gff_df.drop([0,1,5,6,7,9], axis=1, inplace=True)
    gff_df.columns = ['type','seq_start','seq_end','notes']

    return gff_df

class UniProt():
    """UniProt class containing information on a UniProt ID
    """

    def __init__(self, id=None, seq=None, reviewed=None):
        self.id = id
        self.seq = seq
        self.reviewed = self.uniprot_reviewed_checker(self.id)

    def uniprot_valid_id(self, id):
        """Checks if an ID is a valid UniProt ID.
        id: any string
        returns: Boolean of validity
        """
        # regex from: http://www.uniprot.org/help/accession_numbers
        valid_id = re.compile(
            "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
        if valid_id.match(id):
            return True
        else:
            return False

    def uniprot_reviewed_checker(self, uniprot_id):
        """A UniProt ID can be be "unreviewed" or "reviewed". This checks the
        status and returns a boolean indicating this status

        ids: a string UniProt ID

        """
        if not uniprot_id:
            return None

        if not self.uniprot_valid_id(uniprot_id):
            raise Exception(
                'Invalid UniProt ID found. Will not check {}.'.format(uniprot_id))
            return None

        query_string = 'id:' + uniprot_id

        uni_rev_raw = StringIO(bsup.search(
            query_string, columns='id,reviewed', frmt='tab'))
        uni_rev = pd.read_table(uni_rev_raw, sep='\t',
                                index_col=0)['Status'][0]

        if not isinstance(uni_rev, str):
            raise Exception(
                'UniProt ID not found.')
            return None

        if uni_rev.lower() == 'reviewed':
            return True
        else:
            return False

    def parse_uniprot_txt_file(self, cache_txt):
        """
        from: boscoh/uniprot github
        Parses the text of metadata retrieved from uniprot.org.

        Only a few fields have been parsed, but this provides a
        template for the other fields.

        A single description is generated from joining alternative
        descriptions.

        Returns a dictionary with the main UNIPROT ACC as keys.
        """

        tag = None
        uniprot_id = None
        metadata_by_seqid = {}
        for l in cache_txt.splitlines():
            test_tag = l[:5].strip()
            if test_tag and test_tag != tag:
                tag = test_tag
            line = l[5:].strip()
            words = line.split()
            if tag == "ID":
                uniprot_id = words[0]
                is_reviewed = words[1].startswith('Reviewed')
                length = int(words[2])
                metadata_by_seqid[uniprot_id] = {
                    'id': uniprot_id,
                    'is_reviewed': is_reviewed,
                    'length': length,
                    'sequence': '',
                    'accs': [],
                }
                entry = metadata_by_seqid[uniprot_id]
            if tag == "DT":
                # DT   01-OCT-1996, integrated into UniProtKB/Swiss-Prot.
                # DT   17-OCT-2006, sequence version 3.
                # DT   22-JUL-2015, entry version 166.
                comma_split = line.split(',')
                if 'sequence version' in comma_split[1]:
                    # print 'sequence_version', comma_split[0],
                    # dateparse(comma_split[0]).date()
                    entry['sequence_version'] = str(
                        dateparse(comma_split[0]).date())
                elif 'entry version' in comma_split[1]:
                    # print 'entry_version', comma_split[0],
                    # dateparse(comma_split[0]).date()
                    entry['entry_version'] = str(
                        dateparse(comma_split[0]).date())
            if tag == "SQ":
                if words[0] != "SEQUENCE":
                    entry['sequence'] += ''.join(words)
            if tag == "AC":
                accs = [w.replace(";", "") for w in words]
                entry['accs'].extend(accs)
            if tag == "DR":
                if 'PDB' in words[0]:
                    if 'pdb' not in entry:
                        entry['pdb'] = words[1][:-1]
                    if 'pdbs' not in entry:
                        entry['pdbs'] = []
                    entry['pdbs'].append(words[1][:-1])
                if 'RefSeq' in words[0]:
                    if 'refseq' not in entry:
                        entry['refseq'] = []
                    ids = [w[:-1] for w in words[1:]]
                    entry['refseq'].extend(ids)
                if 'KEGG' in words[0]:
                    if 'kegg' not in entry:
                        entry['kegg'] = []
                    ids = [w[:-1] for w in words[1:]]
                    ids = filter(lambda w: len(w) > 1, ids)
                    entry['kegg'].extend(ids)
                if 'GO' in words[0]:
                    if 'go' not in entry:
                        entry['go'] = []
                    entry['go'].append(' '.join(words[1:]))
                if 'Pfam' in words[0]:
                    if 'pfam' not in entry:
                        entry['pfam'] = []
                    entry['pfam'].append(words[1][:-1])
            if tag == "GN":
                if 'gene' not in entry and len(words) > 0:
                    pieces = words[0].split("=")
                    if len(pieces) > 1 and 'name' in pieces[0].lower():
                        entry['gene'] = pieces[1].replace(
                            ';', '').replace(',', '')
            if tag == "OS":
                # OS   Homo sapiens (Human).
                if 'organism' not in entry:
                    entry['organism'] = ""
                entry['organism'] += line
            if tag == "DE":
                if 'descriptions' not in entry:
                    entry['descriptions'] = []
                entry['descriptions'].append(line)
            if tag == "CC":
                if 'comment' not in entry:
                    entry['comment'] = ''
                entry['comment'] += line + '\n'

        for entry in metadata_by_seqid.values():
            descriptions = entry['descriptions']
            for i in reversed(range(len(descriptions))):
                description = descriptions[i]
                if 'Short' in description or 'Full' in description or 'EC=' in description:
                    if 'Short' in description or 'Full' in description:
                        j = description.find('=')
                        descriptions[i] = description[j + 1:].replace(';', '')
                        if 'description' not in entry:
                            entry['description'] = []
                        entry['description'].append(descriptions[i])
                    if 'EC=' in description:
                        j = description.find('=')
                        descriptions[i] = description[j + 1:].replace(';', '')
                        if '{' in descriptions[i]:
                            descriptions[i] = descriptions[i].split(' {')[0]
                        if 'ec' not in entry:
                            entry['ec'] = []
                        entry['ec'].append(descriptions[i])
                else:
                    del descriptions[i]

        return metadata_by_seqid

    def uniprot_metadata(self, uniprot_ids):
        '''
        Input: UniProt ID or IDs
        Output: dictionary of metadata associated with the UniProt IDs
        '''
        counter = 1

        if isinstance(uniprot_ids, str):
            uniprot_ids = [uniprot_ids]
            single = True
        else:
            single = False

        uniprot_metadata_raw = bsup.retrieve(uniprot_ids, frmt='txt')
        uniprot_metadata_final = {}

        for uniprot_id in tqdm(uniprot_ids):

            uniprot_metadata_dict = {}

            if single:
                metadata = self.parse_uniprot_txt_file(uniprot_metadata_raw)
            else:
                metadata = self.parse_uniprot_txt_file(
                    uniprot_metadata_raw[uniprot_ids.index(uniprot_id)])

            metadata_key = list(metadata.keys())[0]

            uniprot_metadata_dict['u_uniprot_acc'] = uniprot_id
            uniprot_metadata_dict['u_seq'] = metadata[metadata_key]['sequence']
            uniprot_metadata_dict['u_seq_len'] = len(
                str(metadata[metadata_key]['sequence']))
            uniprot_metadata_dict['u_reviewed'] = metadata[
                metadata_key]['is_reviewed']
            uniprot_metadata_dict['u_seq_version'] = metadata[
                metadata_key]['sequence_version']
            uniprot_metadata_dict['u_entry_version'] = metadata[
                metadata_key]['entry_version']
            if 'gene' in metadata[metadata_key]:
                uniprot_metadata_dict['u_gene_name'] = metadata[
                    metadata_key]['gene']
            if 'description' in metadata[metadata_key]:
                uniprot_metadata_dict['u_description'] = metadata[
                    metadata_key]['description']
            if 'refseq' in metadata[metadata_key]:
                uniprot_metadata_dict['u_refseq'] = metadata[
                    metadata_key]['refseq']
            if 'kegg' in metadata[metadata_key]:
                uniprot_metadata_dict['u_kegg_id'] = metadata[
                    metadata_key]['kegg']
            if 'go' in metadata[metadata_key]:
                uniprot_metadata_dict['u_go'] = metadata[metadata_key]['go']
            if 'ec' in metadata[metadata_key]:
                uniprot_metadata_dict['u_ec_number'] = metadata[
                    metadata_key]['ec']
            if 'pfam' in metadata[metadata_key]:
                uniprot_metadata_dict['u_pfam'] = metadata[
                    metadata_key]['pfam']

            uniprot_metadata_final[uniprot_id] = uniprot_metadata_dict
        return uniprot_metadata_final

if __name__ == '__main__':
    u = UniProt('P00011')
    print(u.seq)
    print(u.reviewed)
    print(uniprot_reviewed_checker_batch(['P00011', 'P29233']))
    print(u.uniprot_metadata('P29233'))
