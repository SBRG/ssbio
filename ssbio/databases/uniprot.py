import pandas as pd
usimport ssbio.utils
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest

import bioservices
bsup = bioservices.uniprot.UniProt()

from dateutil.parser import parse as dateparse
import warnings
import re
# import cachetools
import requests
import os.path as op

SEVEN_DAYS = 60 * 60 * 24 * 7


def is_valid_uniprot_id(instring):
    """Check if a string is a valid UniProt ID.

    See regex from: http://www.uniprot.org/help/accession_numbers

    Args:
        instring: any string identifier

    Returns: True if the string is a valid UniProt ID

    """
    valid_id = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    if valid_id.match(str(instring)):
        return True
    else:
        return False


# @cachetools.func.ttl_cache(maxsize=800, ttl=SEVEN_DAYS)
def get_fasta(uniprot_id):
    """Get the protein sequence for a UniProt ID as a string.

    Args:
        ident: a UniProt ID

    Returns: string of the protein (amino acid) sequence

    """
    if not is_valid_uniprot_id(uniprot_id):
        raise ValueError("Invalid UniProt ID!")

    # silencing the "Will be moved to Biokit" message
    with ssbio.utils.suppress_stdout():
        return bsup.get_fasta_sequence(uniprot_id)


# @cachetools.func.ttl_cache(maxsize=128, ttl=SEVEN_DAYS)
def uniprot_reviewed_checker(uniprot_id):
    """Check if a single UniProt ID is reviewed or not.

    Args:
        uniprot_id:

    Returns:
        True or False
    """
    if not is_valid_uniprot_id(uniprot_id):
        raise ValueError("Invalid UniProt ID!")

    query_string = 'id:' + uniprot_id

    uni_rev_raw = StringIO(bsup.search(query_string, columns='id,reviewed', frmt='tab'))
    uni_rev_df = pd.read_table(uni_rev_raw, sep='\t', index_col=0)
    uni_rev_df = uni_rev_df.fillna(False)
    uni_rev_df = uni_rev_df[pd.notnull(uni_rev_df.Status)]

    uni_rev_df = uni_rev_df.replace(to_replace="reviewed", value=True)
    uni_rev_df = uni_rev_df.replace(to_replace="unreviewed", value=False)
    uni_rev_dict_adder = uni_rev_df.to_dict()['Status']

    return uni_rev_dict_adder[uniprot_id]


# @cachetools.func.ttl_cache(maxsize=128, ttl=SEVEN_DAYS)
def uniprot_reviewed_checker_batch(uniprot_ids):
    """Batch check if uniprot IDs are reviewed or not

    Args:
        uniprot_ids:

    Returns:
        A dictionary of {UniProtID: Boolean}
    """
    uniprot_ids = ssbio.utils.force_list(uniprot_ids)

    invalid_ids = [i for i in uniprot_ids if not is_valid_uniprot_id(i)]
    uniprot_ids = [i for i in uniprot_ids if is_valid_uniprot_id(i)]

    if invalid_ids:
        warnings.warn("Invalid UniProt IDs {} will be ignored".format(invalid_ids))

    # splitting query up into managable sizes (200 IDs each)
    Nmax = 200
    N, rest = divmod(len(uniprot_ids), Nmax)

    uni_rev_dict = {}

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


def uniprot_ec(uniprot_id):
    r = requests.post('http://www.uniprot.org/uniprot/?query=%s&columns=ec&format=tab' % uniprot_id)

    ec = r.content.decode('utf-8').splitlines()[1]

    if len(ec) == 0:
        ec = None

    return ec


def uniprot_sites(uniprot_id):
    r = requests.post('http://www.uniprot.org/uniprot/%s.gff' % uniprot_id)
    gff = StringIO(r.content.decode('utf-8'))

    try:
        gff_df = pd.read_table(gff, sep='\t', skiprows=2, header=None)
    except ValueError:
        return pd.DataFrame()

    gff_df.drop([0, 1, 5, 6, 7, 9], axis=1, inplace=True)
    gff_df.columns = ['type', 'seq_start', 'seq_end', 'notes']

    return gff_df


def download_uniprot_file(uniprot_id, filetype, outdir='', force_rerun=False):
    """Download any UniProt file for a UniProt ID/ACC

    Args:
        uniprot_id:
        filetype: txt, fasta, xml, rdf, or gff
        outdir:

    Returns:

    """
    url = 'http://www.uniprot.org/uniprot/{}.{}'.format(uniprot_id, filetype)
    outfile = op.join(outdir, '{}.{}'.format(uniprot_id, filetype))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        urlrequest.urlretrieve(url, outfile)

    return outfile


def old_parse_uniprot_txt_file(infile):
    """
    from: boscoh/uniprot github
    Parses the text of metadata retrieved from uniprot.org.

    Only a few fields have been parsed, but this provides a
    template for the other fields.

    A single description is generated from joining alternative
    descriptions.

    Returns a dictionary with the main UNIPROT ACC as keys.
    """
    with open(infile, 'r') as txt:
        cache_txt = txt.read()
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


def parse_uniprot_txt_file(infile):
    uniprot_metadata_dict = {}

    metadata = old_parse_uniprot_txt_file(infile)
    metadata_key = list(metadata.keys())[0]

    # uniprot_metadata_dict['u_seq'] = metadata[metadata_key]['sequence']
    uniprot_metadata_dict['seq_len'] = len(
        str(metadata[metadata_key]['sequence']))
    uniprot_metadata_dict['reviewed'] = metadata[metadata_key]['is_reviewed']
    uniprot_metadata_dict['seq_version'] = metadata[metadata_key]['sequence_version']
    uniprot_metadata_dict['entry_version'] = metadata[metadata_key]['entry_version']
    if 'gene' in metadata[metadata_key]:
        uniprot_metadata_dict['gene_name'] = metadata[metadata_key]['gene']
    if 'description' in metadata[metadata_key]:
        uniprot_metadata_dict['description'] = metadata[metadata_key]['description']
    if 'refseq' in metadata[metadata_key]:
        uniprot_metadata_dict['refseq'] = metadata[metadata_key]['refseq']
    if 'kegg' in metadata[metadata_key]:
        uniprot_metadata_dict['kegg_id'] = metadata[metadata_key]['kegg']
    if 'ec' in metadata[metadata_key]:
        uniprot_metadata_dict['ec_number'] = metadata[metadata_key]['ec']
    if 'pfam' in metadata[metadata_key]:
        uniprot_metadata_dict['pfam'] = metadata[metadata_key]['pfam']
    if 'pdbs' in metadata[metadata_key]:
        uniprot_metadata_dict['pdbs'] = list(set(metadata[metadata_key]['pdbs']))
    return uniprot_metadata_dict


# def uniprot_metadata_batch(uniprot_ids, outdir=None):
#     '''
#     Input: UniProt ID or IDs
#     Output: dictionary of metadata associated with the UniProt IDs
#     '''
#     # TODO: not working as of 2016-05-08 - seems to be a limit on the url requests
#     counter = 1
#
#     # if isinstance(uniprot_ids, str):
#     #     uniprot_ids = [uniprot_ids]
#     #     single = True
#     # else:
#     #     single = False
#     uniprot_ids = utils.force_list(uniprot_ids)
#
#     uniprot_metadata_raw = bs_unip.retrieve(uniprot_ids, frmt='txt')
#     uniprot_metadata_final = {}
#
#     for uniprot_id in tqdm(uniprot_ids):
#
#         uniprot_metadata_dict = {}
#
#         if single:
#             metadata = parse_uniprot_txt_file(uniprot_metadata_raw)
#         else:
#             metadata = parse_uniprot_txt_file(
#                 uniprot_metadata_raw[uniprot_ids.index(uniprot_id)])
#
#         metadata_key = list(metadata.keys())[0]
#
#         uniprot_metadata_dict['u_uniprot_acc'] = uniprot_id
#         uniprot_metadata_dict['u_seq'] = metadata[metadata_key]['sequence']
#         uniprot_metadata_dict['u_seq_len'] = len(
#             str(metadata[metadata_key]['sequence']))
#         uniprot_metadata_dict['u_reviewed'] = metadata[
#             metadata_key]['is_reviewed']
#         uniprot_metadata_dict['u_seq_version'] = metadata[
#             metadata_key]['sequence_version']
#         uniprot_metadata_dict['u_entry_version'] = metadata[
#             metadata_key]['entry_version']
#         if 'gene' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_gene_name'] = metadata[
#                 metadata_key]['gene']
#         if 'description' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_description'] = metadata[
#                 metadata_key]['description']
#         if 'refseq' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_refseq'] = metadata[
#                 metadata_key]['refseq']
#         if 'kegg' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_kegg_id'] = metadata[
#                 metadata_key]['kegg']
#         if 'go' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_go'] = metadata[metadata_key]['go']
#         if 'ec' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_ec_number'] = metadata[
#                 metadata_key]['ec']
#         if 'pfam' in metadata[metadata_key]:
#             uniprot_metadata_dict['u_pfam'] = metadata[
#                 metadata_key]['pfam']
#
#         uniprot_metadata_final[uniprot_id] = uniprot_metadata_dict
#     return uniprot_metadata_final

if __name__ == '__main__':
    pass
