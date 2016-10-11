# TODO: 160501 - copied from ssbio_new - please ORGANIZE!

import gzip
import io
from collections import defaultdict
from io import BytesIO

import cachetools
import pandas as pd
import requests
import xmltodict

from ssbio import utils

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

SEVEN_DAYS = 60 * 60 * 24 * 7

from collections import defaultdict
from Bio.PDB import PDBList
from ssbio.structure.pdbioext import PDBIOExt
from ssbio.structure.cleanpdb import CleanPDB

def download_and_load_pdb(pdb_id, output_dir, file_type='pdb'):
    """Download a PDB ID to a specified output directory, and of a specified file format (currently only .pdb supported)

    Args:
        pdb_id: valid PDB ID
        output_dir: path to output directory
        file_type: pdb, xml, mmcif - TODO: add support

    Returns: PDBIOExt object

    """
    pdbl = PDBList()
    # TODO: currently saves as .ent file
    file_loc = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir)
    return file_loc

def download_and_clean_pdb(pdb_id, chain_id, output_dir, file_type='pdb', output_name=''):
    """Download the original PDB ID as well as a "cleaned" version of it including only the chain of interest.

    Args:
        pdb_id:
        chain_id:
        output_dir:
        file_type:

    Returns:

    """
    my_pdb_file = download_and_load_pdb(pdb_id, output_dir=output_dir)
    my_pdb = PDBIOExt(my_pdb_file)
    my_model = my_pdb.first_model

    # clean pdb and save a file with only these chains of interest
    # suffix appended with chains (1fat_A_B_cleaned.pdb)
    print('Cleaning PDB structure {} and keeping chain {}...'.format(pdb_id, chain_id))
    my_cleaner = CleanPDB(keep_chains=[chain_id])

    if output_name:
        my_new_pdb_id = output_name
    else:
        my_new_pdb_id = '{}_{}_cleaned'.format(pdb_id, chain_id)

    my_clean_pdb = my_pdb.write_pdb(custom_name=my_new_pdb_id, custom_ext='.pdb', out_suffix='',
                                    out_dir=output_dir, custom_selection=my_cleaner)

    return my_clean_pdb


from lxml import etree
from io import StringIO, BytesIO


def map_uniprot_resnum_to_pdb(uniprot_resnum, chain_id, sifts_file):
    """Map a UniProt residue number to its corresponding PDB residue number.

    This function requires that the SIFTS file be downloaded,
    and also a chain ID (as different chains may have different mappings.

    Args:
        uniprot_resnum: integer of the residue number you'd like to map
        chain_id: string of the PDB chain to map to
        sifts_file: path to the SIFTS XML file

    Returns:
        (mapped_resnum, is_observed) tuple - is_observed indicates if the 3D structure actually shows the residue

    """
    # load the xml with lxml
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(sifts_file, parser)
    root = tree.getroot()

    my_pdb_resnum = None

    # TODO: "Engineered_Mutation is also a possible annotation, need to figure out what to do with that
    my_pdb_annotation = False

    # find the right chain (entities in the xml doc)
    for chain in root.findall('.//{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'):
        if chain.attrib['entityId'] == chain_id:
            # need to find the "crossRefDb" tag that has the attributes dbSource="UniProt" and  dbResNum="your_resnum_here"
            # then match it to the crossRefDb dbResNum that has the attribute dbSource="PDBresnum"

            # check if uniprot + resnum even exists in the sifts file (it won't if the pdb doesn't contain the residue)
            my_uniprot_residue = chain.findall(
                './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]' % uniprot_resnum)
            if len(my_uniprot_residue) == 1:
                # get crossRefDb dbSource="PDB"
                parent = my_uniprot_residue[0].getparent()
                my_pdb_residue = parent.findall(
                    './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="PDB"]')
                my_pdb_resnum = int(my_pdb_residue[0].attrib['dbResNum'])

                # get <residueDetail dbSource="PDBe" property="Annotation">
                # will be Not_Observed if it is not seen in the PDB
                my_pdb_annotation = parent.findall(
                    './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail[@dbSource="PDBe"][@property="Annotation"]')
                if len(my_pdb_annotation) == 1:
                    my_pdb_annotation = my_pdb_annotation[0].text
                    if my_pdb_annotation == 'Not_Observed':
                        my_pdb_annotation = False
                else:
                    my_pdb_annotation = True
            else:
                return (None, False)

    return (my_pdb_resnum, my_pdb_annotation)


@cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
def _theoretical_pdbs():
    """Get the list of theoretical PDBs directly from the wwPDB.

    Caches this list for up to seven days for quick access.

    Returns: list of theoretical PDBs

    """
    theoretical_pdbs_link = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/models/index/titles.idx"
    response = urllib2.urlopen(theoretical_pdbs_link)
    theoretical_pdbs_raw = BytesIO(response.read())
    theoretical_pdbs_df = pd.read_table(theoretical_pdbs_raw, sep='\t', header=None)
    list_of_theoretical_pdbs = [x.lower() for x in theoretical_pdbs_df[0].tolist()]

    return list_of_theoretical_pdbs


def is_theoretical_pdb(pdb_id):
    """Check if a PDB ID is theoretical or not

    Args:
        pdb_id (str): PDB ID

    Returns:
        True if the PDB is a theoretical model
        False if not

    """
    pdb_id = pdb_id.lower()
    return pdb_id in _theoretical_pdbs()


@cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
def _obsolete_pdb_mapping():
    """Get the mapping of obsolete PDBs directly from the wwPDB.

    Caches this list for up to seven days for quick access.

    Returns: a dictionary mapping obsolete PDBs to a list of their replacement(kegg).
    If there is no replacement, it is marked as "obsolete"
    Example:
        {
        '1HHB': ['2HHB', '3HHB', '4HHB'],
        '1VSW': ['4V7G'],
        '1FQH': 'obsolete'
        }

    """
    pdb_obsolete_mapping = {}

    req = urllib2.Request('ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat')
    response = urllib2.urlopen(req)
    output = response.read().decode('utf-8').strip().split('\n')
    for line in output:
        entry = line.split()
        if entry[0] == 'LIST':
            continue
        if len(entry) == 3:
            pdb_obsolete_mapping[entry[2].lower()] = 'obsolete'
        if len(entry) >= 4:
            new = [y.lower() for y in entry[3:]]
            pdb_obsolete_mapping[entry[2].lower()] = new

    return pdb_obsolete_mapping


def is_obsolete_pdb(pdb_id):
    """Check if a PDB ID is obsolete or not

    Args:
        pdb_id (str): PDB ID

    Returns:
        True if the PDB is a obsolete model
        False if not

    """
    pdb_id = pdb_id.lower()
    return pdb_id in _obsolete_pdb_mapping().keys()


def get_replacement_pdbs(pdb_id):
    """Get the list of PDB IDs superseding an obsolete PDB ID

    Args:
        pdb_id (str): PDB ID

    Returns: a list of replacement PDB IDs,
        None if the PDB is not obsolete or there is no replacement

    """
    pdb_id = pdb_id.lower()
    if not is_obsolete_pdb(pdb_id):
        return None
    if _obsolete_pdb_mapping()[pdb_id] == 'obsolete':
        return None
    else:
        return _obsolete_pdb_mapping()[pdb_id]


def pdb_current_checker(pdb_ids):
    """Batch check if a list of PDBs is current, theoretical, obsolete. If obsolete, provide new mapping.
    # status can be "theoretical", "current", or a PDB ID giving the non-obsolete entry
    """
    pdb_ids = utils.force_list(pdb_ids)
    pdb_ids = [x.lower() for x in pdb_ids]

    pdb_status = {}

    theoretical = list(set(_theoretical_pdbs()).intersection(pdb_ids))
    obsolete = list(set(_obsolete_pdb_mapping().keys()).intersection(pdb_ids))
    current = list(set(pdb_ids).difference(theoretical, obsolete))

    for t in theoretical:
        pdb_status[t] = 'theoretical'
    for o in obsolete:
        pdb_status[o] = get_replacement_pdbs(o)
    for c in current:
        pdb_status[c] = 'current'

    return pdb_status


@cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
def _sifts_mapping():
    baseURL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"
    filename = "pdb_chain_uniprot.csv.gz"

    response = urllib2.urlopen(baseURL + filename)
    compressed_file = io.BytesIO(response.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file)

    SIFTS = pd.read_csv(decompressed_file, skiprows=1, index_col=['PDB', 'CHAIN'])
    # clean up some columns (some numbers have alternate IDs in them (ie. 15A)
    SIFTS['PDB_BEG_INT'] = SIFTS['PDB_BEG'].replace(to_replace=r'[^\d-]+', value='', regex=True)
    SIFTS['PDB_BEG_INT'] = SIFTS['PDB_BEG_INT'].astype(int)
    SIFTS['SP_BEG_INT'] = SIFTS['SP_BEG'].replace(to_replace=r'[^\d-]+', value='', regex=True)
    SIFTS['SP_BEG_INT'] = SIFTS['SP_BEG_INT'].astype(int)
    SIFTS['OFFSET'] = SIFTS['SP_BEG_INT'] - SIFTS['PDB_BEG_INT']

    return SIFTS.sort_index()


def SIFTS():
    return _sifts_mapping()


def sifts_pdb_chain_to_uniprot(pdb, chain):
    return _sifts_mapping().ix[(pdb.lower(), chain.upper())]['SP_PRIMARY'].unique().tolist()


@cachetools.func.ttl_cache(maxsize=128)
def blast_pdb(seq, evalue=0.001, link=False):
    """
    Returns a DataFrame of BLAST hits to available structures in the PDB

    :param string seq: your sequence, in string format
    :param float evalue: cutoff for the E-value - filters for significant hits. 0.001 is liberal, 0.0001 is stringent
    :return Pandas DataFrame:
    """

    if len(seq) < 12:
        print('ERROR: Sequence must be at least 12 residues long')
        return None
    if link:
        print('Click here for the PDB results page: http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={}&eCutOff={}&maskLowComplexity=yes&matrix=BLOSUM62&outputFormat=HTML'.format(seq, evalue))
    req = urllib2.Request(
        'http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={}&eCutOff={}&maskLowComplexity=yes&matrix=BLOSUM62&outputFormat=XML'.format(
            seq, evalue))
    response = urllib2.urlopen(req)

    if response.getcode() == 200:
        info_dict = defaultdict(dict)
        len_orig = len(seq)

        raw = BytesIO(response.read())
        doc = xmltodict.parse(raw)
        if 'Iteration_hits' not in doc['BlastOutput']['BlastOutput_iterations']['Iteration'].keys():
            # print('ERROR: No hits found')
            return pd.DataFrame()

        else:
            hit_list = utils.force_list(
                doc['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit'])
            for i, hit in enumerate(hit_list):
                hit_pdb = hit['Hit_def'].split('|')[0].split(':')[0]
                hit_pdb_chain = hit['Hit_def'].split('|')[0].split(':')[2]

                hit_dict = hit['Hit_hsps']['Hsp']
                if isinstance(hit_dict, list):
                    hit_dict = hit_dict[0]

                info_dict[i]['hit_pdb'] = hit_pdb
                info_dict[i]['hit_pdb_chain'] = hit_pdb_chain
                info_dict[i]['hit_seq'] = hit_dict['Hsp_hseq']
                info_dict[i]['hit_num_ident'] = hit_dict['Hsp_identity']
                info_dict[i]['hit_percent_ident'] = float(hit_dict['Hsp_identity']) / float(len_orig)
                info_dict[i]['hit_evalue'] = hit_dict['Hsp_evalue']
                info_dict[i]['hit_score'] = hit_dict['Hsp_score']

                # TODO: what about sequence similarity?

            info_df = pd.DataFrame.from_dict(info_dict, orient='index')
            reorg = ['hit_pdb', 'hit_pdb_chain', 'hit_evalue', 'hit_score', 'hit_num_ident', 'hit_percent_ident']
            # print("Number of BLAST hits: {}".format(len(info_df)))
            return info_df[reorg]

    else:
        # print('ERROR: BLAST request timed out')
        return pd.DataFrame()


def top_pdb_blast_hit(seq, evalue=0.001):
    blast_df = blast_pdb(seq, evalue)
    if blast_df.empty:
        return None
    return blast_df.loc[0].to_dict()

@cachetools.func.ttl_cache(maxsize=500)
def best_structures(uniprot_id):
    """Utilize the PDBe REST API to return the rank-ordered list of PDB "best structures"

    The URL to access the results is formatted like so:
    https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/:accession

    Args:
        uniprot_id (str): UniProt Accession ID

    Returns:
        Rank-ordered list of dictionaries representing chain-specific PDB entries

    """

    response = requests.get('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{}'.format(uniprot_id),
                            data={'key': 'value'})
    if response.status_code == 404:
        return {}
    return dict(response.json())[uniprot_id]

#
#
# class Structure():
#     """
#     Class for defining Structure within an SBML model.
#     A structure is mainly:
#     1. a pointer to the file of the protein structure cartesian coordinates
#     2. a representation of a gene or genes within a GEM
#     """
#
#     def __init__(self, id=None, filename=None, seqs=None, current=None):
#         """
#         An object for defining protein structures and how they represent a
#         gene or genes within a SBML model.
#         """
#         self.id = id
#         self.filename = filename
#         self.current = self.pdb_current_checker(self.id)
#
#         # all genes which are transcribed to make up the protein structure
#         self.genes = set()
#
#         # all metabolites/ligands within the structure
#         self.chemicals = set()
#
#         # all reactions which utilize this protein
#         self.reactions = set()
#
#         # dictionary of chain IDs and their sequences as strings
#         self.sequences = seqs
#

#
#
#     def pdb_current_checker(self, pdb_ids):
#         """
#         # status can be "theoretical", "current", or a PDB ID giving the non-obsolete entry
#         """
#         if not pdb_ids:
#             return None
#         if isinstance(pdb_ids, str):
#             pdb_ids = [pdb_ids]
#
#         pdb_status = {}
#
#         theoretical = list(set(PDB_THEORETICAL).intersection(pdb_ids))
#         obsolete = list(set(PDB_OBSOLETE_MAPPING.keys()).intersection(pdb_ids))
#         current = list(set(pdb_ids).difference(theoretical, obsolete))
#
#         for t in theoretical:
#             pdb_status[t] = 'theoretical'
#         for o in obsolete:
#             pdb_status[o] = PDB_OBSOLETE_MAPPING[o]
#         for c in current:
#             pdb_status[c] = 'current'
#
#         return pdb_status
#
#     def pdb_metadata_and_download(pdbids):
#         final_dict = {}
#         counter = 1
#
#         if isinstance(pdbids, str):
#             pdbids = [pdbids]
#
#         for pdbid in pdbids:
#
#             clear_output(wait=True)
#             print('***PROGRESS: PARSED & DOWNLOADED %d/%d PDB IDS***\n' % (counter, len(pdbids)))
#             counter += 1
#             sys.stdout.flush()
#
#             try:
#                 header = pr.parsePDBHeader(pdbid)
#
#             # if PDB isn't able to be downloaded
#             except IOError:
#                 try:
#                     header = pr.parsePDBHeader(pr.fetchPDBviaFTP(pdbid,format='cif'))
#                 except IOError:
#                     continue
#
#             appender = defaultdict(list)
#             for prop,val in header.iteritems():
#                 if isinstance(val, pr.proteins.header.Chemical):
#                     appender['p_chemicals'].append(val.resname)
#                 elif isinstance(val, pr.proteins.header.Polymer):
#                     appender['p_chains'].append(val.chid)
#                     if val.ec:
#                         appender['p_ec_numbers'].append((val.chid,val.ec))
#                 elif prop == 'reference':
#                     if 'doi' in val:
#                         appender['p_doi'] = val['doi']
#                     if 'pmid' in val:
#                         appender['p_pmid'] = val['pmid']
#                 elif prop in ['resolution','space_group','experiment','deposition_date']:
#                     appender['p_' + prop] = val
#
#             tmp = {}
#             for chain in appender['p_chains']:
#                 try:
#                     uniprot_mapping = sifts_pdb_chain_to_uniprot(pdbid, chain)
#                 except KeyError:
#                     uniprot_mapping = ['PDB-'+chain]
#                 tmp[chain] = uniprot_mapping
#             appender['p_chain_uniprot_map'] = tmp
#
#             final_dict[pdbid] = dict(appender)
#
#         return de_unicodeify(final_dict)
#
# if __name__ == '__main__':
#     p = Structure()
#     p.pdb_current_checker('2223')
