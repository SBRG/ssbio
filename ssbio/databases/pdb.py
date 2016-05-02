# TODO: 160501 - copied from ssbio_new - please ORGANIZE!


from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Polypeptide

import pandas as pd
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from io import BytesIO

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2


## OBSOLETE and THEORETICAL STUFF
theoretical_pdbs_link = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/models/index/titles.idx"
response = urllib2.urlopen(theoretical_pdbs_link)
theoretical_pdbs_raw = BytesIO(response.read())
theoretical_pdbs_df = pd.read_table(
    theoretical_pdbs_raw, sep='\t', header=None)

PDB_THEORETICAL = theoretical_pdbs_df[0].tolist()

PDB_OBSOLETE_MAPPING = {}

req = urllib2.Request('ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat')
response = urllib2.urlopen(req)
for line in response:
    entry = line.split()
    if entry[0] == 'LIST':
        continue
    if len(entry) == 3:
        PDB_OBSOLETE_MAPPING[entry[2].upper()] = 'obsolete'
    if len(entry) >= 4:
        new = [y.upper() for y in entry[3:]]
        PDB_OBSOLETE_MAPPING[entry[2].upper()] = new
############


## SIFTS STUFF
now = time.time()
twodays_ago = now - 60*60*24*2
baseURL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"
filename = "pdb_chain_uniprot.csv.gz"

final_filename = tempfile.gettempdir() + '/' + filename.split('.gz')[0]

if os.path.exists(final_filename):
    final_filename_creation_time = os.path.getmtime(final_filename)

    if final_filename_creation_time < twodays_ago:
        response = urllib2.urlopen(baseURL + filename)
        compressedFile = sio.StringIO(response.read())
        decompressedFile = gzip.GzipFile(fileobj=compressedFile)

        with open(final_filename, 'w') as outfile:
            outfile.write(decompressedFile.read())

        SIFTS = pd.read_csv(final_filename, skiprows=1, index_col=['PDB','CHAIN'])

    else:
        SIFTS = pd.read_csv(final_filename, skiprows=1, index_col=['PDB','CHAIN'])
else:
    response = urllib2.urlopen(baseURL + filename)
    compressedFile = sio.StringIO(response.read())
    decompressedFile = gzip.GzipFile(fileobj=compressedFile)

    with open(final_filename, 'w') as outfile:
        outfile.write(decompressedFile.read())

    SIFTS = pd.read_csv(final_filename, skiprows=1, index_col=['PDB','CHAIN'])
SIFTS['PDB_BEG_INT'] = SIFTS['PDB_BEG'].replace(to_replace=r'[^\d-]+', value='', regex=True)
SIFTS['SP_BEG_INT'] = SIFTS['SP_BEG'].replace(to_replace=r'[^\d-]+', value='', regex=True)
SIFTS['OFFSET'] = SIFTS['SP_BEG_INT'].astype(int) - SIFTS['PDB_BEG_INT'].astype(int)

def sifts_pdb_chain_to_uniprot(pdb,chain):
    return SIFTS.ix[(pdb.lower(),chain.upper())]['SP_PRIMARY'].unique().tolist()
####



class Structure():
    """
    Class for defining Structure within an SBML model.
    A structure is mainly:
    1. a pointer to the file of the protein structure cartesian coordinates
    2. a representation of a gene or genes within a GEM
    """

    def __init__(self, id=None, filename=None, seqs=None, current=None):
        """
        An object for defining protein structures and how they represent a
        gene or genes within a SBML model.
        """
        self.id = id
        self.filename = filename
        self.current = self.pdb_current_checker(self.id)

        # all genes which are transcribed to make up the protein structure
        self.genes = set()

        # all metabolites/ligands within the structure
        self.chemicals = set()

        # all reactions which utilize this protein
        self.reactions = set()

        # dictionary of chain IDs and their sequences as strings
        self.sequences = seqs

    def get_pdb_seq(structure):
        '''
        Takes in a Biopython structure object and returns a list of the
        structure's sequences
        :param structure: Biopython structure object
        :return: Dictionary of sequence strings with chain IDs as the key
        '''

        structure_seqs = {}

        # loop over each chain of the PDB
        for chain in structure[0]:

            chain_it = iter(chain)

            chain_seq = ''
            tracker = 0

            # loop over the residues
            for res in chain.get_residues():
                # NOTE: you can get the residue number too
                res_num = res.id[1]

                # double check if the residue name is a standard residue
                # if it is not a standard residue (ie. selenomethionine),
                # it will be filled in with an X on the next iteration)
                if Polypeptide.is_aa(res, standard=True):
                    full_id = res.get_full_id()
                    end_tracker = full_id[3][1]
                    i_code = full_id[3][2]
                    aa = Polypeptide.three_to_one(res.get_resname())

                    # tracker to fill in X's
                    if end_tracker != (tracker + 1):   # and first == False:
                        if i_code != ' ':
                            chain_seq += aa
                            tracker = end_tracker + 1
                            continue
                        else:
                            chain_seq += 'X' * (end_tracker - tracker - 1)

                    chain_seq += aa
                    tracker = end_tracker

                else:
                    continue

            structure_seqs[chain.get_id()] = chain_seq

        return structure_seqs


    def pdb_current_checker(self, pdb_ids):
        """
        # status can be "theoretical", "current", or a PDB ID giving the non-obsolete entry
        """
        if not pdb_ids:
            return None
        if isinstance(pdb_ids, str):
            pdb_ids = [pdb_ids]

        pdb_status = {}

        theoretical = list(set(PDB_THEORETICAL).intersection(pdb_ids))
        obsolete = list(set(PDB_OBSOLETE_MAPPING.keys()).intersection(pdb_ids))
        current = list(set(pdb_ids).difference(theoretical, obsolete))

        for t in theoretical:
            pdb_status[t] = 'theoretical'
        for o in obsolete:
            pdb_status[o] = PDB_OBSOLETE_MAPPING[o]
        for c in current:
            pdb_status[c] = 'current'

        return pdb_status

    def pdb_metadata_and_download(pdbids):
        final_dict = {}
        counter = 1

        if isinstance(pdbids, str):
            pdbids = [pdbids]

        for pdbid in pdbids:

            clear_output(wait=True)
            print('***PROGRESS: PARSED & DOWNLOADED %d/%d PDB IDS***\n' % (counter, len(pdbids)))
            counter += 1
            sys.stdout.flush()

            try:
                header = pr.parsePDBHeader(pdbid)

            # if PDB isn't able to be downloaded
            except IOError:
                try:
                    header = pr.parsePDBHeader(pr.fetchPDBviaFTP(pdbid,format='cif'))
                except IOError:
                    continue

            appender = defaultdict(list)
            for prop,val in header.iteritems():
                if isinstance(val, pr.proteins.header.Chemical):
                    appender['p_chemicals'].append(val.resname)
                elif isinstance(val, pr.proteins.header.Polymer):
                    appender['p_chains'].append(val.chid)
                    if val.ec:
                        appender['p_ec_numbers'].append((val.chid,val.ec))
                elif prop == 'reference':
                    if 'doi' in val:
                        appender['p_doi'] = val['doi']
                    if 'pmid' in val:
                        appender['p_pmid'] = val['pmid']
                elif prop in ['resolution','space_group','experiment','deposition_date']:
                    appender['p_' + prop] = val

            tmp = {}
            for chain in appender['p_chains']:
                try:
                    uniprot_mapping = sifts_pdb_chain_to_uniprot(pdbid, chain)
                except KeyError:
                    uniprot_mapping = ['PDB-'+chain]
                tmp[chain] = uniprot_mapping
            appender['p_chain_uniprot_map'] = tmp

            final_dict[pdbid] = dict(appender)

        return de_unicodeify(final_dict)

if __name__ == '__main__':
    p = Structure()
    p.pdb_current_checker('2223')
