import json
import io
import requests
import pandas as pd
import os.path as op
from io import BytesIO
import cachetools
import gzip

from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from lxml import etree
from ssbio import utils
from ssbio.structure.cleanpdb import CleanPDB
from ssbio.structure.pdbioext import PDBIOExt

import logging
log = logging.getLogger(__name__)

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import gzip


SEVEN_DAYS = 60 * 60 * 24 * 7


def download_structure(pdb_id, file_type, outdir='', outfile='', header=False, force_rerun=False):
    """Download a structure from the RCSB PDB by ID. Specify the file type desired.

    Args:
        pdb_id: PDB ID
        file_type: pdb, pdb.gz, cif, cif.gz, xml.gz
        outdir: optional output directory
        outfile: optional output name

    Returns:

    """
    pdb_id = pdb_id.lower()
    file_type = file_type.lower()
    file_types = ['pdb', 'pdb.gz', 'cif', 'cif.gz', 'xml.gz']
    if file_type not in file_types:
        raise ValueError('Invalid file type, must be either: pdb, pdb.gz, cif, cif.gz, xml.gz')

    if header:
        folder = 'header'
        if outfile:
            outfile = op.join(outdir, outfile)
        else:
            outfile = op.join(outdir, '{}.header.{}'.format(pdb_id, file_type))
    else:
        folder = 'download'
        if outfile:
            outfile = op.join(outdir, outfile)
        else:
            outfile = op.join(outdir, '{}.{}'.format(pdb_id, file_type))

    if not op.exists(outfile) and not force_rerun:
        download_link = 'https://files.rcsb.org/{}/{}.{}'.format(folder, pdb_id, file_type)
        req = requests.get(download_link)

        # Raise error if request fails
        req.raise_for_status()

        with open(outfile, 'w') as f:
            f.write(req.text)

        log.debug('{}: Saved structure file'.format(outfile))
    else:
        log.debug('{}: Structure file already saved'.format(outfile))

    return outfile


def parse_pdb_header(infile):
    """Parse a couple important fields from the mmCIF file format with some manual curation of ligands.

    If you want full access to the mmCIF file just use the MMCIF2Dict class in Biopython.

    Args:
        infile: Path to mmCIF file

    Returns:
        dict: Dictionary of parsed header

    """
    pass


# @cachetools.func.ttl_cache(maxsize=500, ttl=SEVEN_DAYS)
def parse_mmcif_header(infile):
    """Parse a couple important fields from the mmCIF file format with some manual curation of ligands.

    If you want full access to the mmCIF file just use the MMCIF2Dict class in Biopython.

    Args:
        infile: Path to mmCIF file

    Returns:
        dict: Dictionary of parsed header

    """
    newdict = {}
    mmdict = MMCIF2Dict(infile)

    chemical_ids_exclude = ['HOH']
    chemical_types_exclude = ['l-peptide linking','peptide linking']

    if '_exptl.method' in mmdict:
        newdict['experimental_method'] = mmdict['_exptl.method']
    else:
        log.debug('{}: No experimental method field'.format(infile))

    if '_refine.ls_d_res_high' in mmdict:
        newdict['resolution'] = float(mmdict['_refine.ls_d_res_high'])
    else:
        log.debug('{}: No resolution field'.format(infile))

    if '_chem_comp.id' in mmdict:
        chemicals_filtered = utils.filter_list_by_indices(mmdict['_chem_comp.id'],
                                                            utils.not_find(mmdict['_chem_comp.type'],
                                                                           chemical_types_exclude,
                                                                           case_sensitive=False))
        chemicals_fitered = utils.filter_list(chemicals_filtered, chemical_ids_exclude, case_sensitive=True)
        newdict['chemicals'] = chemicals_fitered
    else:
        log.debug('{}: No chemical composition field'.format(infile))

    if '_entity_src_gen.pdbx_gene_src_scientific_name' in mmdict:
        newdict['taxonomy_name'] = mmdict['_entity_src_gen.pdbx_gene_src_scientific_name']
    else:
        log.debug('{}: No organism field'.format(infile))

    return newdict


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


# TODO: check this for python 2
def download_sifts_xml(pdb_id, outdir='', outfile=''):
    """Download the SIFTS file for a PDB ID.

    Args:
        pdb_id:
        outdir:
        outfile:

    Returns:

    """
    baseURL = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    filename = '{}.xml.gz'.format(pdb_id)

    if outfile:
        outfile = op.join(outdir, outfile)
    else:
        outfile = op.join(outdir, filename.split('.')[0] + '.sifts.xml')

    if not op.exists(outfile):
        response = urllib2.urlopen(baseURL + filename)
        with open(outfile, 'wb') as outfile:
            outfile.write(gzip.decompress(response.read()))

    return outfile


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
    # Load the xml with lxml
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(sifts_file, parser)
    root = tree.getroot()

    my_pdb_resnum = None

    # TODO: "Engineered_Mutation is also a possible annotation, need to figure out what to do with that
    my_pdb_annotation = False

    # Find the right chain (entities in the xml doc)
    ent = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'
    for chain in root.findall(ent):
        if chain.attrib['entityId'] == chain_id:
            # Find the "crossRefDb" tag that has the attributes dbSource="UniProt" and  dbResNum="your_resnum_here"
            # Then match it to the crossRefDb dbResNum that has the attribute dbSource="PDBresnum"

            # Check if uniprot + resnum even exists in the sifts file (it won't if the pdb doesn't contain the residue)
            ures = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]' % uniprot_resnum
            my_uniprot_residue = chain.findall(ures)
            if len(my_uniprot_residue) == 1:
                # Get crossRefDb dbSource="PDB"
                parent = my_uniprot_residue[0].getparent()
                pres = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="PDB"]'
                my_pdb_residue = parent.findall(pres)
                my_pdb_resnum = int(my_pdb_residue[0].attrib['dbResNum'])

                # Get <residueDetail dbSource="PDBe" property="Annotation">
                # Will be Not_Observed if it is not seen in the PDB
                anno = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail[@dbSource="PDBe"][@property="Annotation"]'
                my_pdb_annotation = parent.findall(anno)
                if len(my_pdb_annotation) == 1:
                    my_pdb_annotation = my_pdb_annotation[0].text
                    if my_pdb_annotation == 'Not_Observed':
                        my_pdb_annotation = False
                else:
                    my_pdb_annotation = True
            else:
                return None, False

    return my_pdb_resnum, my_pdb_annotation


# @cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
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


# @cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
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

    Returns:
        A list of replacement PDB IDs,
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

    Args:
        pdb_ids (list): List of PDB IDs

    Returns:
        dict: Dictionary of {pdb_id: <status>}, where status can be
            "theoretical", "current", or a list of PDB IDs giving the non-obsolete entry
    """
    pdb_ids = utils.force_lower_list(pdb_ids)

    log.debug('Checking list of PDBs: {}'.format(pdb_ids))

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


def update_pdb_list(pdb_ids):
    """Filter a list of PDBs to remove obsolete and theoretical IDs. Replace obsolete ones with the updated IDs

    Args:
        pdb_ids: List of PDB IDs

    Returns:
        Updated list of PDB IDs

    """
    pdb_ids = utils.force_lower_list(pdb_ids)
    checked = pdb_current_checker(pdb_ids)

    new_pdb_ids = []

    for p in pdb_ids:
        if not checked[p] == 'current':
            if checked[p] == 'theoretical':
                continue
            elif checked[p]:
                new_pdb_ids.extend(p)
        else:
            new_pdb_ids.append(p)

    return list(set(new_pdb_ids))


def best_structures(uniprot_id, outfile='', outdir='', seq_ident_cutoff=0, force_rerun=False):
    """Use the PDBe REST service to query for the best PDB structures for a UniProt ID.

        More information found here: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
        Link used to retrieve results: https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/:accession
        The list of PDB structures mapping to a UniProt accession sorted by coverage of the protein and, if the same, resolution.

        Here is the ranking algorithm described by the PDB paper:
        https://nar.oxfordjournals.org/content/44/D1/D385.full

        "Finally, a single quality indicator is also calculated for each entry by taking the harmonic average
        of all the percentile scores representing model and model-data-fit quality measures and then subtracting
        10 times the numerical value of the resolution (in Angstrom) of the entry to ensure that resolution plays
        a role in characterising the quality of a structure. This single empirical 'quality measure' value is used
        by the PDBe query system to sort results and identify the 'best' structure in a given context. At present,
        entries determined by methods other than X-ray crystallography do not have similar data quality information
        available and are not considered as 'best structures'."

        Args:
            uniprot_id (str): UniProt Accession ID

        Returns:
            list: Rank-ordered list of dictionaries representing chain-specific PDB entries. Keys are:
            pdb_id: the PDB ID which maps to the UniProt ID
            chain_id: the specific chain of the PDB which maps to the UniProt ID
            coverage: the percent coverage of the entire UniProt sequence
            resolution: the resolution of the structure
            start: the structure residue number which maps to the start of the mapped sequence
            end: the structure residue number which maps to the end of the mapped sequence
            unp_start: the sequence residue number which maps to the structure start
            unp_end: the sequence residue number which maps to the structure end
            experimental_method: type of experiment used to determine structure
            tax_id: taxonomic ID of the protein's original organism

    """
    outfile = op.join(outdir, outfile)

    # Load the existing json file
    if outfile and op.exists(outfile) and not force_rerun:
        with open(outfile, 'r') as f:
            raw_data = json.load(f)
        log.debug('{}: Loaded existing json file'.format(uniprot_id))

    # Otherwise run the web request
    else:
        response = requests.get('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{}'.format(uniprot_id),
                                data={'key': 'value'})
        if response.status_code == 404:
            log.debug('{}: 404 returned, probably no structures available.'.format(uniprot_id))
            raw_data = {uniprot_id: {}}
        else:
            log.debug('{}: Obtained best structures'.format(uniprot_id))
            raw_data = response.json()

        # Write the json file if specified
        if outfile:
            with open(outfile, 'w') as f:
                json.dump(raw_data, f)
            log.debug('{}: Saved json file of best structures'.format(uniprot_id))

    data = dict(raw_data)[uniprot_id]

    # Filter for sequence identity percentage
    if seq_ident_cutoff != 0:
        for result in data:
            if result['coverage'] < seq_ident_cutoff:
                data.remove(result)

    return data


# @cachetools.func.ttl_cache(maxsize=1024)
def blast_pdb(seq, outfile='', outdir='', evalue=0.0001, seq_ident_cutoff=0, link=False, force_rerun=False):
    """Returns a list of BLAST hits of a sequence to available structures in the PDB.

    Args:
        seq (str): Your sequence, in string format
        outfile (str): Name of output file
        outdir (str, optional): Path to output directory. Default is the current directory.
        force_rerun (bool, optional): If existing BLAST results should not be used, set to True. Default is False
        evalue (float, optional): Cutoff for the E-value - filters for significant hits. 0.001 is liberal, 0.0001 is stringent (default).
        link (bool, optional): Set to True if a link to the HTML results should be displated

    Returns:
        list: Rank ordered list of BLAST hits in dictionaries.

    """

    if len(seq) < 12:
        raise ValueError('Sequence must be at least 12 residues long.')
    if link:
        page = 'PDB results page: http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={}&eCutOff={}&maskLowComplexity=yes&matrix=BLOSUM62&outputFormat=HTML'.format(seq, evalue)
        print(page)

    parser = etree.XMLParser(ns_clean=True)

    outfile = op.join(outdir, outfile)
    if outfile and op.exists(outfile) and not force_rerun:
        # Parse the existing XML file
        tree = etree.parse(outfile, parser)
        log.debug('{}: Loaded existing BLAST XML results'.format(outfile))
    else:
        # Load the BLAST XML results if force_rerun=True
        page = 'http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={}&eCutOff={}&maskLowComplexity=yes&matrix=BLOSUM62&outputFormat=XML'.format(
                        seq, evalue)
        req = requests.get(page)
        if req.status_code == 200:
            response = req.text

            # Save the XML file
            if outfile:
                with open(outfile, 'w') as f:
                    f.write(response)

            # Parse the XML string
            tree = etree.ElementTree(etree.fromstring(response, parser))
            log.debug('Loaded BLAST results from REST server')
        else:
            log.error('BLAST request timed out')
            return []

    # Get length of original sequence to calculate percentages
    len_orig = float(len(seq))

    root = tree.getroot()
    hit_list = []

    for hit in root.findall('BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
        info = {}

        hitdef = hit.find('Hit_def')
        if hitdef is not None:
            info['hit_pdb'] = hitdef.text.split('|')[0].split(':')[0]
            info['hit_pdb_chains'] = hitdef.text.split('|')[0].split(':')[2].split(',')

        # One PDB can align to different parts of the sequence
        # Will just choose the top hit for this single PDB
        hsp = hit.findall('Hit_hsps/Hsp')[0]

        # Number of identical residues
        hspi = hsp.find('Hsp_identity')
        if hspi is not None:
            info['hit_num_ident'] = int(hspi.text)
            info['hit_percent_ident'] = int(hspi.text)/len_orig

            if int(hspi.text)/len_orig < seq_ident_cutoff:
                log.debug('{}: does not meet sequence identity cutoff'.format(hitdef.text.split('|')[0].split(':')[0]))
                continue

        # Number of similar residues (positive hits)
        hspp = hsp.find('Hsp_positive')
        if hspp is not None:
            info['hit_num_similar'] = int(hspp.text)
            info['hit_percent_similar'] = int(hspp.text) / len_orig

        # Total number of gaps (unable to align in either query or subject)
        hspg = hsp.find('Hsp_gaps')
        if hspg is not None:
            info['hit_num_gaps'] = int(hspg.text)
            info['hit_percent_gaps'] = int(hspg.text) / len_orig

        # E-value of BLAST
        hspe = hsp.find('Hsp_evalue')
        if hspe is not None:
            info['hit_evalue'] = float(hspe.text)

        # Score of BLAST
        hsps = hsp.find('Hsp_score')
        if hsps is not None:
            info['hit_score'] = float(hsps.text)

        hit_list.append(info)

    log.debug("{}: Number of BLAST hits".format(len(hit_list)))
    return hit_list


def blast_pdb_df(seq, xml_outfile='', xml_outdir='', force_rerun=False, evalue=0.001, seq_ident_cutoff=0, link=False):
    """Make a dataframe of BLAST results

    Args: see blast_pdb

    Returns:
        Tuple of (original results, parsed results in Pandas DataFrame)

    """
    blast_results = blast_pdb(seq=seq,
                              outfile=xml_outfile,
                              outdir=xml_outdir,
                              force_rerun=force_rerun,
                              evalue=evalue,
                              seq_ident_cutoff=seq_ident_cutoff,
                              link=link)
    cols = ['hit_pdb', 'hit_pdb_chains', 'hit_evalue', 'hit_score', 'hit_num_ident', 'hit_percent_ident',
            'hit_num_similar', 'hit_percent_similar', 'hit_num_gaps', 'hit_percent_gaps']
    return pd.DataFrame.from_records(blast_results, columns=cols)


# @cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
def _property_table():
    """Download the PDB -> resolution table directly from the RCSB PDB REST service.

    See the other fields that you can get here: http://www.rcsb.org/pdb/results/reportField.do

    Returns:
        Pandas DataFrame: table of structureId as the index, resolution and experimentalTechnique as the columns

    """
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=structureId,resolution,experimentalTechnique,releaseDate&service=wsfile&format=csv'
    r = requests.get(url)
    p = pd.read_csv(StringIO(r.text)).set_index('structureId')
    return p


def get_resolution(pdb_id):
    """Quick way to get the resolution of a PDB ID using the table of results from the REST service

    Returns infinity if the resolution is not available.

    Returns:
        float: resolution of a PDB ID in Angstroms

    TODO:
        - Unit test

    """

    pdb_id = pdb_id.upper()
    if pdb_id not in _property_table().index:
        raise ValueError('PDB ID not in property table')
    else:
        resolution = _property_table().ix[pdb_id, 'resolution']
        if pd.isnull(resolution):
            log.debug('{}: no resolution available, probably not an X-ray crystal structure')
            resolution = float('inf')

    return resolution


# def get_taxonomy(pdb_id):
#     """Quick way to get the taxonomy of a PDB ID using the table of results from the REST service
#
#     Returns None if the taxonomy is not available.
# TODO: taxonomy adds chain ID in the table, need to deal with that
#
#     Returns:
#         str: Organism of a PDB ID
#
#     """
#
#     pdb_id = pdb_id.upper()
#     if pdb_id not in _property_table().index:
#         raise ValueError('PDB ID not in property table')
#     else:
#         taxonomy = _property_table().ix[pdb_id, 'taxonomy']
#         if pd.isnull(taxonomy):
#             log.debug('{}: no taxonomy available')
#             taxonomy = None
#
#     return taxonomy


def get_release_date(pdb_id):
    """Quick way to get the release date of a PDB ID using the table of results from the REST service

    Returns None if the release date is not available.

    Returns:
        str: Organism of a PDB ID

    """

    pdb_id = pdb_id.upper()
    if pdb_id not in _property_table().index:
        raise ValueError('PDB ID not in property table')
    else:
        release_date = _property_table().ix[pdb_id, 'releaseDate']
        if pd.isnull(release_date):
            log.debug('{}: no taxonomy available')
            release_date = None

    return release_date


# @cachetools.func.ttl_cache(maxsize=1, ttl=SEVEN_DAYS)
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


if __name__ == '__main__':
    seq = 'VLSPADKTNVKAAWGVKALSPADKTNVKAALTAVAHVDDMPNAL'
    print(blast_pdb(seq, evalue=1))
    print(len(blast_pdb(seq, evalue=1)))
    print(len(blast_pdb(seq, evalue=1, seq_ident_cutoff=.83)))