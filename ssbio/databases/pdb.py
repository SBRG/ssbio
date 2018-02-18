"""
PDBProp
=======
"""

import gzip
import json
import logging
import os.path as op
import mmtf
import os
from cobra.core import DictList
import pandas as pd
import requests
import deprecation
from Bio.PDB import PDBList
from lxml import etree
from six.moves.urllib_error import URLError
from six.moves.urllib.request import urlopen, urlretrieve

import ssbio.databases.pisa as pisa
import ssbio.utils
from ssbio.protein.structure.structprop import StructProp

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

log = logging.getLogger(__name__)


class PDBProp(StructProp):

    """Store information about a protein structure from the Protein Data Bank.

    Extends the :class:`~ssbio.protein.structure.structprop.StructProp` class to allow initialization of the structure
    by its PDB ID, and then enabling downloads of the structure file as well as parsing its metadata.

    Args:
        ident (str):
        description (str):
        chains (str):
        mapped_chains (str):
        structure_path (str):
        file_type (str): ``pdb``, ``mmCif``, ``xml``, ``mmtf`` - file type for files downloaded from the PDB

    """

    def __init__(self, ident, description=None, chains=None, mapped_chains=None, structure_path=None, file_type=None):
        StructProp.__init__(self, ident, description=description, chains=chains, mapped_chains=mapped_chains,
                            is_experimental=True, structure_path=structure_path, file_type=file_type)
        self.experimental_method = None
        self.resolution = None
        self.date = None
        self.taxonomy_name = None
        self.biological_assemblies = DictList()
        """DictList: A list for storing Bioassembly objects related to this PDB ID"""

    def download_structure_file(self, outdir, file_type=None, load_header_metadata=True, force_rerun=False):
        """Download a structure file from the PDB, specifying an output directory and a file type. Optionally download
        the mmCIF header file and parse data from it to store within this object.

        Args:
            outdir (str): Path to output directory
            file_type (str): ``pdb``, ``mmCif``, ``xml``, ``mmtf`` - file type for files downloaded from the PDB
            load_header_metadata (bool): If header metadata should be loaded into this object, fastest with mmtf files
            force_rerun (bool): If structure file should be downloaded even if it already exists

        """
        ssbio.utils.double_check_attribute(object=self, setter=file_type, backup_attribute='file_type',
                                           custom_error_text='Please set file type to be downloaded from the PDB: '
                                                             'pdb, mmCif, xml, or mmtf')

        # XTODO: check if outfile exists using ssbio.utils.force_rerun, pdblist seems to take long if it exists
        p = PDBList()
        with ssbio.utils.suppress_stdout():
            structure_file = p.retrieve_pdb_file(pdb_code=self.id, pdir=outdir, file_format=file_type, overwrite=force_rerun)
        if not op.exists(structure_file):
            log.debug('{}: {} file not available'.format(self.id, file_type))
            raise URLError('{}.{}: file not available to download'.format(self.id, file_type))
        else:
            log.debug('{}: {} file saved'.format(self.id, file_type))
            self.load_structure_path(structure_file, file_type)
            if load_header_metadata and file_type == 'mmtf':
                self.update(parse_mmtf_header(structure_file))
            if load_header_metadata and file_type != 'mmtf':
                self.update(parse_mmcif_header(download_mmcif_header(pdb_id=self.id, outdir=outdir, force_rerun=force_rerun)))

    def get_pisa_complex_predictions(self, outdir, existing_pisa_multimer_xml=None):
        if not existing_pisa_multimer_xml:
            pisa_xmls = pisa.download_pisa_multimers_xml(pdb_ids=self.id, outdir=outdir,
                                                         save_single_xml_files=True)
        else:
            pisa_xmls = {}
            pisa_xmls[self.id] = existing_pisa_multimer_xml
        pisa_dict = pisa.parse_pisa_multimers_xml(pisa_xmls[self.id], download_structures=True,
                                                  outdir=outdir)

    def __json_encode__(self):
        # TODO: investigate why saving with # does not work!
        to_return = {}
        for x in self.__dict__.keys():
            if x == 'pdb_title' or x == 'description':
                sanitized = ssbio.utils.force_string(getattr(self, x)).replace('#', '-')
            else:
                to_return.update({x: getattr(self, x)})
        return to_return


def parse_mmtf_header(infile):
    """Parse an MMTF file and return basic header-like information.

    Args:
        infile (str): Path to MMTF file

    Returns:
        dict: Dictionary of parsed header

    Todo:
        - Can this be sped up by not parsing the 3D coordinate info somehow?
        - OR just store the sequences when this happens since it is already being parsed.

    """
    infodict = {}

    mmtf_decoder = mmtf.parse(infile)
    infodict['date'] = mmtf_decoder.deposition_date
    infodict['release_date'] = mmtf_decoder.release_date
    infodict['experimental_method'] = [x.decode() for x in mmtf_decoder.experimental_methods]
    infodict['resolution'] = mmtf_decoder.resolution
    infodict['description'] = mmtf_decoder.title

    group_name_exclude = ['HOH']
    chem_comp_type_exclude = ['l-peptide linking', 'peptide linking']
    chemicals = list(set([mmtf_decoder.group_list[idx]['groupName'] for idx in mmtf_decoder.group_type_list if mmtf_decoder.group_list[idx]['chemCompType'].lower() not in chem_comp_type_exclude and mmtf_decoder.group_list[idx]['groupName'] not in group_name_exclude]))
    infodict['chemicals'] = chemicals
    return infodict

def download_mmcif_header(pdb_id, outdir='', force_rerun=False):
    """Download a mmCIF header file from the RCSB PDB by ID.

    Args:
        pdb_id: PDB ID
        outdir: Optional output directory, default is current working directory
        force_rerun: If the file should be downloaded again even if it exists

    Returns:
        str: Path to outfile

    """
    # method in biopython. extra file types have not been added to biopython download yet

    pdb_id = pdb_id.lower()
    file_type = 'cif'
    folder = 'header'
    outfile = op.join(outdir, '{}.header.{}'.format(pdb_id, file_type))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        download_link = 'http://files.rcsb.org/{}/{}.{}'.format(folder, pdb_id, file_type)
        urlretrieve(download_link, outfile)
        log.debug('{}: saved header file'.format(outfile))
    else:
        log.debug('{}: header file already saved'.format(outfile))

    return outfile


def parse_mmcif_header(infile):
    """Parse a couple important fields from the mmCIF file format with some manual curation of ligands.

    If you want full access to the mmCIF file just use the MMCIF2Dict class in Biopython.

    Args:
        infile: Path to mmCIF file

    Returns:
        dict: Dictionary of parsed header

    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    newdict = {}
    try:
        mmdict = MMCIF2Dict(infile)
    except ValueError as e:
        log.exception(e)
        return newdict

    chemical_ids_exclude = ['HOH']
    chemical_types_exclude = ['l-peptide linking','peptide linking']

    if '_struct.title' in mmdict:
        newdict['pdb_title'] = mmdict['_struct.title']
    else:
        log.debug('{}: No title field'.format(infile))

    if '_struct.pdbx_descriptor' in mmdict:
        newdict['description'] = mmdict['_struct.pdbx_descriptor']
    else:
        log.debug('{}: no description field'.format(infile))

    if '_pdbx_database_status.recvd_initial_deposition_date' in mmdict:
        newdict['date'] = mmdict['_pdbx_database_status.recvd_initial_deposition_date']
    elif '_database_PDB_rev.date' in mmdict:
        newdict['date'] = mmdict['_database_PDB_rev.date']
    else:
        log.debug('{}: no date field'.format(infile))

    if '_exptl.method' in mmdict:
        newdict['experimental_method'] = mmdict['_exptl.method']
    else:
        log.debug('{}: no experimental method field'.format(infile))

    # TODO: refactor how to get resolutions based on experimental method
    if '_refine.ls_d_res_high' in mmdict:
        try:
            if isinstance(mmdict['_refine.ls_d_res_high'], list):
                newdict['resolution'] = [float(x) for x in mmdict['_refine.ls_d_res_high']]
            else:
                newdict['resolution'] = float(mmdict['_refine.ls_d_res_high'])
        except:
            try:
                newdict['resolution'] = float(mmdict['_em_3d_reconstruction.resolution'])
            except:
                log.debug('{}: no resolution field'.format(infile))
    else:
        log.debug('{}: no resolution field'.format(infile))

    if '_chem_comp.id' in mmdict:
        chemicals_filtered = ssbio.utils.filter_list_by_indices(mmdict['_chem_comp.id'],
                                                            ssbio.utils.not_find(mmdict['_chem_comp.type'],
                                                                           chemical_types_exclude,
                                                                           case_sensitive=False))
        chemicals_fitered = ssbio.utils.filter_list(chemicals_filtered, chemical_ids_exclude, case_sensitive=True)
        newdict['chemicals'] = chemicals_fitered
    else:
        log.debug('{}: no chemical composition field'.format(infile))

    if '_entity_src_gen.pdbx_gene_src_scientific_name' in mmdict:
        newdict['taxonomy_name'] = mmdict['_entity_src_gen.pdbx_gene_src_scientific_name']
    else:
        log.debug('{}: no organism field'.format(infile))

    return newdict


def download_sifts_xml(pdb_id, outdir='', force_rerun=False):
    """Download the SIFTS file for a PDB ID.

    Args:
        pdb_id (str): PDB ID
        outdir (str): Output directory, current working directory if not specified.
        force_rerun (bool): If the file should be downloaded again even if it exists

    Returns:
        str: Path to downloaded file

    """
    baseURL = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    filename = '{}.xml.gz'.format(pdb_id.lower())

    outfile = op.join(outdir, filename.split('.')[0] + '.sifts.xml')

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        response = urlopen(baseURL + filename)
        with open(outfile, 'wb') as f:
            f.write(gzip.decompress(response.read()))

    return outfile


def map_uniprot_resnum_to_pdb(uniprot_resnum, chain_id, sifts_file):
    """Map a UniProt residue number to its corresponding PDB residue number.

    This function requires that the SIFTS file be downloaded,
    and also a chain ID (as different chains may have different mappings).

    Args:
        uniprot_resnum (int): integer of the residue number you'd like to map
        chain_id (str): string of the PDB chain to map to
        sifts_file (str): Path to the SIFTS XML file

    Returns:
        (tuple): tuple containing:

            mapped_resnum (int): Mapped residue number
            is_observed (bool): Indicates if the 3D structure actually shows the residue

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
        # TODO: IMPORTANT - entityId is not the chain ID!!! it is just in alphabetical order!
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


def best_structures(uniprot_id, outname=None, outdir=None, seq_ident_cutoff=0.0, force_rerun=False):
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
        outname (str): Basename of the output file of JSON results
        outdir (str): Path to output directory of JSON results
        seq_ident_cutoff (float): Cutoff results based on percent coverage (in decimal form)
        force_rerun (bool): Obtain best structures mapping ignoring previously downloaded results

    Returns:
        list: Rank-ordered list of dictionaries representing chain-specific PDB entries. Keys are:
            * pdb_id: the PDB ID which maps to the UniProt ID
            * chain_id: the specific chain of the PDB which maps to the UniProt ID
            * coverage: the percent coverage of the entire UniProt sequence
            * resolution: the resolution of the structure
            * start: the structure residue number which maps to the start of the mapped sequence
            * end: the structure residue number which maps to the end of the mapped sequence
            * unp_start: the sequence residue number which maps to the structure start
            * unp_end: the sequence residue number which maps to the structure end
            * experimental_method: type of experiment used to determine structure
            * tax_id: taxonomic ID of the protein's original organism

    """
    outfile = ''

    if not outdir:
        outdir = ''

    # if output dir is specified but not outname, use the uniprot
    if not outname and outdir:
        outname = uniprot_id

    if outname:
        outname = op.join(outdir, outname)
        outfile = '{}.json'.format(outname)

    # Load a possibly existing json file
    if not ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        with open(outfile, 'r') as f:
            raw_data = json.load(f)
        log.debug('{}: loaded existing json file'.format(uniprot_id))

    # Otherwise run the web request
    else:
        # TODO: add a checker for a cached file of uniprot -> PDBs - can be generated within gempro pipeline and stored
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


def blast_pdb(seq, outfile='', outdir='', evalue=0.0001, seq_ident_cutoff=0.0, link=False, force_rerun=False):
    """Returns a list of BLAST hits of a sequence to available structures in the PDB.

    Args:
        seq (str): Your sequence, in string format
        outfile (str): Name of output file
        outdir (str, optional): Path to output directory. Default is the current directory.
        evalue (float, optional): Cutoff for the E-value - filters for significant hits. 0.001 is liberal, 0.0001 is stringent (default).
        seq_ident_cutoff (float, optional): Cutoff results based on percent coverage (in decimal form)
        link (bool, optional): Set to True if a link to the HTML results should be displayed
        force_rerun (bool, optional): If existing BLAST results should not be used, set to True. Default is False

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
    if ssbio.utils.force_rerun(force_rerun, outfile):
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
    else:
        tree = etree.parse(outfile, parser)
        log.debug('{}: Loaded existing BLAST XML results'.format(outfile))

    # Get length of original sequence to calculate percentages
    len_orig = float(len(seq))

    root = tree.getroot()
    hit_list = []

    for hit in root.findall('BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
        info = {}

        hitdef = hit.find('Hit_def')
        if hitdef is not None:
            info['hit_pdb'] = hitdef.text.split('|')[0].split(':')[0].lower()
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


def blast_pdb_df(blast_results):
    """Make a dataframe of BLAST results"""
    cols = ['hit_pdb', 'hit_pdb_chains', 'hit_evalue', 'hit_score', 'hit_num_ident', 'hit_percent_ident',
            'hit_num_similar', 'hit_percent_similar', 'hit_num_gaps', 'hit_percent_gaps']
    return pd.DataFrame.from_records(blast_results, columns=cols)


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
            log.debug('{}: no release date available')
            release_date = None

    return release_date


def get_num_bioassemblies(pdb_id, cache=False, outdir=None, force_rerun=False):
    """Check if there are bioassemblies using the PDB REST API, and if there are, get the number of bioassemblies
    available.

    See: https://www.rcsb.org/pages/webservices/rest, section 'List biological assemblies'

    Not all PDB entries have biological assemblies available and some have multiple. Details that are necessary to
    recreate a biological assembly from the asymmetric unit can be accessed from the following requests.

    - Number of biological assemblies associated with a PDB entry
    - Access the transformation information needed to generate a biological assembly (nr=0 will return information
      for the asymmetric unit, nr=1 will return information for the first assembly, etc.)

    A query of https://www.rcsb.org/pdb/rest/bioassembly/nrbioassemblies?structureId=1hv4 returns this::

        <nrBioAssemblies structureId="1HV4" hasAssemblies="true" count="2"/>

    Args:
        pdb_id (str): PDB ID
        cache (bool): If the XML file should be downloaded
        outdir (str): If cache, then specify the output directory
        force_rerun (bool): If cache, and if file exists, specify if API should be queried again

    """
    parser = etree.XMLParser(ns_clean=True)

    if not outdir:
        outdir = os.getcwd()
    outfile = op.join(outdir, '{}_nrbiomols.xml'.format(pdb_id))

    if ssbio.utils.force_rerun(force_rerun, outfile):
        page = 'https://www.rcsb.org/pdb/rest/bioassembly/nrbioassemblies?structureId={}'.format(pdb_id)
        req = requests.get(page)

        if req.status_code == 200:
            response = req.text

            # Save the XML file
            if cache:
                with open(outfile, 'w') as f:
                    f.write(response)

            # Parse the XML string
            tree = etree.ElementTree(etree.fromstring(response, parser))
            log.debug('Loaded bioassembly information from REST server')
        else:
            log.error('Request timed out')
            req.raise_for_status()
    else:
        tree = etree.parse(outfile, parser)
        log.debug('{}: Loaded existing XML results'.format(outfile))

    r = tree.getroot()
    has_biomols = r.get('hasAssemblies')
    if has_biomols == 'true':
        has_biomols = True
    else:
        has_biomols = False

    if has_biomols:
        num_biomols = r.get('count')
    else:
        num_biomols = 0

    num_biomols = int(num_biomols)
    return num_biomols


def get_bioassembly_info(pdb_id, biomol_num, cache=False, outdir=None, force_rerun=False):
    """Get metadata about a bioassembly from the RCSB PDB's REST API.

    See: https://www.rcsb.org/pdb/rest/bioassembly/bioassembly?structureId=1hv4&nr=1
    The API returns an XML file containing the information on a biological assembly that looks like this::

        <bioassembly structureId="1HV4" assemblyNr="1" method="PISA" desc="author_and_software_defined_assembly">
            <transformations operator="1" chainIds="A,B,C,D">
                <transformation index="1">
                    <matrix m11="1.00000000" m12="0.00000000" m13="0.00000000" m21="0.00000000" m22="1.00000000" m23="0.00000000" m31="0.00000000" m32="0.00000000" m33="1.00000000"/>
                    <shift v1="0.00000000" v2="0.00000000" v3="0.00000000"/>
                </transformation>
            </transformations>
        </bioassembly>

    Args:
        pdb_id (str): PDB ID
        biomol_num (int): Biological assembly number you are interested in
        cache (bool): If the XML file should be downloaded
        outdir (str): If cache, then specify the output directory
        force_rerun (bool): If cache, and if file exists, specify if API should be queried again

    """
    parser = etree.XMLParser(ns_clean=True)
    #
    # if not outdir:
    #     outdir = os.getcwd()
    # outfile = op.join(outdir, '{}.xml'.format(self.id))
    #
    # if ssbio.utils.force_rerun(force_rerun, outfile):
    #     page = 'https://www.rcsb.org/pdb/rest/bioassembly/bioassembly?structureId={}&nr={}'.format(
    #         self.original_pdb_id, biomol_num)
    #     req = requests.get(page)
    #
    #     if req.status_code == 200:
    #         response = req.text
    #
    #         # Save the XML file
    #         if cache:
    #             with open(outfile, 'w') as f:
    #                 f.write(response)
    #
    #         # Parse the XML string
    #         r = xmltodict.parse(response)
    #         log.debug('Loaded bioassembly information from REST server')
    #     else:
    #         log.error('Request timed out')
    #         req.raise_for_status()
    # else:
    #     with open(outfile, 'r') as f:
    #         r = xmltodict.parse(f.read())
    #     log.debug('{}: Loaded existing XML results'.format(outfile))
    #
    # self.biomol_to_chain_dict[biomol_num] = {'chains': r['bioassembly']['transformations']['@chainIds'],
    #                                          'multiplier': len(r['bioassembly']['transformations']['transformation'])}
    # # TODO: figure out how to store matrices etc.
    #
    # log.info('{}_{}: ')


def download_biomol(pdb_id, biomol_num, outdir, file_type='pdb', force_rerun=False):
    import zlib
    from six.moves.urllib_error import URLError
    from six.moves.urllib.request import urlopen, urlretrieve
    import contextlib

    ssbio.utils.make_dir(outdir)
    server_folder = pdb_id[1:3]

    if file_type == 'pdb':
        # server = 'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/divided/{}/'.format(server_folder)
        server = 'https://files.rcsb.org/download/'
        server_filename = pdb_id + '.pdb%i.gz' % biomol_num
        local_filename = pdb_id + '_bio%i.pdb' % biomol_num
        outfile = op.join(outdir, local_filename)

    elif file_type.lower() == 'mmcif' or file_type.lower() == 'cif':
        server = 'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/mmCIF/divided/{}/'.format(server_folder)
        server_filename = pdb_id + '-assembly%i.cif.gz' % biomol_num
        local_filename = pdb_id + '_bio%i.cif' % biomol_num
        outfile = op.join(outdir, local_filename)

    else:
        raise ValueError('Biological assembly only available in PDB or mmCIF file types.')

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        download_link = op.join(server, server_filename)
        try:
            with contextlib.closing(urlopen(download_link)) as f:
                decompressed_data = zlib.decompress(f.read(), 16 + zlib.MAX_WBITS)
                with open(op.join(outdir, local_filename), 'wb') as f:
                    f.write(decompressed_data)
        except URLError as e:
            print(e)
            return None

    return outfile


########################################################################################################################
########################################################################################################################
# DEPRECATED FUNCTIONS
########################################################################################################################
########################################################################################################################


@deprecation.deprecated(deprecated_in="1.0", removed_in="2.0",
                        details="Use Biopython's PDBList.retrieve_pdb_file function instead")
def download_structure(pdb_id, file_type, outdir='', only_header=False, force_rerun=False):
    """Download a structure from the RCSB PDB by ID. Specify the file type desired.

    Args:
        pdb_id: PDB ID
        file_type: pdb, pdb.gz, mmcif, cif, cif.gz, xml.gz, mmtf, mmtf.gz
        outdir: Optional output directory
        only_header: If only the header file should be downloaded
        force_rerun: If the file should be downloaded again even if it exists

    Returns:
        str: Path to outfile

    """
    # method in biopython. extra file types have not been added to biopython download yet

    pdb_id = pdb_id.lower()
    file_type = file_type.lower()
    file_types = ['pdb', 'pdb.gz', 'mmcif', 'cif', 'cif.gz', 'xml.gz', 'mmtf', 'mmtf.gz']
    if file_type not in file_types:
        raise ValueError('Invalid file type, must be either: pdb, pdb.gz, cif, cif.gz, xml.gz, mmtf, mmtf.gz')

    if file_type == 'mmtf':
        file_type = 'mmtf.gz'

    if file_type.endswith('.gz'):
        gzipped = True
    else:
        gzipped = False

    if file_type == 'mmcif':
        file_type = 'cif'

    if only_header:
        folder = 'header'
        outfile = op.join(outdir, '{}.header.{}'.format(pdb_id, file_type))
    else:
        folder = 'download'
        outfile = op.join(outdir, '{}.{}'.format(pdb_id, file_type))

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        if file_type == 'mmtf.gz' or file_type == 'mmtf':
            mmtf_api = '1.0'
            download_link = 'http://mmtf.rcsb.org/v{}/full/{}.mmtf.gz'.format(mmtf_api, pdb_id)
        else:
            download_link = 'http://files.rcsb.org/{}/{}.{}'.format(folder, pdb_id, file_type)

        urlretrieve(download_link, outfile)

        if gzipped:
            outfile = ssbio.utils.gunzip_file(infile=outfile,
                                              outfile=outfile.strip('.gz'),
                                              outdir=outdir,
                                              delete_original=False,
                                              force_rerun_flag=force_rerun)

        log.debug('{}: saved structure file'.format(outfile))
    else:
        if file_type == 'mmtf.gz':
            outfile = op.join(outdir, '{}.{}'.format(pdb_id, 'mmtf'))
        log.debug('{}: structure file already saved'.format(outfile))

    return outfile
