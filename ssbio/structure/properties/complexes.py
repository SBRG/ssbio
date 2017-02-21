import requests
from collections import defaultdict
from copy import deepcopy
import ssbio.utils
import os
import os.path as op
from lxml import etree
import logging

log = logging.getLogger(__name__)


def download_pisa_multimers_xml(pdb_ids, outdir=None, force_rerun=False):
    """Download the PISA XML file for multimers.

    See: http://www.ebi.ac.uk/pdbe/pisa/pi_download.html for more info

    XML description of macromolecular assemblies:
        http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?pdbcodelist
        where "pdbcodelist" is a comma-separated (strictly no spaces) list of PDB codes. The resulting file contain XML
        output of assembly data, equivalent to that displayed in PISA assembly pages, for each of the specified PDB
        entries.   NOTE: If a mass-download is intended, please minimize the number of retrievals by specifying as many
        PDB codes in the URL as feasible (20-50 is a good range), and never send another URL request until the previous
        one has been completed (meaning that the multimers.pisa file has been downloaded). Excessive requests will
        silently die in the server queue.

    Args:
        pdb_ids (str, list): PDB ID or list of IDs
        outdir (str): Directory to output PISA XML files
        force_rerun (bool): Redownload files if they already exist

    Returns:
        list: of files downloaded

    """
    if not outdir:
        outdir = os.getcwd()

    pdb_ids = ssbio.utils.force_lower_list(sorted(pdb_ids))
    split_list = ssbio.utils.split_list_by_n(pdb_ids, 50)

    files = []

    for l in split_list:
        pdbs = ','.join(l)
        filename = op.join(outdir, '{}.xml'.format(pdbs))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=filename):
            # Request the xml file of the results
            r = requests.get('http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?{}'.format(pdbs))

            # Save results to a file
            with open(filename, 'w') as f:
                f.write(r.text)

            log.debug('Downloaded PISA results')
        else:
            log.debug('PISA results already downloaded')

        files.append(filename)

    return files


def parse_pisa_multimers_xml(pisa_multimers_xml, download_structures=False, outdir=None, force_rerun=False):
    """Retrieve PISA information from an XML results file

    See: http://www.ebi.ac.uk/pdbe/pisa/pi_download.html for more info

    XML description of macromolecular assemblies:
        http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?pdbcodelist
        where "pdbcodelist" is a comma-separated (strictly no spaces) list of PDB codes. The resulting file contain XML
        output of assembly data, equivalent to that displayed in PISA assembly pages, for each of the specified PDB
        entries.   NOTE: If a mass-download is intended, please minimize the number of retrievals by specifying as many
        PDB codes in the URL as feasible (20-50 is a good range), and never send another URL request until the previous
        one has been completed (meaning that the multimers.pisa file has been downloaded). Excessive requests will
        silently die in the server queue.

    Args:
        pdb_id: 4 character PDB ID

    Returns:

    """
    if not outdir:
        outdir = os.getcwd()

    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(pisa_multimers_xml, parser)
    root = tree.getroot()

    pisa = defaultdict(dict)

    for pdb in root.findall('pdb_entry'):
        # Get the PDB ID
        pdb_id = pdb.find('pdb_code').text

        # Check the assembly status
        status = pdb.find('status').text
        errors = ['Entry not found', 'Overlapping structures', 'No symmetry operations']
        if status in errors:
            pisa[pdb_id]['status'] = status
            continue

        # Check monomer status
        num_complexes = int(pdb.find('total_asm').text)
        if num_complexes == 0:
            pisa[pdb_id]['status'] = 'MONOMER'
            continue

        elif num_complexes > 0:
            # All "assembly sets" (see PISA sets for more info)
            sets = pdb.findall('asm_set')

            for s in sets:
                set_id = int(s.find('ser_no').text)

                # All assemblies
                assemblies = s.findall('assembly')

                for cplx in assemblies:

                    ############################################################################################
                    # This part tells you the actual composition of the predicted complex (chains and ligands)
                    parts = cplx.findall('molecule')

                    chains = defaultdict(int)

                    for part in parts:
                        part_id = part.find('chain_id').text
                        if part_id.startswith('['):
                            part_id = 'LIG_' + part_id.split(']')[0].strip('[')
                        chains[str(part_id)] += 1

                    ligands = {}

                    for key in deepcopy(chains).keys():
                        if key.startswith('LIG_'):
                            ligands[str(key.split('_')[1])] = chains.pop(key)
                    ############################################################################################

                    adder = {}

                    cplx_id = int(cplx.find('id').text)
                    cplx_composition = str(cplx.find('composition').text)
                    d_g_diss = float(cplx.find('diss_energy').text)
                    d_g_int = float(cplx.find('int_energy').text)
                    pdb_biomol = int(cplx.find('R350').text)

                    if d_g_diss >= 0:
                        stable = True
                    else:
                        stable = False

                    adder['cplx_composition'] = cplx_composition.strip()
                    adder['cplx_chains'] = chains
                    adder['cplx_ligands'] = ligands
                    adder['stable'] = stable
                    adder['d_g_diss'] = d_g_diss
                    adder['d_g_int'] = d_g_int
                    adder['pdb_biomol'] = pdb_biomol

                    pisa[pdb_id][(set_id, cplx_id)] = adder

                    if download_structures:
                        ident = '{}:{},{}'.format(pdb_id, set_id, cplx_id)
                        filename = op.join(outdir, ident + '.pdb')

                        if ssbio.utils.force_rerun(flag=force_rerun, outfile=filename):
                            download_structure_link = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimer.pdb?{}'.format(
                                ident)
                            r = requests.get(download_structure_link)

                            with open(filename, 'w') as f:
                                f.write(r.text)

                            log.debug('{}: downloaded structure file'.format(ident))
                        else:
                            log.debug('{}: structure file already downloaded'.format(ident))

                        pisa[pdb_id][(set_id, cplx_id)]['structure_file'] = filename

    return pisa


def pdb_chain_stoichiometry_biomolone(pdbid):
    """Get the stoichiometry of the chains in biological assembly 1 as a dictionary.

    Steps taken are:
    1) Download PDB and parse header, make biomolecule if provided
    2) Count how many times each chain appears in biomolecule #1
    3) Convert chain id to uniprot id
    4) Return final dictionary

    Args:
        pdbid (str): 4 character PDB ID

    Returns:
        dict: {(ChainID,UniProtID): # occurences}
    """
    pass