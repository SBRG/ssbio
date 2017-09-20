import requests
from collections import defaultdict
from copy import deepcopy
import ssbio.utils
import os
import os.path as op
from lxml import etree
import logging
import glob

log = logging.getLogger(__name__)


def download_pisa_multimers_xml(pdb_ids, save_single_xml_files=True, outdir=None, force_rerun=False):
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
        save_single_xml_files (bool): If single XML files should be saved per PDB ID. If False, if multiple PDB IDs are
            provided, then a single, combined XML output file is downloaded
        outdir (str): Directory to output PISA XML files
        force_rerun (bool): Redownload files if they already exist

    Returns:
        list: of files downloaded

    """
    if not outdir:
        outdir = os.getcwd()

    files = {}
    pdb_ids = ssbio.utils.force_lower_list(sorted(pdb_ids))

    # If we want to save single PISA XML files per PDB ID...
    if save_single_xml_files:
        # Check for existing PISA XML files
        if not force_rerun:
            existing_files = [op.basename(x) for x in glob.glob(op.join(outdir, '*_multimers.pisa.xml'))]
            # Store the paths to these files to return
            files = {v.split('_')[0]: op.join(outdir, v) for v in existing_files}
            log.debug('Already downloaded PISA files for {}'.format(list(files.keys())))
        else:
            existing_files = []

        # Filter PDB IDs based on existing file
        pdb_ids = [x for x in pdb_ids if '{}_multimers.pisa.xml'.format(x) not in existing_files]

        # Split the list into 50 to limit requests
        split_list = ssbio.utils.split_list_by_n(pdb_ids, 40)

        # Download PISA files
        for l in split_list:
            pdbs = ','.join(l)
            all_pisa_link = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?{}'.format(pdbs)
            r = requests.get(all_pisa_link)

            # Parse PISA file and save individual XML files
            parser = etree.XMLParser(ns_clean=True)
            tree = etree.fromstring(r.text, parser)
            for pdb in tree.findall('pdb_entry'):
                filename = op.join(outdir, '{}_multimers.pisa.xml'.format(pdb.find('pdb_code').text))
                add_root = etree.Element('pisa_multimers')
                add_root.append(pdb)
                with open(filename, 'wb') as f:
                    f.write(etree.tostring(add_root))
                files[pdb.find('pdb_code').text] = filename
                log.debug('{}: downloaded PISA results'.format(pdb))

    else:
        split_list = ssbio.utils.split_list_by_n(pdb_ids, 40)

        for l in split_list:
            pdbs = ','.join(l)
            all_pisa_link = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?{}'.format(pdbs)

            filename = op.join(outdir, '{}_multimers.pisa.xml'.format(pdbs))

            if ssbio.utils.force_rerun(flag=force_rerun, outfile=filename):
                r = requests.get(all_pisa_link)
                with open(filename, 'w') as f:
                    f.write(r.text)
                log.debug('Downloaded PISA results')
            else:
                log.debug('PISA results already downloaded')

            for x in l:
                files[x] = filename

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
        pisa_multimers_xml (str): Path to PISA XML output file
        download_structures (bool): If assembly files should be downloaded
        outdir (str): Directory to output assembly files
        force_rerun (bool): Redownload files if they already exist

    Returns:
        dict: of parsed PISA information

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