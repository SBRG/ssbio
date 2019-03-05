from os import path as op

import requests
from lxml import etree

import ssbio.utils
import logging
log = logging.getLogger(__name__)

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