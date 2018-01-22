import requests
from lxml import etree
import ssbio.utils


def cctop_submit(seq_str):
    """Submit a protein sequence string to CCTOP and return the job ID.

    Args:
        seq_str (str): Protein sequence as a string

    Returns:
        dict: Job ID on the CCTOP server

    """
    url = 'http://cctop.enzim.ttk.mta.hu/php/submit.php?sequence={}&tmFilter&signalPred'.format(seq_str)
    r = requests.post(url)
    jobid = r.text.split('ID: ')[1]

    return jobid


def cctop_check_status(jobid):
    """Check the status of a CCTOP job ID.

    Args:
        jobid (str): Job ID obtained when job was submitted

    Returns:
        str: 'Finished' if the job is finished and results ready to be downloaded, 'Running' if still in progress,
        'Invalid' for any errors.

    """
    status = 'http://cctop.enzim.ttk.mta.hu/php/poll.php?jobId={}'.format(jobid)
    status_text = requests.post(status)
    return status_text.text


def cctop_save_xml(jobid, outpath):
    """Save the CCTOP results file in XML format.

    Args:
        jobid (str): Job ID obtained when job was submitted
        outpath (str): Path to output filename

    Returns:
        str: Path to output filename

    """
    status = cctop_check_status(jobid=jobid)
    if status == 'Finished':
        result = 'http://cctop.enzim.ttk.mta.hu/php/result.php?jobId={}'.format(jobid)
        result_text = requests.post(result)
        with open(outpath, 'w') as f:
            f.write(result_text.text)
        return outpath
    else:
        raise ConnectionRefusedError('CCTOP job incomplete, status is "{}"'.format(status))


def parse_cctop_full(infile):
    """Parse a CCTOP XML results file and return a list of the consensus TM domains in the format::

            [(1, inside_outside_or_tm),
             (2, inside_outside_or_tm),
             ...]

    Where the first value of a tuple is the sequence residue number, and the second is the predicted location with the
    values 'I' (inside), 'O' (outside), or 'M' (membrane).

    Args:
        infile (str): Path to CCTOP XML file

    Returns:
        list: List of tuples in the format described above

    """
    parser = etree.XMLParser(ns_clean=True)
    with open(infile, 'r') as f:
        tree = etree.fromstring(f.read(), parser)

    all_info = []

    if tree.find('Topology') is not None:
        for r in tree.find('Topology').findall('Region'):
            region_start = int(r.attrib['from'])
            region_end = int(r.attrib['to'])
            region = r.attrib['loc']
            for i in range(region_start, region_end + 1):
                all_info.append((i, region))

    return all_info


def parse_cctop_regions(infile):
    """Parse a CCTOP XML results file and return a dictionary of labeled TM domains in the format::

            {'I1': [31, 32, 33, 34, 35, 36],
             'I2': [73, 74, 75, 76, 77, 78],
             'I3': [133, ...],
             ...}

    Where the regions are labeled as they are sequentially found along with a list of the residue numbers that make up
    the region.

    Args:
        infile (str): Path to CCTOP XML file

    Returns:
        dict: Dictionary of lists in the format described above

    """
    parsed = parse_cctop_full(infile)
    if not parsed:
        return {}
    return ssbio.utils.label_sequential_regions(parsed)