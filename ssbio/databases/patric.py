import ftplib
import os.path as op
import logging
import os
log = logging.getLogger(__name__)


def download_coding_sequences(patric_id, seqtype, outdir='', outfile='', force_rerun=False):
    """Download the entire set of DNA or protein sequences from protein-encoding genes in a genome from NCBI.

    Saves a FASTA file in the optional directory specified.

    Args:
        genome_accession_or_id (str): PATRIC ID
        seqtype (str): "dna" or "protein" - if you want the coding sequences in DNA or amino acid formats.
        outdir (str): optional output directory (default is the current directory)
        outfile (str): optional custom name for file
        force_rerun (bool): if you want to redownload existing files

    Returns:
        Path to downloaded FASTA file.

    """
    if seqtype == 'dna':
        extension = 'ffn'
    elif seqtype == 'protein':
        extension = 'faa'
    else:
        raise ValueError('seqtype must be "dna" or "protein"')

    # TODO: use utils functions here
    # path and filename parsing
    if outfile:
        outfile = op.join(outdir, '{}.{}'.format(outfile, extension))
    else:
        # if no outfile is specified, default is "$GI.PATRIC.faa"
        outfile = op.join(outdir, '{}.PATRIC.{}'.format(patric_id, extension))

    if not force_rerun:
        # return the path to the file if it was already downloaded
        if op.exists(outfile) and os.stat(outfile).st_size != 0:
            log.debug('FASTA file already exists at {}'.format(outfile))
            return outfile

    try:
        ftp = ftplib.FTP('ftp.patricbrc.org')
        ftp.login()
        ftp.cwd("/patric2/patric3/genomes/{0}/".format(patric_id))
        with open(outfile, "wb") as gFile:
            ftp.retrbinary('RETR {0}.PATRIC.{1}'.format(patric_id, extension), gFile.write)
        ftp.quit()
    # TODO: check exceptions
    except:
        return None

    return outfile