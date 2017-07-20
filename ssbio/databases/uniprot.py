import os.path as op
import re
import warnings
import bioservices
import pandas as pd
import requests
from dateutil.parser import parse as dateparse
import ssbio.utils
from ssbio.protein.sequence.seqprop import SeqProp

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from six.moves.urllib.request import urlretrieve
from Bio import SeqIO

import logging
log = logging.getLogger(__name__)
bsup = bioservices.uniprot.UniProt()


class UniProtProp(SeqProp):
    """Class to parse through UniProt metadata all at once
    """

    def __init__(self, uniprot_acc, sequence_path=None, xml_path=None, txt_path=None):
        if not is_valid_uniprot_id(uniprot_acc):
            raise ValueError("{}: invalid UniProt ID!".format(uniprot_acc))

        if xml_path:
            metadata_path = xml_path
            self.file_type = 'xml'
        elif txt_path:
            metadata_path = txt_path
            self.file_type = 'txt'
        else:
            metadata_path = None
            self.file_type = None

        SeqProp.__init__(self, ident=uniprot_acc, sequence_path=sequence_path, metadata_path=metadata_path)
        self.uniprot = uniprot_acc
        self.ec_number = None
        self.entry_version = None
        self.pfam = None
        self.reviewed = False
        self.seq_version = None
        self.reviewed = False

    @property
    def seq_record(self):
        # The SeqRecord object is not kept in memory unless explicitly set, so when saving as a json we only save a
        # pointer to the file
        if self.sequence_path:
            seq_record = SeqIO.read(open(self.sequence_path), 'fasta')
        elif self.metadata_path:
            if self.file_type == 'xml':
                with open(self.metadata_path) as handle:
                    seq_record = SeqIO.read(handle, "uniprot-xml")
            elif self.file_type == 'txt':
                with open(self.metadata_path) as handle:
                    seq_record = SeqIO.read(handle, "swiss")
        else:
            seq_record = self._seq_record

        if seq_record:
            seq_record.letter_annotations = self.letter_annotations
            seq_record.annotations = self.annotations
            seq_record.features = self.features

        return seq_record

    def download_seq_file(self, outdir, force_rerun=False):
        """Download and load the UniProt sequence file"""
        uniprot_seq_file = download_uniprot_file(uniprot_id=self.id,
                                                 filetype='fasta',
                                                 outdir=outdir,
                                                 force_rerun=force_rerun)

        self.load_sequence_path(uniprot_seq_file)

    def download_metadata_file(self, outdir, file_type='xml', simple_parse=False, force_rerun=False):
        """Download and load the UniProt metadata file"""
        uniprot_metadata_file = download_uniprot_file(uniprot_id=self.id,
                                                      filetype=file_type,
                                                      outdir=outdir,
                                                      force_rerun=force_rerun)
        self.file_type = file_type
        self.load_metadata_path(metadata_path=uniprot_metadata_file, simple_parse=simple_parse)
        self.sequence_dir = self.metadata_dir

    def load_metadata_path(self, metadata_path, simple_parse=False):
        """Load a metadata file and provide pointers to its location

        Args:
            metadata_path: Path to metadata file

        """
        if not op.dirname(metadata_path):
            self.metadata_dir = '.'
        else:
            self.metadata_dir = op.dirname(metadata_path)
        self.metadata_file = op.basename(metadata_path)

        if simple_parse:
            self.update(parse_uniprot_txt_file(metadata_path), overwrite=True,
                        only_keys=['description', 'kegg', 'refseq', 'ec_number', 'entry_version', 'gene_name',
                                   'pfam', 'pdbs', 'reviewed', 'seq_version'])

        else:
            # TODO: need to rethink the flow here...
            # Additionally load in metadata from xml or txt files
            if self.file_type == 'xml':
                with open(self.metadata_path) as handle:
                    seq_record = SeqIO.read(handle, "uniprot-xml")
            elif self.file_type == 'txt':
                with open(self.metadata_path) as handle:
                    seq_record = SeqIO.read(handle, "swiss")
            self.update(seq_record.__dict__, overwrite=True, only_keys=['description'])
            # TODO: need to parse the dbxrefs for Pfam, refseq, kegg, pdb, and what about "reviewed"?
            self.update(seq_record.annotations, overwrite=True,
                        only_keys=['gene_name_primary', 'modified', 'sequence_modified'])
            self.annotations = seq_record.annotations
            self.letter_annotations = seq_record.letter_annotations
            self.features = seq_record.features

    def ranking_score(self):
        """Provide a score for this UniProt ID based on reviewed (True=1, False=0) + number of PDBs

        Returns:
            int: Scoring for this ID

        """
        return self.reviewed + self.num_pdbs

    def __json_encode__(self):
        # TODO: investigate why saving with # does not work!
        # Recon3D problem genes: testers = ['1588.2', '6819.1', '27233.1', '2264.1']
        to_return = {}
        for x in self.__dict__.keys():
            if x == 'description':
                sanitized = ssbio.utils.force_string(getattr(self, x)).replace('#', '-')
            else:
                to_return.update({x: getattr(self, x)})
        return to_return


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


# TODO: method to blast UniProt to find a 100% sequence match
def blast_uniprot(seq_str, seq_ident=1, evalue=0.0001, reviewed_only=True, organism=None):
    """BLAST the UniProt db to find what IDs match the sequence input

    Args:
        seq_str: Sequence string
        seq_ident: Percent identity to match
        evalue: E-value of BLAST hit

    Returns:


    """
    pass


def get_fasta(uniprot_id):
    """Get the protein sequence for a UniProt ID as a string.

    Args:
        uniprot_id: Valid UniProt ID

    Returns:
        str: String of the protein (amino acid) sequence

    """
    # Silencing the "Will be moved to Biokit" message
    with ssbio.utils.suppress_stdout():
        return bsup.get_fasta_sequence(uniprot_id)


def uniprot_reviewed_checker(uniprot_id):
    """Check if a single UniProt ID is reviewed or not.

    Args:
        uniprot_id:

    Returns:
        bool: If the entry is reviewed

    """

    query_string = 'id:' + uniprot_id

    uni_rev_raw = StringIO(bsup.search(query_string, columns='id,reviewed', frmt='tab'))
    uni_rev_df = pd.read_table(uni_rev_raw, sep='\t', index_col=0)
    uni_rev_df = uni_rev_df.fillna(False)
    uni_rev_df = uni_rev_df[pd.notnull(uni_rev_df.Status)]

    uni_rev_df = uni_rev_df.replace(to_replace="reviewed", value=True)
    uni_rev_df = uni_rev_df.replace(to_replace="unreviewed", value=False)
    uni_rev_dict_adder = uni_rev_df.to_dict()['Status']

    return uni_rev_dict_adder[uniprot_id]


def uniprot_reviewed_checker_batch(uniprot_ids):
    """Batch check if uniprot IDs are reviewed or not

    Args:
        uniprot_ids: UniProt ID or list of UniProt IDs

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
    """Retrieve the EC number annotation for a UniProt ID.

    Args:
        uniprot_id: Valid UniProt ID

    Returns:

    """
    r = requests.post('http://www.uniprot.org/uniprot/?query=%s&columns=ec&format=tab' % uniprot_id)
    ec = r.content.decode('utf-8').splitlines()[1]

    if len(ec) == 0:
        ec = None

    return ec


def uniprot_sites(uniprot_id):
    """Retrive a dictionary of UniProt sites

    Sites are defined here: http://www.uniprot.org/help/site and here: http://www.uniprot.org/help/function_section

    Args:
        uniprot_id: Valid UniProt ID

    Returns:

    """

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
    """Download a UniProt file for a UniProt ID/ACC

    Args:
        uniprot_id: Valid UniProt ID
        filetype: txt, fasta, xml, rdf, or gff
        outdir: Directory to download the file

    Returns:
        str: Absolute path to file

    """

    my_file = '{}.{}'.format(uniprot_id, filetype)
    url = 'http://www.uniprot.org/uniprot/{}'.format(my_file)
    outfile = op.join(outdir, my_file)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        urlretrieve(url, outfile)

    return outfile


def parse_uniprot_txt_file(infile):
    """Parse a raw UniProt metadata file and return a dictionary.

    Args:
        infile: Path to metadata file

    Returns:
        dict: Metadata dictionaruy

    """
    uniprot_metadata_dict = {}

    metadata = old_parse_uniprot_txt_file(infile)
    metadata_keys = list(metadata.keys())

    if metadata_keys:
        metadata_key = metadata_keys[0]
    else:
        return uniprot_metadata_dict

    uniprot_metadata_dict['seq_len'] = len(str(metadata[metadata_key]['sequence']))
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
        uniprot_metadata_dict['kegg'] = metadata[metadata_key]['kegg']
    if 'ec' in metadata[metadata_key]:
        uniprot_metadata_dict['ec_number'] = metadata[metadata_key]['ec']
    if 'pfam' in metadata[metadata_key]:
        uniprot_metadata_dict['pfam'] = metadata[metadata_key]['pfam']
    if 'pdbs' in metadata[metadata_key]:
        uniprot_metadata_dict['pdbs'] = list(set(metadata[metadata_key]['pdbs']))
    return uniprot_metadata_dict


def old_parse_uniprot_txt_file(infile):
    """From: boscoh/uniprot github
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
                if len(words) > 0:
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
