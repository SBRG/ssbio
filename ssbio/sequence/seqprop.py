from ssbio.core.object import Object
from os import path as op
from Bio import SeqIO
import ssbio.utils
import ssbio.sequence.utils
import ssbio.sequence.utils.fasta
import ssbio.sequence.properties.residues
import ssbio.databases.pdb
from cobra.core import DictList
from copy import deepcopy
from collections import defaultdict
import logging
from slugify import Slugify
import requests
custom_slugify = Slugify(safe_chars='-_.')
log = logging.getLogger(__name__)


class SeqProp(Object):
    def __init__(self, ident, sequence_path=None, metadata_path=None, seq=None, description="<unknown description>",
                 write_fasta_file=False, outfile=None, force_rewrite=False):
        """Store basic protein sequence properties.

        Takes as input a sequence FASTA file, sequence string, Biopython Seq or SeqRecord. You can also provide the
        path to the metadata file if available (see UniProtProp and KEGGProp classes for parsers). If provided as a
        string/Seq/SeqRecord object, you can also write them out as FASTA files.

        Args:
            ident (str): Identifier of the sequence
            sequence_path (str): Path to FASTA file of the sequence
            metadata_path (str): Path to metadata file of the sequence
            seq (str, Seq, SeqRecord): Sequence string, Biopython Seq or SeqRecord object
            description (str): Optional description of the sequence
            write_fasta_file (bool): If a FASTA file of the sequence should be written
            outname:
            outdir:
            force_rewrite:
        """

        Object.__init__(self, id=ident, description=description)

        # Database identifiers
        self.bigg = None
        self.kegg = None
        self.refseq = None
        self.uniprot = None
        self.gene_name = None
        self.pdbs = None

        # Sequence information
        self._seq_record = None
        if seq:
            # Load sequence if not provided as a file
            self.seq_record = ssbio.sequence.utils.cast_to_seq_record(obj=seq, id=ident, description=description)

        # File information
        self._sequence_dir = None
        self.sequence_file = None
        if sequence_path:
            self.load_sequence_path(sequence_path)
        self._metadata_dir = None
        self.metadata_file = None
        if metadata_path:
            self.load_metadata_path(metadata_path)
        if write_fasta_file:
            if not outfile:
                raise ValueError('Output path must be specified if you want to write a FASTA file.')
            self.write_fasta_file(outfile=outfile, force_rerun=force_rewrite)

        # Store AlignIO objects of this sequence to others in alignments
        self.sequence_alignments = DictList()
        self.structure_alignments = DictList()

        # Copy SeqRecord annotations and letter annotations for JSON saving capabilities
        self.annotations = {}
        self.letter_annotations = {}

    @property
    def sequence_dir(self):
        return self._sequence_dir

    @sequence_dir.setter
    def sequence_dir(self, path):
        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        self._sequence_dir = path

    @property
    def sequence_path(self):
        if self.sequence_dir and self.sequence_file:
            path = op.join(self.sequence_dir, self.sequence_file)
            if not op.exists(path):
                raise ValueError('{}: file does not exist'.format(path))
            return path
        else:
            if not self.sequence_dir:
                log.debug('{}: sequence directory not set'.format(self.id))
            if not self.sequence_file:
                log.debug('{}: sequence file not available'.format(self.id))
            return None

    @property
    def seq_record(self):
        # The SeqRecord object is not kept in memory unless explicitly set, so when saving as a json we only save a
        # pointer to the file
        if not self.sequence_path:
            return self._seq_record
            # log.error('{}: sequence file not available'.format(self.id))
            # return None


        return SeqIO.read(open(self.sequence_path), 'fasta')

    @seq_record.setter
    def seq_record(self, sr):
        # Copy SeqRecord annotations and letter annotations for JSON saving capabilities
        self.annotations = sr.annotations
        self.letter_annotations = sr.letter_annotations
        self._seq_record = sr

    @property
    def metadata_dir(self):
        return self._metadata_dir

    @metadata_dir.setter
    def metadata_dir(self, path):
        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        self._metadata_dir = path

    @property
    def metadata_path(self):
        if self.metadata_dir and self.metadata_file:
            path = op.join(self.metadata_dir, self.metadata_file)
            if not op.exists(path):
                raise ValueError('{}: file does not exist'.format(path))
            return path
        else:
            if not self.metadata_dir:
                log.debug('{}: metadata directory not set'.format(self.id))
            if not self.metadata_file:
                log.debug('{}: metadata file not available'.format(self.id))
            return None

    @property
    def seq_str(self):
        """Get the sequence formatted as a string"""
        if not self.seq_record:
            return None
        return ssbio.sequence.utils.cast_to_str(self.seq_record)

    @property
    def seq_len(self):
        """Get the sequence length"""
        if not self.seq_record:
            return 0
        return len(self.seq_record.seq)

    @property
    def num_pdbs(self):
        """Report the number of PDB IDs mapped"""
        if not self.pdbs:
            return 0
        else:
            return len(self.pdbs)

    def load_sequence_path(self, sequence_path):
        """Load a sequence file and provide pointers to its location

        Args:
            sequence_path: Path to sequence file

        """
        self.sequence_dir = op.dirname(sequence_path)
        self.sequence_file = op.basename(sequence_path)

    def load_metadata_path(self, metadata_path):
        """Load a metadata file and provide pointers to its location

        Args:
            metadata_path: Path to metadata file

        """
        self.metadata_dir = op.dirname(metadata_path)
        self.metadata_file = op.basename(metadata_path)

    def write_fasta_file(self, outfile, force_rerun=False):
        """Write a FASTA file for the protein sequence.

        Args:
            outfile (str): Path to new FASTA file to be written to
            force_rerun (bool): If an existing file should be overwritten

        """
        outdir, outname, outext = ssbio.utils.split_folder_and_path(outfile)
        sequence_path = ssbio.sequence.utils.fasta.write_fasta_file(self.seq_record,
                                                                    outname=outname,
                                                                    outdir=outdir,
                                                                    outext=outext,
                                                                    force_rerun=force_rerun)
        # Unset the SeqRecord as it will now be dynamically loaded from the file
        self.seq_record = None
        self.load_sequence_path(sequence_path)

    def equal_to(self, seq_prop):
        """Test if the sequence is equal to another SeqProp object's sequence

        Args:
            seq_prop: SeqProp object

        Returns:
            bool: If the sequences are the same

        """
        if not self.seq_str:
            return False

        if not seq_prop:
            return False

        if not seq_prop.seq_str:
            return False

        return self.seq_str == seq_prop.seq_str

    def equal_to_fasta(self, seq_file):
        """Test if this sequence is equal to another sequence file.

        Args:
            seq_file: Path to another sequence file

        Returns:
            bool: If the sequences are the same

        """
        if not self.sequence_path:
            log.error('Sequence file not available')
            return False
        return ssbio.sequence.utils.fasta.fasta_files_equal(self.sequence_path, seq_file)

    def load_letter_annotations(self, annotation_name, letter_annotation):
        """Add per-residue annotations to the seq_record.letter_annotations attribute

        Args:
            annotation_name: Name of the key that will be used to access this annotation
            letter_annotation (str, list, tuple): Sequence that contains information on the sequence, with the length
                being equal to the protein sequence length

        """
        self.letter_annotations[annotation_name] = letter_annotation
        log.debug('{}: loaded letter_annotations and saved as "{}"'.format(self.id, annotation_name))

    def get_biopython_pepstats(self):
        """Run Biopython's built in ProteinAnalysis module.

        Stores statistics in the seq_record.annotations attribute.
        """
        try:
            pepstats = ssbio.sequence.properties.residues.biopython_protein_analysis(self.seq_str)
        except KeyError as e:
            log.error('{}: unable to run ProteinAnalysis module, unknown amino acid {}'.format(self.id, e))
            return
        self.annotations.update(pepstats)

    def get_emboss_pepstats(self):
        """Run the EMBOSS pepstats program on the protein sequence.

        Stores statistics in the seq_record.annotations attribute.
        Saves a .pepstats file of the results where the sequence file is located.
        """
        outfile = ssbio.sequence.properties.residues.emboss_pepstats_on_fasta(infile=self.sequence_path)
        pepstats = ssbio.sequence.properties.residues.emboss_pepstats_parser(outfile)
        self.annotations.update(pepstats)

    def blast_pdb(self, seq_ident_cutoff=0, evalue=0.0001, display_link=False,
                  outdir=None, force_rerun=False):
        """BLAST this sequence to the PDB"""
        if not outdir:
            outdir = self.sequence_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        if not self.seq_str:
            log.error('{}: no sequence loaded'.format(self.id))
            return None

        try:
            blast_results = ssbio.databases.pdb.blast_pdb(self.seq_str,
                                                          outfile='{}_blast_pdb.xml'.format(custom_slugify(self.id)),
                                                          outdir=outdir,
                                                          force_rerun=force_rerun,
                                                          evalue=evalue,
                                                          seq_ident_cutoff=seq_ident_cutoff,
                                                          link=display_link)
        except requests.ConnectionError as e:
            print(e)
            log.error('{}: BLAST request timed out'.format(seq_prop.id))
            return None

        return blast_results

    def summary(self, as_df):
        """Summarize this sequence.

        Returns:
            dict:

        """
        pass

    def sequence_mutation_summary(self):
        """Summarize all mutations found in the sequence_alignments attribute.

        Returns 2 dictionaries, single_counter and fingerprint_counter.

        single_counter:
            Dictionary of {point mutation: list of genes/strains}
            Example: {('A', 24, 'V'): ['Strain1', 'Strain2', 'Strain4'],
                      ('R', 33, 'T'): ['Strain2']}
            Here, we report which genes/strains have the single point mutation.

        fingerprint_counter:
            Dictionary of {mutation group: list of genes/strains}
            Example: {(('A', 24, 'V'), ('R', 33, 'T')): ['Strain2'],
                      (('A', 24, 'V')): ['Strain1', 'Strain4']}
            Here, we report which genes/strains have the specific combinations (or "fingerprints") of point mutations

        Returns:
            dict, dict: single_counter, fingerprint_counter

        """
        if len(self.sequence_alignments) == 0:
            log.error('{}: no sequence alignments'.format(self.id))
            return {}, {}

        fingerprint_counter = defaultdict(list)
        single_counter = defaultdict(list)

        for alignment in self.sequence_alignments:
            other_sequence = alignment.annotations['b_seq']
            mutations = alignment.annotations['mutations']

            if mutations:
                # Turn this list of mutations into a tuple so it can be a dictionary key
                mutations = tuple(tuple(x) for x in mutations)
                fingerprint_counter[mutations].append(other_sequence)

                for m in mutations:
                    single_counter[m].append(other_sequence)

        return dict(single_counter), dict(fingerprint_counter)

    def get_dict(self, only_keys=None, exclude_attributes=None, df_format=False):
        """Get a copy of all attributes as a dictionary, including object properties

        Args:
            only_keys:
            exclude_attributes:
            df_format:

        Returns:

        """
        my_dict = Object.get_dict(self, only_keys=only_keys, exclude_attributes=exclude_attributes, df_format=df_format)

        additional_keys = ['seq_str', 'seq_len', 'sequence_file', 'metadata_file', 'num_pdbs']
        for k in additional_keys:
            my_dict[k] = deepcopy(getattr(self, k))

        # Choose attributes to return, return everything in the object if a list is not specified
        if not only_keys:
            keys = list(my_dict.keys())
        else:
            keys = ssbio.utils.force_list(only_keys)

        # Remove keys you don't want returned
        if exclude_attributes:
            exclude_attributes = ssbio.utils.force_list(exclude_attributes)
            for x in exclude_attributes:
                if x in keys:
                    keys.remove(x)

        # Copy attributes into a new dictionary
        for k, v in my_dict.items():
            if k in keys:
                if df_format:
                    if v and not isinstance(v, str) and not isinstance(v, int) and not isinstance(v, float) and not isinstance(v, bool):
                        try:
                            my_dict[k] = ssbio.utils.force_string(v)
                        except TypeError:
                            log.warning('{}: excluding attribute from dict, cannot transform into string'.format(k))
                    else:
                        my_dict[k] = v
                else:
                    my_dict[k] = v
        return my_dict

    def __json_decode__(self, **attrs):
        for k, v in attrs.items():
            if k == 'sequence_alignments' or k == 'structure_alignments':
                setattr(self, k, DictList(v))
            else:
                setattr(self, k, v)