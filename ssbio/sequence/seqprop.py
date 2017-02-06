from ssbio.core.object import Object
from os import path as op
from Bio import SeqIO
import ssbio.utils
import ssbio.sequence.utils.fasta
import ssbio.sequence.properties.residues
from cobra.core import DictList
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from copy import deepcopy
from collections import defaultdict
import logging

log = logging.getLogger(__name__)


class SeqProp(Object):
    """Class for protein sequence properties.

    SeqProp provides:
    1. Attributes for mapping to other database identifiers
        - .bigg: BiGG Models ID
        - .kegg: KEGG ID
        - .uniprot: UniProt ACC/ID
        - .gene_name: Human readable gene name
        - .pdbs: Mapped PDB IDs

    2. Basic sequence attributes
        - .seq_len: Length of the amino acid sequence

    3. Pointers to filepaths
        - .sequence_file: Basename of sequence file in FASTA format
        - .sequence_path: Full path to sequence file
        - .metadata_file: Basename of metadata file in any format
        - .metadata_path: Full path to metadata file

    4. Pointers to


    """

    def __init__(self, ident, description=None, sequence_file=None, metadata_file=None, seq_record=None, seq_str=None,
                 write_fasta_file=False, outname=None, outdir=None, force_rewrite=False):
        """Parses basic sequence properties like identifiers and sequence length.
        Provides paths to sequence and metadata files.

        Args:
            ident: Identifier of the sequence (can be any string)
            sequence_file: Path to FASTA file of the sequence
            metadata_file: Path to metadata file of the sequence
            seq_record: Biopython SeqRecord object
            seq_str: Sequence formatted as a string

        """
        Object.__init__(self, id=ident, description=description)

        self.bigg = None
        self.kegg = None
        self.refseq = None
        self.uniprot = None
        self.gene_name = None
        self.pdbs = None
        self.seq_record = seq_record

        self.sequence_path = sequence_file
        self.metadata_path = metadata_file

        if seq_str:
            self.seq_record = ssbio.sequence.utils.fasta.load_seq_string_as_seqrecord(seq_str, self.id)
        if sequence_file:
            self.load_seq_file(sequence_file)

        if write_fasta_file:
            self.write_fasta_file(outname=outname, outdir=outdir, force_rerun=force_rewrite)

        # Store AlignIO objects of this sequence to others in alignments
        self.sequence_alignments = DictList()
        self.structure_alignments = DictList()

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
    def sequence_file(self):
        """Get the name of the sequence file"""
        if not self.sequence_path:
            return None
        return op.basename(self.sequence_path)

    @property
    def metadata_file(self):
        """Get the name of the metadata file"""
        if not self.metadata_path:
            return None
        return op.basename(self.metadata_path)

    @property
    def num_pdbs(self):
        """Report the number of PDB IDs mapped"""
        if not self.pdbs:
            return 0
        else:
            return len(self.pdbs)

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
            else:
                log.debug('{}: not copying attribute'.format(k))
        return my_dict

    def load_seq_file(self, sequence_file):
        """Load a sequence file and provide pointers to its location

        Args:
            sequence_file: Path to sequence file

        """
        self.sequence_path = sequence_file
        self.seq_record = SeqIO.read(open(sequence_file), 'fasta')

    def load_metadata_file(self, metadata_file):
        """Load a metadata file and provide pointers to its location

        Args:
            metadata_file: Path to metadata file

        """
        self.metadata_path = metadata_file

    def write_fasta_file(self, outname=None, outdir=None, force_rerun=False):
        """Write a FASTA file for the sequence

        Args:
            outname: Name of the FASTA file (without extension)
            outdir: Path to directory to output file to
            force_rerun: If file should be overwritten if existing

        """
        if not outname:
            outname = self.id
        self.sequence_path = ssbio.sequence.utils.fasta.write_fasta_file(self.seq_record,
                                                                         outname=outname,
                                                                         outdir=outdir,
                                                                         force_rerun=force_rerun)

    def equal_to(self, seq_prop):
        """Test if the sequence is equal to another SeqProp object's sequence

        Args:
            seq_prop: SeqProp object

        Returns:
            bool: If the sequences are the same

        """
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
        self.seq_record.letter_annotations[annotation_name] = letter_annotation
        log.debug('{}: loaded letter_annotations and saved as "{}"'.format(self.id, annotation_name))

    def get_biopython_pepstats(self):
        """Run Biopython's built in ProteinAnalysis module.

        Stores statistics in the seq_record.annotations attribute.
        """
        pepstats = ssbio.sequence.properties.residues.biopython_protein_analysis(self.seq_str)
        self.seq_record.annotations.update(pepstats)

    def get_emboss_pepstats(self):
        """Run the EMBOSS pepstats program on the protein sequence.

        Stores statistics in the seq_record.annotations attribute.
        Saves a .pepstats file of the results where the sequence file is located.
        """
        outfile = ssbio.sequence.properties.residues.emboss_pepstats_on_fasta(infile=self.sequence_path)
        pepstats = ssbio.sequence.properties.residues.emboss_pepstats_parser(outfile)
        self.seq_record.annotations.update(pepstats)

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
            return

        fingerprint_counter = defaultdict(list)
        single_counter = defaultdict(list)

        for alignment in self.sequence_alignments:
            other_sequence = alignment.annotations['b_seq']
            mutations = alignment.annotations['mutations']

            if mutations:
                # Turn this list of mutations into a tuple so it can be a dictionary key
                mutations = tuple(mutations)
                fingerprint_counter[mutations].append(other_sequence)

                for m in mutations:
                    single_counter[m].append(other_sequence)

        return dict(single_counter), dict(fingerprint_counter)