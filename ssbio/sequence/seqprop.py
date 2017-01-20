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
from collections import defaultdict
import logging
log = logging.getLogger(__name__)


class SeqProp(Object):
    """Class for protein sequence properties"""

    def __init__(self, ident, description=None, sequence_file=None, metadata_file=None, seq_record=None, seq_str=None):
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

        self.sequence_len = 0
        self.sequence_file = None
        self.sequence_path = sequence_file
        if sequence_file:
            self.load_seq_file(sequence_file)

        self.metadata_file = None
        self.metadata_path = metadata_file
        if metadata_file:
            self.load_metadata_file(metadata_file)

        self.seq_record = None
        self.seq_str = None
        if seq_record:
            self.load_seq_record(seq_record)
        if seq_str:
            self.load_seq_str(seq_str)

        # Store AlignIO objects of this sequence to others in alignments
        self.sequence_alignments = DictList()
        self.structure_alignments = DictList()

    def load_seq_record(self, seq_record):
        """Load a SeqRecord

        Args:
            seq_record: Sequence string

        """
        seq_record.id = self.id

        self.sequence_len = len(seq_record.seq)
        self.seq_record = seq_record
        self.seq_str = str(seq_record.seq)

    def load_seq_str(self, seq_str):
        """Load a sequence string

        Args:
            seq_str: Sequence string

        """
        self.sequence_len = len(seq_str)
        self.seq_record = SeqRecord(Seq(seq_str, IUPAC.extended_protein), id=self.id)
        self.seq_str = str(self.seq_record.seq)

    def load_seq_file(self, sequence_file):
        """Load a sequence file and provide pointers to its location

        Args:
            sequence_file: Path to sequence file

        """
        self.sequence_len = len(SeqIO.read(open(sequence_file), "fasta"))
        self.sequence_file = op.basename(sequence_file)
        self.sequence_path = sequence_file

    def load_metadata_file(self, metadata_file):
        """Load a metadata file and provide pointers to its location

        Args:
            metadata_file: Path to metadata file

        """
        self.metadata_file = op.basename(metadata_file)
        self.metadata_path = metadata_file

    def get_seq_str(self):
        """Get the sequence formatted as a string

        Returns:
            str: Amino acid sequence in string format

        """
        if self.seq_str:
            return self.seq_str
        elif self.seq_record:
            self.seq_str = str(self.seq_record.seq)
            return self.seq_str
        elif self.get_seq_record():
            return self.seq_str
        else:
            log.debug('{}: no sequence loaded, cannot get sequence string'.format(self.id))
            return None

    def get_seq_record(self):
        """Get the sequence as a SeqRecord object

        Returns:
            SeqRecord: Amino acid sequence as a SeqRecord object

        """
        if self.seq_record:
            return self.seq_record
        elif self.sequence_path:
            self.seq_record = SeqIO.read(open(self.sequence_path), "fasta")
            self.seq_str = str(self.seq_record.seq)
            return self.seq_record
        else:
            log.debug('{}: no sequence loaded, cannot get seqrecord'.format(self.id))
            return None

    def equal_to(self, seq_prop):
        """Test if the sequence is equal to another SeqProp object's sequence

        Args:
            seq_prop: SeqProp object

        Returns:
            bool: If the sequences are the same

        """
        if not self.sequence_path or not seq_prop.sequence_path:
            log.error('Sequence file not available')
            return False
        return ssbio.sequence.utils.fasta.sequences_equal(self.sequence_path, seq_prop.sequence_path)

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
        return ssbio.sequence.utils.fasta.sequences_equal(self.sequence_path, seq_file)

    def num_pdbs(self):
        """Report the number of PDBs mapped
        """

        if not self.pdbs:
            return 0
        else:
            return len(self.pdbs)

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
        pepstats = ssbio.sequence.properties.residues.biopython_protein_analysis(self.get_seq_str())
        self.seq_record.annotations.update(pepstats)

    def get_emboss_pepstats(self):
        """Run the EMBOSS pepstats program on the protein sequence.

        Stores statistics in the seq_record.annotations attribute.
        Saves a .pepstats file of the results where the sequence file is located.
        """
        outfile = ssbio.sequence.properties.residues.emboss_pepstats_on_fasta(infile=self.sequence_path)
        pepstats = ssbio.sequence.properties.residues.emboss_pepstats_parser(outfile)
        self.seq_record.annotations.update(pepstats)

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