import logging
from collections import defaultdict
from copy import deepcopy
from os import path as op
import requests
from Bio import SeqIO
from slugify import Slugify
import ssbio.protein.sequence.properties.residues
import ssbio.protein.sequence.utils
import ssbio.protein.sequence.utils.fasta
import ssbio.databases.pdb
import ssbio.utils
from ssbio.core.object import Object
from BCBio import GFF

custom_slugify = Slugify(safe_chars='-_.')
log = logging.getLogger(__name__)

"""seqprop.py

Todo:
    * Include methods to read and write GFF files so features don't need to be stored in memory
    * load_json for SeqProp objects needs to load letter_annotations as a Bio.SeqRecord._RestrictedDict object, 
    otherwise newly stored annotations can be of any length
    
"""


class SeqProp(Object):
    """Generic class to store information on a protein sequence.

    The main utilities of this class are to:
    
    #. Provide database identifier mappings as top level attributes
    #. Manipulate a sequence as a Biopython SeqRecord
    #. Calculate, store, and access sequence features and annotations
    #. File I/O (sequence and feature files)
    
    Attributes:
        bigg (str): BiGG ID for this protein
        kegg (str): KEGG ID for this protein
        refseq (str): RefSeq ID for this protein
        uniprot (str): UniProt ID for this protein
        gene_name (str): Gene name encoding this protein
        pdbs (list): List of PDB IDs mapped to this protein
        go (list): List of GO terms
        pfam (list): List of PFAMs
        ec_number: EC numbers for this protein
        seq_record (SeqRecord): Biopython ``SeqRecord`` representation of sequence
        sequence_file (str): Path to FASTA file
        metadata_file (str): Path to generic metadata file
        annotations (dict): Freeform dictionary of annotations, copied from any SeqRecord ``annotations``
        letter_annotations (dict): Per-residue annotations, copied from any SeqRecord ``letter_annotations``
        features (list): Sequence features, copied from any SeqRecord ``features``
    
    """
    def __init__(self, ident, description='<unknown description>', seq=None,
                 sequence_path=None, metadata_path=None, feature_path=None,
                 sync_description=True, sync_annotations=True, sync_letter_annotations=True, sync_features=True,
                 write_fasta_file=False, outfile=None, force_rewrite=False):
        """Store basic protein sequence properties.

        Takes as input a sequence FASTA file, sequence string, Biopython Seq or SeqRecord. You can also provide the
        path to the metadata file if available (see UniProtProp and KEGGProp classes for parsers). If provided as a
        string/Seq/SeqRecord object, you can also set write_fasta_file to write them out as FASTA files.

        Args:
            ident (str): Identifier of the sequence
            sequence_path (str): Path to FASTA file of the sequence
            metadata_path (str): Path to metadata file of the sequence
            seq (str, Seq, SeqRecord): Sequence string, Biopython Seq or SeqRecord object
            description (str): Optional description of the sequence
            sync_description (bool): If the ``description`` field from a SeqRecord object should be copied
            sync_annotations (bool): If the ``annotations`` field from a SeqRecord object should be copied
            sync_letter_annotations (bool): If the ``letter_annotations`` field from a SeqRecord object should be copied
            sync_features (bool): If the ``features`` field from a SeqRecord object should be copied
            write_fasta_file (bool): If a FASTA file of the sequence should be written
            outfile (str): Path to output file if write_fasta_file is True
            force_rewrite (bool): If existing outfile should be overwritten
        """

        Object.__init__(self, id=ident, description=description)

        # Sync options
        self.sync_description = sync_description
        self.sync_annotations = sync_annotations
        self.sync_letter_annotations = sync_letter_annotations
        self.sync_features = sync_features

        # Database identifiers
        self.bigg = None
        self.kegg = None
        self.refseq = None
        self.uniprot = None
        self.gene_name = None
        self.pdbs = None
        self.go = None
        self.pfam = None
        self.ec_number = None

        # Annotations and features
        # If metadata file is provided, any existing letter annotations and features are copied
        self.annotations = {}
        self.letter_annotations = {}
        self.features = []

        # Sequence information
        self._seq_record = None
        self._sequence_dir = None
        self.sequence_file = None
        self._metadata_dir = None
        self.metadata_file = None

        # SeqRecord is read directly from a provided file and not kept in memory
        if seq:
            # Load sequence if not provided as a file
            self.seq_record = seq
            self.load_from_seq_record_attributes(self.seq_record)

        # Sequence file information
        if sequence_path:
            self.load_sequence_path(sequence_path)

        # Metadata file information
        if metadata_path:
            self.load_metadata_path(metadata_path)

        # Feature file information
        if feature_path:
            self.load_feature_path(feature_path)

        # Write a standalone fasta file if specified
        if write_fasta_file:
            if not outfile:
                raise ValueError('Output path must be specified if you want to write a FASTA file.')
            self.write_fasta_file(outfile=outfile, force_rerun=force_rewrite)

    @property
    def sequence_dir(self):
        if not self._sequence_dir:
            raise OSError('No sequence folder set')
        return self._sequence_dir

    @sequence_dir.setter
    def sequence_dir(self, path):
        if not op.exists(path):
            raise OSError('{}: folder does not exist'.format(path))

        self._sequence_dir = path

    @property
    def sequence_path(self):
        if not self.sequence_file:
            raise OSError('Sequence file not loaded')

        path = op.join(self.sequence_dir, self.sequence_file)
        if not op.exists(path):
            raise OSError('{}: file does not exist'.format(path))
        return path

    @property
    def metadata_dir(self):
        if not self._metadata_dir:
            raise OSError('No metadata folder set')
        return self._metadata_dir

    @metadata_dir.setter
    def metadata_dir(self, path):
        if not op.exists(path):
            raise OSError('{}: folder does not exist'.format(path))

        self._metadata_dir = path

    @property
    def metadata_path(self):
        if not self.metadata_file:
            raise OSError('Metadata file not loaded')

        path = op.join(self.metadata_dir, self.metadata_file)
        if not op.exists(path):
            raise OSError('{}: file does not exist'.format(path))
        return path

    @property
    def feature_dir(self):
        if not self._feature_dir:
            raise OSError('No feature folder set')
        return self._feature_dir

    @feature_dir.setter
    def feature_dir(self, path):
        if not op.exists(path):
            raise OSError('{}: folder does not exist'.format(path))

        self._feature_dir = path

    @property
    def feature_path(self):
        path = op.join(self.feature_dir, self.feature_file)
        if not op.exists(path):
            raise OSError('{}: file does not exist'.format(path))
        return path

    @property
    def seq_record(self):
        """SeqRecord: Dynamically loaded SeqRecord object from the sequence or metadata file"""
        # If there was no sequence file provided, return the provided sequence, if any
        if not self.sequence_file:
            sr = self._seq_record
        # Otherwise read directly from the FASTA file
        else:
            sr = SeqIO.read(self.sequence_path, 'fasta')

        if sr:
            # And associate attributes to the SeqRecord
            self.load_to_seq_record_attributes(sr)

        return sr

    @seq_record.setter
    def seq_record(self, seq):
        # Only used when a sequence is manually provided as a string, Seq, or SeqRecord
        if seq:
            self._seq_record = ssbio.protein.sequence.utils.cast_to_seq_record(obj=seq,
                                                                               id=self.id,
                                                                               description=self.description)
        else:
            self._seq_record = None

    @property
    def seq_str(self):
        """str: Get the sequence formatted as a string"""
        if not self.seq_record:
            return None
        return ssbio.protein.sequence.utils.cast_to_str(self.seq_record)

    @property
    def seq_len(self):
        """int: Get the sequence length"""
        if not self.seq_record:
            return 0
        return len(self.seq_record.seq)

    @property
    def num_pdbs(self):
        """int: Report the number of PDB IDs stored in the ``pdbs`` attribute"""
        if not self.pdbs:
            return 0
        else:
            return len(self.pdbs)

    def load_to_seq_record_attributes(self, sr):
        """Load specific attributes from the SeqProp to a SeqRecord object
        """
        if self.sync_description:
            if self.description != '<unknown description>':
                sr.description = self.description
        if self.sync_annotations:
            sr.annotations = self.annotations
        if self.sync_letter_annotations:
            sr.letter_annotations = self.letter_annotations
        if self.sync_features:
            sr.features = self.features

    def load_from_seq_record_attributes(self, sr):
        """Load specific attributes from a SeqRecord object into a SeqProp.
        """
        if self.sync_description:
            if sr.description != '<unknown description>':
                self.description = sr.description
        if self.sync_annotations:
            self.annotations = sr.annotations
        if self.sync_letter_annotations:
            self.letter_annotations = sr.letter_annotations
        if self.sync_features:
            self.features = sr.features

    def load_sequence_path(self, sequence_path):
        """Provide pointers to the paths of the sequence file

        Args:
            sequence_path: Path to sequence file

        """
        if not op.exists(sequence_path):
            raise OSError('{}: file does not exist'.format(sequence_path))

        if not op.dirname(sequence_path):
            self.sequence_dir = '.'
        else:
            self.sequence_dir = op.dirname(sequence_path)
        self.sequence_file = op.basename(sequence_path)

        # Load the SeqRecord so attributes are copied
        tmp_loader = self.seq_record
        self.load_from_seq_record_attributes(tmp_loader)

    def load_metadata_path(self, metadata_path):
        """Provide pointers to the paths of the metadata file

        Args:
            metadata_path: Path to metadata file

        """
        if not op.exists(metadata_path):
            raise OSError('{}: file does not exist!'.format(metadata_path))

        if not op.dirname(metadata_path):
            self.metadata_dir = '.'
        else:
            self.metadata_dir = op.dirname(metadata_path)
        self.metadata_file = op.basename(metadata_path)

    def load_feature_path(self, gff_path):
        """Load a GFF file with information on a single sequence and store features in the ``features`` attribute
        
        Args:
            gff_path: Path to GFF file.

        """
        if not op.exists(gff_path):
            raise OSError('{}: file does not exist!'.format(gff_path))

        if not op.dirname(gff_path):
            self.feature_dir = '.'
        else:
            self.feature_dir = op.dirname(gff_path)
        self.feature_file = op.basename(gff_path)

        with open(self.feature_path) as handle:
            feats = list(GFF.parse(handle))
            if len(feats) > 1:
                log.warning('Too many sequences in GFF')
            else:
                self.features = feats[0].features

    def write_fasta_file(self, outfile, force_rerun=False):
        """Write a FASTA file for the protein sequence

        If sequence was initialized in the object, it is cleared and ``seq_record`` will instead load from this file.

        Args:
            outfile (str): Path to new FASTA file to be written to
            force_rerun (bool): If an existing file should be overwritten

        """
        outdir, outname, outext = ssbio.utils.split_folder_and_path(outfile)
        sequence_path = ssbio.protein.sequence.utils.fasta.write_fasta_file(self.seq_record,
                                                                            outname=outname,
                                                                            outdir=outdir,
                                                                            outext=outext,
                                                                            force_rerun=force_rerun)
        # Unset the SeqRecord as it will now be dynamically loaded from the file
        self.seq_record = None
        self.load_sequence_path(sequence_path)

    def equal_to(self, seq_prop):
        """Test if the sequence is equal to another SeqProp's sequence

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
        return ssbio.protein.sequence.utils.fasta.fasta_files_equal(self.sequence_path, seq_file)

    def get_biopython_pepstats(self):
        """Run Biopython's built in ProteinAnalysis module.

        Stores statistics in the ``annotations`` attribute.
        """
        try:
            pepstats = ssbio.protein.sequence.properties.residues.biopython_protein_analysis(self.seq_str)
        except KeyError as e:
            log.error('{}: unable to run ProteinAnalysis module, unknown amino acid {}'.format(self.id, e))
            return
        self.annotations.update(pepstats)

    def get_emboss_pepstats(self):
        """Run the EMBOSS pepstats program on the protein sequence.

        Stores statistics in the ``annotations`` attribute.
        Saves a ``.pepstats`` file of the results where the sequence file is located.
        """
        outfile = ssbio.protein.sequence.properties.residues.emboss_pepstats_on_fasta(infile=self.sequence_path)
        pepstats = ssbio.protein.sequence.properties.residues.emboss_pepstats_parser(outfile)
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
            log.error('{}: BLAST request timed out'.format(self.id))
            print(e)
            return None

        return blast_results

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