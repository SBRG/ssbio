"""
SeqProp
=======
"""

import os.path as op
import pandas as pd
import requests
import logging
from copy import copy, deepcopy
from slugify import Slugify

from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

from ssbio.core.object import Object
import ssbio.utils
import ssbio.databases.pdb
import ssbio.protein.sequence.utils
import ssbio.protein.sequence.utils.fasta
import ssbio.protein.sequence.properties.residues

custom_slugify = Slugify(safe_chars='-_.')
log = logging.getLogger(__name__)


class SeqProp(SeqRecord):

    """Generic class to represent information for a protein sequence.

    Extends the Biopython SeqRecord class. The main functionality added is the ability to set and load directly from
    sequence, metadata, and feature files. Additionally, methods are provided to calculate and store sequence
    properties in the ``annotations`` and ``letter_annotations`` field of a SeqProp. These can then be accessed
    for a range of residue numbers.

    Attributes:
        id (str): Unique identifier for this protein sequence
        seq (Seq): Protein sequence as a Biopython Seq object
        name (str): Optional name for this sequence
        description (str): Optional description for this sequence
        bigg (str, list): BiGG IDs mapped to this sequence
        kegg (str, list): KEGG IDs mapped to this sequence
        refseq (str, list): RefSeq IDs mapped to this sequence
        uniprot (str, list): UniProt IDs mapped to this sequence
        gene_name (str, list): Gene names mapped to this sequence
        pdbs (list): PDB IDs mapped to this sequence
        go (str, list): GO terms mapped to this sequence
        pfam (str, list): PFAMs mapped to this sequence
        ec_number (str, list): EC numbers mapped to this sequence
        sequence_file (str): FASTA file for this sequence
        metadata_file (str): Metadata file (any format) for this sequence
        feature_file (str): GFF file for this sequence
        features (list): List of protein sequence features, which define regions of the protein
        annotations (dict): Annotations of this protein sequence, which summarize global properties
        letter_annotations (RestrictedDict): Residue-level annotations, which describe single residue properties

    Todo:
        - Properly inherit methods from the Object class...

    """

    def __init__(self, seq, id, name='<unknown name>', description='<unknown description>',
                 sequence_path=None, metadata_path=None, feature_path=None):
        """Initialize a SeqProp object with the standard SeqRecord inputs, along with optional paths to files.

        Args:
            seq (str, Seq, SeqRecord): Protein sequence to load
            id (str): Unique identifier for this protein sequence
            name (str): Optional name for this sequence
            description (str): Optional description for this sequence
            sequence_path (str): Absolute or relative path to sequence (FASTA) file
            metadata_path (str): Absolute or relative path to metadata file
            feature_path (str): Absolute or relative path to feature (GFF) file

        """

        # Top level database identifiers
        self.id = id
        self.description = description
        self.notes = {}

        self.bigg = None
        self.kegg = None
        self.refseq = None
        self.uniprot = None
        self.gene_name = None
        self.pdbs = None
        self.go = None
        self.pfam = None
        self.ec_number = None

        # Files, sequences, features
        self._sequence_dir = None
        self.sequence_file = None
        self._metadata_dir = None
        self.metadata_file = None
        self._feature_dir = None
        self.feature_file = None
        self._features = None

        # Stuff from SeqRecord
        self.name = None
        self.description = None
        self.dbxrefs = None
        self.features = None
        self.annotations = None
        self._per_letter_annotations = None

        self._seq = None
        self.seq = seq

        super(SeqProp, self).__init__(seq=self.seq, id=id, name=name, description=description)

        if sequence_path:
            self.sequence_path = sequence_path
        if metadata_path:
            self.metadata_path = metadata_path
        if feature_path:
            self.feature_path = feature_path

    @property
    def seq(self):
        """Seq: Dynamically loaded Seq object from the sequence file"""

        if self.sequence_file:
            file_to_load = copy(self.sequence_path)
            log.debug('{}: reading sequence from sequence file {}'.format(self.id, file_to_load))
            tmp_sr = SeqIO.read(file_to_load, 'fasta')
            return tmp_sr.seq

        else:
            if not self._seq:
                log.debug('{}: no sequence stored in memory'.format(self.id))
            else:
                log.debug('{}: reading sequence from memory'.format(self.id))

            return self._seq

    @seq.setter
    def seq(self, s):
        if not s:
            self._seq = None

        elif self.sequence_file:
            raise ValueError('{}: unable to set sequence, sequence file is associated with this object'.format(self.id))

        elif type(s) == str or type(s) == Seq:
            self._seq = ssbio.protein.sequence.utils.cast_to_seq(obj=s)

        # If a SeqRecord, copy all attributes
        elif type(s) == SeqRecord:
            self._seq = s.seq
            if self.name == '<unknown name>':
                self.name = s.name
            if self.description == '<unknown description>':
                self.description = s.description
            if not self.dbxrefs:
                self.dbxrefs = s.dbxrefs
            if not self.features:
                self.features = s.features
            if not self.annotations:
                self.annotations = s.annotations
            if not self.letter_annotations:
                self.letter_annotations = s.letter_annotations

    @property
    def features(self):
        """list: Get the features stored in memory or in the GFF file"""

        if self.feature_file:
            log.debug('{}: reading features from feature file {}'.format(self.id, self.feature_path))
            with open(self.feature_path) as handle:
                feats = list(GFF.parse(handle))
                if len(feats) > 1:
                    log.warning('Too many sequences in GFF')
                else:
                    return feats[0].features

        else:
            return self._features

    @features.setter
    def features(self, feats):
        if self.feature_file:
            raise ValueError('{}: unable to set features, feature file is associated with this object'.format(self.id))
        elif feats:
            self._features = feats
        else:
            self._features = []

    @property
    def seq_str(self):
        """str: Get the sequence formatted as a string"""

        if not self.seq:
            return None

        return ssbio.protein.sequence.utils.cast_to_str(self.seq)

    @property
    def seq_len(self):
        """int: Get the sequence length"""

        if not self.seq:
            return 0

        return len(self.seq)

    @property
    def num_pdbs(self):
        """int: Report the number of PDB IDs stored in the ``pdbs`` attribute"""

        if not self.pdbs:
            return 0

        else:
            return len(self.pdbs)

    @property
    def sequence_dir(self):
        if not self._sequence_dir:
            raise OSError('No sequence folder set')
        return self._sequence_dir

    @sequence_dir.setter
    def sequence_dir(self, path):
        if path and not op.exists(path):
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

    @sequence_path.setter
    def sequence_path(self, fasta_path):
        """Provide pointers to the paths of the FASTA file

        Args:
            fasta_path: Path to FASTA file

        """
        if not fasta_path:
            self.sequence_dir = None
            self.sequence_file = None

        else:
            if not op.exists(fasta_path):
                raise OSError('{}: file does not exist'.format(fasta_path))

            if not op.dirname(fasta_path):
                self.sequence_dir = '.'
            else:
                self.sequence_dir = op.dirname(fasta_path)
            self.sequence_file = op.basename(fasta_path)

            tmp_sr = SeqIO.read(fasta_path, 'fasta')
            if self.name == '<unknown name>':
                self.name = tmp_sr.name
            if self.description == '<unknown description>':
                self.description = tmp_sr.description
            if not self.dbxrefs:
                self.dbxrefs = tmp_sr.dbxrefs
            if not self.features:
                self.features = tmp_sr.features
            if not self.annotations:
                self.annotations = tmp_sr.annotations
            if not self.letter_annotations:
                self.letter_annotations = tmp_sr.letter_annotations

    @property
    def metadata_dir(self):
        if not self._metadata_dir:
            raise OSError('No metadata folder set')
        return self._metadata_dir

    @metadata_dir.setter
    def metadata_dir(self, path):
        if path and not op.exists(path):
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

    @metadata_path.setter
    def metadata_path(self, m_path):
        """Provide pointers to the paths of the metadata file

        Args:
            m_path: Path to metadata file

        """
        if not m_path:
            self.metadata_dir = None
            self.metadata_file = None

        else:
            if not op.exists(m_path):
                raise OSError('{}: file does not exist!'.format(m_path))

            if not op.dirname(m_path):
                self.metadata_dir = '.'
            else:
                self.metadata_dir = op.dirname(m_path)
            self.metadata_file = op.basename(m_path)

    @property
    def feature_dir(self):
        if not self._feature_dir:
            raise OSError('No feature folder set')
        return self._feature_dir

    @feature_dir.setter
    def feature_dir(self, path):
        if path and not op.exists(path):
            raise OSError('{}: folder does not exist'.format(path))

        self._feature_dir = path

    @property
    def feature_path(self):
        path = op.join(self.feature_dir, self.feature_file)
        if not op.exists(path):
            raise OSError('{}: file does not exist'.format(path))
        return path

    @feature_path.setter
    def feature_path(self, gff_path):
        """Load a GFF file with information on a single sequence and store features in the ``features`` attribute

        Args:
            gff_path: Path to GFF file.

        """
        if not gff_path:
            self.feature_dir = None
            self.feature_file = None

        else:
            if not op.exists(gff_path):
                raise OSError('{}: file does not exist!'.format(gff_path))

            if not op.dirname(gff_path):
                self.feature_dir = '.'
            else:
                self.feature_dir = op.dirname(gff_path)
            self.feature_file = op.basename(gff_path)

    def feature_path_unset(self):
        """Copy features to memory and remove the association of the feature file."""
        if not self.feature_file:
            raise IOError('No feature file to unset')

        with open(self.feature_path) as handle:
            feats = list(GFF.parse(handle))
            if len(feats) > 1:
                log.warning('Too many sequences in GFF')
            else:
                tmp = feats[0].features

        self.feature_dir = None
        self.feature_file = None
        self.features = tmp

    def __str__(self):
        return Object.__str__(self)

    def __repr__(self):
        return Object.__repr__(self)

    def update(self, newdata, overwrite=False, only_keys=None):
        Object.update(self, newdata=newdata, overwrite=overwrite, only_keys=only_keys)

    def get_dict(self, only_attributes=None, exclude_attributes=None, df_format=False):
        """Get a dictionary of this object's attributes. Optional format for storage in a Pandas DataFrame.

        Args:
            only_attributes (str, list): Attributes that should be returned. If not provided, all are returned.
            exclude_attributes (str, list): Attributes that should be excluded.
            df_format (bool): If dictionary values should be formatted for a dataframe
                (everything possible is transformed into strings, int, or float -
                if something can't be transformed it is excluded)

        Returns:
            dict: Dictionary of attributes

        """

        # Choose attributes to return, return everything in the object if a list is not specified
        if not only_attributes:
            keys = list(self.__dict__.keys())
        else:
            keys = ssbio.utils.force_list(only_attributes)

        # Remove keys you don't want returned
        if exclude_attributes:
            exclude_attributes = ssbio.utils.force_list(exclude_attributes)
            for x in exclude_attributes:
                if x in keys:
                    keys.remove(x)

        # Copy attributes into a new dictionary
        df_dict = {}
        for k, orig_v in self.__dict__.items():
            if k in keys:
                v = deepcopy(orig_v)
                if df_format:
                    if v and not isinstance(v, str) and not isinstance(v, int) and not isinstance(v,
                                                                                                  float) and not isinstance(
                            v, bool):
                        try:
                            df_dict[k] = ssbio.utils.force_string(deepcopy(v))
                        except TypeError:
                            log.warning('{}: excluding attribute from dict, cannot transform into string'.format(k))
                    elif not v and not isinstance(v, int) and not isinstance(v, float):
                        df_dict[k] = None
                    else:
                        df_dict[k] = deepcopy(v)
                else:
                    df_dict[k] = deepcopy(v)
        return df_dict

    def save_dataframes(self, outdir, prefix='df_'):
        Object.save_dataframes(self, outdir=outdir, prefix=prefix)

    def save_pickle(self, outfile, protocol=2):
        Object.save_pickle(self, outfile=outfile, protocol=protocol)

    def __json_encode__(self):
        Object.__json_encode__(self)

    def save_json(self, outfile, compression=False):
        Object.save_json(self, outfile=outfile, compression=compression)

    def equal_to(self, seq_prop):
        """Test if the sequence is equal to another SeqProp's sequence

        Args:
            seq_prop: SeqProp object

        Returns:
            bool: If the sequences are the same

        """
        if not self.seq or not seq_prop or not seq_prop.seq:
            return False

        return self.seq == seq_prop.seq

    def write_fasta_file(self, outfile, force_rerun=False):
        """Write a FASTA file for the protein sequence, ``seq`` will now load directly from this file.

        Args:
            outfile (str): Path to new FASTA file to be written to
            force_rerun (bool): If an existing file should be overwritten

        """
        if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
            SeqIO.write(self, outfile, "fasta")

        # The Seq as it will now be dynamically loaded from the file
        self.sequence_path = outfile

    def write_gff_file(self, outfile, force_rerun=False):
        """Write a GFF file for the protein features, ``features`` will now load directly from this file.

        Args:
            outfile (str): Path to new FASTA file to be written to
            force_rerun (bool): If an existing file should be overwritten

        """
        if ssbio.utils.force_rerun(outfile=outfile, flag=force_rerun):
            with open(outfile, "w") as out_handle:
                GFF.write([self], out_handle)

        self.feature_path = outfile

    def add_point_feature(self, resnum, feat_type=None, feat_id=None):
        """Add a feature to the features list describing a single residue.

        Args:
            resnum (int): Protein sequence residue number
            feat_type (str, optional): Optional description of the feature type (ie. 'catalytic residue')
            feat_id (str, optional): Optional ID of the feature type (ie. 'TM1')

        """
        if self.feature_file:
            raise ValueError('Feature file associated with sequence, please remove file association to append '
                             'additional features.')

        if not feat_type:
            feat_type = 'Manually added protein sequence single residue feature'
        newfeat = SeqFeature(location=FeatureLocation(ExactPosition(resnum-1), ExactPosition(resnum)),
                             type=feat_type,
                             id=feat_id)

        self.features.append(newfeat)

    def add_region_feature(self, start_resnum, end_resnum, feat_type=None, feat_id=None):
        """Add a feature to the features list describing a region of the protein sequence.

        Args:
            start_resnum (int): Start residue number of the protein sequence feature
            end_resnum (int): End residue number of the protein sequence feature
            feat_type (str, optional): Optional description of the feature type (ie. 'binding domain')
            feat_id (str, optional): Optional ID of the feature type (ie. 'TM1')

        """
        if self.feature_file:
            raise ValueError('Feature file associated with sequence, please remove file association to append '
                             'additional features.')

        if not feat_type:
            feat_type = 'Manually added protein sequence region feature'
        newfeat = SeqFeature(location=FeatureLocation(start_resnum-1, end_resnum),
                             type=feat_type,
                             id=feat_id)

        self.features.append(newfeat)

    def get_biopython_pepstats(self):
        """Run Biopython's built in ProteinAnalysis module and store statistics in the ``annotations`` attribute."""

        if self.seq:
            try:
                pepstats = ssbio.protein.sequence.properties.residues.biopython_protein_analysis(self.seq)
            except KeyError as e:
                log.error('{}: unable to run ProteinAnalysis module, unknown amino acid {}'.format(self.id, e))
                return
            except ValueError as e:
                log.error('{}: unable to run ProteinAnalysis module, {}'.format(self.id, e))
                return
            self.annotations.update(pepstats)
        else:
            raise ValueError('{}: no sequence available, unable to run ProteinAnalysis'.format(self.id))

    def get_emboss_pepstats(self):
        """Run the EMBOSS pepstats program on the protein sequence.

        Stores statistics in the ``annotations`` attribute.
        Saves a ``.pepstats`` file of the results where the sequence file is located.
        """
        if not self.sequence_file:
            raise IOError('FASTA file needs to be written for EMBOSS pepstats to be run')
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

    def get_residue_annotations(self, start_resnum, end_resnum=None):
        """Retrieve letter annotations for a residue or a range of residues

        Args:
            start_resnum (int): Residue number
            end_resnum (int): Optional residue number, specify if a range is desired

        Returns:
            dict: Letter annotations for this residue or residues

        """
        if not end_resnum:
            end_resnum = start_resnum

        # Create a new SeqFeature
        f = SeqFeature(FeatureLocation(start_resnum - 1, end_resnum))

        # Get sequence properties
        return f.extract(self).letter_annotations

    def get_aggregation_propensity(self, email, password, cutoff_v=5, cutoff_n=5, run_amylmuts=False, outdir=None):
        """Run the AMYLPRED2 web server to calculate the aggregation propensity of this protein sequence, which is
        the number of aggregation-prone segments on the unfolded protein sequence.

        Stores statistics in the ``annotations`` attribute, under the key `aggprop-amylpred`.

        See :mod:`ssbio.protein.sequence.properties.aggregation_propensity` for instructions and details.

        """
        if not outdir:
            outdir = self.sequence_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        import ssbio.protein.sequence.properties.aggregation_propensity as agg

        agg_predictions = agg.AMYLPRED(email=email, password=password)
        result = agg_predictions.get_aggregation_propensity(seq=self, outdir=outdir,
                                                            cutoff_v=cutoff_v, cutoff_n=cutoff_n,
                                                            run_amylmuts=run_amylmuts)
        self.annotations['aggprop-amylpred'] = result

    def get_kinetic_folding_rate(self, secstruct, at_temp=None):
        """Run the FOLD-RATE web server to calculate the kinetic folding rate given an amino acid sequence and its
        structural classficiation (alpha/beta/mixed)

        Stores statistics in the ``annotations`` attribute, under the key `kinetic_folding_rate_<TEMP>-foldrate`.

        See :func:`ssbio.protein.sequence.properties.kinetic_folding_rate.get_foldrate` for instructions and details.

        """

        import ssbio.protein.sequence.properties.kinetic_folding_rate as kfr

        rate = kfr.get_foldrate(seq=self, secstruct=secstruct)
        self.annotations['kinetic_folding_rate_37.0_C-foldrate'] = rate

        if at_temp:
            new_rate = kfr.get_foldrate_at_temp(ref_rate=rate, new_temp=at_temp)
            self.annotations['kinetic_folding_rate_{}_C-foldrate'.format(at_temp)] = new_rate

    def get_thermostability(self, at_temp):
        """Run the thermostability calculator using either the Dill or Oobatake methods.

        Stores calculated (dG, Keq) tuple in the ``annotations`` attribute, under the key
        `thermostability_<TEMP>-<METHOD_USED>`.

        See :func:`ssbio.protein.sequence.properties.thermostability.get_dG_at_T` for instructions and details.

        """

        import ssbio.protein.sequence.properties.thermostability as ts

        dG = ts.get_dG_at_T(seq=self, temp=at_temp)
        self.annotations['thermostability_{}_C-{}'.format(at_temp, dG[2].lower())] = (dG[0], dG[1])