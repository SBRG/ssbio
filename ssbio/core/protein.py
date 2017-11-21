import logging
import requests
import shutil
import pandas as pd
import numpy as np
import os.path as op
import seaborn as sns
from slugify import Slugify
from collections import defaultdict
from six.moves.urllib.error import URLError

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.PDB.PDBExceptions import PDBException, PDBConstructionException
from msgpack.exceptions import ExtraData

from cobra.core import DictList
import ssbio.utils
import ssbio.databases.pdb
import ssbio.protein.sequence.utils.alignment
import ssbio.protein.sequence.utils.fasta
import ssbio.protein.structure.properties.quality
from ssbio.core.object import Object
from ssbio.protein.sequence.seqprop import SeqProp
from ssbio.databases.kegg import KEGGProp
from ssbio.databases.uniprot import UniProtProp
from ssbio.protein.structure.structprop import StructProp
from ssbio.databases.pdb import PDBProp
from ssbio.protein.structure.homology.itasser.itasserprop import ITASSERProp
from ssbio.protein.structure.homology.itasser.itasserprep import ITASSERPrep


custom_slugify = Slugify(safe_chars='-_.')
log = logging.getLogger(__name__)


class Protein(Object):

    """Store information on a protein that represents the translated unit of a gene.

    The main utilities of this class are to:
    
    #. Load, parse, and store multiple sources of the same or similar (ie. from different strains) \
        protein sequences as ``SeqProp`` objects in the ``sequences`` attribute
    #. Load, parse, and store multiple experimental or predicted protein structures as ``StructProp`` \
        objects in the ``structures`` attribute
    #. Calculate, store, and access sequence alignments to stored sequences or structures
    #. Provide summaries of alignments and mutations seen
    #. Map between residue numbers of sequences and structures

    Args:
        ident (str): Unique identifier for this protein
        description (str): Optional description for this protein
        root_dir (str): Path to where the folder named by this protein's ID will be created.
            Default is current working directory.
        pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` -
            choose a file type for files downloaded from the PDB

    Todo:
        - Implement structural alignment objects
    
    """

    __representative_sequence_attributes = ['gene', 'uniprot', 'kegg', 'pdbs',
                                            'sequence_path', 'metadata_path']
    __representative_structure_attributes = ['is_experimental', 'reference_seq_top_coverage', 'date', 'description',
                                             'resolution','taxonomy_name']

    def __init__(self, ident, description=None, root_dir=None, pdb_file_type='mmtf'):
        """Initialize a Protein object.

        A Protein contains sequences, structures, and a single representative sequence and structure.

        """
        Object.__init__(self, id=ident, description=description)

        self.pdb_file_type = pdb_file_type
        """str: ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` - choose a file 
        type for files downloaded from the PDB"""

        # Create directories
        self._root_dir = None
        if root_dir:
            self.root_dir = root_dir

        # Sequences
        self.sequences = DictList()
        """DictList: Stored amino acids which are related to this protein"""
        self.representative_sequence = None
        """SeqProp: Sequence set to represent this protein"""

        # Structures
        self.structures = DictList()
        """DictList: Stored protein structures which are related to this protein"""
        self.representative_structure = None
        """StructProp: Structure set to represent this protein, optionally in monomeric form"""
        self.representative_chain = None
        """str: Chain ID in the representative structure which best represents a sequence"""
        self.representative_chain_seq_coverage = 0
        """float: Percent identity of sequence coverage for the representative chain"""

        # Alignments
        self.sequence_alignments = DictList()
        """DictList: Pairwise or multiple sequence alignments stored as ``Bio.Align.MultipleSeqAlignment`` objects"""
        self.structure_alignments = DictList()
        """DictList: Pairwise or multiple structure alignments - currently a placeholder"""

    @property
    def root_dir(self):
        """str: Path to where the folder named by this protein's ID will be created. Default is current working
        directory."""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.debug('Changing root directory of Protein "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)) and not op.exists(op.join(path, self.id) + '_protein'):
                raise IOError('Protein "{}" does not exist in folder {}'.format(self.id, path))

        self._root_dir = path

        for d in [self.protein_dir, self.sequence_dir, self.structure_dir]:
            ssbio.utils.make_dir(d)

        # Propagate changes to SeqProps
        if hasattr(self, 'sequences'):
            for seq in self.sequences:
                seq.sequence_dir = self.sequence_dir
                seq.metadata_dir = self.sequence_dir
        if hasattr(self, 'representative_sequence'):
            if self.representative_sequence:
                self.representative_sequence.sequence_dir = self.sequence_dir
                self.representative_sequence.metadata_dir = self.sequence_dir

        # Propagate changes to StructProps
        if hasattr(self, 'structures'):
            for struct in self.structures:
                struct.structure_dir = self.structure_dir
        if hasattr(self, 'representative_structure'):
            if self.representative_structure:
                self.representative_structure.structure_dir = self.structure_dir

    @property
    def protein_dir(self):
        """str: Protein folder"""
        if self.root_dir:
            # Add a _protein suffix to the folder if it has the same name as its root folder
            folder_name = self.id
            if folder_name == op.basename(self.root_dir):
                folder_name = self.id + '_protein'
            return op.join(self.root_dir, folder_name)
        else:
            log.warning('Root directory not set')
            return None

    @property
    def sequence_dir(self):
        """str: Directory where sequence related files are stored"""
        if self.root_dir:
            return op.join(self.protein_dir, 'sequences')
        else:
            log.debug('Root directory not set')
            return None

    @property
    def structure_dir(self):
        """str: Directory where structure related files are stored"""
        if self.root_dir:
            return op.join(self.protein_dir, 'structures')
        else:
            log.debug('Root directory not set')
            return None

    @property
    def protein_statistics(self):
        """Get a dictionary of basic statistics describing this protein"""
        d = {}

        d['id'] = self.id
        d['sequences'] = [x.id for x in self.sequences]
        d['num_sequences'] = self.num_sequences
        if self.representative_sequence:
            d['representative_sequence'] = self.representative_sequence.id
        d['num_structures'] = self.num_structures
        d['experimental_structures'] = [x.id for x in self.get_experimental_structures()]
        d['num_experimental_structures'] = self.num_structures_experimental
        d['homology_models'] = [x.id for x in self.get_homology_models()]
        d['num_homology_models'] = self.num_structures_homology
        if self.representative_structure:
            d['representative_structure'] = self.representative_structure.id
            d['representative_chain'] = self.representative_chain
            d['representative_chain_seq_coverage'] = self.representative_chain_seq_coverage
        d['num_sequence_alignments'] = len(self.sequence_alignments)
        d['num_structure_alignments'] = len(self.structure_alignments)

        return d

    @property
    def num_sequences(self):
        """int: Return the total number of sequences"""
        return len(self.sequences)

    @property
    def num_structures(self):
        """int: Return the total number of structures"""
        return len(self.structures)

    @property
    def num_structures_experimental(self):
        """int: Return the total number of experimental structures"""
        return len(self.get_experimental_structures())

    @property
    def num_structures_homology(self):
        """int: Return the total number of homology models"""
        return len(self.get_homology_models())

    def get_experimental_structures(self):
        """DictList: Return a DictList of all experimental structures in self.structures"""
        return DictList(x for x in self.structures if x.is_experimental)

    def get_homology_models(self):
        """DictList: Return a DictList of all homology models in self.structures"""
        return DictList(x for x in self.structures if not x.is_experimental)

    def filter_sequences(self, seq_type):
        """Return a DictList of only specified types in the sequences attribute.

        Args:
            seq_type (SeqProp): Object type

        Returns:
            DictList: A filtered DictList of specified object type only

        """
        return DictList(x for x in self.sequences if isinstance(x, seq_type))

    def load_kegg(self, kegg_id, kegg_organism_code=None, kegg_seq_file=None, kegg_metadata_file=None,
                  set_as_representative=False, download=False, outdir=None, force_rerun=False):
        """Load a KEGG ID, sequence, and metadata files into the sequences attribute.

        Args:
            kegg_id (str): KEGG ID
            kegg_organism_code (str): KEGG organism code to prepend to the kegg_id if not part of it already.
                Example: ``eco:b1244``, eco is the organism code
            kegg_seq_file (str): Path to KEGG FASTA file
            kegg_metadata_file (str): Path to KEGG metadata file (raw KEGG format)
            set_as_representative (bool): If this KEGG ID should be set as the representative sequence
            download (bool): If the KEGG sequence and metadata files should be downloaded if not provided
            outdir (str): Where the sequence and metadata files should be downloaded to
            force_rerun (bool): If ID should be reloaded and files redownloaded

        Returns:
            KEGGProp: object contained in the sequences attribute

        """
        if download:
            if not outdir:
                outdir = self.sequence_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')

        if kegg_organism_code:
            kegg_id = kegg_organism_code + ':' + kegg_id

        # If we have already loaded the KEGG ID
        if self.sequences.has_id(kegg_id):
            # Remove it if we want to force rerun things
            if force_rerun:
                existing = self.sequences.get_by_id(kegg_id)
                self.sequences.remove(existing)
            # Otherwise just get that KEGG object
            else:
                log.debug('{}: KEGG ID already present in list of sequences'.format(kegg_id))
                kegg_prop = self.sequences.get_by_id(kegg_id)

        # Check again (instead of else) in case we removed it if force rerun
        if not self.sequences.has_id(kegg_id):
            kegg_prop = KEGGProp(id=kegg_id, seq=None, fasta_path=kegg_seq_file, txt_path=kegg_metadata_file)
            if download:
                kegg_prop.download_seq_file(outdir, force_rerun)
                kegg_prop.download_metadata_file(outdir, force_rerun)

            # Check if KEGG sequence matches a potentially set representative sequence
            # Do not add any info if a UniProt ID was already mapped though, we want to use that
            if self.representative_sequence:
                if not self.representative_sequence.uniprot:
                    if kegg_prop.equal_to(self.representative_sequence):
                        # Update the representative sequence field with KEGG metadata
                        self.representative_sequence.update(kegg_prop.get_dict(), only_keys=['sequence_path',
                                                                                             'metadata_path',
                                                                                             'kegg',
                                                                                             'description',
                                                                                             'taxonomy',
                                                                                             'id',
                                                                                             'pdbs',
                                                                                             'uniprot',
                                                                                             'seq_record',
                                                                                             'gene_name',
                                                                                             'refseq'])
                    else:
                        # TODO: add option to use manual or kegg sequence if things do not match
                        log.warning('{}: representative sequence does not match mapped KEGG sequence.'.format(self.id))

            self.sequences.append(kegg_prop)

        if set_as_representative:
            self.representative_sequence = kegg_prop

        return self.sequences.get_by_id(kegg_id)

    def load_uniprot(self, uniprot_id, uniprot_seq_file=None, uniprot_xml_file=None, download=False, outdir=None,
                     set_as_representative=False, force_rerun=False):
        """Load a UniProt ID and associated sequence/metadata files into the sequences attribute.
        
        Sequence and metadata files can be provided, or alternatively downloaded with the download flag set to True.
            Metadata files will be downloaded as XML files.
        
        Args:
            uniprot_id (str): UniProt ID/ACC
            uniprot_seq_file (str): Path to FASTA file
            uniprot_xml_file (str): Path to UniProt XML file
            download (bool): If sequence and metadata files should be downloaded
            outdir (str): Output directory for sequence and metadata files
            set_as_representative (bool): If this sequence should be set as the representative one
            force_rerun (bool): If files should be redownloaded and metadata reloaded

        Returns:
            UniProtProp: sequence that was loaded into the sequences attribute

        """
        if download:
            if not outdir:
                outdir = self.sequence_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')

        # If we have already loaded the UniProt ID
        if self.sequences.has_id(uniprot_id):
            # Remove it if we want to force rerun things
            if force_rerun:
                existing = self.sequences.get_by_id(uniprot_id)
                self.sequences.remove(existing)
            # Otherwise just get that UniProt object
            else:
                log.debug('{}: UniProt ID already present in list of sequences'.format(uniprot_id))
                uniprot_prop = self.sequences.get_by_id(uniprot_id)

        if not self.sequences.has_id(uniprot_id):
            uniprot_prop = UniProtProp(id=uniprot_id,
                                       seq=None,
                                       fasta_path=uniprot_seq_file,
                                       xml_path=uniprot_xml_file)
            if download:
                uniprot_prop.download_metadata_file(outdir=outdir,
                                                    force_rerun=force_rerun)
                uniprot_prop.download_seq_file(outdir=outdir, force_rerun=force_rerun)

            # Also check if UniProt sequence matches a potentially set representative sequence
            if self.representative_sequence:
                # Test equality
                if uniprot_prop.equal_to(self.representative_sequence):
                    # Update the representative sequence field with UniProt metadata if it is the same
                    self.representative_sequence.update(uniprot_prop.get_dict(), only_keys=['sequence_path',
                                                                                            'metadata_path',
                                                                                            'kegg',
                                                                                            'description',
                                                                                            'taxonomy',
                                                                                            'id',
                                                                                            'pdbs',
                                                                                            'uniprot',
                                                                                            'seq_record',
                                                                                            'gene_name',
                                                                                            'refseq'])
                else:
                    # TODO: add option to use manual or uniprot sequence if things do not match
                    log.warning('{}: representative sequence does not match mapped UniProt sequence'.format(self.id))
            self.sequences.append(uniprot_prop)

        if set_as_representative:
            self.representative_sequence = uniprot_prop

        return self.sequences.get_by_id(uniprot_id)

    def load_manual_sequence_file(self, ident, seq_file, copy_file=False, outdir=None, set_as_representative=False):
        """Load a manual sequence given as a FASTA file and optionally set it as the representative sequence.
        Also store it in the sequences attribute.

        Args:
            ident:
            seq_file:
            copy_file:
            outdir:
            set_as_representative:

        """
        if copy_file:
            if not outdir:
                outdir = self.sequence_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')
            shutil.copy(seq_file, outdir)
            seq_file = op.join(outdir, seq_file)

        manual_sequence = SeqProp(id=ident, sequence_path=seq_file, seq=None)
        self.sequences.append(manual_sequence)

        if set_as_representative:
            self.representative_sequence = manual_sequence

        return self.sequences.get_by_id(ident)

    def load_manual_sequence(self, seq, ident=None, write_fasta_file=False, outname=None, outdir=None,
                             set_as_representative=False, force_rewrite=False):
        """Load a manual sequence given as a string and optionally set it as the representative sequence.
        Also store it in the sequences attribute.

        Args:
            seq (str, Seq, SeqRecord): Sequence string, Biopython Seq or SeqRecord object
            ident (str): Optional identifier for the sequence, required if seq is a string. Also will override existing
                IDs in Seq or SeqRecord objects if set.
            write_fasta_file:
            outname:
            outdir:
            set_as_representative:
            force_rewrite:

        Returns:

        """

        if not outname:
            outname = ident

        if write_fasta_file:
            if not outdir:
                outdir = self.sequence_dir
                if not outdir:
                    raise ValueError('Output directory must be specified')
            outfile = op.join(outdir, '{}.faa'.format(outname))
        else:
            outfile = None

        if isinstance(seq, str) or isinstance(seq, Seq):
            if not ident:
                raise ValueError('ID must be specified if sequence is a string or Seq object')
        else:
            if not ident:
                # Use ID from SeqRecord ID if new ID not provided
                ident = seq.id
            else:
                # Overwrite SeqRecord ID with new ID if provided
                seq.id = ident

        manual_sequence = SeqProp(id=ident, seq=seq)
        if write_fasta_file:
            manual_sequence.write_fasta_file(outfile=outfile, force_rerun=force_rewrite)
        self.sequences.append(manual_sequence)

        if set_as_representative:
            self.representative_sequence = manual_sequence

        return self.sequences.get_by_id(ident)

    def set_representative_sequence(self, force_rerun=False):
        """Automatically consolidate loaded sequences (manual, UniProt, or KEGG) and set a single representative sequence.

        Manually set representative sequences override all existing mappings. UniProt mappings override KEGG mappings
        except when KEGG mappings have PDBs associated with them and UniProt doesn't.

        Args:
            force_rerun (bool): Set to True to recheck stored sequences

        Returns:
            SeqProp: Which sequence was set as representative

        """

        if len(self.sequences) == 0:
            log.error('{}: no sequences mapped'.format(self.id))
            return self.representative_sequence

        kegg_mappings = self.filter_sequences(KEGGProp)
        if len(kegg_mappings) > 0:
            kegg_to_use = kegg_mappings[0]
            if len(kegg_mappings) > 1:
                log.warning('{}: multiple KEGG mappings found, using the first entry {}'.format(self.id, kegg_to_use.id))

        uniprot_mappings = self.filter_sequences(UniProtProp)

        # If a representative sequence has already been set, nothing needs to be done
        if self.representative_sequence and not force_rerun:
            log.debug('{}: representative sequence already set'.format(self.id))

        # If there is a KEGG annotation and no UniProt annotations, set KEGG as representative
        elif len(kegg_mappings) > 0 and len(uniprot_mappings) == 0:
            self.representative_sequence = kegg_to_use
            log.debug('{}: representative sequence set from KEGG ID {}'.format(self.id, kegg_to_use.id))

        # If there are UniProt annotations and no KEGG annotations, set UniProt as representative
        elif len(kegg_mappings) == 0 and len(uniprot_mappings) > 0:
            # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
            # This way, UniProts with PDBs get ranked to the top, or if no PDBs, reviewed entries
            u_ranker = []
            for u in uniprot_mappings:
                u_ranker.append((u.id, u.ranking_score()))
            sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
            best_u_id = sorted_by_second[0][0]

            best_u = uniprot_mappings.get_by_id(best_u_id)
            self.representative_sequence = best_u
            log.debug('{}: representative sequence set from UniProt ID {}'.format(self.id, best_u_id))

        # If there are both UniProt and KEGG annotations...
        elif len(kegg_mappings) > 0 and len(uniprot_mappings) > 0:
            # Use KEGG if the mapped UniProt is unique, and it has PDBs
            if kegg_to_use.num_pdbs > 0 and not uniprot_mappings.has_id(kegg_to_use.uniprot):
                self.representative_sequence = kegg_to_use
                log.debug('{}: representative sequence set from KEGG ID {}'.format(self.id, kegg_to_use.id))
            else:
                # If there are multiple uniprots rank them by the sum of reviewed (bool) + num_pdbs
                u_ranker = []
                for u in uniprot_mappings:
                    u_ranker.append((u.id, u.ranking_score()))
                sorted_by_second = sorted(u_ranker, key=lambda tup: tup[1], reverse=True)
                best_u_id = sorted_by_second[0][0]

                best_u = uniprot_mappings.get_by_id(best_u_id)
                self.representative_sequence = best_u
                log.debug('{}: representative sequence set from UniProt ID {}'.format(self.id, best_u_id))

        return self.representative_sequence

    def pairwise_align_sequences_to_representative(self, gapopen=10, gapextend=0.5, outdir=None,
                                                   engine='needle', parse=True, force_rerun=False):
        """Align all sequences in the sequences attribute to the representative sequence.

        Stores the alignments the ``sequence_alignments`` DictList

        Args:
            gapopen (int): Only for `needle` - Gap open penalty is the score taken away when a gap is created
            gapextend (float): Only for `needle` - Gap extension penalty is added to the standard gap penalty for each
                base or residue in the gap
            outdir (str): Only for `needle` - Path to output directory. Default is the protein sequence directory.
            engine (str): `biopython` or `needle` - which pairwise alignment program to use
            parse (bool): Store locations of mutations, insertions, and deletions in the alignment object (as an
                annotation)
            force_rerun (bool): Only for `needle` - Default False, set to True if you want to rerun the alignment
                if outfile exists.

        """

        if not outdir:
            outdir = self.sequence_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        for seq in self.sequences:
            aln_id = '{}_{}'.format(self.id, seq.id)
            outfile = '{}.needle'.format(aln_id)

            if not self.representative_sequence:
                log.error('{}: no representative sequence set, skipping gene'.format(self.id))
                continue

            if self.sequence_alignments.has_id(aln_id):
                log.debug('{}: alignment already completed'.format(seq.id))
                continue

            if not seq.seq_str:
                log.error('{}: no sequence stored, skipping alignment'.format(seq.id))
                continue

            # Don't need to compare sequence to itself
            if seq.id == self.representative_sequence.id:
                continue

            aln = ssbio.protein.sequence.utils.alignment.pairwise_sequence_alignment(a_seq=self.representative_sequence.seq_str,
                                                                                     a_seq_id=self.id,
                                                                                     b_seq=seq.seq_str,
                                                                                     b_seq_id=seq.id,
                                                                                     gapopen=gapopen, gapextend=gapextend,
                                                                                     engine=engine,
                                                                                     outdir=outdir,
                                                                                     outfile=outfile,
                                                                                     force_rerun=force_rerun)
            # Add an identifier to the MultipleSeqAlignment object for storage in a DictList
            aln.id = aln_id
            aln.annotations['a_seq'] = self.representative_sequence.id
            aln.annotations['b_seq'] = seq.id

            if parse:
                aln_df = ssbio.protein.sequence.utils.alignment.get_alignment_df(a_aln_seq=str(list(aln)[0].seq),
                                                                                 b_aln_seq=str(list(aln)[1].seq))
                aln.annotations['ssbio_type'] = 'seqalign'
                aln.annotations['mutations'] = ssbio.protein.sequence.utils.alignment.get_mutations(aln_df)
                aln.annotations['deletions'] = ssbio.protein.sequence.utils.alignment.get_deletions(aln_df)
                aln.annotations['insertions'] = ssbio.protein.sequence.utils.alignment.get_insertions(aln_df)

            self.sequence_alignments.append(aln)

    def get_sequence_properties(self, representative_only=True):
        """Run Biopython ProteinAnalysis and EMBOSS pepstats to summarize basic statistics of the protein sequences.
        Annotations are stored in the protein's respective SeqProp objects at ``.annotations``
        
        Args:
            representative_only (bool): If analysis should only be run on the representative sequence

        """
        if representative_only:
            # Check if a representative sequence was set
            if not self.representative_sequence:
                log.warning('{}: no representative sequence set, cannot get sequence properties'.format(self.id))
                return

            # Also need to check if a sequence has been stored
            if not self.representative_sequence.seq:
                log.warning('{}: representative sequence {} set, but no sequence stored. '
                            'Cannot get sequence properties.'.format(self.id, self.representative_sequence.id))
                return

            self.representative_sequence.get_biopython_pepstats()
            self.representative_sequence.get_emboss_pepstats()

        if not representative_only:
            for s in self.sequences:
                # Need to check if a sequence has been stored
                if not s.seq:
                    log.warning('{}: no sequence stored. '
                                'Cannot get sequence properties.'.format(s.id))
                    continue

                else:
                    s.get_biopython_pepstats()
                    s.get_emboss_pepstats()

    def prep_itasser_modeling(self, itasser_installation, itlib_folder, runtype, create_in_dir=None,
                            execute_from_dir=None, print_exec=False, **kwargs):
        """Prepare to run I-TASSER homology modeling for a sequence.

        Args:
            itasser_installation (str): Path to I-TASSER folder, i.e. ``~/software/I-TASSER4.4``
            itlib_folder (str): Path to ITLIB folder, i.e. ``~/software/ITLIB``
            runtype: How you will be running I-TASSER - local, slurm, or torque
            create_in_dir (str): Local directory where folders will be created
            execute_from_dir (str): Optional path to execution directory - use this if you are copying the homology
                models to another location such as a supercomputer for running
            all_genes (bool): If all genes should be prepped, or only those without any mapped structures
            print_exec (bool): If the execution statement should be printed to run modelling
            **kwargs:

        Todo:
            * kwargs - extra options for SLURM or Torque execution

        """
        # TODO: kwargs for slurm/torque options

        if not create_in_dir:
            if not self.structure_dir:
                raise ValueError('Output directory must be specified')
            self.homology_models_dir = op.join(self.structure_dir, 'homology_models')
        else:
            self.homology_models_dir = create_in_dir

        ssbio.utils.make_dir(self.homology_models_dir)

        if not execute_from_dir:
            execute_from_dir = self.homology_models_dir

        repseq = self.representative_sequence

        itasser_kwargs = {'light': True,
                          'java_home': None,
                          'binding_site_pred': False,
                          'ec_pred': False,
                          'go_pred': False,
                          'job_scheduler_header': None,
                          'additional_options': None}

        if kwargs:
            itasser_kwargs.update(kwargs)

        ITASSERPrep(ident=self.id, seq_str=repseq.seq_str, root_dir=self.homology_models_dir,
                    itasser_path=itasser_installation, itlib_path=itlib_folder,
                    runtype=runtype, print_exec=print_exec, execute_dir=execute_from_dir,
                    java_home=itasser_kwargs['java_home'],
                    light=itasser_kwargs['light'],
                    binding_site_pred=itasser_kwargs['binding_site_pred'],
                    ec_pred=itasser_kwargs['ec_pred'],
                    go_pred=itasser_kwargs['go_pred'],
                    job_scheduler_header=itasser_kwargs['job_scheduler_header'],
                    additional_options=itasser_kwargs['additional_options'])

        log.debug('Prepared I-TASSER modeling folder {}'.format(self.homology_models_dir))

    def blast_representative_sequence_to_pdb(self, seq_ident_cutoff=0, evalue=0.0001, display_link=False,
                                             outdir=None, force_rerun=False):
        """BLAST the representative protein sequence to the PDB. Saves a raw BLAST result file (XML file).

        Args:
            seq_ident_cutoff (float, optional): Cutoff results based on percent coverage (in decimal form)
            evalue (float, optional): Cutoff for the E-value - filters for significant hits. 0.001 is liberal,
                0.0001 is stringent (default).
            display_link (bool, optional): Set to True if links to the HTML results should be displayed
            outdir (str): Path to output directory of downloaded XML files, must be set if protein directory
                was not initialized
            force_rerun (bool, optional): If existing BLAST results should not be used, set to True. Default is False

        Returns:
            list: List of new ``PDBProp``s added to the ``structures`` attribute

        """
        # Check if a representative sequence was set
        if not self.representative_sequence:
            log.warning('{}: no representative sequence set, cannot BLAST'.format(self.id))
            return None

        # Also need to check if a sequence has been stored
        if not self.representative_sequence.seq:
            log.warning('{}: no representative sequence loaded, cannot BLAST'.format(self.id))
            return None

        # BLAST the sequence to the PDB
        blast_results = self.representative_sequence.blast_pdb(seq_ident_cutoff=seq_ident_cutoff,
                                                               evalue=evalue,
                                                               display_link=display_link,
                                                               outdir=outdir,
                                                               force_rerun=force_rerun)

        new_pdbs = []

        # Add BLAST results to the list of structures
        if blast_results:

            # Filter for new BLASTed PDBs
            pdbs = [x['hit_pdb'].lower() for x in blast_results]
            new_pdbs = [y for y in pdbs if not self.structures.has_id(y)]
            if new_pdbs:
                log.debug('{}: adding {} PDBs from BLAST results'.format(self.id, len(new_pdbs)))
            blast_results = [z for z in blast_results if z['hit_pdb'].lower() in new_pdbs]

            for blast_result in blast_results:
                pdb = blast_result['hit_pdb'].lower()
                chains = blast_result['hit_pdb_chains']

                for chain in chains:
                    # load_pdb will append this protein to the list of structures
                    new_pdb = self.load_pdb(pdb_id=pdb, mapped_chains=chain)
                    new_pdb.add_chain_ids(chain)
                    new_chain = new_pdb.chains.get_by_id(chain)

                    # Store BLAST results within the chain
                    new_chain.blast_results = blast_result

        return new_pdbs

    @property
    def df_pdb_blast(self):
        """Get a dataframe of PDB BLAST results"""

        blast_results_pre_df = []

        for p in self.get_experimental_structures():
            for c in p.chains:
                if hasattr(c, 'blast_results'):
                    # Summary dataframe
                    infodict = p.get_dict_with_chain(chain=c.id)['blast_results']
                    infodict['pdb_id'] = p.id
                    infodict['pdb_chain_id'] = c.id
                    blast_results_pre_df.append(infodict)

        cols = ['pdb_id', 'pdb_chain_id', 'hit_score', 'hit_evalue', 'hit_percent_similar',
                'hit_percent_ident', 'hit_percent_gaps', 'hit_num_ident', 'hit_num_similar', 'hit_num_gaps']
        df = pd.DataFrame.from_records(blast_results_pre_df, columns=cols).set_index('pdb_id')
        return ssbio.utils.clean_df(df)

    def map_uniprot_to_pdb(self, seq_ident_cutoff=0.0, outdir=None, force_rerun=False):
        """Map the representative sequence's UniProt ID to PDB IDs using the PDBe "Best Structures" API.

        Args:
            seq_ident_cutoff (float): Sequence identity cutoff in decimal form
            outdir (str): Output directory to cache JSON results of search
            force_rerun (bool): Force re-downloading of JSON results if they already exist

        Returns:
            list: A rank-ordered list of PDB IDs that map to the UniProt ID

        """
        if not self.representative_sequence:
            log.error('{}: no representative sequence set, cannot use best structures API'.format(self.id))
            return None

        # Check if a UniProt ID is attached to the representative sequence
        uniprot_id = self.representative_sequence.uniprot
        if not uniprot_id:
            log.error('{}: no representative UniProt ID set, cannot use best structures API'.format(self.id))
            return None

        if '-' in uniprot_id:
            log.debug('{}: "-" detected in UniProt ID, isoform specific sequences are ignored with best structures API'.format(self.id))
            uniprot_id = uniprot_id.split('-')[0]

        if not outdir:
            outdir = self.sequence_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        best_structures = ssbio.databases.pdb.best_structures(uniprot_id,
                                                              outname='{}_best_structures'.format(custom_slugify(uniprot_id)),
                                                              outdir=outdir,
                                                              seq_ident_cutoff=seq_ident_cutoff,
                                                              force_rerun=force_rerun)

        new_pdbs = []
        if best_structures:
            rank = 1
            for best_structure in best_structures:
                currpdb = str(best_structure['pdb_id'].lower())
                new_pdbs.append(currpdb)
                currchain = str(best_structure['chain_id'])

                # load_pdb will append this protein to the list
                new_pdb = self.load_pdb(pdb_id=currpdb, mapped_chains=currchain)

                # Also add this chain to the chains attribute so we can save the
                # info we get from best_structures
                new_pdb.add_chain_ids(currchain)

                pdb_specific_keys = ['experimental_method', 'resolution']
                chain_specific_keys = ['coverage', 'start', 'end', 'unp_start', 'unp_end']

                new_pdb.update(best_structure, only_keys=pdb_specific_keys)
                new_chain = new_pdb.chains.get_by_id(currchain)
                new_chain.update(best_structure, only_keys=chain_specific_keys)
                new_chain.update({'rank': rank})

                rank += 1

            log.debug('{}, {}: {} PDB/chain pairs mapped'.format(self.id, uniprot_id, len(best_structures)))
        else:
            log.debug('{}, {}: no PDB/chain pairs mapped'.format(self.id, uniprot_id))

        return new_pdbs

    @property
    def df_pdb_ranking(self):
        """Get a dataframe of UniProt -> best structure in PDB results"""

        best_structures_pre_df = []

        chain_specific_keys = ['coverage', 'start', 'end', 'unp_start', 'unp_end', 'rank']

        for p in self.get_experimental_structures():
            for c in p.chains:
                if hasattr(c, 'rank'):
                    # Summary dataframe
                    infodict = p.get_dict_with_chain(chain=c.id, df_format=True, chain_keys=chain_specific_keys)
                    infodict['pdb_id'] = p.id
                    infodict['pdb_chain_id'] = c.id
                    infodict['uniprot'] = self.representative_sequence.uniprot
                    best_structures_pre_df.append(infodict)

        cols = ['uniprot', 'pdb_id', 'pdb_chain_id', 'experimental_method', 'resolution', 'coverage',
                'taxonomy_name', 'start', 'end', 'unp_start', 'unp_end', 'rank']
        df = pd.DataFrame.from_records(best_structures_pre_df, columns=cols).set_index(['pdb_id', 'pdb_chain_id'])
        return ssbio.utils.clean_df(df)

    def load_pdb(self, pdb_id, mapped_chains=None, pdb_file=None, file_type=None, is_experimental=True,
                 set_as_representative=False, representative_chain=None, force_rerun=False):
        """Load a structure ID and optional structure file into the structures attribute.

        Args:
            pdb_id (str): PDB ID
            mapped_chains (str, list): Chain ID or list of IDs which you are interested in
            pdb_file (str): Path to PDB file
            file_type (str): Type of PDB file
            is_experimental (bool): If this structure file is experimental
            set_as_representative (bool): If this structure should be set as the representative structure
            representative_chain (str): If ``set_as_representative`` is ``True``, provide the representative chain ID
            force_rerun (bool): If the PDB should be reloaded if it is already in the list of structures

        Returns:
            PDBProp: The object that is now contained in the structures attribute

        """

        if self.structures.has_id(pdb_id):
            # Remove the structure if set to force rerun
            if force_rerun:
                existing = self.structures.get_by_id(pdb_id)
                self.structures.remove(existing)
            # Otherwise just retrieve it
            else:
                log.debug('{}: PDB ID already present in list of structures'.format(pdb_id))
                pdb = self.structures.get_by_id(pdb_id)
                if pdb_file:
                    pdb.load_structure_path(pdb_file, file_type)
                if mapped_chains:
                    pdb.add_mapped_chain_ids(mapped_chains)

        # Create a new StructProp entry
        if not self.structures.has_id(pdb_id):
            if is_experimental:
                pdb = PDBProp(ident=pdb_id, mapped_chains=mapped_chains, structure_path=pdb_file, file_type=file_type)
            else:
                pdb = StructProp(ident=pdb_id, mapped_chains=mapped_chains, structure_path=pdb_file, file_type=file_type)
            self.structures.append(pdb)

        if set_as_representative:
            # Parse structure so chains are stored before setting representative
            pdb.parse_structure()
            self._representative_structure_setter(structprop=pdb, keep_chain=representative_chain, force_rerun=force_rerun)

        return self.structures.get_by_id(pdb_id)

    def load_itasser_folder(self, ident, itasser_folder, organize=False, outdir=None, organize_name=None,
                            set_as_representative=False, representative_chain='X', force_rerun=False):
        """Load the results folder from an I-TASSER run (local, not from the server), copy structure files over, and create summary dataframes.

        Args:
            ident: I-TASSER ID
            itasser_folder: Path to results folder
            organize (bool): If select files from modeling should be copied to the Protein directory
            outdir (str): Path to directory where files will be copied and organized to
            organize_name (str): Basename of files to rename results to. If not provided, will use id attribute.
            set_as_representative: If this structure should be set as the representative structure
            representative_chain (str): If ``set_as_representative`` is ``True``, provide the representative chain ID
            force_rerun (bool): If the PDB should be reloaded if it is already in the list of structures

        Returns:
            ITASSERProp: The object that is now contained in the structures attribute

        """
        if organize:
            if not outdir:
                outdir = self.structure_dir
                if not outdir:
                    raise ValueError('Directory to copy results to must be specified')

        if self.structures.has_id(ident):
            if force_rerun:
                existing = self.structures.get_by_id(ident)
                self.structures.remove(existing)
            else:
                log.debug('{}: already present in list of structures'.format(ident))
                itasser = self.structures.get_by_id(ident)

        if not self.structures.has_id(ident):
            itasser = ITASSERProp(ident, itasser_folder)
            self.structures.append(itasser)

        if set_as_representative:
            self._representative_structure_setter(structprop=itasser, keep_chain=representative_chain, force_rerun=force_rerun)

        if organize:
            if itasser.structure_file:
                # The name of the actual pdb file will be $GENEID_model1.pdb
                if not organize_name:
                    new_itasser_name = self.id + '_model1'
                else:
                    new_itasser_name = organize_name

                # Additional results will be stored in a subdirectory
                dest_itasser_extra_dir = op.join(outdir, '{}_itasser'.format(new_itasser_name))
                ssbio.utils.make_dir(dest_itasser_extra_dir)

                # Copy the model1.pdb and also create summary dataframes
                itasser.copy_results(copy_to_dir=outdir, rename_model_to=new_itasser_name, force_rerun=force_rerun)

        return self.structures.get_by_id(ident)

    @property
    def df_homology_models(self):
        """Get a dataframe of homology models"""

        itasser_pre_df = []

        df_cols = ['id', 'structure_file', 'model_date', 'difficulty',
                   'top_template_pdb', 'top_template_chain', 'c_score',
                   'tm_score', 'tm_score_err', 'rmsd', 'rmsd_err',
                   'top_bsite_site_num', 'top_bsite_c_score', 'top_bsite_cluster_size', 'top_bsite_algorithm',
                   'top_bsite_pdb_template_id', 'top_bsite_pdb_template_chain', 'top_bsite_pdb_ligand',
                   'top_bsite_binding_location_coords', 'top_bsite_c_score_method', 'top_bsite_binding_residues',
                   'top_bsite_ligand_cluster_counts',
                   'top_ec_pdb_template_id', 'top_ec_pdb_template_chain', 'top_ec_tm_score', 'top_ec_rmsd',
                   'top_ec_seq_ident', 'top_ec_seq_coverage', 'top_ec_c_score', 'top_ec_ec_number',
                   'top_ec_binding_residues',
                   'top_go_mf_go_id', 'top_go_mf_go_term', 'top_go_mf_c_score', 'top_go_bp_go_id', 'top_go_bp_go_term',
                   'top_go_bp_c_score', 'top_go_cc_go_id', 'top_go_cc_go_term', 'top_go_cc_c_score']

        for p in self.get_homology_models():
            # Summary dataframe
            new_homology_dict = p.get_dict(df_format=True, only_attributes=df_cols)
            itasser_pre_df.append(new_homology_dict)

        df = pd.DataFrame.from_records(itasser_pre_df, columns=df_cols).set_index('id')
        return ssbio.utils.clean_df(df)

    def pdb_downloader_and_metadata(self, outdir=None, pdb_file_type=None, force_rerun=False):
        """Download experimental structures"""
        if not outdir:
            outdir = self.structure_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        if not pdb_file_type:
            pdb_file_type = self.pdb_file_type

        # Check if we have any PDBs
        if self.num_structures_experimental == 0:
            log.debug('{}: no structures available - nothing will be downloaded'.format(self.id))
            return

        downloaded_pdb_ids = []
        # Download the PDBs
        for s in self.get_experimental_structures():
            log.debug('{}: downloading structure file from the PDB...'.format(s.id))
            s.download_structure_file(outdir=outdir, file_type=pdb_file_type, force_rerun=force_rerun)
            # Download the mmCIF header file to get additional information
            if 'cif' not in pdb_file_type:
                s.download_cif_header_file(outdir=outdir, force_rerun=force_rerun)
            downloaded_pdb_ids.append(s.id)

        return downloaded_pdb_ids

    @property
    def df_pdb_metadata(self):
        """Get a dataframe of PDB metadata (PDBs have to be downloaded first)"""
        if self.num_structures == 0:
            log.error('No experimental PDB structures have been mapped to protein')
            return pd.DataFrame()

        pdb_pre_df = []

        for p in self.get_experimental_structures():
            if not p.structure_file:
                log.error('{}: PDB file has not been downloaded, run "pdb_downloader_and_metadata" \
                    to download all structures'.format(p.id))
                continue

            infodict = p.get_dict(df_format=True)
            infodict['pdb_id'] = p.id
            pdb_pre_df.append(infodict)

        cols = ['pdb_id', 'pdb_title', 'description', 'experimental_method', 'mapped_chains',
                'resolution', 'chemicals', 'date', 'taxonomy_name',
                'structure_file']
        df = pd.DataFrame.from_records(pdb_pre_df, columns=cols).set_index('pdb_id')
        return ssbio.utils.clean_df(df)

    def align_seqprop_to_structprop(self, seqprop, structprop, chains=None, outdir=None,
                                    engine='needle', parse=True, force_rerun=False,
                                    **kwargs):
        """Run and store alignments of a SeqProp to chains in the ``mapped_chains`` attribute of a StructProp.

        Alignments are stored in the sequence_alignments attribute, with the IDs formatted
            as <SeqProp_ID>_<StructProp_ID>-<Chain_ID>. Although it is more intuitive to align to individual 
            ChainProps, StructProps should be loaded as little as possible to reduce run times.

        Args:
            seqprop (SeqProp): SeqProp object with a loaded sequence
            structprop (StructProp): StructProp object with a loaded structure
            chains (str, list): Chain ID or IDs to map to. If not specified, ``mapped_chains`` attribute is inspected
                for chains. If no chains there, all chains will be aligned to.
            outdir (str): Directory to output sequence alignment files (only if running with needle)
            engine (str): Which pairwise alignment tool to use ("needle" or "biopython")
            parse (bool): Store locations of mutations, insertions, and deletions in the alignment object (as an annotation)
            force_rerun (bool): If alignments should be rerun
            **kwargs: Other alignment options
            
        """
        # TODO: **kwargs for alignment options

        if not outdir:
            outdir = self.sequence_dir

        # Parse the structure so chain sequences are stored
        structprop.parse_structure()

        if chains:
            chains_to_align_to = ssbio.utils.force_list(chains)
        elif structprop.mapped_chains:
            chains_to_align_to = structprop.mapped_chains
        else:
            log.warning('{}-{}: no chains specified in structure to align to, all chains will be aligned to'.format(seqprop.id,
                                                                                                                    structprop.id))
            chains_to_align_to = structprop.chains.list_attr('id')

        for chain_id in chains_to_align_to:
            full_structure_id = '{}-{}'.format(structprop.id, chain_id)
            aln_id = '{}_{}'.format(seqprop.id, full_structure_id)
            outfile = '{}.needle'.format(aln_id)

            if self.sequence_alignments.has_id(aln_id) and not force_rerun:
                log.debug('{}: alignment already completed, skipping'.format(aln_id))
                continue

            log.debug('{}: running pairwise alignment to structure {}, chain {}'.format(seqprop.id,
                                                                                        structprop.id,
                                                                                        chain_id))

            # Check if the chain ID actually exists
            if structprop.chains.has_id(chain_id):
                chain_prop = structprop.chains.get_by_id(chain_id)
                chain_seq_record = chain_prop.seq_record
            else:
                log.error('{}: chain not present in structure file!'.format(full_structure_id))
                continue

            # Check if the chain sequence was parsed
            if not chain_seq_record:
                log.error('{}: chain sequence not available, was structure parsed?'.format(full_structure_id))
                continue

            # Run the pairwise alignment
            try:
                aln = ssbio.protein.sequence.utils.alignment.pairwise_sequence_alignment(a_seq=seqprop,
                                                                                         a_seq_id=seqprop.id,
                                                                                         b_seq=chain_seq_record,
                                                                                         b_seq_id=full_structure_id,
                                                                                         engine=engine,
                                                                                         outdir=outdir,
                                                                                         outfile=outfile,
                                                                                         force_rerun=force_rerun)
            except ValueError:
                log.error('{}: alignment failed to run, unable to check structure'.format(full_structure_id))
                continue

            # Add an identifier to the MultipleSeqAlignment object for storage in a DictList
            aln.id = aln_id
            # Add annotations to keep track of what was aligned
            aln.annotations['a_seq'] = seqprop.id
            aln.annotations['b_seq'] = full_structure_id
            aln.annotations['structure_id'] = structprop.id
            aln.annotations['chain_id'] = chain_id
            aln.annotations['ssbio_type'] = 'structalign'

            # Optionally parse for locations of mutations, deletions, and insertions
            # Store mapping to chain index as letter annotations in the sequence
            # Store locations in the alignment's annotations
            if parse:
                aln_df = ssbio.protein.sequence.utils.alignment.get_alignment_df(a_aln_seq=str(list(aln)[0].seq),
                                                                                 b_aln_seq=str(list(aln)[1].seq))

                chain_indices = aln_df[pd.notnull(aln_df.id_a_pos)].id_b_pos.tolist()
                seqprop.letter_annotations['{}_chain_index'.format(aln_id)] = chain_indices

                aln.annotations['mutations'] = ssbio.protein.sequence.utils.alignment.get_mutations(aln_df)
                aln.annotations['deletions'] = ssbio.protein.sequence.utils.alignment.get_deletions(aln_df)
                aln.annotations['insertions'] = ssbio.protein.sequence.utils.alignment.get_insertions(aln_df)

            if force_rerun and self.sequence_alignments.has_id(aln.id):
                self.sequence_alignments.remove(aln.id)
            self.sequence_alignments.append(aln)

    def find_representative_chain(self, seqprop, structprop, chains_to_check=None,
                                  seq_ident_cutoff=0.5, allow_missing_on_termini=0.2,
                                  allow_mutants=True, allow_deletions=False,
                                  allow_insertions=False, allow_unresolved=True):
        """Set and return the representative chain based on sequence quality checks to a reference sequence.

        Args:
            seqprop (SeqProp): SeqProp object to compare to chain sequences
            structprop (StructProp): StructProp object with chains to compare to in the ``mapped_chains`` attribute. If
                there are none present, ``chains_to_check`` can be specified, otherwise all chains are checked.
            chains_to_check (str, list): Chain ID or IDs to check for sequence coverage quality
            seq_ident_cutoff (float): Percent sequence identity cutoff, in decimal form
            allow_missing_on_termini (float): Percentage of the total length of the reference sequence which will be ignored
                when checking for modifications. Example: if 0.1, and reference sequence is 100 AA, then only residues
                5 to 95 will be checked for modifications.
            allow_mutants (bool): If mutations should be allowed or checked for
            allow_deletions (bool): If deletions should be allowed or checked for
            allow_insertions (bool): If insertions should be allowed or checked for
            allow_unresolved (bool): If unresolved residues should be allowed or checked for

        Returns:
            str: the best chain ID, if any

        """
        if chains_to_check:
            chains_to_check = ssbio.utils.force_list(chains_to_check)
        elif structprop.mapped_chains:
            chains_to_check = structprop.mapped_chains
        else:
            log.warning('{}-{}: no chains specified in structure to align to, all chains will be checked'.format(seqprop.id,
                                                                                                                 structprop.id))
            chains_to_check = structprop.chains.list_attr('id')

        for chain_id in chains_to_check:
            full_structure_id = '{}-{}'.format(structprop.id, chain_id)
            aln_id = '{}_{}'.format(seqprop.id, full_structure_id)

            if self.sequence_alignments.has_id(aln_id):
                alignment = self.sequence_alignments.get_by_id(aln_id)
            else:
                log.error('{}: structure alignment not found, please run the alignment first'.format(aln_id))
                continue

            # Compare sequence to structure's sequence using the alignment
            found_good_chain = ssbio.protein.structure.properties.quality.sequence_checker(reference_seq_aln=alignment[0],
                                                                                           structure_seq_aln=alignment[1],
                                                                                           seq_ident_cutoff=seq_ident_cutoff,
                                                                                           allow_missing_on_termini=allow_missing_on_termini,
                                                                                           allow_mutants=allow_mutants,
                                                                                           allow_deletions=allow_deletions,
                                                                                           allow_insertions=allow_insertions,
                                                                                           allow_unresolved=allow_unresolved)

            # If found_good_chain = True, return chain ID
            # If not, move on to the next potential chain
            if found_good_chain:
                self.representative_chain = chain_id
                self.representative_chain_seq_coverage = alignment.annotations['percent_identity']
                return chain_id
        else:
            log.debug('{}: no chains meet quality checks'.format(structprop.id))
            return None

    def map_seqprop_resnums_to_structprop_chain_index(self, resnums, seqprop=None, structprop=None, chain_id=None,
                                                      use_representatives=False):
        """Map a residue number in the seqprop to the mapping index in the structprop/chain_id

        Use this to get the indices of the chain to then get structure residue number.

        Args:
            resnums (int, list): Residue numbers in the sequence
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): Chain ID to map to index
            use_representatives (bool): If representative sequence/structure/chain should be used in mapping

        Returns:
            dict: Mapping of resnums to indices

        """
        resnums = ssbio.utils.force_list(resnums)

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain
            full_structure_id = '{}-{}'.format(structprop.id, chain_id).replace('REP-', '')
            aln_id = '{}_{}'.format(seqprop.id, full_structure_id)
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Please specify sequence, structure, and chain ID')
            full_structure_id = '{}-{}'.format(structprop.id, chain_id)
            aln_id = '{}_{}'.format(seqprop.id, full_structure_id)

        access_key = '{}_chain_index'.format(aln_id)
        if access_key not in seqprop.letter_annotations:
            raise KeyError('{}: structure mapping {} not available in sequence letter annotations. Was alignment parsed? '
                           'Run ``align_seqprop_to_structprop`` with ``parse=True``.'.format(access_key, aln_id))
        chain_index_mapping = seqprop.letter_annotations[access_key]

        resnum_to_chain_index = {}
        for x in resnums:
            ix = chain_index_mapping[x - 1] - 1

            if np.isnan(ix):
                log.warning('{}-{}, {}: no equivalent residue found in structure sequence'.format(structprop.id,
                                                                                                  chain_id,
                                                                                                  x))
            else:
                resnum_to_chain_index[int(x)] = int(ix)

        return resnum_to_chain_index

    def map_seqprop_resnums_to_structprop_resnums(self, resnums, seqprop=None, structprop=None, chain_id=None,
                                                  use_representatives=False):
        """Map a residue number in the seqprop to the structure's residue number for a specified chain.

        Args:
            resnums (int, list): Residue numbers in the sequence
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): Chain ID to map to
            use_representatives (bool): If the representative sequence and structure should be used. If True, seqprop,
                structprop, and chain_id do not need to be defined.

        Returns:
            dict: Mapping of sequence residue numbers to structure residue numbers

        """
        resnums = ssbio.utils.force_list(resnums)

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Please specify sequence, structure, and chain ID')

        mapping_to_repchain_index = self.map_seqprop_resnums_to_structprop_chain_index(resnums=resnums,
                                                                                       seqprop=seqprop,
                                                                                       structprop=structprop,
                                                                                       chain_id=chain_id,
                                                                                       use_representatives=use_representatives)

        chain = structprop.chains.get_by_id(chain_id)
        chain_structure_resnum_mapping = chain.seq_record.letter_annotations['structure_resnums']

        final_mapping = {}
        for k, v in mapping_to_repchain_index.items():
            k = int(k)
            rn = chain_structure_resnum_mapping[v]

            if rn == float('Inf'):
                log.warning('{}-{}, {}: structure file does not contain coordinates for this residue'.format(structprop.id,
                                                                                                             chain_id,
                                                                                                             k))
            else:
                rn = int(rn)
                final_mapping[k] = rn

                # Additionally report if residues are the same - they could be different in the structure though
                format_data = {'seqprop_id'       : seqprop.id,
                               'seqprop_resid'    : seqprop[k - 1],
                               'seqprop_resnum'   : k,
                               'structprop_id'    : structprop.id,
                               'structprop_chid'  : chain_id,
                               'structprop_resid' : chain.seq_record[rn - 1],
                               'structprop_resnum': rn}

                if seqprop[k-1] != chain.seq_record[rn-1]:
                    log.warning('Sequence {seqprop_id} residue {seqprop_resid}{seqprop_resnum} does not match to '
                                'structure {structprop_id}-{structprop_chid} residue '
                                '{structprop_resid}{structprop_resnum}. NOTE: this may be due to '
                                'structural differences'.format(**format_data))
                else:
                    log.debug('Sequence {seqprop_id} residue {seqprop_resid}{seqprop_resnum} is mapped to '
                              'structure {structprop_id}-{structprop_chid} residue '
                              '{structprop_resid}{structprop_resnum}'.format(**format_data))

        return final_mapping

    def map_structprop_resnums_to_seqprop_resnums(self, resnums, structprop=None, chain_id=None, seqprop=None,
                                                  use_representatives=False):
        """Map a residue number in the structprop+chain to the seqprop's residue number.

        Args:
            resnums (int, list): Residue numbers in the structure
            structprop (StructProp): StructProp object
            chain_id (str): Chain ID to map from
            seqprop (SeqProp): SeqProp object
            use_representatives (bool): If the representative sequence and structure should be used. If True, seqprop,
                structprop, and chain_id do not need to be defined.

        Returns:
            dict: Mapping of structure residue numbers to sequence residue numbers

        """
        pass
        # resnums = ssbio.utils.force_list(resnums)
        #
        # if use_representatives:
        #     seqprop = self.representative_sequence
        #     structprop = self.representative_structure
        #     chain_id = self.representative_chain
        # else:
        #     if not seqprop or not structprop or not chain_id:
        #         raise ValueError('Please specify sequence, structure, and chain ID')
        #
        # mapping_to_repchain_index = self.map_seqprop_resnums_to_structprop_chain_index(resnums=resnums,
        #                                                                                seqprop=seqprop,
        #                                                                                structprop=structprop,
        #                                                                                chain_id=chain_id,
        #                                                                                use_representatives=use_representatives)
        #
        # chain = structprop.chains.get_by_id(chain_id)
        # chain_structure_resnum_mapping = chain.seq_record.letter_annotations['structure_resnums']
        #
        # final_mapping = {}
        # for k, v in mapping_to_repchain_index.items():
        #     k = int(k)
        #     rn = chain_structure_resnum_mapping[v]
        #
        #     if rn == float('Inf'):
        #         log.warning('{}-{}, {}: structure file does not contain coordinates for this residue'.format(structprop.id,
        #                                                                                                      chain_id,
        #                                                                                                      k))
        #     else:
        #         rn = int(rn)
        #         final_mapping[k] = rn
        #
        #         # Additionally report if residues are the same - they could be different in the structure though
        #         format_data = {'seqprop_id'       : seqprop.id,
        #                        'seqprop_resid'    : seqprop[k - 1],
        #                        'seqprop_resnum'   : k,
        #                        'structprop_id'    : structprop.id,
        #                        'structprop_chid'  : chain_id,
        #                        'structprop_resid' : chain.seq_record[rn - 1],
        #                        'structprop_resnum': rn}
        #
        #         if seqprop[k-1] != chain.seq_record[rn-1]:
        #             log.warning('Sequence {seqprop_id} residue {seqprop_resid}{seqprop_resnum} does not match to '
        #                         'structure {structprop_id}-{structprop_chid} residue '
        #                         '{structprop_resid}{structprop_resnum}. NOTE: this may be due to '
        #                         'structural differences'.format(**format_data))
        #         else:
        #             log.debug('Sequence {seqprop_id} residue {seqprop_resid}{seqprop_resnum} is mapped to '
        #                       'structure {structprop_id}-{structprop_chid} residue '
        #                       '{structprop_resid}{structprop_resnum}'.format(**format_data))
        #
        # return final_mapping


    def _representative_structure_setter(self, structprop, keep_chain, clean=True, keep_chemicals=None,
                                         out_suffix='_clean', outdir=None, force_rerun=False):
        """Set the representative structure by 1) cleaning it and 2) copying over attributes of the original structure.

        The structure is copied because the chains stored may change, and cleaning it makes a new PDB file.

        Args:
            structprop (StructProp): StructProp object to set as representative
            keep_chain (str): Chain ID to keep
            clean (bool): If the PDB file should be cleaned (see ssbio.structure.utils.cleanpdb)
            keep_chemicals (str, list): Keep specified chemical names
            out_suffix (str): Suffix to append to clean PDB file
            outdir (str): Path to output directory

        Returns:
            StructProp: representative structure

        """

        # Set output directory for cleaned PDB file
        if not outdir:
            outdir = self.structure_dir
            if not outdir:
                raise ValueError('Output directory must be specified')

        # Create new ID for this representative structure, it cannot be the same as the original one
        new_id = 'REP-{}'.format(structprop.id)

        # Remove the previously set representative structure if set to force rerun
        if self.structures.has_id(new_id):
            if force_rerun:
                existing = self.structures.get_by_id(new_id)
                self.structures.remove(existing)

        # If the structure is to be cleaned, and which chain to keep
        if clean:
            final_pdb = structprop.clean_structure(outdir=outdir, out_suffix=out_suffix,
                                                   keep_chemicals=keep_chemicals, keep_chains=keep_chain,
                                                   force_rerun=force_rerun)
            log.debug('{}: cleaned structure and saved new file at {}'.format(structprop.id, final_pdb))
        else:
            final_pdb = structprop.structure_path

        self.representative_structure = StructProp(ident=new_id, chains=keep_chain, mapped_chains=keep_chain,
                                                   structure_path=final_pdb, file_type='pdb')
        self.representative_chain = keep_chain

        self.representative_structure.update(structprop.get_dict_with_chain(chain=keep_chain),
                                             only_keys=self.__representative_structure_attributes,
                                             overwrite=True)

        # Save the original structure ID as an extra attribute
        self.representative_structure.original_structure_id = structprop.id

        # Also need to parse the clean structure and save its sequence..
        self.representative_structure.parse_structure()

        # And finally add it to the list of structures
        self.structures.append(self.representative_structure)

    def set_representative_structure(self, seq_outdir=None, struct_outdir=None, pdb_file_type=None,
                                     engine='needle', always_use_homology=False, rez_cutoff=0.0,
                                     seq_ident_cutoff=0.5, allow_missing_on_termini=0.2,
                                     allow_mutants=True, allow_deletions=False,
                                     allow_insertions=False, allow_unresolved=True, clean=True,
                                     keep_chemicals=None,
                                     force_rerun=False):
        """Set a representative structure from a structure in self.structures

        Args:
            seq_outdir (str): Path to output directory of sequence alignment files
            struct_outdir (str): Path to output directory of structure files
            pdb_file_type (str): pdb, pdb.gz, mmcif, cif, cif.gz, xml.gz, mmtf, mmtf.gz - PDB structure file type that
                should be downloaded
            engine (str): "needle" or "biopython" - which pairwise sequence alignment engine should be used
                needle is the standard EMBOSS tool to run pairwise alignments
                biopython is Biopython's implementation of needle. Results can differ!
            always_use_homology (bool): If homology models should always be set as the representative structure
            rez_cutoff (float): Resolution cutoff, in Angstroms (only if experimental structure)
            seq_ident_cutoff (float): Percent sequence identity cutoff, in decimal form
            allow_missing_on_termini (float): Percentage of the total length of the reference sequence which will be ignored
                when checking for modifications. Example: if 0.1, and reference sequence is 100 AA, then only residues
                5 to 95 will be checked for modifications.
            allow_mutants (bool): If mutations should be allowed or checked for
            allow_deletions (bool): If deletions should be allowed or checked for
            allow_insertions (bool): If insertions should be allowed or checked for
            allow_unresolved (bool): If unresolved residues should be allowed or checked for
            clean (bool): If structure should be cleaned
            keep_chemicals (str, list): Keep specified chemical names if structure is to be cleaned
            force_rerun (bool): If sequence to structure alignment should be rerun

        Returns:
            StructProp: Representative structure from the list of structures.

            This is a not a map to the original structure, it is copied from its reference.

        """
        log.debug('{}: setting representative structure'.format(self.id))

        if len(self.structures) == 0:
            log.debug('{}: no structures available'.format(self.id))
            return None

        if not self.representative_sequence:
            log.error('{}: no representative sequence to compare structures to'.format(self.id))
            return None

        if self.representative_structure and not force_rerun:
            log.debug('{}: representative structure already set'.format(self.id))
            return self.representative_structure

        if not pdb_file_type:
            pdb_file_type = self.pdb_file_type

        if not seq_outdir:
            seq_outdir = self.sequence_dir
            if not seq_outdir:
                raise ValueError('Sequence output directory must be specified')

        if not struct_outdir:
            struct_outdir = self.structure_dir
            if not struct_outdir:
                raise ValueError('Structure output directory must be specified')

        has_homology = False
        has_pdb = False
        use_homology = False
        use_pdb = False

        if self.num_structures_homology > 0:
            has_homology = True
        if self.num_structures_experimental > 0:
            has_pdb = True

        # If we mark to always use homology, use it if it exists
        if always_use_homology:
            if has_homology:
                use_homology = True
            elif has_pdb:
                use_pdb = True
        # If we don't always want to use homology, use PDB if it exists
        else:
            if has_homology and has_pdb:
                use_pdb = True
                use_homology = True
            elif has_homology and not has_pdb:
                use_homology = True
            elif has_pdb and not has_homology:
                use_pdb = True

        if use_pdb:
            # Put PDBs through QC/QA
            all_pdbs = self.get_experimental_structures()
            log.debug('{}: checking quality of {} experimental structures'.format(self.id, len(all_pdbs)))

            for pdb in all_pdbs:
                # Download the structure and parse it
                # This will add all chains to the mapped_chains attribute if there are none
                try:
                    pdb.download_structure_file(outdir=struct_outdir, file_type=pdb_file_type, force_rerun=force_rerun)
                    # Download the mmCIF header file to get additional information
                    if 'cif' not in pdb_file_type:
                        pdb.download_cif_header_file(outdir=struct_outdir, force_rerun=force_rerun)
                except (requests.exceptions.HTTPError, URLError):
                    log.error('{}: structure file could not be downloaded in {} format'.format(pdb, pdb_file_type))
                    continue
                    # TODO: add try/except to download cif file as fallback like below?

                if rez_cutoff and pdb.resolution:
                    if pdb.resolution > rez_cutoff:
                        log.debug('{}: structure does not meet experimental resolution cutoff'.format(pdb, pdb_file_type))
                        continue
                # TODO: clean up these try/except things
                try:
                    self.align_seqprop_to_structprop(seqprop=self.representative_sequence,
                                                     structprop=pdb,
                                                     outdir=seq_outdir,
                                                     engine=engine,
                                                     parse=True,
                                                     force_rerun=force_rerun)
                except (PDBConstructionException, ExtraData) as e:
                    log.error('{}: unable to parse structure file as {}. Falling back to mmCIF format.'.format(pdb, pdb_file_type))
                    print(e)
                    # Fall back to using mmCIF file if structure cannot be parsed
                    try:
                        pdb.download_structure_file(outdir=struct_outdir, file_type='cif',
                                                    force_rerun=force_rerun)
                        # Download the mmCIF header file to get additional information
                        if 'cif' not in pdb_file_type:
                            pdb.download_cif_header_file(outdir=struct_outdir, force_rerun=force_rerun)
                    except (requests.exceptions.HTTPError, URLError):
                        log.error('{}: structure file could not be downloaded'.format(pdb))
                        continue
                    self.align_seqprop_to_structprop(seqprop=self.representative_sequence,
                                                     structprop=pdb,
                                                     outdir=seq_outdir,
                                                     engine=engine,
                                                     parse=True,
                                                     force_rerun=force_rerun)

                best_chain = self.find_representative_chain(seqprop=self.representative_sequence,
                                                            structprop=pdb,
                                                            seq_ident_cutoff=seq_ident_cutoff,
                                                            allow_missing_on_termini=allow_missing_on_termini,
                                                            allow_mutants=allow_mutants, allow_deletions=allow_deletions,
                                                            allow_insertions=allow_insertions, allow_unresolved=allow_unresolved)

                if best_chain:
                    try:
                        self._representative_structure_setter(structprop=pdb,
                                                              clean=clean,
                                                              out_suffix='-{}_clean'.format(best_chain),
                                                              keep_chain=best_chain,
                                                              keep_chemicals=keep_chemicals,
                                                              outdir=struct_outdir,
                                                              force_rerun=force_rerun)
                    except:
                        # TODO: inspect causes of these errors - most common is Biopython PDBParser error
                        logging.exception("Unknown error with PDB ID {}".format(pdb.id))
                        continue
                    log.debug('{}-{}: set as representative structure'.format(pdb.id, best_chain))
                    pdb.reset_chain_seq_records()
                    return self.representative_structure
                else:
                    pdb.reset_chain_seq_records()
            else:
                log.debug('{}: no experimental structures meet cutoffs'.format(self.id))

        # If we are to use homology, save its information in the representative structure field
        if use_homology:
            log.debug('{}: checking quality of homology models'.format(self.id))
            all_models = self.get_homology_models()

            # TODO: homology models are not ordered in any other way other than how they are loaded,
            # rethink this for multiple homology models
            for homology in all_models:
                if not homology.structure_file:
                    log.debug('{}: no homology structure file'.format(self.id))
                    continue

                self.align_seqprop_to_structprop(seqprop=self.representative_sequence,
                                                 structprop=homology,
                                                 outdir=seq_outdir,
                                                 engine=engine,
                                                 parse=True,
                                                 force_rerun=force_rerun)
                best_chain = self.find_representative_chain(seqprop=self.representative_sequence,
                                                            structprop=homology,
                                                            seq_ident_cutoff=seq_ident_cutoff,
                                                            allow_missing_on_termini=allow_missing_on_termini,
                                                            allow_mutants=allow_mutants,
                                                            allow_deletions=allow_deletions,
                                                            allow_insertions=allow_insertions,
                                                            allow_unresolved=allow_unresolved)

                if best_chain:
                    # If chain ID is empty (some homology models are like that), use ID "X"
                    if not best_chain.strip():
                        best_chain = 'X'
                    try:
                        self._representative_structure_setter(structprop=homology,
                                                              # new_id='{}-{}'.format(homology.id, best_chain), # 170906 Deprecated use of new_id
                                                              clean=True,
                                                              out_suffix='-{}_clean'.format(best_chain),
                                                              keep_chain=best_chain,
                                                              outdir=struct_outdir,
                                                              force_rerun=force_rerun)
                    except:
                        # TODO: inspect causes of these errors - most common is Biopython PDBParser error
                        logging.exception("Unknown error with homology model {}".format(homology.id))
                        continue
                    log.debug('{}-{}: set as representative structure'.format(homology.id, best_chain))
                    homology.reset_chain_seq_records()
                    return self.representative_structure
                else:
                    homology.reset_chain_seq_records()

        log.warning('{}: no structures meet quality checks'.format(self.id))
        return None

    def get_dssp_annotations(self, representative_only=True, force_rerun=False):
        """Run DSSP on structures and store calculations.
        Annotations are stored in the protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-dssp']``

        Todo:
            * Some errors arise from storing annotations for nonstandard amino acids, need to run DSSP separately for those
            
        Args:
            representative_only (bool): If analysis should only be run on the representative structure
            force_rerun (bool): If calculations should be rerun even if an output file exists

        """
        if representative_only:
            if self.representative_structure:
                try:
                    self.representative_structure.get_dssp_annotations(outdir=self.structure_dir, force_rerun=force_rerun)
                except PDBException as e:
                    log.error('{}: Biopython error, issue matching sequences with {}'.format(self.id, self.representative_structure))
                    print(e)
                except TypeError as e:
                    log.error('{}: Biopython error, DSSP SeqRecord length mismatch with {}'.format(self.id, self.representative_structure))
                    print(e)
                except Exception as e:
                    log.error('{}: DSSP failed to run on {}'.format(self.id, self.representative_structure))
                    print(e)
            else:
                log.warning('{}: no representative structure set, cannot run DSSP'.format(self.id))
        else:
            for s in self.structures:
                try:
                    s.get_dssp_annotations(outdir=self.structure_dir)
                except PDBException as e:
                    log.error('{}: Biopython error, issue matching sequences with {}'.format(self.id, s.id))
                    print(e)
                except TypeError as e:
                    log.error('{}: Biopython error, DSSP SeqRecord length mismatch with {}'.format(self.id, s.id))
                    print(e)
                except Exception as e:
                    log.error('{}: DSSP failed to run on {}'.format(self.id, s.id))
                    print(e)

    def get_msms_annotations(self, representative_only=True, force_rerun=False):
        """Run MSMS on structures and store calculations.
        Annotations are stored in the protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-msms']``
        
        Args:
            representative_only (bool): If analysis should only be run on the representative structure
            force_rerun (bool): If calculations should be rerun even if an output file exists

        """
        if representative_only:
            if self.representative_structure:
                try:
                    self.representative_structure.get_residue_depths(outdir=self.structure_dir, force_rerun=force_rerun)
                except TypeError:
                    log.error('{}: MSMS SeqRecord length mismatch with {}'.format(self.id, self.representative_structure))
                except:
                    log.error('{}: unknown MSMS error with {}'.format(self.id, self.representative_structure))
            else:
                log.warning('{}: no representative structure set, cannot run MSMS'.format(self.id))
        else:
            for s in self.structures:
                try:
                    s.get_residue_depths(outdir=self.structure_dir)
                except TypeError:
                    log.error('{}: MSMS SeqRecord length mismatch with {}'.format(self.id, s.id))
                except Exception as e:
                    log.error('{}: unknown MSMS error with {}'.format(self.id, s.id))
                    print(e)

    def get_freesasa_annotations(self, include_hetatms=False, representative_only=True, force_rerun=False):
        """Run freesasa on structures and store calculations.
        Annotations are stored in the protein structure's chain sequence at:
        ``seq_record.letter_annotations['*-freesasa']``

        Args:
            include_hetatms (bool): If HETATMs should be included in calculations. Defaults to ``False``.
            representative_only (bool): If analysis should only be run on the representative structure
            force_rerun (bool): If calculations should be rerun even if an output file exists

        """
        if representative_only:
            if self.representative_structure:
                try:
                    self.representative_structure.get_freesasa_annotations(outdir=self.structure_dir,
                                                                           include_hetatms=include_hetatms,
                                                                           force_rerun=force_rerun)
                except TypeError:
                    log.error('{}: freesasa SeqRecord length mismatch with {}'.format(self.id, self.representative_structure))
                except:
                    log.error('{}: unknown freesasa error with {}'.format(self.id, self.representative_structure))
            else:
                log.warning('{}: no representative structure set, cannot run freesasa'.format(self.id))
        else:
            for s in self.structures:
                try:
                    s.get_freesasa_annotations(outdir=self.structure_dir, include_hetatms=include_hetatms)
                except TypeError:
                    log.error('{}: freesasa SeqRecord length mismatch with {}'.format(self.id, s.id))
                except Exception as e:
                    log.error('{}: unknown freesasa error with {}'.format(self.id, s.id))
                    print(e)

    def find_disulfide_bridges(self, representative_only=True):
        """Run Biopython's disulfide bridge finder and store found bridges.
        Annotations are stored in the protein structure's chain sequence at:
        ``seq_record.annotations['SSBOND-biopython']``
        
        Args:
            representative_only (bool): If analysis should only be run on the representative structure

        """
        if representative_only:
            if self.representative_structure:
                try:
                    self.representative_structure.find_disulfide_bridges()
                except KeyError:
                    log.error('{}: unable to run disulfide bridge finder on {}'.format(self.id, self.representative_structure))
            else:
                log.warning('{}: no representative structure set, cannot run disulfide bridge finder'.format(self.id))
        else:
            for s in self.structures:
                try:
                    s.find_disulfide_bridges()
                except KeyError:
                    log.error('{}: unable to run disulfide bridge finder on {}'.format(self.id, s.id))

    def get_residue_annotations(self, seq_resnum, seqprop=None, structprop=None, chain_id=None,
                                use_representatives=False):
        """Get all residue-level annotations stored in the SeqProp letter_annotations field for a given residue number

        Uses the representative sequence, structure, and chain ID stored by default. If other properties from other
            structures are desired, input the proper IDs. An alignment for the given sequence to the structure must
            be present in the sequence_alignments list.

        Args:
            seq_resnum (int): Residue number in the sequence
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): ID of the structure's chain to get annotation from
            use_representatives (bool): If the representative sequence/structure/chain IDs should be used

        Returns:
            dict: All available letter_annotations for this residue number

        """

        if use_representatives:
            if seqprop and structprop and chain_id:
                raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                                 'Set use_representatives to False if custom IDs are to be used.')
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                                 'to True.')

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain

        chain = structprop.chains.get_by_id(chain_id)

        log.debug('Using sequence: {}, structure: {}, chain: {}'.format(seqprop.id, structprop.id, chain_id))

        # Create a new SeqFeature
        f = SeqFeature(FeatureLocation(seq_resnum-1, seq_resnum))

        # Get sequence properties
        seq_features = f.extract(seqprop)

        # Store in dictionary to return, clean it up
        all_info = ssbio.utils.clean_single_dict(indict=seq_features.letter_annotations,
                                                 prepend_to_keys='seq_',
                                                 remove_keys_containing='_chain_index')
        all_info['seq_resnum'] = seq_resnum
        all_info['seq_residue'] = str(seq_features.seq)

        if structprop:
            # Get structure properties
            mapping_to_structure_resnum = self.map_seqprop_resnums_to_structprop_resnums(resnums=seq_resnum,
                                                                                         seqprop=seqprop,
                                                                                         structprop=structprop,
                                                                                         chain_id=chain_id,
                                                                                         use_representatives=use_representatives)

            # Try finding the residue in the structure
            if f.location.end.position in mapping_to_structure_resnum:
                struct_resnum = mapping_to_structure_resnum[f.location.end.position]
                struct_f = SeqFeature(FeatureLocation(struct_resnum-1, struct_resnum))

                struct_seq_features = struct_f.extract(chain.seq_record)
                struct_info = ssbio.utils.clean_single_dict(indict=struct_seq_features.letter_annotations,
                                                            prepend_to_keys='struct_',
                                                            remove_keys_containing='structure_resnums')
                struct_info['struct_resnum'] = struct_resnum
                struct_info['struct_residue'] = str(struct_seq_features.seq)
                all_info.update(struct_info)

                # Warn if residue differs from sequence
                if seq_features.seq != struct_seq_features.seq:
                    log.warning('Sequence residue ({}{}) does not match structure residue ({}{}). '
                                'This may simply be due to differences in the structure'.format(seq_features.seq,
                                                                                                seq_resnum,
                                                                                                struct_seq_features.seq,
                                                                                                struct_resnum))

        return all_info

    def sequence_mutation_summary(self, alignment_ids=None, alignment_type=None):
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
            
        Args:
            alignment_ids (str, list): Specified alignment ID or IDs to use
            alignment_type (str): Specified alignment type contained in the ``annotation`` field of an alignment object,
                ``seqalign`` or ``structalign`` are the current types.

        Returns:
            dict, dict: single_counter, fingerprint_counter

        """
        if alignment_ids:
            ssbio.utils.force_list(alignment_ids)

        if len(self.sequence_alignments) == 0:
            log.error('{}: no sequence alignments'.format(self.id))
            return {}, {}

        fingerprint_counter = defaultdict(list)
        single_counter = defaultdict(list)

        for alignment in self.sequence_alignments:
            # Ignore alignments if a list of identifiers is provided
            if alignment_ids:
                if alignment.id not in alignment_ids:
                    continue
            # Ignore alignments if type is specified
            if alignment_type:
                if alignment.annotations['ssbio_type'] != alignment_type:
                    continue

            other_sequence = alignment.annotations['b_seq']
            mutations = alignment.annotations['mutations']

            if mutations:
                # Turn this list of mutations into a tuple so it can be a dictionary key
                mutations = tuple(tuple(x) for x in mutations)
                fingerprint_counter[mutations].append(other_sequence)

                for m in mutations:
                    single_counter[m].append(other_sequence)

        return dict(single_counter), dict(fingerprint_counter)

    def add_features_to_nglview(self, view, seqprop=None, structprop=None, chain_id=None, use_representatives=False):
        """Add select features from the selected SeqProp object to an NGLWidget view object.

        Currently parsing for:
        * Single residue features (ie. metal binding sites)
        * Disulfide bonds

        Args:
            view (NGLWidget): NGLWidget view object
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): ID of the structure's chain to get annotation from
            use_representatives (bool): If the representative sequence/structure/chain IDs should be used

        """
        if use_representatives:
            if seqprop and structprop and chain_id:
                raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                                 'Set use_representatives to False if custom IDs are to be used.')
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                                 'to True.')

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain

        # Parse and store chain seq if not already stored
        if not structprop.chains.has_id(chain_id):
            structprop.parse_structure()
            if not structprop.chains.has_id(chain_id):
                raise ValueError('Chain {} not present in structure {}'.format(chain_id, structprop.id))

        if not seqprop.features:
            log.warning('{}: no stored features'.format(seqprop.id))


        # Loop through any stored features
        for f in seqprop.features:

            # Display disulfide bonds
            if f.type.lower() == 'disulfide bond':
                # TODO: double check if .start or .start + 1
                disulfide = self.map_seqprop_resnums_to_structprop_resnums(resnums=[f.location.start + 1, f.location.end],
                                                                           seqprop=seqprop,
                                                                           structprop=structprop,
                                                                           chain_id=chain_id,
                                                                           use_representatives=use_representatives)
                to_view = [str(x)+'.CA' for x in list(disulfide.values())]
                view.add_distance(atom_pair=[to_view], color='black')
                log.info('Disulfide bridge at residues {} & {}'.format(f.location.start + 1, f.location.end))

            # Display DNA-binding regions
            if f.type.lower() == 'dna-binding region' or f.type.lower() == 'nucleotide phosphate-binding region':
                impres = self.map_seqprop_resnums_to_structprop_resnums(resnums=[f.location.start + 1,
                                                                                 f.location.end],
                                                                        seqprop=seqprop,
                                                                        structprop=structprop,
                                                                        chain_id=chain_id,
                                                                        use_representatives=use_representatives)

                # TODO: need to check if f.location.start was mapped and if not, try incrementing. or  input the list
                # of resnums, not just the start and end
                if f.location.start + 1 in impres and f.location.end in impres:
                    mapped_start = impres[f.location.start + 1]
                    mapped_end = impres[f.location.end]
                    view.add_ball_and_stick(selection=':{} and ( {}-{} )'.format(chain_id,
                                                                                 mapped_start,
                                                                                 mapped_end), color='black')
                    log.info('{} at sequence region {}-{}, structure residues {}-{}'.format(f.type,
                                                                                            f.location.start,
                                                                                            f.location.end,
                                                                                            mapped_start,
                                                                                            mapped_end))

            # Display other single residues
            if f.location.end - 1 == f.location.start:
                if f.type.lower() == 'sequence variant' or f.type.lower() == 'mutagenesis site':
                    continue
                impres = self.map_seqprop_resnums_to_structprop_resnums(resnums=f.location.end,
                                                                        seqprop=seqprop,
                                                                        structprop=structprop,
                                                                        chain_id=chain_id,
                                                                        use_representatives=use_representatives)
                if f.location.end in impres:
                    impres_mapped = impres[f.location.end]
                    view.add_ball_and_stick(selection=str(impres_mapped), color='black')
                    view.add_label(selection=':{} and {}'.format(chain_id, impres_mapped), label_type='res', color='black')
                    log.info('{} at sequence residue {}, structure residue {}'.format(f.type, f.location.end, impres_mapped))

    def add_mutations_to_nglview(self, view, alignment_type='seqalign', alignment_ids=None,
                                 seqprop=None, structprop=None, chain_id=None, use_representatives=False,
                                 grouped=False, color='red', unique_colors=True,
                                 opacity_range=(0.8,1), scale_range=(1,5)):
        """Add representations to an NGLWidget view object for residues that are mutated in the
            ``sequence_alignments`` attribute.

        Args:
            view (NGLWidget): NGLWidget view object
            alignment_type (str): Specified alignment type contained in the ``annotation`` field of an alignment object,
                ``seqalign`` or ``structalign`` are the current types.
            alignment_ids (str, list): Specified alignment ID or IDs to use
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): ID of the structure's chain to get annotation from
            use_representatives (bool): If the representative sequence/structure/chain IDs should be used
            grouped (bool): If groups of mutations should be colored and sized together
            color (str): Color of the mutations (overridden if unique_colors=True)
            unique_colors (bool): If each mutation/mutation group should be colored uniquely
            opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
            scale_range (tuple): Min/max size values (mutations that show up more will be bigger)

        """
        if use_representatives:
            if seqprop and structprop and chain_id:
                raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                                 'Set use_representatives to False if custom IDs are to be used.')
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                                 'to True.')

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain

        log.debug('Using sequence: {}, structure: {}, chain: {}'.format(seqprop.id, structprop.id, chain_id))

        # Get the summary of mutations
        single, fingerprint = self.sequence_mutation_summary(alignment_type=alignment_type, alignment_ids=alignment_ids)

        # Map residues from sequence to structure
        if not grouped:
            single_lens = {k: len(v) for k, v in single.items()}
            single_map_to_structure = {}

            for k, v in single_lens.items():
                resnum = int(k[1])
                resnum_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=resnum,
                                                                                     seqprop=seqprop,
                                                                                     structprop=structprop,
                                                                                     chain_id=chain_id,
                                                                                     use_representatives=use_representatives)
                if resnum not in resnum_to_structure:
                    log.warning('{}: residue is not available in structure {}'.format(resnum, structprop.id))
                    continue
                new_key = resnum_to_structure[resnum]
                single_map_to_structure[new_key] = v

            structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                                structure_resnums=single_map_to_structure,
                                                                chain=chain_id,
                                                                color=color,
                                                                unique_colors=unique_colors,
                                                                opacity_range=opacity_range,
                                                                scale_range=scale_range)
        else:
            log.warning('Viewing mutation groups is currently in beta -- groups may overwrite each other')
            fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
            fingerprint_map_to_structure = {}

            for k, v in fingerprint_lens.items():
                k_list = [int(x[1]) for x in k]
                resnums_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=k_list,
                                                                                      seqprop=seqprop,
                                                                                      structprop=structprop,
                                                                                      chain_id=chain_id,
                                                                                      use_representatives=use_representatives)
                new_key = tuple(y for y in resnums_to_structure.values())
                fingerprint_map_to_structure[new_key] = v

            structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                                structure_resnums=fingerprint_map_to_structure,
                                                                chain=chain_id,
                                                                color=color,
                                                                unique_colors=unique_colors,
                                                                opacity_range=opacity_range,
                                                                scale_range=scale_range)

    def add_fingerprint_to_nglview(self, view, fingerprint,
                                   seqprop=None, structprop=None, chain_id=None, use_representatives=False,
                                   color='red', opacity_range=(0.8, 1), scale_range=(1, 5)):
        """Add representations to an NGLWidget view object for residues that are mutated in the
            ``sequence_alignments`` attribute.

        Args:
            view (NGLWidget): NGLWidget view object
            fingerprint (dict): Single mutation group from the ``sequence_mutation_summary`` function
            seqprop (SeqProp): SeqProp object
            structprop (StructProp): StructProp object
            chain_id (str): ID of the structure's chain to get annotation from
            use_representatives (bool): If the representative sequence/structure/chain IDs should be used
            color (str): Color of the mutations (overridden if unique_colors=True)
            opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
            scale_range (tuple): Min/max size values (mutations that show up more will be bigger)

        """
        if use_representatives:
            if seqprop and structprop and chain_id:
                raise ValueError('Overriding sequence, structure, and chain IDs with representatives. '
                                 'Set use_representatives to False if custom IDs are to be used.')
        else:
            if not seqprop or not structprop or not chain_id:
                raise ValueError('Input sequence, structure, and chain to map between, or set use_representatives '
                                 'to True.')

        if use_representatives:
            seqprop = self.representative_sequence
            structprop = self.representative_structure
            chain_id = self.representative_chain

        log.debug('Using sequence: {}, structure: {}, chain: {}'.format(seqprop.id, structprop.id, chain_id))

        fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
        fingerprint_map_to_structure = {}
        for k, v in fingerprint_lens.items():
            k_list = [int(x[1]) for x in k]
            resnums_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=k_list,
                                                                                  seqprop=seqprop,
                                                                                  structprop=structprop,
                                                                                  chain_id=chain_id,
                                                                                  use_representatives=use_representatives)
            new_key = tuple(y for y in resnums_to_structure.values())
            fingerprint_map_to_structure[new_key] = v

        structprop.add_scaled_residues_highlight_to_nglview(view=view,
                                                            structure_resnums=fingerprint_map_to_structure,
                                                            chain=chain_id,
                                                            color=color,
                                                            opacity_range=opacity_range,
                                                            scale_range=scale_range)

    # def view_all_mutations(self, alignment_type, grouped=False, color='red', unique_colors=True, structure_opacity=0.5,
    #                        opacity_range=(0.8,1), scale_range=(1,5), gui=False):
    #     """Map all sequence alignment mutations to the structure.
    #
    #     Args:
    #         alignment_type (str): Specified alignment type contained in the ``annotation`` field of an alignment object,
    #             ``seqalign`` or ``structalign`` are the current types.
    #         grouped (bool): If groups of mutations should be colored and sized together
    #         color (str): Color of the mutations (overridden if unique_colors=True)
    #         unique_colors (bool): If each mutation/mutation group should be colored uniquely
    #         structure_opacity (float): Opacity of the protein structure cartoon representation
    #         opacity_range (tuple): Min/max opacity values (mutations that show up more will be opaque)
    #         scale_range (tuple): Min/max size values (mutations that show up more will be bigger)
    #         gui (bool): If the NGLview GUI should show up
    #
    #     Returns:
    #         NGLviewer object
    #
    #     """
    #     single, fingerprint = self.sequence_mutation_summary(alignment_type=alignment_type)
    #
    #     single_lens = {k: len(v) for k, v in single.items()}
    #     single_map_to_structure = {}
    #     for k, v in single_lens.items():
    #         resnum = int(k[1])
    #         resnum_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=resnum, use_representatives=True)
    #         if resnum not in resnum_to_structure:
    #             log.warning('{}: residue is not available in structure {}'.format(resnum, self.representative_structure.id))
    #             continue
    #         new_key = resnum_to_structure[resnum]
    #         single_map_to_structure[new_key] = v
    #
    #     if not grouped:
    #         view = self.representative_structure.view_structure_and_highlight_residues_scaled(single_map_to_structure,
    #                                                                                    color=color, unique_colors=unique_colors,
    #                                                                                    structure_opacity=structure_opacity,
    #                                                                                    opacity_range=opacity_range,
    #                                                                                    scale_range=scale_range,
    #                                                                                    gui=gui)
    #         return view
    #
    #     else:
    #         fingerprint_lens = {k: len(v) for k, v in fingerprint.items()}
    #         fingerprint_map_to_structure = {}
    #         for k, v in fingerprint_lens.items():
    #             k_list = [int(x[1]) for x in k]
    #             resnums_to_structure = self.map_seqprop_resnums_to_structprop_resnums(resnums=k_list, use_representatives=True)
    #             new_key = tuple(y for y in resnums_to_structure.values())
    #             fingerprint_map_to_structure[new_key] = v
    #
    #         view = self.representative_structure.view_structure_and_highlight_residues_scaled(fingerprint_map_to_structure,
    #                                                                                    color=color, unique_colors=unique_colors,
    #                                                                                    opacity_range=opacity_range,
    #                                                                                    scale_range=scale_range,
    #                                                                                    gui=gui)
    #         return view

    # def summarize_protein(self):
    #     """Gather all possible attributes in the sequences and structures and summarize everything.
    #
    #     Returns:
    #         dict:
    #
    #     """
    #     d = OrderedDict()
    #     repseq = self.representative_sequence
    #     if not self.representative_structure:
    #         repstruct = StructProp(self.id)
    #         repchain = repstruct.representative_chain
    #     else:
    #         repstruct = self.representative_structure
    #         repchain = self.representative_structure.representative_chain
    #     single, fingerprint = self.representative_sequence.sequence_mutation_summary()
    #     numstrains = len(self.sequences) - 1
    #
    #     d['Gene ID'] = self.id
    #     d['Number of sequences'] = len(self.sequences)
    #     d['Number of structures (total)'] = self.num_structures
    #     d['Number of structures (experimental)'] = self.num_structures_experimental
    #     d['Number of structures (homology models)'] = self.num_structures_homology
    #
    #     # d['------REPRESENTATIVE SEQUENCE PROPERTIES------')
    #     #     d['Sequence ID'] = repseq.id
    #     d['Sequence length'] = repseq.sequence_len
    #     d['Predicted number of transmembrane helices'] = repseq.seq_record.annotations['num_tm_helix-tmhmm']
    #
    #     # d['------REPRESENTATIVE STRUCTURE PROPERTIES------')
    #     d['Structure ID'] = repstruct.id
    #     #     d['Structure representative chain'] = format(repchain.id)))
    #     d['Structure is experimental'] = repstruct.is_experimental
    #     # d['Structure origin'] = repstruct.taxonomy_name))
    #     # d['Structure description'] = repstruct.description))
    #     if repstruct.reference_seq_top_coverage:
    #         d['Structure coverage of sequence'] = str(repstruct.reference_seq_top_coverage) + '%'
    #
    #     # d['------ALIGNMENTS SUMMARY------')
    #     #     d['Number of sequence alignments'] = len(repseq.sequence_alignments)))
    #     #     d['Number of structure alignments'] = len(repseq.structure_alignments)))
    #
    #     singles = []
    #     for k, v in single.items():
    #         k = [str(x) for x in k]
    #         if len(v) / numstrains >= 0.01:
    #             singles.append(''.join(k))  # len(v) is the number of strains
    #     d['Mutations that show up in more than 10% of strains'] = ';'.join(singles)
    #
    #     allfingerprints = []
    #     for k, v in fingerprint.items():
    #         if len(v) / numstrains >= 0.01:
    #             fingerprints = []
    #             for m in k:
    #                 y = [str(x) for x in m]
    #                 fingerprints.append(''.join(y))
    #             allfingerprints.append('-'.join(fingerprints))
    #     d['Mutation groups that show up in more than 10% of strains'] = ';'.join(allfingerprints)
    #
    #     return d

    def __json_decode__(self, **attrs):
        for k, v in attrs.items():
            if k == 'sequences' or k == 'structures' or k == 'sequence_alignments' or k == 'structure_alignments':
                setattr(self, k, DictList(v))
            else:
                setattr(self, k, v)