"""
StructProp
==========
"""

import logging
import os.path as op

import numpy as np
import pandas as pd
from cobra.core import DictList
from collections import defaultdict
import ssbio.protein.sequence.utils.alignment
import ssbio.protein.structure.properties.dssp
import ssbio.protein.structure.properties.msms
import ssbio.protein.structure.properties.residues
import ssbio.protein.structure.properties.quality
import ssbio.protein.structure.properties.freesasa as fs
import ssbio.utils
from ssbio.core.object import Object
from ssbio.protein.structure.chainprop import ChainProp
from ssbio.protein.structure.utils.structureio import StructureIO

import seaborn as sns
import ssbio.protein.structure.utils.cleanpdb
log = logging.getLogger(__name__)


class StructProp(Object):

    """Generic class to represent information for a protein structure.

    Provides access to the 3D coordinates using a Biopython Structure object through the method ``parse_structure``.
    The main functionality added is the ability to set and load directly from any supported structure and metadata
    file. Additionally, the ``mapped_chains`` attribute allows for analysis of a subset of chains, which will map
    to a gene of interest. Also provides methods through ``nglview`` to view the structure in a Jupyter notebook.

    Attributes:
        id (str): Unique identifier for this protein structure
        name (str): Optional name for this structure
        description (str): Optional description for this structure
        is_experimental (bool): Flag to note if this structure is an experimental model or a homology model
        chains (DictList): A DictList of chains have their sequence stored in them, along with residue-specific
            annotations
        mapped_chains (list): A simple list of chain IDs (strings) that will be used to subset analyses
        file_type (str): Type of structure file
        structure_file (str): Name of the structure file

    """

    def __init__(self, ident, description=None, chains=None, mapped_chains=None,
                 is_experimental=False, structure_path=None, file_type=None):
        """Initialize a StructProp object.
        
        Args:
            ident (str): Unique identifier for this structure 
            description (str): Optional human-readable description
            chains (str, list): Chain ID or list of IDs
            mapped_chains (str, list): A chain ID or IDs to indicate what chains should be analyzed
            is_experimental (bool): Flag to indicate if structure is an experimental or computational model
            structure_path (str): Path to structure file 
            file_type (str): Type of structure file - ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``,
                ``xml.gz``, ``mmtf``, ``mmtf.gz``

        """
        Object.__init__(self, id=ident, description=description)

        self.is_experimental = is_experimental

        # Chain information
        # chains is a DictList of ChainProp objects
        # If you run self.parse_structure(), all chains will be parsed and stored here
        # Use mapped_chains below to keep track of chains you are interested in
        self.chains = DictList()
        if chains:
            self.add_chain_ids(chains)
        # mapped_chains is an ordered list of mapped chain IDs which would come from BLAST or the best_structures API
        self.mapped_chains = []
        if mapped_chains:
            self.add_mapped_chain_ids(mapped_chains)

        # Simple flag to track if this structure has had its structure + chain sequences parsed
        self.parsed = False

        # File information
        self.file_type = file_type
        self._structure_dir = None
        self.structure_file = None
        if structure_path:
            self.load_structure_path(structure_path, file_type)

    @property
    def structure_dir(self):
        if not self._structure_dir:
            raise OSError('No sequence folder set')
        return self._structure_dir

    @structure_dir.setter
    def structure_dir(self, path):
        if path and not op.exists(path):
            raise OSError('{}: folder does not exist'.format(path))

        self._structure_dir = path

    @property
    def structure_path(self):
        if not self.structure_file:
            raise OSError('{}: structure file not available'.format(self.id))

        path = op.join(self.structure_dir, self.structure_file)
        if not op.exists(path):
            raise ValueError('{}: file does not exist'.format(path))
        return path

    def load_structure_path(self, structure_path, file_type):
        """Load a structure file and provide pointers to its location

        Args:
            structure_path (str): Path to structure file
            file_type (str): Type of structure file

        """

        if not file_type:
            raise ValueError('File type must be specified')

        self.file_type = file_type
        self.structure_dir = op.dirname(structure_path)
        self.structure_file = op.basename(structure_path)

    def parse_structure(self):
        """Read the 3D coordinates of a structure file and return it as a Biopython Structure object

        Also create ChainProp objects in the chains attribute

        Returns:
            Structure: Biopython Structure object

        """
        # TODO: perhaps add option to parse into ProDy object?
        if not self.structure_path:
            log.error('{}: no structure file, unable to parse'.format(self.id))
            return None
        else:
            # Add Biopython structure object
            structure = StructureIO(self.structure_path, self.file_type)

            # Add all chains to self.chains as ChainProp objects
            structure_chains = [x.id for x in structure.first_model.child_list]
            self.add_chain_ids(structure_chains)
            self.get_structure_seqs(structure.first_model)

            # Also add all chains to self.mapped_chains ONLY if there are none specified
            if not self.mapped_chains:
                self.add_mapped_chain_ids(structure_chains)

            self.parsed = True
            return structure

    def clean_structure(self, out_suffix='_clean', outdir=None, force_rerun=False,
                        remove_atom_alt=True, keep_atom_alt_id='A',remove_atom_hydrogen=True,  add_atom_occ=True,
                        remove_res_hetero=True, keep_chemicals=None, keep_res_only=None,
                        add_chain_id_if_empty='X', keep_chains=None):
        """Clean the structure file associated with this structure, and save it as a new file. Returns the file path.

        Args:
            out_suffix (str): Suffix to append to original filename
            outdir (str): Path to output directory
            force_rerun (bool): If structure should be re-cleaned if a clean file exists already
            remove_atom_alt (bool): Remove alternate positions
            keep_atom_alt_id (str): If removing alternate positions, which alternate ID to keep
            remove_atom_hydrogen (bool): Remove hydrogen atoms
            add_atom_occ (bool): Add atom occupancy fields if not present
            remove_res_hetero (bool): Remove all HETATMs
            keep_chemicals (str, list): If removing HETATMs, keep specified chemical names
            keep_res_only (str, list): Keep ONLY specified resnames, deletes everything else!
            add_chain_id_if_empty (str): Add a chain ID if not present
            keep_chains (str, list): Keep only these chains

        Returns:
            str: Path to cleaned PDB file

        """

        if not self.structure_path:
            log.error('{}: no structure file, unable to clean'.format(self.id))
            return None

        clean_pdb_file = ssbio.protein.structure.utils.cleanpdb.clean_pdb(self.structure_path, out_suffix=out_suffix,
                                                                          outdir=outdir, force_rerun=force_rerun,
                                                                          remove_atom_alt=remove_atom_alt,
                                                                          remove_atom_hydrogen=remove_atom_hydrogen,
                                                                          keep_atom_alt_id=keep_atom_alt_id,
                                                                          add_atom_occ=add_atom_occ,
                                                                          remove_res_hetero=remove_res_hetero,
                                                                          keep_chemicals=keep_chemicals,
                                                                          keep_res_only=keep_res_only,
                                                                          add_chain_id_if_empty=add_chain_id_if_empty,
                                                                          keep_chains=keep_chains)

        return clean_pdb_file

    def add_mapped_chain_ids(self, mapped_chains):
        """Add chains by ID into the mapped_chains attribute

        Args:
            mapped_chains (str, list): Chain ID or list of IDs

        """
        mapped_chains = ssbio.utils.force_list(mapped_chains)

        for c in mapped_chains:
            if c not in self.mapped_chains:
                self.mapped_chains.append(c)
                log.debug('{}: added to list of mapped chains'.format(c))
            else:
                log.debug('{}: chain already in list of mapped chains, not adding'.format(c))

    def add_chain_ids(self, chains):
        """Add chains by ID into the chains attribute

        Args:
            chains (str, list): Chain ID or list of IDs

        """
        chains = ssbio.utils.force_list(chains)

        for c in chains:
            if self.chains.has_id(c):
                log.debug('{}: chain already present'.format(c))
            else:
                chain_prop = ChainProp(ident=c, pdb_parent=self.id)
                self.chains.append(chain_prop)
                log.debug('{}: added to chains list'.format(c))

    def get_structure_seqs(self, model):
        """Gather chain sequences and store in their corresponding ``ChainProp`` objects in the ``chains`` attribute.

        Args:
            model (Model): Biopython Model object of the structure you would like to parse

        """

        # Don't overwrite existing ChainProp objects
        dont_overwrite = []
        chains = list(model.get_chains())
        for x in chains:
            if self.chains.has_id(x.id):
                if self.chains.get_by_id(x.id).seq_record:
                    dont_overwrite.append(x.id)
        if len(dont_overwrite) == len(chains):
            log.debug('Not writing structure sequences, already stored')
            return

        # Returns the structures sequences with Xs added
        structure_seqs = ssbio.protein.structure.properties.residues.get_structure_seqrecords(model)
        log.debug('{}: gathered chain sequences'.format(self.id))

        # Associate with ChainProps
        for seq_record in structure_seqs:
            log.debug('{}: adding chain sequence to ChainProp'.format(seq_record.id))
            my_chain = self.chains.get_by_id(seq_record.id)
            my_chain.seq_record = seq_record

    def reset_chain_seq_records(self):
        for x in self.chains:
            x.reset_seq_record()

    def get_dict_with_chain(self, chain, only_keys=None, chain_keys=None, exclude_attributes=None, df_format=False):
        """get_dict method which incorporates attributes found in a specific chain. Does not overwrite any attributes
        in the original StructProp.

        Args:
            chain:
            only_keys:
            chain_keys:
            exclude_attributes:
            df_format:

        Returns:
            dict: attributes of StructProp + the chain specified

        """

        # Choose attributes to return, return everything in the object if a list is not specified
        if not only_keys:
            keys = list(self.__dict__.keys())
        else:
            keys = ssbio.utils.force_list(only_keys)

        # Remove keys you don't want returned
        if exclude_attributes:
            exclude_attributes = ssbio.utils.force_list(exclude_attributes)
            for x in exclude_attributes:
                if x in keys:
                    keys.remove(x)
        else:
            exclude_attributes = []

        exclude_attributes.extend(['mapped_chains', 'chains'])

        final_dict = {k: v for k, v in Object.get_dict(self, only_attributes=keys, exclude_attributes=exclude_attributes,
                                                       df_format=df_format).items()}


        chain_prop = self.chains.get_by_id(chain)
        # Filter out keys that show up in StructProp
        if not chain_keys:
            chain_keys = [x for x in chain_prop.get_dict().keys() if x not in final_dict]

        chain_dict = chain_prop.get_dict(only_attributes=chain_keys, df_format=df_format)
        final_dict.update(chain_dict)

        return final_dict

    def find_disulfide_bridges(self, threshold=3.0):
        """Run Biopython's search_ss_bonds to find potential disulfide bridges for each chain and store in ChainProp."""
        parsed = self.parse_structure()
        if not parsed:
            log.error('{}: unable to open structure to find S-S bridges'.format(self.id))
            return

        disulfide_bridges = ssbio.protein.structure.properties.residues.search_ss_bonds(parsed.first_model,
                                                                                        threshold=threshold)
        if not disulfide_bridges:
            log.debug('{}: no disulfide bridges found'.format(self.id))

        for chain, bridges in disulfide_bridges.items():
            self.chains.get_by_id(chain).seq_record.annotations['SSBOND-biopython'] = disulfide_bridges[chain]
            log.debug('{}: found {} disulfide bridges'.format(chain, len(bridges)))
            log.debug('{}: stored disulfide bridges in the chain\'s seq_record letter_annotations'.format(chain))

    def get_dssp_annotations(self, outdir, force_rerun=False):
        """Run DSSP on this structure and store the DSSP annotations in the corresponding ChainProp SeqRecords

        Calculations are stored in the ChainProp's ``letter_annotations`` at the following keys:

            * ``SS-dssp``
            * ``RSA-dssp``
            * ``ASA-dssp``
            * ``PHI-dssp``
            * ``PSI-dssp``

        Args:
            outdir (str): Path to where DSSP dataframe will be stored.
            force_rerun (bool): If DSSP results should be recalculated

        TODO:
            * Also parse global properties, like total accessible surface area. Don't think Biopython parses those?

        """
        parsed = self.parse_structure()
        if not parsed:
            log.error('{}: unable to open structure to run DSSP'.format(self.id))
            return

        log.debug('{}: running DSSP'.format(self.id))
        dssp_results = ssbio.protein.structure.properties.dssp.get_dssp_df(model=parsed.first_model,
                                                                           pdb_file=self.structure_path,
                                                                           outdir=outdir,
                                                                           force_rerun=force_rerun)

        if dssp_results.empty:
            log.error('{}: unable to run DSSP'.format(self.id))
            return

        chains = dssp_results.chain.unique()
        dssp_summary = ssbio.protein.structure.properties.dssp.secondary_structure_summary(dssp_results)

        for chain in chains:
            ss = dssp_results[dssp_results.chain == chain].ss.tolist()
            exposure_rsa = dssp_results[dssp_results.chain == chain].exposure_rsa.tolist()
            exposure_asa = dssp_results[dssp_results.chain == chain].exposure_asa.tolist()
            phi = dssp_results[dssp_results.chain == chain].phi.tolist()
            psi = dssp_results[dssp_results.chain == chain].psi.tolist()

            chain_prop = self.chains.get_by_id(chain)
            chain_seq = chain_prop.seq_record

            # Making sure the X's are filled in
            ss = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                      new_seq=ss,
                                                                                      fill_with='-')
            exposure_rsa = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                                new_seq=exposure_rsa,
                                                                                                fill_with=float('Inf'))
            exposure_asa = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                                new_seq=exposure_asa,
                                                                                                fill_with=float('Inf'))
            phi = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                       new_seq=phi,
                                                                                       fill_with=float('Inf'))
            psi = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                       new_seq=psi,
                                                                                       fill_with=float('Inf'))

            chain_prop.seq_record.annotations.update(dssp_summary[chain])

            chain_prop.seq_record.letter_annotations['SS-dssp'] = ss
            chain_prop.seq_record.letter_annotations['RSA-dssp'] = exposure_rsa
            chain_prop.seq_record.letter_annotations['ASA-dssp'] = exposure_asa
            chain_prop.seq_record.letter_annotations['PHI-dssp'] = phi
            chain_prop.seq_record.letter_annotations['PSI-dssp'] = psi
            log.debug('{}: stored DSSP annotations in chain seq_record letter_annotations'.format(chain))

    def get_residue_depths(self, outdir, force_rerun=False):
        """Run MSMS on this structure and store the residue depths/ca depths in the corresponding ChainProp SeqRecords
        """
        # TODO: rename to get_msms_annotations
        if self.file_type != 'pdb':
            raise ValueError('{}: unable to run MSMS with "{}" file type. Please change file type to "pdb"'.format(self.id,
                                                                                                            self.file_type))

        parsed = self.parse_structure()
        if not parsed:
            log.error('{}: unable to open structure to run MSMS'.format(self.id))
            return

        log.debug('{}: running MSMS'.format(self.id))
        msms_results = ssbio.protein.structure.properties.msms.get_msms_df(model=parsed.first_model,
                                                                           pdb_file=self.structure_path,
                                                                           outdir=outdir, force_rerun=force_rerun)
        if msms_results.empty:
            log.error('{}: unable to run MSMS'.format(self.id))
            return

        chains = msms_results.chain.unique()

        for chain in chains:
            res_depths = msms_results[msms_results.chain == chain].res_depth.tolist()
            ca_depths = msms_results[msms_results.chain == chain].ca_depth.tolist()

            chain_prop = self.chains.get_by_id(chain)
            chain_seq = chain_prop.seq_record

            # Making sure the X's are filled in
            res_depths = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                              new_seq=res_depths,
                                                                                              fill_with=float('Inf'))

            ca_depths = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain_seq,
                                                                                             new_seq=ca_depths,
                                                                                             fill_with=float('Inf'))

            chain_prop.seq_record.letter_annotations['RES_DEPTH-msms'] = res_depths
            chain_prop.seq_record.letter_annotations['CA_DEPTH-msms'] = ca_depths
            log.debug('{}: stored residue depths in chain seq_record letter_annotations'.format(chain))

    def get_freesasa_annotations(self, outdir, include_hetatms=False, force_rerun=False):
        """Run ``freesasa`` on this structure and store the calculated properties in the corresponding ChainProps
        """
        if self.file_type != 'pdb':
            log.error('{}: unable to run freesasa with "{}" file type. Please change file type to "pdb"'.format(self.id,
                                                                                                                self.file_type))
            return

        # Parse the structure to store chain sequences
        parsed = self.parse_structure()
        if not parsed:
            log.error('{}: unable to open structure to run freesasa'.format(self.id))
            return

        # Set outfile name
        log.debug('{}: running freesasa'.format(self.id))
        if include_hetatms:
            outfile = '{}.freesasa_het.rsa'.format(self.id)
        else:
            outfile = '{}.freesasa_nohet.rsa'.format(self.id)

        # Run freesasa
        result = fs.run_freesasa(infile=self.structure_path,
                                 outfile=outfile,
                                 include_hetatms=include_hetatms,
                                 outdir=outdir,
                                 force_rerun=force_rerun)

        # Parse results
        result_parsed = fs.parse_rsa_data(result)
        prop_dict = defaultdict(lambda: defaultdict(list))
        for k, v in result_parsed.items():
            chain = k[0]
            for prop, calc in v.items():
                prop_dict[chain][prop].append(calc)

        # Reorganize and store results
        all_props = ['all_atoms_abs', 'all_atoms_rel', 'side_chain_abs', 'side_chain_rel', 'main_chain_abs',
                     'main_chain_rel', 'non_polar_abs', 'non_polar_rel', 'all_polar_abs', 'all_polar_rel']
        all_props_renamed = {'all_atoms_abs' : 'ASA_ALL-freesasa',
                             'all_atoms_rel' : 'RSA_ALL-freesasa',
                             'all_polar_abs' : 'ASA_POLAR-freesasa',
                             'all_polar_rel' : 'RSA_POLAR-freesasa',
                             'main_chain_abs': 'ASA_BACKBONE-freesasa',
                             'main_chain_rel': 'RSA_BACKBONE-freesasa',
                             'non_polar_abs' : 'ASA_NONPOLAR-freesasa',
                             'non_polar_rel' : 'RSA_NONPOLAR-freesasa',
                             'side_chain_abs': 'ASA_RESIDUE-freesasa',
                             'side_chain_rel': 'RSA_RESIDUE-freesasa'}

        ## Rename dictionary keys based on if HETATMs were included
        if include_hetatms:
            suffix = '_het'
        else:
            suffix = '_nohet'

        for k, v in all_props_renamed.items():
            all_props_renamed[k] = v + suffix

        for chain in self.chains:
            for prop in all_props:
                prop_list = ssbio.protein.structure.properties.residues.match_structure_sequence(orig_seq=chain.seq_record,
                                                                                                 new_seq=prop_dict[chain.id][prop],
                                                                                                 fill_with=float('Inf'),
                                                                                                 ignore_excess=True)
                chain.seq_record.letter_annotations[all_props_renamed[prop]] = prop_list
            log.debug('{}: stored freesasa calculations in chain seq_record letter_annotations'.format(chain))

    def view_structure(self, only_chains=None, opacity=1.0, recolor=False, gui=False):
        """Use NGLviewer to display a structure in a Jupyter notebook

        Args:
            only_chains (str, list): Chain ID or IDs to display
            opacity (float): Opacity of the structure
            recolor (bool): If structure should be cleaned and recolored to silver
            gui (bool): If the NGLview GUI should show up

        Returns:
            NGLviewer object

        """
        # TODO: show_structure_file does not work for MMTF files - need to check for that and load accordingly

        if ssbio.utils.is_ipynb():
            import nglview as nv
        else:
            raise EnvironmentError('Unable to display structure - not running in a Jupyter notebook environment')

        if not self.structure_path:
            raise ValueError("Structure file not loaded")

        only_chains = ssbio.utils.force_list(only_chains)
        to_show_chains = '( '
        for c in only_chains:
            to_show_chains += ':{} or'.format(c)
        to_show_chains = to_show_chains.strip(' or ')
        to_show_chains += ' )'

        if self.file_type == 'mmtf' or self.file_type == 'mmtf.gz':
            view = nv.NGLWidget()
            view.add_component(self.structure_path)
        else:
            view = nv.show_structure_file(self.structure_path, gui=gui)

        if recolor:
            view.clear_representations()
            if only_chains:
                view.add_cartoon(selection='{} and (not hydrogen)'.format(to_show_chains), color='silver', opacity=opacity)
            else:
                view.add_cartoon(selection='protein', color='silver', opacity=opacity)
        elif only_chains:
            view.clear_representations()
            view.add_cartoon(selection='{} and (not hydrogen)'.format(to_show_chains), color='silver', opacity=opacity)

        return view

    def add_residues_highlight_to_nglview(self, view, structure_resnums, chain=None, res_color='red'):
        """Add a residue number or numbers to an NGLWidget view object.

        Args:
            view (NGLWidget): NGLWidget view object
            structure_resnums (int, list): Residue number(s) to highlight, structure numbering
            chain (str, list): Chain ID or IDs of which residues are a part of. If not provided, all chains in the
                mapped_chains attribute will be used. If that is also empty, and exception is raised.
            res_color (str): Color to highlight residues with

        """
        if not chain:
            chain = self.mapped_chains
            if not chain:
                raise ValueError('Please input chain ID to display residue on')

        if isinstance(structure_resnums, list):
            structure_resnums = list(set(structure_resnums))
        elif isinstance(structure_resnums, int):
            structure_resnums = ssbio.utils.force_list(structure_resnums)
        else:
            raise ValueError('Input must either be a residue number of a list of residue numbers')

        to_show_chains = '( '
        for c in chain:
            to_show_chains += ':{} or'.format(c)
        to_show_chains = to_show_chains.strip(' or ')
        to_show_chains += ' )'

        to_show_res = '( '
        for m in structure_resnums:
            to_show_res += '{} or '.format(m)
        to_show_res = to_show_res.strip(' or ')
        to_show_res += ' )'

        log.info('Selection: {} and not hydrogen and {}'.format(to_show_chains, to_show_res))

        view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res), color=res_color)

    def add_scaled_residues_highlight_to_nglview(self, view, structure_resnums, chain=None, color='red',
                                                 unique_colors=False, opacity_range=(0.5,1), scale_range=(.7, 10)):
        """Add a list of residue numbers (which may contain repeating residues) to a view, or add a dictionary of
            residue numbers to counts. Size and opacity of added residues are scaled by counts.

        Args:
            view (NGLWidget): NGLWidget view object
            structure_resnums (int, list, dict): Residue number(s) to highlight, or a dictionary of residue number to
                frequency count
            chain (str, list): Chain ID or IDs of which residues are a part of. If not provided, all chains in the
                mapped_chains attribute will be used. If that is also empty, and exception is raised.
            color (str): Color to highlight residues with
            unique_colors (bool): If each mutation should be colored uniquely (will override color argument)
            opacity_range (tuple): Min/max opacity values (residues that have higher frequency counts will be opaque)
            scale_range (tuple): Min/max size values (residues that have higher frequency counts will be bigger)

        """
        # TODO: likely to move these functions to a separate nglview/utils folder since they are not coupled to the structure
        # TODO: add color by letter_annotations!
        if not chain:
            chain = self.mapped_chains
            if not chain:
                raise ValueError('Please input chain ID to display residue on')
        else:
            chain = ssbio.utils.force_list(chain)

        if isinstance(structure_resnums, dict):
            opacity_dict = ssbio.utils.scale_calculator(opacity_range[0], structure_resnums, rescale=opacity_range)
            scale_dict = ssbio.utils.scale_calculator(scale_range[0], structure_resnums, rescale=scale_range)
        else:
            opacity_dict = {x: max(opacity_range) for x in ssbio.utils.force_list(structure_resnums)}
            scale_dict = {x: max(scale_range) for x in ssbio.utils.force_list(structure_resnums)}

        if isinstance(structure_resnums, list):
            structure_resnums = list(set(structure_resnums))
        elif isinstance(structure_resnums, dict):
            structure_resnums = list(structure_resnums.keys())
        elif isinstance(structure_resnums, int):
            structure_resnums = ssbio.utils.force_list(structure_resnums)
        else:
            raise ValueError('Input must either be a list of residue numbers or a dictionary of residue numbers '
                             'and their frequency.')

        colors = sns.color_palette("hls", len(structure_resnums)).as_hex()

        to_show_chains = '( '
        for c in chain:
            to_show_chains += ':{} or'.format(c)
        to_show_chains = to_show_chains.strip(' or ')
        to_show_chains += ' )'

        for i, x in enumerate(structure_resnums):
            if isinstance(x, tuple):
                to_show_res = '( '
                for mut in x:
                    to_show_res += '{} or '.format(mut)
                to_show_res = to_show_res.strip(' or ')
                to_show_res += ' )'
            else:
                to_show_res = x

            log.info('Selection: {} and not hydrogen and {}'.format(to_show_chains, to_show_res))

            if unique_colors:
                view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res),
                                        color=colors[i], opacity=opacity_dict[x], scale=scale_dict[x])
            else:
                view.add_ball_and_stick(selection='{} and not hydrogen and {}'.format(to_show_chains, to_show_res),
                                        color=color, opacity=opacity_dict[x], scale=scale_dict[x])

    def __json_decode__(self, **attrs):
        for k, v in attrs.items():
            if k == 'chains':
                setattr(self, k, DictList(v))
            else:
                setattr(self, k, v)