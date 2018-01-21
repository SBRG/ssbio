"""
Complex
=======
"""

import os.path as op
import logging
import ssbio.utils
from ssbio.core.object import Object
from ssbio.core.protein import Protein
from cobra.core import DictList

log = logging.getLogger(__name__)


class Complex(Object):

    """Store information about a protein complex, a generic representation of a 3D oligomeric complex composed of
    individual protein subunits.

    The main utilities of this class are to:

    * Allow as input a name for the complex and a dictionary of its subunit composition
    * Map each single subunit to its available experimental structures and homology models using methods in the
      :class:`~ssbio.core.protein.Protein` class
    * Map experimental structures and homology models to their available oliogmeric states
    * Select a single :attr:`~ssbio.core.complex.Complex.representative_complex` to which best represents the 3D
      oligomeric structure that best matches the defined subunit composition
    * Calculate, store, and access properties related to this complex
    * Provide summaries of available structures and selection details for the representative complex

    Args:
        ident (str): Unique identifier for this protein
        subunits (dict): Subunit composition defined as ``{protein_subunit_id: number_of_times_used_in_complex}``
        description (str): Optional description for this complex
        root_dir (str): Path to where the folder named by this complex's ID will be created.
            Default is current working directory.
        pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` -
            choose a file type for files downloaded from the PDB

    """

    def __init__(self, ident, subunits, description=None, root_dir=None, pdb_file_type='mmtf'):
        Object.__init__(self, id=ident, description=description)

        self.subunit_dict = subunits
        """dict: Subunit composition defined as ``{protein_subunit_id: number_of_times_used_in_complex}``"""

        self._subunits = None
        self._oligomeric_state = None

        self.pdb_file_type = pdb_file_type
        """str: ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` - choose a file 
        type for files downloaded from the PDB"""

        # Create directories
        self._root_dir = None
        if root_dir:
            self.root_dir = root_dir

    @property
    def root_dir(self):
        """str: Path to where the folder named by this complex's ID will be created. Default is current working
        directory."""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.debug('Changing root directory of Complex "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)) and not op.exists(op.join(path, self.id) + '_complex'):
                raise IOError('Complex "{}" does not exist in folder {}'.format(self.id, path))

        self._root_dir = path

        for d in [self.complex_dir]:
            ssbio.utils.make_dir(d)

    @property
    def complex_dir(self):
        """str: Complex folder"""
        if self.root_dir:
            # Add a _complex suffix to the folder if it has the same name as its root folder
            folder_name = self.id
            if folder_name == op.basename(self.root_dir):
                folder_name = self.id + '_complex'
            return op.join(self.root_dir, folder_name)
        else:
            log.warning('Root directory not set')
            return None

    @property
    def oligomeric_state(self):
        """Return the oligomeric state based on the contents of the :attr:`~ssbio.core.complex.Complex.subunit_dict`
        dictionary.

        Returns:
            str: Oligomeric state of the complex, currently can be ``monomer``, ``homomer``, or ``heteromer``

        """
        # TODO: [VizRecon]
        # Check the number of keys in the self.subunit_dict dictionary
        # Set as monomer if there is one key, and the value == 1
        # Set as homomer if there is one key and value > 1
        # Set as heteromer otherwise
        return None

    @property
    def subunits(self):
        """DictList: Subunits represented as a DictList of Protein objects"""

        # TODO: [VizRecon]
        # TODO: will need to adapt this to allow for input of previously created Protein objects

        subunits = DictList()

        for s in self.subunit_dict:
            subunits.append(Protein(ident=s, description='Subunit of complex {}'.format(self.id),
                                    root_dir=self.complex_dir, pdb_file_type=self.pdb_file_type))

        return subunits

    def map_subunit_to_sequence_and_structures(self, subunit_id):
        """Run the sequence and structure mapping code for a specified protein subunit.

        This stores the mapping information directly inside the Protein subunit object itself

        Args:
            subunit_id (str): ID of protein subunit to run mapping code for

        """
        # TODO: Nathan

    def set_representative_complex(self):
        """Set the representative 3D structure for this complex based on coverage of subunits.

        Args:

        Returns:

        """
        # TODO: [VizRecon]

        if self.oligomeric_state == 'monomer':
            pass

        elif self.oligomeric_state == 'homomer':
            pass

        elif self.oligomeric_state == 'heteromer':
            pass


def get_structure_stoichiometry(structure_file):
    """Parse a structure file and return a chain stoichiometry dictionary.

    Args:
        structure_file (str): Path to protein structure file (experimental or homology model

    Returns:
        dict: Dictionary of ``{chain_id: number_of_times_used}``

    """
    # TODO: [VizRecon]
    # TODO: most likely to be moved to another module, need to brainstorm
