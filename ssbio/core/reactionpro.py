from cobra.core import Reaction
from ssbio.core.complex import Complex
import os.path as op
import ssbio.utils
import logging
log = logging.getLogger(__name__)

class ReactionPro(Reaction):

    """Extends the COBRAPy Reaction object to add GEM-PRO annotations in the form of a Complex attribute.

    Additionally adds optional directory paths for files.

    """

    def __init__(self, id=None, name='', subsystem='', lower_bound=0., upper_bound=1000., objective_coefficient=0.,
                 root_dir=None, pdb_file_type='mmtf'):
        Reaction.__init__(self, id=id, name=name, subsystem=subsystem, lower_bound=lower_bound, upper_bound=upper_bound,
                          objective_coefficient=objective_coefficient)

        self.pdb_file_type = pdb_file_type

        # Create directories
        self._root_dir = None
        if root_dir:
            self.root_dir = root_dir

        self.complex = None

    @property
    def root_dir(self):
        """Directory where Reaction folder is located"""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.debug('Changing root directory of Reaction "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)):
                raise IOError('Reaction "{}" does not exist in folder {}'.format(self.id, path))

        self._root_dir = path

        ssbio.utils.make_dir(self.reaction_dir)

        # Propagate changes to complex
        if hasattr(self, 'complex'):
            self.complex.root_dir = self.reaction_dir

    @property
    def reaction_dir(self):
        """Reaction folder"""
        if self.root_dir:
            return op.join(self.root_dir, self.id)
        else:
            return None

    def set_complex(self):
        self.complex = Complex(self.id, subunits=self.gene_reaction_rule,
                               root_dir=self.reaction_dir, pdb_file_type=self.pdb_file_type)
