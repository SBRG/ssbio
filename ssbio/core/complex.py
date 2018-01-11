"""
Complex
=======
"""

import os.path as op
import logging
import ssbio.utils
from ssbio.core.object import Object

log = logging.getLogger(__name__)


class Complex(Object):

    """Store information about a protein complex, which is composed of individual protein subunits
    """

    def __init__(self, ident, subunits, description=None, root_dir=None, pdb_file_type='mmtf'):
        Object.__init__(self, id=ident, description=description)

        self.subunits = subunits
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