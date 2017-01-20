from ssbio.core.object import Object
import logging
log = logging.getLogger(__name__)


class ChainProp(Object):
    """Class for protein structural properties of a specific chain"""

    def __init__(self, ident, pdb_parent, seq_record=None, description=None):
        Object.__init__(self, id=ident, description=description)

        self.pdb_parent = pdb_parent
        self.seq_record = seq_record