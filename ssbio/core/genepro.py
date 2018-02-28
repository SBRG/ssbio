from cobra.core import Gene
from ssbio.core.protein import Protein
import os.path as op
import ssbio.utils
import cobra.io.dict
import logging
log = logging.getLogger(__name__)

class GenePro(Gene):

    """Extends the COBRAPy Gene object to add GEM-PRO annotations in the form of a Protein attribute.

    Additionally adds optional directory paths for files.

    Args:
        id (str): ID for the gene
        name (str): Name of the gene
        functional (bool): If this gene is functional in the linked genome-scale model
        root_dir (str): Path to directory to store protein related files
        pdb_file_type (str): Type of PDB structure file to use

    """

    def __init__(self, id, name='', functional=True, root_dir=None, pdb_file_type='mmtf'):
        Gene.__init__(self, id=id, name=name, functional=functional)

        self.pdb_file_type = pdb_file_type

        # Create directories
        self._root_dir = None
        if root_dir:
            self.root_dir = root_dir

        self.protein = Protein(ident=id, root_dir=self.gene_dir, pdb_file_type=self.pdb_file_type)

    @property
    def root_dir(self):
        """Directory where Gene folder is located"""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            raise ValueError('No path specified')

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.debug('Changing root directory of Gene "{}" from {} to {}'.format(self.id, self.root_dir, path))

            if not op.exists(op.join(path, self.id)):
                raise IOError('Gene "{}" does not exist in folder {}'.format(self.id, path))

        self._root_dir = path

        ssbio.utils.make_dir(self.gene_dir)

        # Propagate changes to protein
        if hasattr(self, 'protein'):
            self.protein.root_dir = self.gene_dir

    @property
    def gene_dir(self):
        """Gene folder"""
        if self.root_dir:
            return op.join(self.root_dir, self.id)
        else:
            return None

    def copy_modified_gene(self, modified_gene, ignore_model_attributes=True):
        """Copy attributes of a Gene object over to this Gene, given that the modified gene has the same ID.

        Args:
            modified_gene (Gene, GenePro): Gene with modified attributes that you want to copy over.
            ignore_model_attributes (bool): If you want to ignore copying over attributes related to metabolic models.

        """
        ignore = ['_model', '_reaction', '_functional', 'model', 'reaction', 'functional']
        for attr in filter(lambda a: not a.startswith('__') and not isinstance(getattr(type(self), a, None), property) and not callable(getattr(self, a)),
                           dir(modified_gene)):
            if attr not in ignore and ignore_model_attributes:
                setattr(self, attr, getattr(modified_gene, attr))

    def reset_protein(self):
        self.protein = Protein(self.id, root_dir=self.gene_dir)

    def __json_encode__(self):
        reqd_attribs = cobra.io.dict._REQUIRED_GENE_ATTRIBUTES
        optn_attribs = cobra.io.dict._OPTIONAL_GENE_ATTRIBUTES
        ordered_optn_attribs = cobra.io.dict._ORDERED_OPTIONAL_GENE_KEYS
        optn_attribs.update({'_root_dir': None})
        ordered_optn_attribs.append('_root_dir')
        new_gene = {key: str(getattr(self, key))
                    for key in reqd_attribs}
        cobra.io.dict._update_optional(self, new_gene, optn_attribs, ordered_optn_attribs)

        # Protein
        new_gene['protein'] = self.protein

        return new_gene

    def __json_decode__(self, **attrs):
        Gene.__init__(self, id=attrs['id'])
        for k, v in attrs.items():
            if k not in ['id']:
                setattr(self, k, v)
