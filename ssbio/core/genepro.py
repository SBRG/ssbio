from cobra.core import Gene
from ssbio.core.protein import Protein
import os.path as op
import ssbio.utils
import cobra.io.dict
import logging
log = logging.getLogger(__name__)

class GenePro(Gene):
    """Extends the COBRAPy Gene object to add GEM-PRO annotations in the form of a Protein attribute
    """

    def __init__(self, id, name='', functional=True, root_dir=None, pdb_file_type='cif'):
        Gene.__init__(self, id=id, name=name, functional=functional)
        self.protein = Protein(ident=id, root_dir=root_dir, pdb_file_type=pdb_file_type)

        self.root_dir = root_dir

    @property
    def root_dir(self):
        """Directory where Gene folder is located"""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        self._root_dir = path

        if path:
            if not op.exists(path):
                raise ValueError('{}: folder does not exist'.format(path))

            ssbio.utils.make_dir(self.gene_dir)

            # Additionally set the protein folder path
            self.protein.root_dir = self.gene_dir

    @property
    def gene_dir(self):
        """Gene folder"""
        if self.root_dir:
            return op.join(self.root_dir, self.id)
        else:
            log.warning('Root directory not set')
            return None

    def reset_protein(self):
        self.protein = Protein(self.id)

    def __json_encode__(self):
        reqd_attribs = cobra.io.dict._REQUIRED_GENE_ATTRIBUTES
        optn_attribs = cobra.io.dict._OPTIONAL_GENE_ATTRIBUTES
        optn_attribs.update({'root_dir': None})
        new_gene = {key: str(getattr(self, key))
                    for key in reqd_attribs}
        cobra.io.dict._update_optional(self, new_gene, optn_attribs)

        # Protein
        new_gene['protein'] = self.protein

        return new_gene

    def __json_decode__(self, **attrs):
        Gene.__init__(self, id=attrs['id'])
        self.protein = attrs['protein']
        for k, v in cobra.io.dict._OPTIONAL_GENE_ATTRIBUTES.items():
            if k not in attrs:
                setattr(self, k, v)
        for k, v in attrs.items():
            if k == 'root_dir':
                if op.exists(v):
                    setattr(self, k, v)
                else:
                    log.debug('Directory does not exist, files will not be mapped')
                    continue
            if k != 'protein' or k != 'id':
                setattr(self, k, v)
