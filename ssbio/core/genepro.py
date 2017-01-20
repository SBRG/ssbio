from cobra.core import Gene
from ssbio.core.protein import Protein


class GenePro(Gene):
    """Extends the COBRAPy Gene object to add GEM-PRO annotations
    """

    def __init__(self, id, name='', functional=True):
        Gene.__init__(self, id=id, name=name, functional=functional)
        self.protein = Protein(ident=id)