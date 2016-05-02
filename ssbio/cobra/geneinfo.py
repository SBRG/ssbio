from cobra.core import Gene
import pandas as pd
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from bioservices.uniprot import UniProt
bsup = UniProt()


class GeneInfo(Gene):
    """Extension of the cobra Gene class to provide ID mapping
    """

    def __init__(self, id=None):
        Gene.__init__(self, id)
        pass

if __name__ == '__main__':
    pass
