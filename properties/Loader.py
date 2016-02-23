from Bio import PDB
from Bio import Struct

# NOTE: MMCIFParser was modified according to
# https://github.com/biopython/biopython/issues/523
# this will throw an error with the normal Biopython distribution
cifp = PDB.MMCIFParser(QUIET=True)
pdbp = PDB.PDBParser(PERMISSIVE=True, QUIET=True)


class Loader():
    """
    Class to load either PDB or mmCIF files into a Struct/Protein object
    Please see GSOC 2010 Biopython project:
    http://biopython.org/wiki/Struct
    """

    def structure_reader(self, filename):
        if '.cif' in filename:
            s = cifp.get_structure('mycif', filename)
        else:
            s = pdbp.get_structure('mypdb', filename)

        return s

if __name__ == '__main__':
    l = Loader()
    l.structure_reader("test_structures/1u8f.cif")
