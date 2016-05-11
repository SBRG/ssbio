import os
import os.path as op
from Bio import PDB
from Bio.PDB import PDBIO

# NOTE: MMCIFParser needs to be modified according to
# https://github.com/biopython/biopython/issues/523
# Otherwise this will throw an error with the normal Biopython distribution
cifp = PDB.MMCIFParser(QUIET=True)
pdbp = PDB.PDBParser(PERMISSIVE=True, QUIET=True)


class PDBIOExt(PDBIO):
    """Class to load either PDB or mmCIF files into a Biopython Structure object
    """

    def __init__(self, in_file):
        super().__init__()
        self.in_file = in_file

        if '.cif' in self.in_file:
            structure = cifp.get_structure('mycif', self.in_file)
        else:
            structure = pdbp.get_structure('mypdb', self.in_file)
        self.structure = structure
        self.first_model = structure[0]

    def _output_filepath(self, out_suffix, out_dir=None):
        """Returns an output file path based on the input filename to write a modified file.

        Args:
            out_suffix (str): string to append to the filename of the new PDB file
            out_dir (str): optional working directory where cleaned PDB file should be written to

        Returns:
            out_file (str): file path of the new PDB file
        """

        # Parsing the input filename
        filename_full = op.basename(self.in_file)
        filename, ext = op.splitext(filename_full)

        # Assembling the new output filename
        out_file = '{}_{}{}'.format(filename, out_suffix, ext)

        if out_dir:
            if not op.exists(out_dir):
                os.mkdir(out_dir)
            out_file = op.join(out_dir, out_file)

        return out_file

    def write_pdb(self, out_suffix='new', out_dir=None, custom_selection=None):
        """Write a new PDB file for the first model of this Structure appended with a suffix.

        Set custom_selection to a PDB.Select class for custom SMCRA selections

        Args:
            out_suffix: string to append to new PDB file - default is "_new"
            out_dir: optional directory to output the file
            custom_selection: optional custom selection class

        Returns:
            out_file: filepath of new PDB file

        """

        # Prepare the output file path
        out_file = self._output_filepath(out_suffix=out_suffix, out_dir=out_dir)
        self.set_structure(self.first_model)
        self.save(out_file, custom_selection)

        return out_file


if __name__ == '__main__':
    pass
