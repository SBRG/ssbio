import os
import os.path as op
from Bio import PDB
from Bio.PDB import PDBIO
import logging
log = logging.getLogger(__name__)

# NOTE: MMCIFParser needs to be modified according to
# https://github.com/biopython/biopython/issues/523
# Otherwise this will throw an error with the normal Biopython distribution
cifp = PDB.MMCIFParser(QUIET=True)
# TODO: fix the mmcif stuff and add a new mmcifparser class in ssbio
# see for fixes: https://github.com/nmih/biopython/commit/d03425227fac016f3d4d59ac6d3332fec8ebfead

pdbp = PDB.PDBParser(PERMISSIVE=True, QUIET=True)


class PDBIOExt(PDBIO):
    """Class to load either PDB or mmCIF files into a Biopython Structure object.
        This is needed for the CleanPDB and MutatePDB classes, in order to use the custom_selection flag.
    """

    def __init__(self, in_file):
        super(PDBIO, self).__init__()
        self.in_file = in_file

        log.debug('{}: Loading structure...'.format(in_file))
        # TODO: make infile type explicit instead of parsing the file name
        if '.cif' in self.in_file:
            structure = cifp.get_structure('mycif', self.in_file)
        else:
            structure = pdbp.get_structure('mypdb', self.in_file)

        self.structure = structure
        self.first_model = structure[0]

        # TODO: need to properly learn about extending a class
        self.use_model_flag=0

    # TODO: replace with function in utils module
    def _output_filepath(self, custom_name='', custom_ext='', out_suffix='', out_dir=None):
        """Returns an output file path based on the input filename to write a modified file.

        Args:
            custom_name (str): optional custom name to name the file
            out_suffix (str): optional string to append to the filename of the new PDB file
            out_dir (str): optional working directory where cleaned PDB file should be written to

        Returns:
            out_file (str): file path of the new PDB file
        """

        # Parsing the input filename
        filename_full = op.basename(self.in_file)
        filename, ext = op.splitext(filename_full)

        # TODO: improve logic here
        if custom_name:
            filename = custom_name
        if custom_ext:
            ext = custom_ext

        # Assembling the new output filename
        if out_suffix:
            out_file = '{}_{}{}'.format(filename, out_suffix, ext)
        else:
            out_file = '{}{}'.format(filename, ext)

        if out_dir:
            if not op.exists(out_dir):
                os.mkdir(out_dir)
            out_file = op.join(out_dir, out_file)

        return out_file

    def write_pdb(self, custom_name='', custom_ext='', out_suffix='new', out_dir=None, custom_selection=None):
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
        out_file = self._output_filepath(custom_name=custom_name, custom_ext=custom_ext, out_suffix=out_suffix, out_dir=out_dir)
        self.set_structure(self.first_model)
        self.save(out_file, custom_selection)

        return out_file


if __name__ == '__main__':
    pass
