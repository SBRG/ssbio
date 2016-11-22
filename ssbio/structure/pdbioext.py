import ssbio.utils
from ssbio.structure.bp_mmcifparser import MMCIFParserFix
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

import logging
log = logging.getLogger(__name__)

cifp = MMCIFParserFix(QUIET=True)
pdbp = PDBParser(PERMISSIVE=True, QUIET=True)


class PDBIOExt(PDBIO):
    """Extended class to load a PDB or mmCIF file into a Biopython PDBIO object.

    Loads the first model when there are multiple available.
    Also adds some logging methods.
    """

    def __init__(self, structure_file, file_type='pdb'):
        super(PDBIOExt, self).__init__()

        self.structure_file = structure_file

        # Load the structure
        if file_type.lower() == 'pdb':
            structure = pdbp.get_structure(id='ssbio_pdb', file=structure_file)
        if file_type.lower() == 'mmcif' or file_type.lower() == 'cif':
            structure = cifp.get_structure(structure_id='ssbio_cif', filename=structure_file)

        # If there are multiple models (NMR), use the first model as the structure
        if len(structure) > 1:
            structure = structure[0]
            log.debug('{}: using first model'.format(structure_file))

        if len(structure) == 0:
            log.error('{}: no models in structure!'.format(structure_file))

        # Set this structure as the main one
        self.set_structure(structure)

    def write_pdb(self, custom_name='', custom_ext='', out_suffix='new', out_dir=None, custom_selection=None):
        """Write a new PDB file for the Structure's FIRST MODEL.

        Set custom_selection to a PDB.Select class for custom SMCRA selections.

        Args:
            out_suffix: string to append to new PDB file - default is "_new"
            out_dir: optional directory to output the file
            custom_selection: optional custom selection class

        Returns:
            out_file: filepath of new PDB file

        """

        # Prepare the output file path
        outfile = ssbio.utils.outfile_name_maker(inname=self.structure_file, outext=custom_ext, outfile=custom_name, outdir=out_dir, append_to_name=out_suffix)
        self.save(outfile, custom_selection)

        return outfile


if __name__ == '__main__':
    pass
