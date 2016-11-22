import ssbio.utils
from ssbio.structure.bp_mmcifparser import MMCIFParserFix
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select

import logging
log = logging.getLogger(__name__)

cifp = MMCIFParserFix(QUIET=True)
pdbp = PDBParser(PERMISSIVE=True, QUIET=True)


class PDBIOExt(PDBIO):
    """Extended class to load a PDB or mmCIF file into a Biopython PDBIO object.

    Loads the first model when there are multiple available.
    Also adds some logging methods.
    """

    def __init__(self, structure_file, file_type):
        super(PDBIOExt, self).__init__()

        self.structure_file = structure_file

        # Load the structure
        if file_type.lower() == 'pdb':
            structure = pdbp.get_structure(id='ssbio_pdb', file=structure_file)
        if file_type.lower() == 'mmcif' or file_type.lower() == 'cif':
            structure = cifp.get_structure(structure_id='ssbio_cif', filename=structure_file)

        # If there are multiple models (NMR), use the first model as the representative structure
        if len(structure) > 1:
            self.first_model = structure[0]
            structure = structure[0]
            self.set_structure(structure)
            log.debug('{}: using first model'.format(structure_file))
        elif len(structure) == 0:
            log.error('{}: no models in structure!'.format(structure_file))
        else:
            self.set_structure(structure)
            self.first_model = structure[0]

    def write_pdb(self, custom_name='', out_suffix='_new', out_dir=None, custom_selection=None):
        """Write a new PDB file for the Structure's FIRST MODEL.

        Set custom_selection to a PDB.Select class for custom SMCRA selections.

        Args:
            out_suffix: string to append to new PDB file - default is "_new"
            outdir: optional directory to output the file
            custom_selection: optional custom selection class

        Returns:
            out_file: filepath of new PDB file

        """
        if not custom_selection:
            custom_selection = Select()

        # Prepare the output file path
        outfile = ssbio.utils.outfile_name_maker(inname=self.structure_file,
                                                 outfile=custom_name,
                                                 append_to_name=out_suffix,
                                                 outdir=out_dir,
                                                 outext='.pdb')
        try:
            self.save(outfile, custom_selection)
        except TypeError:
            # If trying to save something that can't be saved as a PDB (example: 5iqr.cif), log an error and return None
            # The error thrown by PDBIO.py is "TypeError: %c requires int or char"
            log.error('{}: unable to save structure in PDB file format'.format(self.structure_file))
            return None

        return outfile


if __name__ == '__main__':
    pass
