"""
Oligomer
========
"""

from ssbio.protein.structure.structprop import StructProp
import os.path as op
from string import ascii_uppercase
from copy import copy
import ssbio.utils


class Oligomer(StructProp):

    """Methods to deal with a generic oligomeric structure, extends StructProp.

    The main utilities of this class are to:

    * Parse a structure file that represents an oligomeric structure
    * Return stoichiometric coefficients of the chains in the bioassembly
    * More to come...

    Args:
        ident (str): Oligomer ID
        description (str): Optional description for this bioassembly
        is_experimental (bool): Flag to indicate if structure is an experimental or computational model
        structure_path (str): Path to structure file
        file_type (str): Type of structure file - ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``,
            ``xml.gz``, ``mmtf``, ``mmtf.gz``

    Todo:
        - Check out http://biopython.org/wiki/Interface_Analysis (Bio.PDB.NACCESS, Interface, InterfaceAnalysis modules)

    """

    def __init__(self, ident, description=None, is_experimental=False, structure_path=None, file_type=None):
        StructProp.__init__(self, ident=ident, description=description, chains=None, mapped_chains=None,
                            is_experimental=is_experimental, structure_path=structure_path, file_type=file_type)
        self.subunits = None
        """dict: Subunit composition defined as ``{chain_id: number_of_times_used_in_bioassembly}``"""

    def merge_models(self, store_in_memory=False, outdir=None, outname=None, force_rerun=False):
        """Merge all existing models into a Structure's first_model attribute.

        This directly modifies the Biopython Structure object. Chains IDs will start from A and increment for each new
        chain (which is a Model that is converted).

        Args:
            store_in_memory (bool): If the modified Biopython Structure object should be stored in the attribute
                ``structure``
            outdir (str): If ``store_in_memory`` is False, the structure file has to be written somewhere so an output
                directory must be specified here
            outname (str): If ``store_in_memory`` is False, the structure file has to be written somewhere so an output
                filename must be specified here (i.e. 4BXI_bio1)
            force_rerun (bool): If merged file should be overwritten if it already exists

        """

        if store_in_memory:
            if self.structure:
                parsed = copy(self.structure)
            else:
                parsed = self.parse_structure()
            self.structure = merge_all_models_into_first_model(parsed)
        else:
            new_structure_path = write_merged_bioassembly(inpath=self.structure_path,
                                                          outdir=outdir, outname=outname,
                                                          force_rerun=force_rerun)
            self.load_structure_path(new_structure_path, file_type='pdb')


def merge_all_models_into_first_model(biop_structure):
    """Merge all existing models into a Structure's first_model attribute.

    This directly modifies the Biopython Structure object. Chains IDs will start from A and increment for each new
    chain (model that is converted).

    Args:
        biop_structure (Structure): Structure with multiple models that should be merged

    """
    from string import ascii_uppercase
    idx = 1
    first_model = biop_structure[0]

    for m in biop_structure.get_models():
        # Don't duplicate the original model
        if first_model.id == m.id:
            continue
        for c in m.get_chains():
            c.id = ascii_uppercase[idx]
            first_model.add(c)
        idx += 1


def write_merged_bioassembly(inpath, outdir, outname, force_rerun=False):
    """Utility to take as input a bioassembly file and merge all its models into multiple chains in a single model.

    Args:
        infile (str): Path to input PDB file with multiple models that represent an oligomeric form of a structure.
        outdir (str): Path to output directory
        outname (str): New filename of structure file
        force_rerun (bool): If a new PDB should be written if the file exists

    Returns:
        str: Path to newly written PDB file.

    """
    outpath = outfile=op.join(outdir, outname + '.pdb')

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=op.join(outdir, outname + '.pdb')):
        s = StructProp('Model merging', structure_path=inpath, file_type='pdb')
        ss = s.parse_structure()
        merge_all_models_into_first_model(ss.structure)
        outpath = ss.write_pdb(custom_name=outname, out_dir=outdir, force_rerun=force_rerun)
    else:
        return outpath