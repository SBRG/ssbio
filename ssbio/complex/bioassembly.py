"""
Bioassembly
===========
"""

from ssbio.core.object import Object
import os.path as op


class Bioassembly(Object):

    """Methods to deal with a PDB biological assembly file.

    The main utilities of this class are to:

    * Specify a PDB file type and download bioassemblies given a PDB ID, or download the original PDB and work with
      information contained in the BIOMT remark field.
    * Combine bioassembly models into one single model so it is viewed as a unit when parsed with Biopython
    * Return stoichiometric coefficients of the chains in the bioassembly

    Args:
        ident (str): PDB ID
        biomol (str): Biological assembly number as defined in the PDB
        subunits (dict): Subunit composition defined as ``{chain_id: number_of_times_used_in_bioassembly}``

        description (str): Optional description for this bioassembly
        root_dir (str): Path to where bioassemblies will be downloaded. Default is current working directory.
        pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` -
            choose a file type for files downloaded from the PDB

    """

    def __init__(self):
        pass


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


def write_merged_bioassembly(inpath, outpath, force_rerun=False):
    """Utility to take as input a bioassembly file and merge all its models into multiple chains in a single model.

    Args:
        infile (str): Path to input PDB file with multiple models that represent an oligomeric form of a structure.
        outfile (str): Path to new output file.

    Returns:
        str: Path to newly written PDB file.

    """
    import ssbio.protein.structure.structprop as sp
    s = sp.StructProp('Model merging', structure_path=inpath, file_type='pdb')
    ss = s.parse_structure()
    merge_all_models_into_first_model(ss.structure)

    tmp = op.split(outpath)
    outdir = tmp[0]
    name = op.splitext(tmp[1])[0]
    ss.write_pdb(custom_name=name, out_dir=outdir, force_rerun=force_rerun)