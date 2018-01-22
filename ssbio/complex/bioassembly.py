"""
Bioassembly
===========
"""

from ssbio.core.object import Object


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