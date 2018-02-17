"""
Bioassembly
===========
"""

from ssbio.complex.oligomer import Oligomer
import os.path as op
import requests
from lxml import etree
import os
import ssbio.utils
import xmltodict
import logging

log = logging.getLogger(__name__)


class Bioassembly(Oligomer):

    """Methods to deal with a PDB biological assembly file, extends Oligomer.

    The main utilities of this class are to:

    * Specify a PDB ID and biomolecule number to download bioassembly from the PDB, or download the original PDB and
      assemble the bioassembly using the transformation matrix
    * Combine bioassembly models into one single model so it is viewed as a unit when parsed with Biopython

    Args:
        ident (str): PDB ID - the full ID of this object will be "<PDB_ID>_biomol" (i.e. 4bxi_biomol)
        biomol (int, optional): Biological assembly number as defined in the PDB - if not specified, first bioassembly
            is assumed.
        description (str): Optional description for this bioassembly
        root_dir (str): Path to where bioassemblies will be downloaded. Default is current working directory.
        pdb_file_type (str): ``pdb``, ``pdb.gz``, ``mmcif``, ``cif``, ``cif.gz``, ``xml.gz``, ``mmtf``, ``mmtf.gz`` -
            choose a file type for files downloaded from the PDB

    Todo:
        - Check out http://biopython.org/wiki/Interface_Analysis (

    """

    def __init__(self, pdb_id, biomol_num, description=None, is_experimental=True, structure_path=None, file_type=None, ):
        full_id = '{}_bio{}'.format(pdb_id, biomol_num)
        super(Bioassembly, self).__init__(ident=full_id, description=description, is_experimental=is_experimental,
                                          structure_path=structure_path, file_type=file_type)

        self.original_pdb_id = pdb_id
        self.num_biomols = 0
        """int: Total number of bioassemblies available for this PDB"""
        self.biomol_to_chain_dict = {}
        """dict: Bioassembly number to a dictionary containing info about the utilized chains and their stoichiometry"""

    def load_bioassembly_info_from_file(self, biomol_num):
        """Load metadata about a bioassembly (such as chains and their transformations) from a structure file.
        """

        # current functionality is to take in a pre-assembled bioassembly file, parse it's MODELs and get info from that
        pass
        from Bio.PDB import PDBParser
        p = PDBParser(PERMISSIVE=True, QUIET=True)

        structure = StructureIO(self.structure_path, self.file_type)

        structure = p.get_structure('tmp', infile)
        bioass_to_chain_stoich = defaultdict(int)
        for model in structure:
            for chain in model:
                bioass_to_chain_stoich[chain.id] += 1
        return dict(bioass_to_chain_stoich)

