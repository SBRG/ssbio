"""
FoldX
=====
"""

from ssbio.core.object import Object


class FoldX(Object):

    """Class to run various commands from the FoldX suite of tools on a protein structure file.

    Args:
        structure_id (str): Name of your structure or this mini FoldX project
        pdb_path (str): Path to PDB file (PDB file format only!)
        rotabase_path (str): Path to the rotabase.txt file
        foldx_path (str, optional): Path to the FoldX executable, if empty, tries to execute ``foldx`` from the shell
        root_dir (str, optional): Path to

    """

    def __init__(self, structure_id, pdb_path, rotabase_path, foldx_path=None, root_dir=None):
        super(FoldX, self).__init__(id=structure_id, description='FoldX project')