# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for dealing with PDB structures.
"""

import os

def read( handle, id=None ):
    """
    Reads a structure via PDBParser.
    Simplifies life..
    """

    from Bio.PDB import PDBParser
    
    if not id:
        id = os.path.basename(handle).split('.')[0] # Get from filename
    
    p = PDBParser()
    s = p.get_structure(id, handle)

    return s

def write( structure, name=None ):
    """
    Writes a Structure in PDB format through PDBIO.
    Simplifies life..
    """
    
    from Bio.PDB import PDBIO
    
    io = PDBIO()
    
    io.set_structure(structure)
    
    if not name:
        s_name = structure.id
    else:
        s_name = name

    name = "%s.pdb" %s_name
    seed = 0
        
    while 1:
        if os.path.exists(name):
            name = "%s_%s.pdb" %(s_name, seed)
            seed +=1
        else:
            break
                
    io.save(name)
        
    return name