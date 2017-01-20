#!/usr/bin/env python
# encoding: utf-8
# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from copy import deepcopy

# For Encad 3P CGing
from Bio.Struct.Geometry import center_of_mass
from Bio.PDB.Atom import Atom

def CA_TRACE(residue):
    """
    Reduces protein residues to the alpha carbon:
    CA trace only.
    """

    cg_residue = deepcopy(residue)

    for atom in cg_residue.child_dict.keys():
        if cg_residue[atom].name != 'CA':
            cg_residue.detach_child(cg_residue[atom].name)
    
    return cg_residue  

def ENCAD_3P(residue):
    """
    Reduces complexity of protein residue to a 3 points coarse grained model:
    CA, O, Bead in specific atom location.
    
    Based on Michael Levitt's ENCAD coarse graining model (2010).
    """
    
    conversion_table = {
                        "LEU": ["CG"],
                        "ALA": ["CB"],  
                        "GLY": [],  
                        "VAL": ["CB"],
                        "GLU": ["CD"],
                        "ASP": ["CG"],
                        "THR": ["CB"],
                        "LYS": ["CD"],
                        "ILE": ["CG1"],
                        "ARG": ["NE"],
                        "ASN": ["CG"],
                        "PRO": ["CG"],
                        "PHE": ["CG"],
                        "GLN": ["CD"],
                        "SER": ["CB", "OG"],
                        "HIS": ["CG"],
                        "MET": ["CG"],
                        "TYR": ["CD1", "CD2"],
                        "TRP": ["CD2"],
                        "CYS": ["CB", "SG"],
                        }
    
    cg_residue = deepcopy(residue)
    sc_atoms = []
    
    residue_name = cg_residue.resname
    
    if residue_name not in conversion_table:
        raise ValueError("Residue %s not recognized" %residue_name)
        
    for atom in cg_residue.child_dict.keys():
        if cg_residue[atom].name in conversion_table[residue_name]:
            sc_atoms.append(cg_residue[atom])
            cg_residue.detach_child(cg_residue[atom].name)            
        elif not cg_residue[atom].name in ["CA", "O"]:
            cg_residue.detach_child(cg_residue[atom].name)    

    # Side chain bead
    if len(sc_atoms):
        dummy = sorted(sc_atoms)[-1]
        if len(sc_atoms) > 1: # Merge - C.o.M        
            bead_coord = center_of_mass(sc_atoms, geometric=True)
        else:
            bead_coord = dummy.coord
            
        sc_bead = Atom( 'CMA', 
                        bead_coord, 
                        dummy.bfactor,
                        dummy.occupancy,
                        dummy.altloc,
                        ' CMA',
                        dummy.serial_number,
                        dummy.element,
                        dummy.hetatm)

        cg_residue.add(sc_bead)
        
    return cg_residue

def MARTINI(residue):
    """
    Reduces complexity of protein residue to the MARTINI coarse grained model:
    CA, O, Bead(s) in specific atom location.
    
    Reference:
    Monticelli et al. The MARTINI coarse-grained force field: extension to proteins. 
    J. Chem. Theory Comput. (2008) vol. 4 (5) pp. 819-834
    """
    
    conversion_table = {
                        "LEU": ["CG"],
                        "ALA": [], # No Side Chain  
                        "GLY": [],  
                        "VAL": ["CB"],
                        "GLU": ["CB"],
                        "ASP": ["CG"],
                        "THR": ["CB"],
                        "LYS": ["CG", "CE"],
                        "ILE": ["CD1"],
                        "ARG": ["CG", "NE"],
                        "ASN": ["CG"],
                        "PRO": ["CG"],
                        "PHE": ["CG", "CE1", "CE2"],
                        "GLN": ["CB"],
                        "SER": ["CB"],
                        "HIS": ["CB", "ND1", "NE2"],
                        "MET": ["CG"],
                        "TYR": ["CG", "CE1", "CE2"],
                        "TRP": ["CB", "CD1", "CD2", "CE2"],
                        "CYS": ["SG"],
                        }
                        
    cg_residue = deepcopy(residue)
    
    residue_name = cg_residue.resname
    
    if residue_name not in conversion_table:
        raise ValueError("Residue %s not recognized" %residue_name)
    
    for atom in cg_residue.child_dict.keys():
        if cg_residue[atom].name in conversion_table[residue_name]:
            # S1, S2, S3, or S4
            name = "S%s" % (int(conversion_table[residue_name].index(cg_residue[atom].name)) + 1)
            cg_residue[atom].name = name
            cg_residue[atom].fullname = " %s " %name
            cg_residue[atom].id = name
            cg_residue[atom].occupancy = int(cg_residue[atom].occupancy)
            cg_residue[atom].bfactor = int(cg_residue[atom].bfactor)
            cg_residue[atom].coord = [ int(xyz) for xyz in cg_residue[atom].coord ]
            
        elif cg_residue[atom].name == "CA":
            cg_residue[atom].name = "BB"
            cg_residue[atom].fullname = " BB "
            cg_residue[atom].id = "BB"
            cg_residue[atom].occupancy = int(cg_residue[atom].occupancy)
            cg_residue[atom].bfactor = int(cg_residue[atom].bfactor)
            cg_residue[atom].coord = [ int(xyz) for xyz in cg_residue[atom].coord ]                        
        else:
            cg_residue.detach_child(cg_residue[atom].name)
            
            
    return cg_residue