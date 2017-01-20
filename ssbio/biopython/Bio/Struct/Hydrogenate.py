# Copyright (C) 2010,  Joao Rodrigues (anaryin@gmail.com)
# This module is heavily based on PyMol's code.
# !! Similarities are not a coincidence. !!
# PyMol: chempy/protein.py chempy/place.py
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# General

from random import randint

# Vector Operations

from cpv import normalize,  scale,  cross_product,  add,  random_sphere,  sub

# BioPython modules
from Bio.PDB.Atom import Atom

# Bond information (Taken from PyMol)
import protein_residues
import protein_amber
import bond_amber

TET_TAN = 1.41
TRI_TAN = 1.732

class Hydrogenate_Protein:

    def __init__(self,  forcefield=protein_amber,  template=protein_residues):
        
        # FF and Tmplt pre-load
        
        self.tmpl = {   'cter': template.c_terminal, 
                        'nter': template.n_terminal, 
                        'nrml': template.normal
                    }
        self.ffld = {   'cter': forcefield.c_terminal, 
                        'nter': forcefield.n_terminal, 
                        'nrml': forcefield.normal
                    }

        # Constants (Taken from PyMol)
        self.N_TERMINAL_ATOMS = set(['HT', 'HT1', 'HT2', 'HT3', 'H1', 'H2', 'H3', 
                                  '1H', '2H', '3H', '1HT', '2HT', '3HT'])

        self.C_TERMINAL_ATOMS = set(['OXT', 'O2', 'OT1', 'OT2'])
        
        # Protonation Methods
        
        self.protonation_methods = {
                                    1: self._add_1, 
                                    2: self._add_2, 
                                    3: self._add_3, 
                                    4: self._add_4
                                    }
    
    def add_hydrogens(self,  u_input,  bondfield=bond_amber):
        
        self.nh_structure = u_input
        self._exclude_ss_bonded_cysteines()
        
        last_count = -1        
        
        while 1:
            n_missing = self._build_bonding_network()
            # If goes around twice and no H is added,  breaks.
            if n_missing == last_count:
                break
            n_hadded = self.place_hydrogens(bondfield)
            last_count = n_missing - n_hadded # If all added is 0
        
        return last_count
    
    
    def _build_bonding_network(self):
        """
        Evaluates atoms per residue for missing and known bonded partners.
        Based on bond_amber.
        A better alternative would be to iterate over the entire list of residues and 
        use NeighborSearch to probe neighbors for atom X in residue i,  i-1 and i+1
        """
        
        self.bonds = {} # [Residue][Atom]: [ [missing],  [bonded] ]
        self.selection = {} # [Residue]: 'nrml'/'nter'/'cter'
        
        missing = 0
        
        for residue in self.nh_structure.get_residues():
            
            bond_dict = self.bonds[residue] = {}
            
            atom_dict = residue.child_dict
            atom_names = set(atom_dict.keys())
            
            # Pre-Populate Dictionary
            
            for name in atom_names:
                bond_dict[name] = [ [],  [] ]
            
            # Define Template
            if atom_names.intersection(self.C_TERMINAL_ATOMS):
                selection = 'cter'
            elif atom_names.intersection(self.N_TERMINAL_ATOMS):
                selection = 'nter'
            else:
                selection = 'nrml'
            
            tmpl = self.tmpl[selection]
            self.selection[residue] = selection # For place_hs
            
            # Iterate Template Bonds and record info
            
            if not tmpl.has_key(residue.resname):
                raise ValueError("Unknown Residue Type: %s" %residue.resname)
            
            template_bonds = tmpl[residue.resname]['bonds']
            
            for bond in template_bonds.keys():
                
                a1,  a2 = bond
                
                if a1 in atom_names and not a2 in atom_names:
                    bond_dict[a1][0].append(a2)
                    missing += 1
                elif a1 not in atom_names and a2 in atom_names:
                    bond_dict[a2][0].append(a1)
                    missing += 1
                else: # 
                    bond_dict[a1][1].append(atom_dict[a2])
                    bond_dict[a2][1].append(atom_dict[a1])   
        
        return missing
        
    def _exclude_ss_bonded_cysteines(self):
        """
        Pre-compute ss bonds to discard cystines for H-adding.
        """
        
        ss_bonds =  self.nh_structure.search_ss_bonds()
        for cys_pair in ss_bonds:
            cys1,  cys2 = cys_pair
            
            cys1.resname = 'CYX'
            cys2.resname = 'CYX'
    
    def _find_secondary_anchors(self,  residue,  heavy_atom,  anchor):
        """
        Searches through the bond network for atoms bound to the anchor.
        Returns a secondary and tertiary anchors.
        Example,  for CA,  returns C and O.
        """
        
        for secondary in self.bonds[residue][anchor.name][1]:
            for tertiary in self.bonds[residue][secondary.name][1]:
                if (tertiary.name != heavy_atom.name 
                    and tertiary.name != anchor.name):
                    
                    return (secondary,  tertiary)
        
        return None
        
    def place_hydrogens(self,  bondfield):

        n_added = 0
        
        # define bondfield
        self.bondfield = bondfield
          
        for residue in sorted(self.bonds.keys(),  key=lambda x: x.get_id()[1]):
            incomplete_atoms = self.bonds[residue]

            for atom in incomplete_atoms:
                missing_atoms = incomplete_atoms[atom][0]
                if len(missing_atoms):                
                
                    if len(missing_atoms) == 1:
                        
                        h_coord = self._add_1(  missing_atoms[0],  
                                                residue.child_dict[atom],  
                                                incomplete_atoms[atom][1])
                                                
                        new_at = Atom(  name=missing_atoms[0],  
                                        coord=h_coord,  bfactor=1.0,  
                                        occupancy=0.0,  altloc=' ',  
                                        fullname=missing_atoms[0],  
                                        serial_number=randint(5000, 9999), 
                                        element='H' )
                                        
                        residue.add(new_at)
                        n_added += 1
                        
                    elif len(missing_atoms) == 2:
                        coordinates = self._add_2(  missing_atoms,  
                                                    residue.child_dict[atom],  
                                                    incomplete_atoms[atom][1])
                        for name in coordinates:
                            new_at = Atom(  name=name,  
                                            coord=coordinates[name],  bfactor=1.0,  
                                            occupancy=0.0,  altloc=' ',  
                                            fullname=name,  
                                            serial_number=randint(5000, 9999),  
                                            element='H')
                                            
                            residue.add(new_at)
                            n_added += 1
                    
                    elif len(missing_atoms) == 3:
                        coordinates = self._add_3(  missing_atoms,  
                                                    residue.child_dict[atom],  
                                                    incomplete_atoms[atom][1])
                                                    
                        for name in coordinates:
                            new_at = Atom(  name=name,  
                                            coord=coordinates[name],  
                                            bfactor=1.0,  occupancy=0.0,  
                                            altloc=' ',  
                                            fullname=name,  
                                            serial_number=randint(5000, 9999),  
                                            element='H')
                                            
                            residue.add(new_at)
                            n_added += 1

                    elif len(missing_atoms) == 4:
                        coordinates = self._add_4(  missing_atoms,  
                                                    residue.child_dict[atom],  
                                                    incomplete_atoms[atom][1])
                                                    
                        new_at = Atom(  name=missing_atoms[0],  
                                        coord=coordinates,  
                                        bfactor=1.0,  occupancy=0.0,  
                                        altloc=' ',  
                                        fullname=missing_atoms[0],  
                                        serial_number=randint(5000, 9999),  
                                        element='H')
                                        
                        residue.add(new_at)
                        n_added += 1
        
        return n_added
                            
    def _add_1(self,  hydrogen_name,  heavy_atom,  bonds):
        """
        Adds the missing proton to single protonated heavy_atoms.
        """
        
        residue = heavy_atom.parent
        ffld = self.ffld[self.selection[residue]]
        bnd_len = self.bondfield.length
        anchor = bonds[0]
                
        # If not linear
        if self.bondfield.nonlinear.has_key(ffld[(  residue.resname,  
                                                    heavy_atom.name)]['type']):
            # returns tuple of two atoms
            bonded = self._find_secondary_anchors(residue,  heavy_atom,  anchor)

            if bonded:
                # Phenolic hydrogens,  etc.
                if self.bondfield.planer.has_key(ffld[( residue.resname,  
                                                        anchor.name)]['type']):
                    secondary_anchor = bonded[0]
                   
                    p0 = heavy_atom.coord - anchor.coord
                    d2 = secondary_anchor.coord - anchor.coord
                    p1 = normalize(cross_product(d2,  p0))
                    p2 = normalize(cross_product(p0,  p1))                     
                    vector = scale(p2,  TRI_TAN)
                    vector = normalize(add(p0,  vector))
                    
                    hydrogen_coord = add(heavy_atom.coord, 
                                        scale(vector,  
                                              bnd_len[(ffld[(residue.resname, 
                                                             heavy_atom.name)]['type'], 
                                              ffld[(residue.resname,  
                                                    hydrogen_name)]['type'])]))
                
                else: # Ser,  Cys,  Thr hydroxyl hydrogens
                    secondary_anchor = bonded[0]
                    vector = anchor.coord - secondary_anchor.coord
                    hydrogen_coord = add(heavy_atom.coord, 
                                        scale(vector,  
                                              bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'],  
                                              ffld[(residue.resname,  hydrogen_name)]['type'])]    ))
                    
            elif len(bonds):
                d2 = [1.0, 0, 0]
                p0 = heavy_atom.coord - anchor.coord
                p1 = normalize(cross_product(d2,  p0))
                vector = scale(p1,  TET_TAN)
                vector = normalize(add(p0,  vector))
                
                hydrogen_coord = add(heavy_atom.coord, 
                                     scale(  vector,  
                                             bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'],  
                                             ffld[(residue.resname,  hydrogen_name)]['type'])]    ))
            else:
                hydrogen_coord = random_sphere(heavy_atom.coord, 
                                               bnd_len[ (ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                         ffld[(residue.resname,  hydrogen_name)]['type']) ])

        elif len(bonds): # linear sum...amide,  tbu,  etc
            vector = [0.0, 0.0, 0.0]
            if heavy_atom.name == 'N': # Fix to get previous atom O from peptide bond. Ugly.
                prev_res = list(residue.get_id())
                prev_res[1] -= 1
                prev_res = tuple(prev_res)
                if residue.parent.child_dict.has_key(prev_res):
                    prev_res = residue.parent.child_dict[prev_res]
                    bonds.append(prev_res.child_dict['O'])
            for b in bonds:
                d = heavy_atom.coord - b.coord
                vector = add(vector,  d)
            vector = normalize(vector)
            hydrogen_coord = add(heavy_atom.coord, 
                                 scale(vector, 
                                 bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                          ffld[(residue.resname,  hydrogen_name)]['type']) ]))
        else:
            hydrogen_coord = random_sphere(heavy_atom.coord, 
                                           bnd_len[ (ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                     ffld[(residue.resname,  hydrogen_name)]['type']) ])
        
        return hydrogen_coord
        
    def _add_2(self,  hydrogens,  heavy_atom,  bonds):
        
        hydrogen_coord = {} # Returns two coordinate sets [Atom]: [x,  y,  z]
        
        residue = heavy_atom.parent
        ffld = self.ffld[self.selection[residue]]
        bnd_len = self.bondfield.length
        anchor = bonds[0]
        first_hydrogen,  second_hydrogen = hydrogens
        
        if self.bondfield.planer.has_key(ffld[(residue.resname,  heavy_atom.name)]['type']): # guanido,  etc
            bonded = self._find_secondary_anchors(residue,  heavy_atom,  bonds[0])
            if bonded: # 1-4 present
                secondary_anchor = bonded[0]
                
                # Add first H
                
                d1 = heavy_atom.coord - anchor.coord
                d2 = secondary_anchor.coord - anchor.coord
                p1 = normalize(cross_product(d2,  d1))
                p2 = normalize(cross_product(d1,  p1))
                vector = scale(p2,  TRI_TAN)
                vector = normalize(add(d1,  vector))
                hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                            scale(vector, 
                                                  bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                           ffld[(residue.resname,  first_hydrogen)]['type'])]))                                                         
                # Add second H

                vector = scale(p2,  -TRI_TAN)
                vector = normalize(add(d1,  vector))
                hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                      scale(vector, 
                                                            bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                     ffld[(residue.resname,  second_hydrogen)]['type'])]))
            
            elif len(bonds): # no 1-4 found

                # H1

                d2 = [1.0, 0, 0]
                d1 = heavy_atom.coord - anchor.coord
                p1 = normalize(cross_product(d2,  d1))
                p2 = normalize(cross_product(d1,  p1))                
                vector = scale(p2,  TRI_TAN)
                vector = normalize(add(d1,  vector))
                
                hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                                      scale(vector, 
                                                            bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                     ffld[(residue.resname,  first_hydrogen)]['type'])]))
                
                # H2

                vector = scale(p2,  -TRI_TAN)
                vector = normalize(add(d1,  vector))
                hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                      scale(vector, 
                                                            bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                     ffld[(residue.resname,  second_hydrogen)]['type'])]))
            else: # If all else fails,  add one random
                hydrogen_coord[first_hydrogen] = random_sphere(heavy_atom.coord, 
                                                                bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                         ffld[(residue.resname,  first_hydrogen)]['type'])])
        
        elif len(bonds)>=2: # simple tetrahedral
            
            secondary_anchor = bonds[1]
            
            # H1
            vector = [0.0, 0.0, 0.0]
            d1 = heavy_atom.coord - anchor.coord
            d2 = heavy_atom.coord - secondary_anchor.coord
            vector = add(d1,  d2)
            p0 = normalize(vector)
            p1 = normalize(cross_product(d2,  p0))
            vector = scale(p1,  TET_TAN)
            vector = normalize(add(p0,  vector))
            
            hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  first_hydrogen)]['type'])]))
            # H2
               
            vector = scale(p1,  -TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  second_hydrogen)]['type'])]))
        else:
            # H1
            if len(bonds): # sulfonamide? 
                
                d2 = [1.0, 0, 0]
                d1 = heavy_atom.coord - anchor.coord              
                p1 = normalize(cross_product(d2,  d1))
                vector = scale(p1,  TET_TAN)
                vector = normalize(add(p0,  vector))
                hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                                      scale(vector, 
                                                            bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                     ffld[(residue.resname,  first_hydrogen)]['type'])]))
            else: # blind
                hydrogen_coord[first_hydrogen] = random_sphere(heavy_atom.coord, 
                                                                bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                         ffld[(residue.resname,  first_hydrogen)]['type'])])
            # H2
            first_hydrogen_coord = hydrogen_coord[first_hydrogen]

            vector = [0.0, 0.0, 0.0]
            d1 = heavy_atom.coord - anchor.coord
            d2 = sub(heavy_atom.coord,  first_hydrogen_coord)
            vector = add(normalize(d1),  normalize(d2))
            p0 = normalize(vector)
            p1 = normalize(cross_product(d2,  p0))
            vector = scale(p1,  TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  second_hydrogen)]['type'])]))
        return hydrogen_coord
        
    def _add_3(self,  hydrogens,  heavy_atom,  bonds):
        
        
        hydrogen_coord = {} # Returns two coordinate sets [Atom]: [x,  y,  z]
        
        residue = heavy_atom.parent
        ffld = self.ffld[self.selection[residue]]
        bnd_len = self.bondfield.length
        
        anchor = bonds[0]
        first_hydrogen,  second_hydrogen,  third_hydrogen = hydrogens

        bonded = self._find_secondary_anchors(residue,  heavy_atom,  bonds[0])
        
        if bonded: # 1-4 present
            
            secondary_anchor = bonded[0]
            
            # H1
            d1 = heavy_atom.coord - anchor.coord
            d2 = secondary_anchor.coord - anchor.coord
            p1 = normalize(cross_product(d2,  d1))
            p2 = normalize(cross_product(d1,  p1))
            vector = scale(p2,  -TET_TAN)
            vector = normalize(add(d1,  vector))
            hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  first_hydrogen)]['type'])]))                                                         
            
            # H2
            first_hydrogen_coord = hydrogen_coord[first_hydrogen]
            
            d1 = heavy_atom.coord - anchor.coord
            d2 = sub(heavy_atom.coord,  first_hydrogen_coord)
            vector = add(d1,  normalize(d2))
            p0 = normalize(vector)
            p1 = normalize(cross_product(d2,  p0))
            vector = scale(p1,  TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  second_hydrogen)]['type'])]))
            
            # H3

            vector = scale(p1,  -TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[third_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  third_hydrogen)]['type'])]))
        elif len(bonds): # fall-back
            # H1
            d2 = [1.0, 0, 0]
            d1 = heavy_atom.coord - anchor.coord              
            p1 = normalize(cross_product(d2,  d1))
            vector = scale(p1,  TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[first_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  first_hydrogen)]['type'])]))
            
            # H2
            first_hydrogen_coord = hydrogen_coord[first_hydrogen]
            vector = [0.0, 0.0, 0.0]
            d1 = heavy_atom.coord - anchor.coord
            d2 = sub(heavy_atom.coord,  first_hydrogen_coord)
            vector = add(d1,  normalize(d2))
            p0 = normalize(vector)
            p1 = normalize(cross_product(d2,  p0))
            vector = scale(p1,  TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[second_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  second_hydrogen)]['type'])]))
            
            # H3

            vector = scale(p1,  -TET_TAN)
            vector = normalize(add(p0,  vector))
            hydrogen_coord[third_hydrogen] = add(heavy_atom.coord, 
                                                  scale(vector, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  third_hydrogen)]['type'])]))
        
        else: # worst case: add one and get rest next time around
            hydrogen_coord[first_hydrogen] = random_sphere(heavy_atom.coord, 
                                                            bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                     ffld[(residue.resname,  first_hydrogen)]['type'])])
        return hydrogen_coord
    
    def _add_4(self,  hydrogens,  heavy_atom,  bonds):
                
        residue = heavy_atom.parent
        ffld = self.ffld[self.selection[residue]]
        bnd_len = self.bondfield.length

        hydrogen_coord = random_sphere(heavy_atom.coord, 
                                                        bnd_len[(ffld[(residue.resname,  heavy_atom.name)]['type'], 
                                                                 ffld[(residue.resname,  hydrogens[0])]['type'])])
        return hydrogen_coord