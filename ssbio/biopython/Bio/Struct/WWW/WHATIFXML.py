# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os, sys
import numpy
import xml.dom.minidom
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

from Bio.PDB.Atom import Atom # for single atom additions

__doc__="Parser for WHATIF server XML Structure format."


# According to WHAT IF, their XML format is similar to PDB XML
# so this might serve as a base for future PDB XML parsers.
# Yet, since WHAT IF already does extensive validation,
# this implementation performs very basic checks.

# ie. not suitable for general wide use of PDB XML parsing..

# Allowed AA List for HETATM distinguishing.
# N_TERMINAL_ATOMS and C_TERMINAL_ATOMS for correct atom attribution

AA_LIST = set([ 'ILE', 'GLN', 'GLY', 'GLP', 'GLU', 
                'CYS', 'ASP', 'LYS', 'PRO', 'CYX', 
                'HID', 'HIE', 'ASH', 'ASN', 'CYM', 
                'HIP', 'VAL', 'THR', 'TRP', 'SER', 
                'PHE', 'ALA', 'MET', 'LEU', 'ARG', 
                'HIS', 'TYR'])

N_TERMINAL_ATOMS = set(['HT','HT1','HT2','HT3','H1','H2','H3',
                          '1H','2H','3H','1HT','2HT','3HT'])

C_TERMINAL_ATOMS = set(['OXT','O2','OT1','OT2'])


class ParserException(Exception):
    pass

class XMLParser:
    
    def __init__(self, structure_builder=None):
        
        if structure_builder!=None:
            self.structure_builder=structure_builder
        else:
            self.structure_builder=StructureBuilder()
        
    def read(self, f, id):

        try:
            self.handle = xml.dom.minidom.parse(f)
        except xml.parsers.expat.ExpatError, e:
            self._handle_parser_exception(e)
        
        self.structure_builder.init_structure(id)
        
        self._parse()
        
        return self.structure_builder.get_structure()
    
    
    def _handle_parser_exception(self, message):

        message = "%s. Error when parsing input file (not XML format or corrupted)" %message
        
        raise ParserException(message)
    
    def _handle_builder_exception(self, message, residue):
        """
        Makes a PDB Construction Error a bit more verbose and informative
        """
        
        message = "%s. Error when parsing residue %s:%s" %(message, residue['number'], residue['name'])
        
        raise PDBConstructionException(message)
    
    def _parse(self):
        """
        Parse atomic data of the XML file.
        """
        
        atom_counter = 0
        
        structure_build = self.structure_builder
        
        residues = self._extract_residues()
        
        cur_model = None
        cur_chain = None
        structure_build.init_seg(' ') # There is never a SEGID present
        
        for r in residues:

            # New model?
            if cur_model != r['model']:
                cur_model = r['model']
                try:
                    structure_build.init_model(cur_model)
                except PDBConstructionException, message:
                    self._handle_builder_exception(message, r) 
            
            # New chain?
            if cur_chain != r['chain']:
                cur_chain = r['chain']
                try:
                    structure_build.init_chain(cur_chain)
                except PDBConstructionException, message:
                    self._handle_builder_exception(message, r)

            # Create residue
            if r['name'] in AA_LIST: # Get residue type crudely since there is no HETATM / ATOM
                hetero_flag = ' '
            elif r['name'] == 'WAT' or r['name'] == 'HOH':
                hetero_flag = 'W'
            else:
                hetero_flag = 'H'
            
            # Some terminal atoms are added at residue 0. This residue has a small number of atoms.
            # Protonated non-terminal glycine has 7 atoms. Any of these residues is smaller.
            # HETATMs have only a couple of atoms (3 for water for example) and they are ok.
            
            if (len(r['atoms']) >= 7) or (hetero_flag != " "):
                try:
                    structure_build.init_residue(r['name'], hetero_flag, r['number'], r['icode'])
                except PDBConstructionException, message:
                    self._handle_builder_exception(message, r)
            
                # Create Atoms
                for atom in r['atoms']:
                    a = self._parse_atom(atom)

                    if not sum(a['coord']): # e.g. HG of metal bound CYS coords are 0,0,0.
                        continue
                        
                    try:
                        atom_counter += 1
                        # fullname = name; altloc is empty;
                        structure_build.init_atom(a['name'], a['coord'], a['bfactor'], a['occupancy'], ' ',
                                                    a['name'], atom_counter, a['element'],  hetero_flag)
                    except PDBConstructionException, message:
                        self._handle_builder_exception(message, r)      
            
            elif len(r['atoms']) < 7: # Terminal Residues
                                
                for atom in r['atoms']:
                    
                    a = self._parse_atom(atom)
                    
                    if not sum(a['coord']): # e.g. HG of metal bound CYS coords are 0,0,0.
                        continue

                    atom_counter += 1                    
                    ter_atom = Atom(a['name'], a['coord'], a['bfactor'], a['occupancy'], ' ',
                                    a['name'], atom_counter, a['element'], hetero_flag)
                    
                    if a['name'] in N_TERMINAL_ATOMS:
                        
                        inc_struct = self.structure_builder.get_structure()
                        for model in inc_struct:
                            for chain in model:
                                if chain.id == r['chain']:
                                    for residue in chain: # Find First residue matching name
                                        if residue.resname == r['name']:
                                            residue.add(ter_atom)
                                            break
                    
                    elif a['name'] in C_TERMINAL_ATOMS:

                        inc_struct = self.structure_builder.get_structure()
                        c_ter = None
                        for model in inc_struct:
                            for chain in model:
                                if chain.id == r['chain']:
                                    for residue in chain: # Find Last residue matching name
                                        if residue.resname == r['name']:
                                            c_ter = residue
                                    if c_ter:
                                        c_ter.add(ter_atom)
                    
                    # Else, discard atom...
                    

    def _extract_residues(self):
        """
        WHAT IF puts terminal atoms in new residues at the end for some reason..
        """        
        
        r_list = self.handle.getElementsByTagName("response")
        
        r_data = {}
        
        for r in r_list:
            data = self._parse_residue(r)
            res_id = (data['model'], data['chain'], data['number']) # (A, 1, 1), (A, 1, 2), ...
            if not r_data.has_key(res_id):
                r_data[res_id] = data
            else: # Some atoms get repeated at the end with TER Hydrogens/oxygens
                r_data[res_id]['atoms'] += data['atoms'] # Append Atoms
        
        for key in sorted(r_data.keys()):
            yield r_data[key]
        
    def _parse_residue(self, residue):
        """
        Extracts Residue Name, Number, Chain, Model, Atoms.
        I/O: xml object <response> / dictionary
        """
        
        # Filter Element Nodes
        childs = [ child for child in residue.childNodes if child.nodeType == child.ELEMENT_NODE ]
        
        # Parse info out
        resi = int(childs[0].firstChild.data.strip())
        resn = childs[1].firstChild.data.strip()
        icode = childs[3].firstChild.data
        chain = childs[4].firstChild.data.strip()
        model = int(childs[5].firstChild.data.strip())
        atoms = childs[6:]
        
        # Output
        
        return {'name':  resn,  'number': resi,
                'icode': icode, 'chain':  chain, 
                'model': model, 'atoms': atoms}
    
    def _parse_atom(self, atom):
        
        # Filter Element Nodes
        childs = [ child.firstChild.data.strip() for child in atom.childNodes if child.nodeType == child.ELEMENT_NODE ]
        
        number = int(childs[0])
        coord = numpy.array([float(childs[1]), float(childs[2]), float(childs[3])], 'f')
        bfactor = float(childs[4])
        #chain = childs[5]
        name = childs[6]
        occupancy = int(childs[7])
        element = childs[8]
        
        return {'name':      name,      'number':   number,
                'coord':     coord,     'bfactor':  bfactor,
                'occupancy': occupancy, 'element':  element}
