# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__doc__ = """
Interface for servers provided at http://swift.cmbi.kun.nl/.

Available Functions:
add_hydrogens        Add Hydrogens to a protein using the PDBasXMLwithSymwithPolarH service.

WARNING:
Some services might be experimental. Examine the results thouroughly!
If you find any problem, please contact us.
"""

import urllib
import os, tempfile, time
from Bio.PDB.Entity import Entity
from Bio.PDB import PDBIO
import xml.dom.minidom

# WHAT IF returns structures in XML format
from WHATIFXML import XMLParser

# Warning message to be outputted when experimental service is requested.

_WARNING =  "!! WARNING !!\nThis service is still experimental. It may yield WRONG results. \
            \nPlease contact the Biopython developers if you find your results suspicious (or just plain wrong)."

class ConnectionError(Exception):
    pass

class WHATIF:
    
    def __init__(self):
        
        retry = 0
        
        while 1:

            if not self.is_alive(): # Problem...
                if retry == 3:
                    raise ConnectionError(  "WHATIF Server appears to be down.. \
                                            \nContact the Biopython developers if the problem persists.")
                else:
                    retry += 1
                    time.sleep(5*retry)
            else:
                break # Up and running!
                
        # Initiate XMLParser
        
        self.parser = XMLParser()
    
    
    # Utility functions
    def _smcra_to_str(self, smcra, temp_dir='/tmp/'):
        """
        WHATIF's input are PDB format files.
        Converts a SMCRA object to a PDB formatted string.
        """
        
        temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
        
        io = PDBIO()
        io.set_structure(smcra)
        io.save(temp_path)
        
        f = open(temp_path, 'r')
        string = f.read()
        f.close()
        
        os.remove(temp_path)
        
        return string

        
    def _str_to_smcra(self, string, temp_dir='/tmp/'):
        
        temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
        f = open(temp_path, 'r')
        f.write(string)
        f.close()
        
        p = PDBParser()
        s = p.get_structure('structure', temp_path)
 
        os.remove(temp_path)
        
        return s
    

    def is_alive(self):
        """
        Test Function to check WHAT IF servers are up and running.
        """
        
        u = urllib.urlopen("http://wiws.cmbi.ru.nl/rest/TestEmpty/id/1crn/")
        x = xml.dom.minidom.parse(u)
        self.alive = len(x.getElementsByTagName("TestEmptyResponse"))

        return self.alive
    
    
    def UploadPDB(self, structure):
        """
        Uploads a structure to the Server.
        Allowed input format: Bio.PDB Structure object, PDB formatted string
        Returns id for future services.
        """
        
        if isinstance(structure, Entity): # SMCRA
            s = self._smcra_to_str(structure)

        elif isinstance(structure, str): # String
            s = structure
        
        else:
            raise ValueError('Unknown format. Use SMCRA object or string.')
            
        u = urllib.urlopen("http://www.cmbi.ru.nl/wiwsd/rest/UploadPDB", s)
        x = xml.dom.minidom.parse(u)
        
        id = x.getElementsByTagName("response")[0].childNodes[0].data

        return id
        
        
    def PDBasXMLwithSymwithPolarH(self, id):
        """
        Adds Hydrogen Atoms to a Structure.
        """
        
        print _WARNING
        # Protonated Structure in XML Format
        h_s_xml = urllib.urlopen("http://www.cmbi.ru.nl/wiwsd/rest/PDBasXMLwithSymwithPolarH/id/" + id)
        self.raw = h_s_xml
        p = self.parser
        h_s_smcra = p.read(h_s_xml, 'WHATIF_Output')
        
        return h_s_smcra