# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from copy import deepcopy

from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import is_aa

class Protein(Structure):
    
    def __init__(self, struct_id):
        
        Structure.__init__(self, struct_id)
        
    @classmethod
    def from_structure(cls, original, filter_residues):
        """
        Loads structure as a protein, exposing
        protein-specific methods.
        """
        P = cls(original.id)
        P.full_id = original.full_id
        
        for child in original.child_dict.values():
            copycat = deepcopy(child)
            P.add(copycat)
        
        # Discriminate non-residues (is_aa function)
        remove_list = []
        if filter_residues:
            for model in P:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] != ' ' or not is_aa(residue):
                            remove_list.append(residue)
            
            for residue in remove_list:
                residue.parent.detach_child(residue.id)
            
            for chain in P.get_chains(): # Remove empty chains
                if not len(chain.child_list):
                    model.detach_child(chain.id)
                    
        P.header = deepcopy(original.header)
        P.xtra = deepcopy(original.xtra)
        
        return P
          
    def search_ss_bonds(self, threshold=3.0):
        """ Searches S-S bonds based on distances
            between atoms in the structure (first model only).
            Average distance is 2.05A. Threshold is 3A default.
            Returns iterator with tuples of residues.
        """
        
        # Taken from http://docs.python.org/library/itertools.html
        # Python 2.4 does not include itertools.combinations

        def combinations(iterable, r):
            # combinations('ABCD', 2) --> AB AC AD BC BD CD
            # combinations(range(4), 3) --> 012 013 023 123
            pool = tuple(iterable)
            n = len(pool)
            if r > n:
                return
            indices = range(r)
            yield tuple(pool[i] for i in indices)
            while True:
                for i in reversed(range(r)):
                    if indices[i] != i + n - r:
                        break
                else:
                    return
                indices[i] += 1
                for j in range(i+1, r):
                    indices[j] = indices[j-1] + 1
                yield tuple(pool[i] for i in indices)

        
        
        model = self.child_list[0]
        cysteines = [r for r in model.get_residues() if r.get_resname() == 'CYS']

        pairs = combinations(cysteines, 2) # Iterator with pairs

        for cys_pair in pairs:
            if cys_pair[0]['SG'] - cys_pair[1]['SG'] < threshold:
                yield cys_pair
    
    def check_missing_atoms(self, template=None, ha_only=True):
        """
        Checks for missing atoms based on a template.
        Default: Searches for missing heavy atoms (not Hydrogen) based on Bio.Struct.protein_residues
        
        Arguments:
          - template, dictionary, keys are residue names, values list of atom names.
          - ha_only, boolean, default True, restrict check to heavy atoms.
          
        Returns a dictionary of tuples with the missing atoms per residue.
        """
        
        missing_atoms = {}
        
        if not template:
            import protein_residues
            template = protein_residues.normal # Don't care for terminal residues here..
            
        for residue in self.get_residues():
            
            if not template.has_key(residue.resname):
                # Maybe add this as a warning instead of exception?
                raise ValueError('Residue name (%s) not in the template' %residue.resname )
            
            if ha_only:
                heavy_atoms = [ atom for atom in template[residue.resname]['atoms'].keys() 
                                if atom[0] != 'H' and not (atom[0].isdigit() and atom[1] == 'H')]
                reference_set = set(heavy_atoms)
            else:
                reference_set = set(template[residue.resname]['atoms'].keys())
            
            structure_set = set(residue.child_dict.keys())
            
            diff = reference_set.difference(structure_set)
            
            if diff:
                residue_uniq_id = (residue.parent.id, residue.resname, residue.get_id()[1]) # Chain, Name, Number
                missing_atoms[residue_uniq_id] = list(diff)
        
        return missing_atoms
        
    
    def coarse_grain(self, cg_type="CA_TRACE"):
        """ 
            Reduces the protein structure complexity to a few (pseudo-)atoms per residue.

            Parameters:
              - cg_type:      CA_TRACE (Ca-only) [Default]
                              ENCAD_3P (CA, O, SC Beads)
                              MARTINI (CA, O, SC Beads)
            
            Returns a new structure object.
        """
        
        # Import CG Types
        
        import CG_Models
        
        CG_Library = {  "CA_TRACE": CG_Models.CA_TRACE,
                        "ENCAD_3P": CG_Models.ENCAD_3P,
                        "MARTINI":  CG_Models.MARTINI   }
        
        CG_Method = CG_Library[cg_type]
        
        # Creates a brand new structure object
        from Bio.PDB.StructureBuilder import StructureBuilder
        
        structure_builder=StructureBuilder()
        
        cg_id = "CG_" + self.id
        structure_builder.init_structure(cg_id)
        structure_builder.init_seg(' ') # Empty SEGID

        for model in self:
            structure_builder.init_model(model.id)
            
            for chain in model:
                structure_builder.init_chain(chain.id)
                cur_chain = structure_builder.chain

                for residue in chain:
                    cg_residue = CG_Method(residue)
                    cur_chain.add(cg_residue)
        
        cg_structure = structure_builder.get_structure()
        
        return cg_structure


    def find_seq_homologues(self, return_raw=False):
        """
        Uses NCBI BLAST to look for structures deposited in the PDB database
        that share __sequence__ homology with the target protein/chain.
        Bridges to Bio.BLAST.NCBIWWW.
        """
        
        
        # Get sequence from structure/chain
        # We could use Bio.PDB.Polypeptide?
        
        from Bio.SCOP.Raf import to_one_letter_code
        
        s = self
        seq_iter = s.get_residues()
        seq_str = ''

        for aa in seq_iter:
            if aa.resname in to_one_letter_code:
                seq_str += to_one_letter_code[aa.resname]
        
        # Use BLAST to find homologous sequences with associated
        # structures in the PDB database.
        # Perhaps include local BLAST?
        
        from Bio.Blast.NCBIWWW import qblast

        # Adapt for short query sequences if needed
        # From http://www.ncbi.nlm.nih.gov/blast/producttable.shtml#shortp
        
        if len(seq_str) < 15:
            word_size = 2
            expect = 20000
            matrix_name = 'PAM30'
            filter = None
        else:
            word_size = 3
            expect = 10.0
            matrix_name = 'BLOSUM62'
            filter = 'SEG'
             
        query_result = qblast(  "blastp", "pdb", seq_str, 
                                word_size, expect, matrix_name, filter)
        
        if return_raw:
           return query_result 
            
        # Parse BLAST result to yield results
        # PDBID : (E-Value, Identity, Positives, Gaps, Alignment)
        
        from Bio.Blast import NCBIXML
        
        blast_records = NCBIXML.read(query_result)
        
        results = []
        
        for alignment in blast_records.alignments:
            for hsps in alignment.hsps:
                id_perc = "%s/%s" %(hsps.identities, alignment.length)
                pos_perc = "%s/%s" %(hsps.positives, alignment.length)
                gaps_perc = "%s/%s" %(hsps.gaps, alignment.length)

                pdb_id = alignment.title.split('|')[3]
                e_value = hsps.expect
                
                results.append( ( pdb_id,
                                "%2.5e" %e_value,
                                id_perc, 
                                pos_perc, 
                                gaps_perc,
                                '\n'.join([hsps.query, hsps.match, hsps.sbjct]) ) )

        return results