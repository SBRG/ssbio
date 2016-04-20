#!/usr/bin/env python

# ##############################
# ##  for standalone testing
# import sys
# import os.path as op
# new_path = op.join(op.expanduser('~'), 'Dropbox/Projects/ssbio')
# if new_path not in sys.path:
#     sys.path.append(new_path)
# ## end for standalone testing
# ##############################

from Bio import PDB

from Bio.PDB.Polypeptide import aa1
from Bio.PDB.Polypeptide import aa3
from Bio.PDB.Polypeptide import one_to_three

class MutatePDB(PDB.Select):
    """Selection rules to mutate a PDB file

    These rules aim to:
    - Mutate a specified residue number to a new amino acid
    """

    keep_atom_list = ['N', 'C', 'O', 'CA']

    def __init__(self, mutation_list):
        """Initialize the parameters which indicate what mutations will occur

        Args:
            chain:
            residue_number:
            mutate_to:
        """
        self.mutation_list = [(i[0], int(i[1]), self.standard_resname(i[2])) for i in mutation_list]
        self.chains_and_residues = [(i[0], int(i[1])) for i in mutation_list]

    def standard_resname(self, res):
        resname3 = res.upper()
        if resname3 not in list(aa3) and resname3 not in list(aa1):
            raise ValueError("Unrecognised residue {}".format(res))
        if len(resname3) == 1:
            resname3 = one_to_three(resname3)

        return resname3

    def accept_residue(self, residue):
        hetfield, resseq, icode = residue.get_id()

        chain = residue.get_parent()
        chain_id = chain.get_id()

        if (chain_id,resseq) in self.chains_and_residues:
            prev_resname = residue.resname
            get_index = self.chains_and_residues.index((chain_id,resseq))
            residue.resname = self.mutation_list[get_index][2]
            print("Mutated {0}.{1}.{2} to {0}.{1}.{3}".format(chain_id, resseq, prev_resname, residue.resname))
        return True

    def accept_atom(self, atom):
        residue = atom.get_parent()
        hetfield, resseq, icode = residue.get_id()

        chain = residue.get_parent()
        chain_id = chain.get_id()

        if (chain_id,resseq) in self.chains_and_residues and atom.get_id() not in self.keep_atom_list:
            # print("Removing atom {}.{}.{}".format(chain_id, resseq, atom.get_id()))
            return False

        return True

def parse_mutation_input(instr):
    init_split = instr.split(',')
    second_split = [tuple(i.split('.')) for i in init_split]
    return second_split

if __name__ == '__main__':
    import tempfile
    from ssbio.tools.pdbioext import PDBIOExt
    from ssbio.tools.cleanpdb import CleanPDB

    # # load inputs from command line
    import argparse
    p = argparse.ArgumentParser(description='Mutates a PDB file')
    p.add_argument('infile', help='PDB file you want to mutate')
    p.add_argument('mutations', help='Mutations in the form of Chain1.ResNum1.Mutation1,Chain2.ResNum2.Mutation2. Example: A.4.TYR,B.4.TYR')
    p.add_argument('--outsuffix', '-o', default='mutated', help='Suffix appended to PDB file')
    p.add_argument('--clean', '-c', action='store_true', help='Clean PDB and keep only chain with mutation')
    args = p.parse_args()

    mutations = parse_mutation_input(args.mutations)

    my_pdb = PDBIOExt(args.infile)
    if args.clean:
        my_cleaner = CleanPDB(keep_chains=[x[0] for x in mutations])
        my_clean_pdb = my_pdb.write_pdb(out_suffix='clean', out_dir=tempfile.gettempdir(), custom_selection=my_cleaner)
        my_pdb = PDBIOExt(my_clean_pdb)

    my_mutation = MutatePDB(mutations)
    my_mutated_pdb = my_pdb.write_pdb(out_suffix=args.outsuffix, out_dir='mutated_pdbs', custom_selection=my_mutation)
    print('Mutated PDB at: {}'.format(my_mutated_pdb))
