#!/usr/bin/env python

##############################
##  for standalone testing
import sys
new_path = '/home/nathan/Dropbox/Projects/ssbio/'
if new_path not in sys.path:
    sys.path.append(new_path)
## end for standalone testing
##############################

from Bio.PDB.Polypeptide import aa1
from Bio.PDB.Polypeptide import aa3
from Bio.PDB.Polypeptide import one_to_three

import ssbio.tools.iotools as iotools

class MutatePDB:
    """Mutate a PDB file
    """

    keep_atom_list = ['N', 'C', 'O', 'CA']

    def __init__(self, in_file):
        """Return a MutatePDB object which is ready for mutating.

        Attributes are simply the input PDB file name and the first model.

        Args:
            in_file: PDB input file path
        """
        l = iotools.IOTools()
        structure = l.structure_reader(in_file)
        self.in_file = in_file
        self.model = structure[0]

    def write_pdb(self, out_suffix='mutated', out_dir=None):
        """Write a new PDB file from a Biopython Model object and a given filename, appended with a suffix.

        Args:
            out_suffix: string to append to new PDB file - default is "_mutated"
            out_dir: optional directory to output the file

        Returns:
            out_file: filepath of new PDB file

        """

        # Prepare the output file path
        out_file = self._output_filepath(out_suffix=out_suffix, out_dir=out_dir)

        # IO object creation
        io = PDBIO()
        io.set_structure(self.model)

        if not_disordered:
            io.save(out_file, select=NotDisordered())
        else:
            io.save(out_file)

        return out_file

    def mutate_chain_residue(self, chain_id, residue_number, mutate_to):
        """Mutates a specified residue number to another residue in the specified chain.

        Args:
            chain_id:
            residue_number:
            mutate_to:

        Returns:

        """
        # Standardizing the residue to mutate to
        mutate_to = mutate_to.upper()
        if mutate_to not in list(aa3) and mutate_to not in list(aa1):
            raise ValueError("Unrecognised residue {}".format(mutate_to))
        if len(mutate_to) == 1:
            mutate_to = one_to_three(mutate_to)

        chain = self.model[chain_id]
        residue = chain[residue_number]

        # Remove all atoms except protein backbone atoms
        to_detach = []
        for atom in residue.child_dict.keys():
            if atom in self.keep_atom_list:
                pass
            else:
                to_detach.append(atom)
        for atom in to_detach:
            residue.detach_child(atom)

        # Mutate the residue
        residue.resname = mutate_to

        # TODO: investigate this action
        # If residue is non-standard, make it a non-heteroatom
        if residue.id[0] != ' ':
            residue.id = (' ', residue.id[1], ' ')

        return model

# def monomer_mutation_function(corrected_df, pdb_file, custom_id='', custom_dir=''):
#
#     if len(custom_id) == 0:
#         custom_id = os.path.basename(pdb_file)[:4]
#         output_file = os.path.basename(pdb_file)[:4] + '_modified.pdb'
#     else:
#         output_file = custom_id + '.pdb'
#
#     # read in structure to parser
#     struct = l.structure_reader(pdb_file)
#     # create PDBIO class for writing
#
#     temp_pdb_file = os.path.join(
#         tempfile.gettempdir(), '%s_no_disorder.pdb' % custom_id)
#     io.set_structure(struct)
#     io.save(temp_pdb_file, select=NotDisordered(alt))
#
#     # load structure (no disorder) into biopython & perform modifications
#     structure_new = setup_pdb_for_amber_monomer(
#         corrected_df, temp_pdb_file)
#     #  write out new structure
#     write_to_pdb(structure_new, os.path.join(custom_dir, output_file))
#
#     return os.path.join(custom_dir, output_file)
#
# def setup_pdb_for_amber_monomer(df, file_name):
#     '''Takes in any PDB file and
#     (1) performs site-directed mutagenesis given a user-defined sequence change
#     (2) removes all hetero atoms (e.g. water and metals/cofactors)
#     (3) removes disordered atoms
#     (4) keep only 1 chain (df needs to only contain the gene of interest)
#     '''
#
#     my_structure = l.structure_reader(file_name)
#     model = my_structure[0]
#
#     # adjust the alignment dataframe to add pdb residue numbers
#     # adding pdb_start and pdb_stop as empty columns
#     df['id_b_start'] = np.nan
#     df['id_b_stop'] = np.nan
#     # now adding in pdb_start and pdb_stop
#     for pdb, chain, pdb_start in get_pdb_res_starts(file_name):
#         subset = df[df.chain == chain]
#         for idx, row in subset.iterrows():
#             if row['type'] != 'deletion':
#                 adder = row['id_a_stop'] - row['id_a_start']
#                 pdb_stop = adder + pdb_start
#                 df.loc[idx, 'id_b_start'] = pdb_start
#                 df.loc[idx, 'id_b_stop'] = pdb_stop
#                 pdb_start = pdb_stop + 1
#
#     # choosing a representative chain (multiple could align to one gene -
#     # we choose the best aligning one)
#     sorted_by_match_and_chain = df[df['type'] == 'match'].sort_values(
#         by=['count', 'chain'], ascending=[False, True]).set_index(['chain', 'count']).index.tolist()
#     chain_id = str(sorted_by_match_and_chain[0][0])
#
#     # get mutations, if any
#     df = df[df.type == 'mutation']
#
#     # operate only on this chain
#     chain = model[chain_id]
#     residue_list = chain.child_list
#
#     for residue in list(residue_list):
#         # mutate residues according to sequence alignment data
#         if len(df) > 0:
#             if residue.id[1] in df[df.id_b.str.endswith(chain_id)].id_b_start.tolist():
#                 res_w = df[df.id_b.str.endswith(chain_id)][
#                     df.id_b_start == residue.id[1]].id_b_aa.values[0]
#                 res_mut = df[df.id_b.str.endswith(chain_id)][
#                     df.id_b_start == residue.id[1]].id_a_aa.values[0]
#
#                 if res_mut not in list(aa1):
#                     warnings.warn(
#                         '***UNKNOWN AMINO ACID IN UNIPROT SEQUENCE. SKIPPING***')
#                     continue
#
#                 # print(residue.id[1], residue.resname.upper(),
#                     #   " | ", res_w, "mutate to:", res_mut)
#
#                 # Remove all atoms except protein backbone atoms:
#                 to_detach = []
#                 for atom in residue.child_dict.keys():
#                     if atom in keep_atom_list:
#                         pass
#                     else:
#                         to_detach.append(atom)
#                 for atom in to_detach:
#                     residue.detach_child(atom)
#
#                 residue.resname = one_to_three(res_mut)
#                 # if residue is non-standard, make it a non-heteroatom
#                 if residue.id[0] != ' ':
#                     resnum = residue.id[1]
#                     residue.id = (' ', resnum, ' ')
#
#     # Remove all hydrogens and heteroatom residues (e.g. WAT or metal) from
#     # structure:
#     for residue in list(chain):
#         id = residue.id
#         if id[0] != ' ':
#             chain.detach_child(id)
#         if len(chain) == 0:
#             model.detach_child(chain.id)
#         for atom in residue.get_list():
#             # print residue.resname, residue.id[1], atom.element
#             if atom.element == 'H':
#                 residue.detach_child(atom.id)
#
#     return chain

if __name__ == '__main__':
    # load inputs from command line
    import argparse

    p = argparse.ArgumentParser(description='Mutates a PDB file')
    p.add_argument('infile', help='PDB file you want to mutate')
    # p.add_argument('chain', help='Chain in which the residue is contained')
    # p.add_argument('resnum', help='Residue number you want to mutate')
    # p.add_argument('resmut', help='Residue you want to mutate to')
    p.add_argument('--clean', '-c', action='store_false', help='Clean the PDB and only keep the chain of interest')
    args = p.parse_args()

    my_pdb = MutatePDB(args.infile)
    my_pdb.mutate_chain_residue('A', 9, 'I')