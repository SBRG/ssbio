import os
import tempfile

import numpy as np
import pandas as pd

from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import aa1
from Bio.PDB.Polypeptide import one_to_three

import sys
sys.path.insert(0, os.path.abspath('..'))
from ssbio.properties.loader import Loader
from ssbio.properties.pdbstart import get_pdb_res_starts

io = PDBIO()
parser = PDBParser()
l = Loader()

# default is to keep atoms with alternate location ID 'A'
alt = 'A'
# default is to not keep "alt" identifier
keep_alt = 0
keep_atom_list = ['N', 'C', 'O', 'CA']


class NotDisordered(PDB.Select):
    """
    class to select non disordered atoms and those with the correct
    alternate location ID
    """

    def __init__(self, alt):
        self.alt = alt

    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        elif atom.get_altloc() == self.alt:
            if not keep_alt:
                atom.set_altloc(' ')
            return True
        else:
            return False


def setup_pdb_for_amber_monomer(df, file_name):
    '''Takes in any PDB file and
    (1) performs site-directed mutagenesis given a user-defined sequence change
    (2) removes all hetero atoms (e.g. water and metals/cofactors)
    (3) removes disordered atoms
    (4) keep only 1 chain (df needs to only contain the gene of interest)
    '''

    my_structure = l.structure_reader(file_name)
    model = my_structure[0]

    # adjust the alignment dataframe to add pdb residue numbers
    # adding pdb_start and pdb_stop as empty columns
    df['id_b_start'] = np.nan
    df['id_b_stop'] = np.nan
    # now adding in pdb_start and pdb_stop
    for pdb, chain, pdb_start in get_pdb_res_starts(file_name):
        subset = df[df.chain == chain]
        for idx, row in subset.iterrows():
            if row['type'] != 'deletion':
                adder = row['id_a_stop'] - row['id_a_start']
                pdb_stop = adder + pdb_start
                df.loc[idx, 'id_b_start'] = pdb_start
                df.loc[idx, 'id_b_stop'] = pdb_stop
                pdb_start = pdb_stop + 1

    # choosing a representative chain (multiple could align to one gene -
    # we choose the best aligning one)
    sorted_by_match_and_chain = df[df['type'] == 'match'].sort_values(
        by=['count', 'chain'], ascending=[False, True]).set_index(['chain', 'count']).index.tolist()
    chain_id = str(sorted_by_match_and_chain[0][0])

    # get mutations, if any
    df = df[df.type == 'mutation']

    # operate only on this chain
    chain = model[chain_id]
    residue_list = chain.child_list

    for residue in list(residue_list):
        # mutate residues according to sequence alignment data
        if len(df) > 0:
            if residue.id[1] in df[df.id_b.str.endswith(chain_id)].id_b_start.tolist():
                res_w = df[df.id_b.str.endswith(chain_id)][
                    df.id_b_start == residue.id[1]].id_b_aa.values[0]
                res_mut = df[df.id_b.str.endswith(chain_id)][
                    df.id_b_start == residue.id[1]].id_a_aa.values[0]

                if res_mut not in list(aa1):
                    warnings.warn(
                        '***UNKNOWN AMINO ACID IN UNIPROT SEQUENCE. SKIPPING***')
                    continue

                print(residue.id[1], residue.resname.upper(),
                      " | ", res_w, "mutate to:", res_mut)

                # Remove all atoms except protein backbone atoms:
                to_detach = []
                for atom in residue.child_dict.keys():
                    if atom in keep_atom_list:
                        pass
                    else:
                        to_detach.append(atom)
                for atom in to_detach:
                    residue.detach_child(atom)

                residue.resname = one_to_three(res_mut)
                # if residue is non-standard, make it a non-heteroatom
                if residue.id[0] != ' ':
                    resnum = residue.id[1]
                    residue.id = (' ', resnum, ' ')

    # Remove all hydrogens and heteroatom residues (e.g. WAT or metal) from
    # structure:
    for residue in list(chain):
        id = residue.id
        if id[0] != ' ':
            chain.detach_child(id)
        if len(chain) == 0:
            model.detach_child(chain.id)
        for atom in residue.get_list():
            # print residue.resname, residue.id[1], atom.element
            if atom.element == 'H':
                residue.detach_child(atom.id)

    return chain


def write_to_pdb(struct, new_file_name):
    io.set_structure(struct)
    io.save(new_file_name)


def monomer_mutation_function(corrected_df, pdb_file, custom_id=''):

    if len(custom_id) == 0:
        custom_id = os.path.basename(pdb_file)[:4]

    output_file = custom_id + '_modified.pdb'

    # read in structure to parser
    struct = l.structure_reader(pdb_file)
    # create PDBIO class for writing

    temp_pdb_file = os.path.join(
        tempfile.gettempdir(), '%s_no_disorder.pdb' % custom_id)
    io.set_structure(struct)
    io.save(temp_pdb_file, select=NotDisordered(alt))

    # load structure (no disorder) into biopython & perform modifications
    structure_new = setup_pdb_for_amber_monomer(
        corrected_df, temp_pdb_file)
    #  write out new structure
    write_to_pdb(structure_new, output_file)

    return output_file


if __name__ == '__main__':
    l.structure_reader('test_files/1kf6.pdb')

    df = pd.read_csv('test_files/tmp.csv')
    mutated_file = monomer_mutation_function(df, 'test_files/1kf6.pdb')
    print(mutated_file)
