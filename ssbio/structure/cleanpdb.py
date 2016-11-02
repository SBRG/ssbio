from Bio import PDB
import argparse
import glob
import os
import os.path as op

from ssbio.structure.pdbioext import PDBIOExt
from tqdm import tqdm


class CleanPDB(PDB.Select):
    """Selection rules to clean a PDB file

    These rules aim to:
    - Add missing chains to a PDB file
    - Select a single chain if noted
    - Remove alternate atom locations
    - Add atom occupancies
    - Add B (temperature) factors (default Biopython behavior)
    """

    def __init__(self, remove_atom_alt=True, remove_atom_hydrogen=True, keep_atom_alt_id='A', add_atom_occ=True,
                 remove_res_hetero=True, add_chain_id_if_empty='X', keep_chains=[]):
        """Initialize the parameters which indicate what cleaning will occur

        Args:
            remove_atom_alt:
            remove_atom_hydrogen:
            keep_atom_alt_id:
            add_atom_occ:
            remove_res_hetero:
            add_chain_id_if_empty:
            keep_chains:
        """
        self.remove_atom_alt = remove_atom_alt
        self.remove_atom_hydrogen = remove_atom_hydrogen
        self.keep_atom_alt_id = keep_atom_alt_id
        self.add_atom_occ = add_atom_occ
        self.remove_res_hetero = remove_res_hetero
        self.add_chain_id_if_empty = add_chain_id_if_empty
        self.keep_chains = keep_chains

    def accept_chain(self, chain):
        # If the chain does not have an ID, add one to it and keep it
        # http://comments.gmane.org/gmane.comp.python.bio.devel/10639
        if self.add_chain_id_if_empty and not chain.id.strip():
            chain.id = self.add_chain_id_if_empty
            return True
        # If a chain is specified and the current chain equals that specified chain, keep it
        elif self.keep_chains and chain.id in self.keep_chains:
            return True
        # If a chain is specified but the current chain does not equal that specified chain, remove it
        elif self.keep_chains and chain.id not in self.keep_chains:
            return False
        # If no chain is specified, keep all chains
        else:
            return True

    def accept_residue(self, residue):
        hetfield, resseq, icode = residue.get_id()
        # If you want to remove residues that are not normal, remove them
        if self.remove_res_hetero and hetfield[0] != ' ':
            return False
        else:
            return True

    def accept_atom(self, atom):
        # If the you want to remove hydrogens and the atom is a H, remove it
        if self.remove_atom_hydrogen and atom.element == 'H':
            return False
        # If you want to remove alternate locations, and the alternate location is not the one you want to keep, remove it
        elif self.remove_atom_alt and atom.is_disordered() and atom.get_altloc() != self.keep_atom_alt_id:
            return False
        else:
            # Add occupancies if there are none and you want to
            # http://comments.gmane.org/gmane.comp.python.bio.general/6289
            if self.add_atom_occ and atom.occupancy is None:
                atom.set_occupancy(1)
            if self.remove_atom_alt:
                atom.set_altloc(' ')
            return True


# def clean_pdb(infile, chain_id, outdir, outfile=''):
#     """Clean a PDB file and keep only the chain of interest.
#
#     Args:
#         infile: Path to PDB file
#         chain_id: Chain to keep
#         outdir: Directory to output the clean PDB file
#
#     Returns:
#
#     """
#     my_pdb = PDBIOExt(infile)
#     my_model = my_pdb.first_model
#
#     # clean pdb and save a file with only these chains of interest
#     # suffix appended with chains (1fat_A_B_cleaned.pdb)
#     log.debug('Cleaning PDB file {} and keeping chain {}...'.format(infile, chain_id))
#     my_cleaner = CleanPDB(keep_chains=[chain_id])
#
#     if output_name:
#         my_new_pdb_id = output_name
#     else:
#         my_new_pdb_id = '{}_{}_cleaned'.format(pdb_id, chain_id)
#
#     my_clean_pdb = my_pdb.write_pdb(custom_name=my_new_pdb_id, custom_ext='.pdb', out_suffix='',
#                                     out_dir=output_dir, custom_selection=my_cleaner)
#
#     return my_clean_pdb


# # TODO: does CleanPDB add the TERs?
# def add_ter_to_pdb(infile, outfile_name=None):
#     '''
#     Adds 'TER' cards to a PDB file when encountering:
#     - a OXT atom - indicating the end of an amino acid chain
#     - a ATOM to HETATM change - indicating a cofactor or ligand
#     - a HETATM change to a new residue - indicating a new cofactor or ligand
#     Input: any PDB file
#     Output: the path to the fixed pdb file
#     '''
#     with open(infile,'r') as pdb_file:
#         lines = pdb_file.readlines()
#
#     # open new file to write to
#     if not outfile_name:
#         outfile = os.path.splitext(os.path.basename(infile))[0] + '_fix.pdb'
#     else:
#         outfile = outfile_name
#
#     with open(outfile,'w') as new_pdb_file:
#
#         for line in lines:
#             # grab residue name to compare with previous residue name
#             resname = line[17:20]
#             # if AMBER added an OXT, that is usually the end of the protein chain
#             if 'OXT' in line:
#                 new_pdb_file.write(line)
#                 new_pdb_file.write('TER\n')
#
#             # TODO: this should be manual input
#             elif 'MG' in line:
#                 new_pdb_file.write(line)
#                 new_pdb_file.write('TER\n')
#
#             # if there is a change from ATOM to HETATM, that usually indicates the presence of a cofactor/ligand
#             # also check if the previous resname was different - could be a start of a new cofactor/ligand
#             elif 'HETATM' in line and ('ATOM' in lines[lines.index(line)-1] or resname != prev_resname):
#                 resname = line[17:20]
#                 new_pdb_file.write('TER\n')
#                 new_pdb_file.write(line)
#             else:
#                 new_pdb_file.write(line)
#             prev_resname = resname
#
#     return outfile


if __name__ == '__main__':
    # load inputs from command line

    p = argparse.ArgumentParser(description='Cleans a PDB file')
    p.add_argument('infile', help='PDB file or folder you want to clean', nargs='+', type=str)
    p.add_argument('--outsuffix', '-o', default='clean', help='Suffix appended to PDB file')
    p.add_argument('--chain', '-c', help='Keep only specified chains')
    p.add_argument('--keephydro', '-hy', action='store_false', help='Keep hydrogen atoms')
    p.add_argument('--keephetero', '-ht', action='store_false', help='Keep hetero atoms')
    # TODO: if this flag is present, the alternate positions seem to switch line positions
    p.add_argument('--keepalt', '-ka', action='store_false', help='Keep alternate positions')
    args = p.parse_args()

    out_dir = 'clean_pdbs'
    if args.chain:
        chains = args.chain.split(',')
    else:
        chains = args.chain

    # TODO: improve arg parsing for files/dirs
    if len(args.infile) == 1 and op.isdir(args.infile[0]):
        os.chdir(args.infile[0])
        pdbs = glob.glob('*.pdb')
    else:
        pdbs = args.infile

    for pdb in tqdm(pdbs):
        if op.isdir(pdb):
            continue
        # print('Cleaning PDB: {}'.format(pdb))
        my_pdb = PDBIOExt(pdb)
        my_cleaner = CleanPDB(remove_atom_alt=args.keepalt, remove_atom_hydrogen=args.keephydro, keep_atom_alt_id='A', add_atom_occ=True,
                              remove_res_hetero=args.keephetero, add_chain_id_if_empty='X', keep_chains=chains)
        my_clean_pdb = my_pdb.write_pdb(out_suffix=args.outsuffix, out_dir=out_dir, custom_selection=my_cleaner)
    print('Clean PDBs at: {}'.format(out_dir))