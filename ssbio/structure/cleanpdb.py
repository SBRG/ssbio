from Bio import PDB
import argparse
import textwrap
import os
import ssbio.utils
import os.path as op

from ssbio.structure.pdbioext import PDBIOExt
from tqdm import tqdm

import logging
log = logging.getLogger(__name__)


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
                 remove_res_hetero=True, add_chain_id_if_empty='X', keep_chains=None):
        """Initialize the parameters which indicate what cleaning will occur

        Args:
            remove_atom_alt: Remove alternate positions
            remove_atom_hydrogen: Remove hydrogen atoms
            keep_atom_alt_id: If removing alternate positions, which alternate ID to keep
            add_atom_occ: Add atom occupancy fields if not present
            remove_res_hetero: Remove all HETATMs
            add_chain_id_if_empty: Add a chain ID if not present
            keep_chains: Keep only these chains
        """
        self.remove_atom_alt = remove_atom_alt
        self.remove_atom_hydrogen = remove_atom_hydrogen
        self.keep_atom_alt_id = keep_atom_alt_id
        self.add_atom_occ = add_atom_occ
        self.remove_res_hetero = remove_res_hetero
        self.add_chain_id_if_empty = add_chain_id_if_empty
        if not keep_chains:
            self.keep_chains = []
        else:
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
        # If you want to remove alternate locations, and the alt location isn't the one you want to keep, remove it
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


if __name__ == '__main__':
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=textwrap.dedent("""\
                                Clean PDB files - cleanpdb.py
                                -----------------------------
                                This script will automatically:

                                * Add missing chains to a PDB file
                                * Select a single chain or chains if noted
                                * Remove alternate atom locations
                                * Add atom occupancies
                                * Add B (temperature) factors (default Biopython behavior)

                                Cleaned PDBs will be in a clean_pdbs folder where the script is executed.
                                Example: script help
                                $ cleanpdb --help

                                Example: clean one PDB file
                                $ cleanpdb 1kf6.pdb

                                Example: clean one PDB file and keep only chains A and B
                                $ cleanpdb 1kf6.pdb --chain A,B

                                Example: clean multiple PDB files
                                $ cleanpdb *.pdb

                                Example: clean a whole directory of PDB
                                $ cleanpdb /path/to/pdb/files
                                """))
    p.add_argument('infile', help='PDB file or folder you want to clean', nargs='+', type=str)
    p.add_argument('--outsuffix', '-os', default='_clean', help='Suffix appended to PDB file')
    p.add_argument('--outdir', '-od', default='clean_pdbs', help='Directory to output clean PDBs')
    p.add_argument('--chain', '-c', default=None, help='Keep only specified chains')
    p.add_argument('--keephydro', '-hy', action='store_false', help='Keep hydrogen atoms (default is to remove)')
    p.add_argument('--keephetero', '-ht', action='store_false', help='Keep hetero atoms (default is to remove)')
    # TODO: if this flag is present, the alternate positions seem to switch line positions
    p.add_argument('--keepalt', '-ka', action='store_false', help='Keep alternate positions (default is to remove)')
    p.add_argument('--force', '-f', action='store_true', help='Force rerunning of cleaning even if the clean PDB exists')
    args = p.parse_args()

    if args.chain:
        args.chain = args.chain.split(',')

    if not op.isdir(args.outdir):
        os.mkdir(args.outdir)

    infiles = ssbio.utils.input_list_parser(args.infile)

    for pdb in tqdm(infiles):

        outfile = ssbio.utils.outfile_maker(inname=pdb,
                                            append_to_name=args.outsuffix,
                                            outdir=args.outdir,
                                            outext='.pdb')

        if ssbio.utils.force_rerun(flag=args.force, outfile=outfile):

            my_pdb = PDBIOExt(pdb, file_type='pdb')
            my_cleaner = CleanPDB(remove_atom_alt=args.keepalt,
                                  remove_atom_hydrogen=args.keephydro,
                                  keep_atom_alt_id='A',
                                  add_atom_occ=True,
                                  remove_res_hetero=args.keephetero,
                                  add_chain_id_if_empty='X',
                                  keep_chains=args.chain)

            my_clean_pdb = my_pdb.write_pdb(out_suffix=args.outsuffix,
                                            out_dir=args.outdir,
                                            custom_selection=my_cleaner,
                                            force_rerun=args.force)

    print('Clean PDBs at: {}'.format(args.outdir))
