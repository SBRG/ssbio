#!/usr/bin/env python

##############################
##  for standalone testing
import sys
new_path = '/home/nathan/Dropbox/Projects/ssbio/'
if new_path not in sys.path:
    sys.path.append(new_path)
## end for standalone testing
##############################

import os.path as op
from Bio import PDB
from Bio.PDB import PDBIO

import ssbio.tools.iotools as iotools


class NotDisordered(PDB.Select):
    """
    Class to select non disordered atoms and those with the correct
    alternate location ID - http://biopython.org/wiki/Remove_PDB_disordered_atoms
    """

    # default is to not keep "alt" identifier
    # default is to keep atoms with alternate location ID 'A'
    def __init__(self, alt='A', keep_alt=False):
        self.alt = alt
        self.keep_alt = keep_alt

    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        elif atom.get_altloc() == self.alt:
            if not self.keep_alt:
                atom.set_altloc(' ')
            return True
        else:
            return False


class CleanPDB:
    """Various tools to clean PDB files.

    These functions aim to:
    - Add missing chains to a PDB file
    - Add atom occupancies
    - Add B (temperature) factors
    - Select only chains of interest
    - Mutate residues (to be run through AMBER leap)
    """
    def __init__(self, in_file):
        """Return a CleanPDB object which is ready for cleaning.

        Attributes are simply the input PDB file name and the first model.

        Args:
            in_file: PDB input file path
        """
        l = iotools.IOTools()
        structure = l.structure_reader(in_file)
        self.in_file = in_file
        self.model = structure[0]

    def _output_filepath(self, out_suffix, out_dir=None):
        """Parses a PDB input filename and returns a output file path to write a modified file.

        Args:
            out_suffix (str): string to append to the filename of the new PDB file
            out_dir (str): optional working directory where cleaned PDB file should be written to

        Returns:
            out_file (str): file path of the new PDB file
        """

        # Parsing the input filename
        filename_full = op.basename(self.in_file)
        filename, ext = op.splitext(filename_full)

        # Assembling the new output filename
        out_file = '{}_{}{}'.format(filename, out_suffix, ext)

        if out_dir:
            out_file = op.join(out_dir, out_file)

        return out_file

    def write_pdb(self, out_suffix='modified', out_dir=None, not_disordered=False):
        """Write a new PDB file from a Biopython Model object and a given filename, appended with a suffix.

        Set not_disordered to True to remove alternate locations of atoms.

        Args:
            out_suffix: string to append to new PDB file - default is "_modified"
            out_dir: optional directory to output the file
            not_disordered: optional flag to remove alternate locations of atoms

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

    def add_chain_id(self, chain_id='X', write_file=False):
        """Add missing chain IDs - default is X

        See: http://comments.gmane.org/gmane.comp.python.bio.devel/10639

        Args:
            chain_id: chain ID to add to empty chains
            write_file: flag to write a PDB file

        Returns:
            model: New Biopython Model object with added chain IDs
        """
        for chain in self.model.get_chains():
            if not chain.id.strip():  # chain could be an empty string ' ' so strip it!
                chain.id = chain_id

        if write_file:
            self.write_pdb(out_suffix='chainAdded')

        return self.model

    def add_occupancies(self, write_outfile=False):
        """Adding occupancies if there are none

        See: http://comments.gmane.org/gmane.comp.python.bio.general/6289

        Args:
            write_outfile: flag to write a PDB file

        Returns:
            model: New Biopython Model object with added chain IDs

        """
        for atom in self.model.get_atoms():
            if atom.occupancy is None:
                atom.set_occupancy(1)

        if write_outfile:
            self.write_pdb(out_suffix='occupancyAdded')

        return self.model

    def strip_hydro_hetero(self, write_outfile=False):
        """Remove all hydrogens and heteroatom residues (e.g. WAT or metal)

        Args:
            write_outfile: flag to write a PDB file

        Returns:
            model: New Biopython Model object with hydrogens and heteroatoms removed
        """

        for chain_object in self.model.get_chains():
            for residue in list(chain_object):
                # print(residue)
                res_id = residue.id
                if res_id[0] != ' ':
                    chain_object.detach_child(res_id)
                if len(chain_object) == 0:
                    self.model.detach_child(chain_object.id)
                for atom in residue.get_list():
                    if atom.element == 'H':
                        residue.detach_child(atom.id)

        if write_outfile:
            self.write_pdb(out_suffix='stripped')

        return self.model

    def clean_pdb(self, add_chain=True, add_occ=True, strip_hydro_hetero=True, not_disordered=True,
                  out_suffix='clean', out_dir=None):
        """All in one function to clean a PDB file and save a new file

        This function does by default these functions to a PDB file:
        1) Adding a chain ID ('X' if there is none)
        2) Adding atom occupancies if there are none (default 1)
        3) Removes all hetero atoms and hydrogens (e.g. water and metals/cofactors)
        4) Save only the "A" alternate locations (removes disordered atoms)
        5) Adds B-factors (default functionality of Bio.PDB)

        Args:
            add_chain: add chain ID if missing
            add_occ: add occupancies if missing
            strip_hydro_hetero: strip hydrogens and hetero atoms
            not_disordered: save only alternate location A
            out_suffix: string appended to new PDB file
            out_dir: directory where new PDB file will be saved to

        Returns:
            cleaned_pdb: file path of new PDB file

        """

        if add_chain:
            self.add_chain_id()
        if add_occ:
            self.add_occupancies()
        if strip_hydro_hetero:
            self.strip_hydro_hetero()

        cleaned_pdb = self.write_pdb(out_suffix=out_suffix, out_dir=out_dir, not_disordered=not_disordered)

        return cleaned_pdb


if __name__ == '__main__':
    # load inputs from command line
    import argparse
    p = argparse.ArgumentParser(description='Cleans a PDB file')
    p.add_argument('infile', help='PDB file you want to clean')
    p.add_argument('--nostrip', '-ns', action='store_false')
    p.add_argument('--keepalt', '-ka', action='store_false')
    args = p.parse_args()

    my_pdb = CleanPDB(args.infile)
    my_pdb.clean_pdb(strip_hydro_hetero=args.nostrip, not_disordered=args.keepalt)
