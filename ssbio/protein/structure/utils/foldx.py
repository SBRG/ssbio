"""
FoldX
=====
"""

import os
import logging
import shutil
import ssbio.utils
import os.path as op
import pandas as pd
from ssbio.core.object import Object
log = logging.getLogger(__name__)


class FoldX(Object):

    """Class to run various commands from the FoldX suite of tools on a protein structure file.

    Args:
        pdb_path (str): Path to PDB file (PDB file format only!)
        project_id (str): Name of your structure or this mini FoldX project
        description (str): Short description of the PDB file being run on, ie. describe the oligomeric state.
        rotabase_path (str): Path to the rotabase.txt file
        foldx_exec (str, optional): Path to the FoldX executable, if empty, tries to execute ``foldx`` from the shell
        root_dir (str, optional): Path to directory where a new directory named after ``project_id`` and
            ``description`` with "_foldx" appended to it will be created.

    Todo:
        - Need checks to ensure only PDB files are being used
        - STDOUT and STDERR logging in a file

    """

    def __init__(self, pdb_path, project_id, description, rotabase_path, foldx_exec=None, root_dir=None):
        super(FoldX, self).__init__(id=project_id, description=description)

        # Create directories
        self._foldx_dirname = '{}_{}_foldx'.format(self.id, self.description)
        self._root_dir = None
        self.root_dir = root_dir

        # FoldX related
        self.foldx_exec = foldx_exec
        if not foldx_exec:
            self.foldx_exec = 'foldx'
        self.rotabase_path = rotabase_path

        # Copy PDB file
        self.pdb_path = shutil.copy2(pdb_path, self.foldx_dir)
        self.pdb_file = op.basename(self.pdb_path)

        # Copy rotabase.txt
        self.rotabase_path = shutil.copy2(self.rotabase_path, self.foldx_dir)

        # Output files
        self.repaired_pdb_outfile = None
        self.mutation_infile = None
        self.mutation_ddG_avg_outfile = None
        self.mutation_ddG_raw_outfile = None

        # To keep track of results
        self.mutation_index_to_group = {}

    @property
    def root_dir(self):
        """str: Path to where the folder named by this protein's ID will be created. Default is current working
        directory."""
        return self._root_dir

    @root_dir.setter
    def root_dir(self, path):
        if not path:
            path = os.getcwd()

        if not op.exists(path):
            raise ValueError('{}: folder does not exist'.format(path))

        if self._root_dir:
            log.debug('Changing root directory of FoldX {} from {} to {}'.format(self._foldx_dirname, self.root_dir, path))

            if not op.exists(op.join(path, self.id)):
                raise IOError('FoldX {} does not exist in folder {}'.format(self._foldx_dirname, path))

        self._root_dir = path

        for d in [self.foldx_dir]:
            ssbio.utils.make_dir(d)

    @property
    def foldx_dir(self):
        """str: FoldX folder"""
        if self.root_dir:
            return op.join(self.root_dir, self._foldx_dirname)
        else:
            log.warning('Root directory not set')
            return None

    @property
    def df_mutation_ddG_avg(self):
        return pd.read_csv(op.join(self.foldx_dir, self.mutation_ddG_avg_outfile), skiprows=8, sep='\t')

    @property
    def df_mutation_ddG_raw(self):
        return pd.read_csv(op.join(self.foldx_dir, self.mutation_ddG_raw_outfile), skiprows=8, sep='\t')

    def run_repair_pdb(self, silent=False, force_rerun=False):
        """Run FoldX RepairPDB on this PDB file.

        Original command::

            foldx --command=RepairPDB --pdb=4bxi.pdb

        Args:
            silent (bool): If FoldX output should be silenced from printing to the shell.
            force_rerun (bool): If FoldX RepairPDB should be rerun even if a repaired file exists.

        """
        # Create RepairPDB command
        foldx_repair_pdb = 'foldx --command=RepairPDB --pdb={}'.format(self.pdb_file)

        # Repaired PDB output file name
        foldx_repair_outfile = '{}_Repair.pdb'.format(op.splitext(self.pdb_file)[0])

        # Run RepairPDB
        ssbio.utils.command_runner(shell_command=foldx_repair_pdb, force_rerun_flag=force_rerun, silent=silent,
                                   outfile_checker=foldx_repair_outfile, cwd=self.foldx_dir)

        # TODO: write stdout/stderr to log file somewhere!

        self.repaired_pdb_outfile = foldx_repair_outfile

    def create_mutation_file(self, list_of_tuples):
        """Create the FoldX file 'individual_list.txt' to run BuildModel upon.

        Args:
            list_of_tuples (list): A list of tuples indicating mutation groups to carry out BuildModel upon. Example::

                [
                    (('N', 'A', 308, 'S'), ('S', 'A', 320, 'T'), ('S', 'A', 321, 'H')),  # Mutation group 1
                    (('S', 'A', 321, 'R'), ('T', 'A', 345, 'S'))  # Mutation group 2
                ]

        """

        self.mutation_infile = op.join(self.foldx_dir, 'individual_list.txt')

        idx = 1

        with open(self.mutation_infile, 'w') as f:
            for mutant_group in list_of_tuples:
                # Write the mutation string to the file
                mutstring = ''.join(list(map(lambda x: '{}{}{}{};'.format(x[0], x[1], x[2], x[3]), mutant_group)))
                f.write(mutstring + '\n')

                # Also keep track of the index being used for this mutation
                self.mutation_index_to_group[idx] = mutant_group
                idx += 1

    def create_random_mutation_file(self, list_of_tuples, original_sequence,
                                    randomize_resnums=False, randomize_resids=False,
                                    skip_resnums=None):
        """Create the FoldX file 'individual_list.txt', but randomize the mutation numbers or residues that were input.

        The randomize combinations can be a little confusing - this is what can happen:

            - randomize_resnums=False, randomize_resids=False: no change, original mutations are carried out
            - randomize_resnums=True, randomize_resids=False: mutations of resid X to resid Y will be carried out,
                but on a different residue number where resid X is found
            - randomize_resnums=False, randomize_resids=True: mutations of residue X# to a random residue will be
                carried out
            - randomize_resnums=True, randomize_resids=True: original mutations will be ignored, random mutation of
                any residue will be carried out

        Args:
            list_of_tuples (list): A list of tuples indicating mutation groups to be randomized.
            original_sequence (str, Seq, SeqRecord): Original amino acid sequence
            randomize_resnums (bool): If residue numbers should be randomized
            randomize_resids (bool): If residues themselves should be randomized
            skip_resnums (list):

        """

        import random

        def find(s, ch):
            return [i for i, ltr in enumerate(s) if ltr == ch]



    def run_build_model(self, num_runs=5, silent=False, force_rerun=False):
        """Run FoldX BuildModel command with a mutant file input.

        Original command::

            foldx --command=BuildModel --pdb=4bxi_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=5

        Args:
            num_runs (int):
            silent (bool): If FoldX output should be silenced from printing to the shell.
            force_rerun (bool): If FoldX BuildModel should be rerun even if the results file exists.

        """
        # BuildModel output files
        self.mutation_ddG_avg_outfile = 'Average_{}.fxout'.format(op.splitext(self.repaired_pdb_outfile)[0])
        self.mutation_ddG_raw_outfile = 'Raw_{}.fxout'.format(op.splitext(self.repaired_pdb_outfile)[0])

        # BuildModel command
        foldx_build_model = 'foldx --command=BuildModel --pdb={} --mutant-file={} --numberOfRuns={}'.format(self.repaired_pdb_outfile,
                                                                                                            op.basename(self.mutation_infile),
                                                                                                            num_runs)

        ssbio.utils.command_runner(shell_command=foldx_build_model, force_rerun_flag=force_rerun, silent=silent,
                                   outfile_checker=self.mutation_ddG_avg_outfile, cwd=self.foldx_dir)

    def get_ddG_results(self):
        """Parse the results from BuildModel and get the delta delta G's.

        A positive ddG means that the mutation(s) is destabilzing, negative means stabilizing.

            - highly stabilising (ΔΔG < −1.84 kcal/mol);
            - stabilising (−1.84 kcal/mol ≤ ΔΔG < −0.92 kcal/mol);
            - slightly stabilising (−0.92 kcal/mol ≤ ΔΔG < −0.46 kcal/mol);
            - neutral (−0.46 kcal/mol < ΔΔG ≤ +0.46 kcal/mol);
            - slightly destabilising (+0.46 kcal/mol < ΔΔG ≤ +0.92 kcal/mol);
            - destabilising (+0.92 kcal/mol < ΔΔG ≤ +1.84 kcal/mol);
            - highly destabilising (ΔΔG > +1.84 kcal/mol).

        Returns:
            dict: Dictionary of mutation group to predicted ddG.

        """
        foldx_avg_df = self.df_mutation_ddG_avg

        foldx_avg_ddG = {}
        results = foldx_avg_df[['Pdb', 'total energy', 'SD']].T.to_dict().values()
        for r in results:
            ident = r['Pdb'].split('_')[-1]
            ddG = r['total energy']
            ddG_sd = r['SD']

            foldx_avg_ddG[self.mutation_index_to_group[int(ident)]] = (ddG, ddG_sd)

        return foldx_avg_ddG