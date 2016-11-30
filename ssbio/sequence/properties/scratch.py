import os.path as op
import ssbio.utils
import ssbio.sequence.fasta

import logging
log = logging.getLogger(__name__)


class SCRATCH():
    """Provide wrappers for running and parsing SCRATCH on a sequence file or sequence string.

    SCRATCH predicts:
        - Secondary structure
            - 3 classes (helix, strand, other) using SSpro
            - 8 classes (standard DSSP definitions) using SSpro8
        - Relative solvent accessibility (RSA, also known as relative accessible surface area)
            - @ 25% exposed RSA cutoff (<25% RSA means it is buried)
            - @ all cutoffs in 5% increments from 0 to 100
    """

    # TODO: also provide summary dataframes

    def __init__(self, project_name, seq_file=None, seq_str=None):
        self.project_name = project_name
        self.seq_file = seq_file
        if seq_str:
            self.seq_file = ssbio.sequence.fasta.load_seq_str_as_temp_file(seq_str)

    def run_scratch(self, path_to_scratch, num_cores=1, outname=None, outdir=None, force_rerun=False):
        """Run SCRATCH on the seq_file that was loaded into the class.

        Args:
            path_to_scratch: Path to the SCRATCH executable, run_SCRATCH-1D_predictors.sh
            outname: Prefix to name the output files
            outdir: Directory to store the output files
            force_rerun: Flag to force rerunning of SCRATCH even if the output files exist

        Returns:

        """
        if not outname:
            outname = self.project_name
        if not outdir:
            outdir = ''

        outname = op.join(outdir, outname)

        self.out_sspro = '{}.ss'.format(outname)
        self.out_sspro8 = '{}.ss8'.format(outname)
        self.out_accpro = '{}.acc'.format(outname)
        self.out_accpro20 = '{}.acc20'.format(outname)

        # TODO: check for multiple output files in command_runner
        ssbio.utils.command_runner(shell_command='{} {} {} {}'.format(path_to_scratch, self.seq_file, outname, num_cores),
                                   outfile='{}.ss'.format(outname),
                                   force_rerun_flag=force_rerun)

    def sspro_summary(self):
        """Parse the SSpro output file and return a summary of secondary structure composition.

        The output file is just a FASTA formatted file, so you can get residue level
            information by parsing it like a normal sequence file.

        Returns:
            dict: Percentage of:
                H: helix
                E: strand
                C: the rest

        """
        summary = {}

        records = ssbio.sequence.fasta.load_fasta_file(self.out_sspro)
        for r in records:
            seq_summary = {}
            seq_summary['H'] = r.seq.count('H')/float(len(r))
            seq_summary['E'] = r.seq.count('E')/float(len(r))
            seq_summary['C'] = r.seq.count('C')/float(len(r))

            summary[r.id] = seq_summary

        return summary

    def sspro8_summary(self):
        """Parse the SSpro8 output file and return a summary of secondary structure composition.

        The output file is just a FASTA formatted file, so you can get residue level
            information by parsing it like a normal sequence file.

        Returns:
            dict: Percentage of:
                H: alpha-helix
                G: 310-helix
                I: pi-helix (extremely rare)
                E: extended strand
                B: beta-bridge
                T: turn
                S: bend
                C: the rest

        """
        summary = {}

        records = ssbio.sequence.fasta.load_fasta_file(self.out_sspro8)
        for r in records:
            seq_summary = {}
            seq_summary['H'] = r.seq.count('H') / float(len(r))
            seq_summary['G'] = r.seq.count('G') / float(len(r))
            seq_summary['I'] = r.seq.count('I') / float(len(r))
            seq_summary['E'] = r.seq.count('E') / float(len(r))
            seq_summary['B'] = r.seq.count('B') / float(len(r))
            seq_summary['T'] = r.seq.count('T') / float(len(r))
            seq_summary['S'] = r.seq.count('S') / float(len(r))
            seq_summary['C'] = r.seq.count('C') / float(len(r))

            summary[r.id] = seq_summary

        return summary

    def accpro_summary(self):
        """Parse the ACCpro output file and return a summary of percent exposed/buried residues.

        The output file is just a FASTA formatted file, so you can get residue level
            information by parsing it like a normal sequence file.

        Returns:
            dict: Percentage of buried and exposed residues

        """
        summary = {}

        records = ssbio.sequence.fasta.load_fasta_file(self.out_accpro)
        for r in records:
            seq_summary = {}
            seq_summary['exposed'] = r.seq.count('e') / float(len(r))
            seq_summary['buried'] = r.seq.count('-') / float(len(r))

            summary[r.id] = seq_summary

        return summary

    def accpro20_summary(self, cutoff):
        """Parse the ACCpro output file and return a summary of percent exposed/buried residues based on a cutoff.

        The output file is just a FASTA formatted file, so you can get residue level
            information by parsing it like a normal sequence file.

        Args:
            cutoff (float): Cutoff for defining a buried or exposed residue.

        Returns:
            dict: Percentage of buried and exposed residues

        """
        summary = {}

        if cutoff < 1:
            cutoff = 1 * 100

        records = read_accpro20(self.out_accpro20)

        for k,v in records.items():
            seq_summary = {}

            seq_acc = [int(x) for x in v.split(' ')]

            exposed = 0
            buried = 0
            for s in seq_acc:
                if s > cutoff:
                    exposed += 1
                else:
                    buried += 1

            seq_summary['exposed'] = exposed / float(len(seq_acc))
            seq_summary['buried'] = buried / float(len(seq_acc))

            summary[k] = seq_summary

        return summary


def read_accpro20(infile):
    """Read the accpro20 output (.acc20) and return the parsed FASTA records.

    Keeps the spaces between the accessibility numbers.

    Args:
        infile: Path to .acc20 file

    Returns:
        dict: Dictionary of accessibilities with keys as the ID

    """
    with open(infile) as f:
        records = f.read().splitlines()

    accpro20_dict = {}
    for i, r in enumerate(records):
        if i % 2 == 0:
            accpro20_dict[records[i][1:]] = records[i + 1]

    return accpro20_dict