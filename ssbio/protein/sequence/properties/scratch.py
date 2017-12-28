import logging
import os.path as op

import ssbio.protein.sequence.utils.fasta
import ssbio.utils

log = logging.getLogger(__name__)


class SCRATCH():
    """Provide wrappers for running and parsing SCRATCH on a sequence file or sequence string.

    To run from the command line:

        ./run_SCRATCH-1D_predictors.sh  input_fasta  output_prefix  [num_threads]

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
            self.seq_file = ssbio.protein.sequence.utils.fasta.write_seq_as_temp_fasta(seq_str)

    def run_scratch(self, path_to_scratch, num_cores=1, outname=None, outdir=None, force_rerun=False):
        """Run SCRATCH on the sequence_file that was loaded into the class.

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

    def sspro_results(self):
        """Parse the SSpro output file and return a dict of secondary structure compositions.

        Returns:
            dict: Keys are sequence IDs, values are the lists of secondary structure predictions.
                H: helix
                E: strand
                C: the rest

        """
        return ssbio.protein.sequence.utils.fasta.load_fasta_file_as_dict_of_seqs(self.out_sspro)

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

        records = ssbio.protein.sequence.utils.fasta.load_fasta_file(self.out_sspro)
        for r in records:
            seq_summary = {}
            seq_summary['percent_H-sspro'] = r.seq.count('H')/float(len(r))
            seq_summary['percent_E-sspro'] = r.seq.count('E')/float(len(r))
            seq_summary['percent_C-sspro'] = r.seq.count('C')/float(len(r))

            summary[r.id] = seq_summary

        return summary

    def sspro8_results(self):
        """Parse the SSpro8 output file and return a dict of secondary structure compositions.
        """
        return ssbio.protein.sequence.utils.fasta.load_fasta_file_as_dict_of_seqs(self.out_sspro8)

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

        records = ssbio.protein.sequence.utils.fasta.load_fasta_file(self.out_sspro8)
        for r in records:
            seq_summary = {}
            seq_summary['percent_H-sspro8'] = r.seq.count('H') / float(len(r))
            seq_summary['percent_G-sspro8'] = r.seq.count('G') / float(len(r))
            seq_summary['percent_I-sspro8'] = r.seq.count('I') / float(len(r))
            seq_summary['percent_E-sspro8'] = r.seq.count('E') / float(len(r))
            seq_summary['percent_B-sspro8'] = r.seq.count('B') / float(len(r))
            seq_summary['percent_T-sspro8'] = r.seq.count('T') / float(len(r))
            seq_summary['percent_S-sspro8'] = r.seq.count('S') / float(len(r))
            seq_summary['percent_C-sspro8'] = r.seq.count('C') / float(len(r))

            summary[r.id] = seq_summary

        return summary

    def accpro_results(self):
        """Parse the ACCpro output file and return a dict of secondary structure compositions.
        """
        return ssbio.protein.sequence.utils.fasta.load_fasta_file_as_dict_of_seqs(self.out_accpro)

    def accpro_summary(self):
        """Parse the ACCpro output file and return a summary of percent exposed/buried residues.

        The output file is just a FASTA formatted file, so you can get residue level
            information by parsing it like a normal sequence file.

        Returns:
            dict: Percentage of buried and exposed residues

        """
        summary = {}

        records = ssbio.protein.sequence.utils.fasta.load_fasta_file(self.out_accpro)
        for r in records:
            seq_summary = {}
            seq_summary['percent_exposed-accpro'] = r.seq.count('e') / float(len(r))
            seq_summary['percent_buried-accpro'] = r.seq.count('-') / float(len(r))

            summary[r.id] = seq_summary

        return summary

    def accpro20_results(self):
        """Parse the ACCpro output file and return a dict of secondary structure compositions"""
        return read_accpro20(self.out_accpro20)
        # return ssbio.sequence.utils.fasta.load_fasta_file_as_dict_of_seqs(self.out_accpro20)

    def accpro20_summary(self, cutoff):
        """Parse the ACCpro output file and return a summary of percent exposed/buried residues based on a cutoff.

        Below the cutoff = buried
        Equal to or greater than cutoff = exposed
        The default cutoff used in accpro is 25%.

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

            exposed = 0
            buried = 0
            for s in v:
                if s > cutoff:
                    exposed += 1
                else:
                    buried += 1

            seq_summary['percent_exposed-accpro20'] = exposed / float(len(v))
            seq_summary['percent_buried-accpro20'] = buried / float(len(v))

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
            # TODO: Double check how to parse FASTA IDs (can they have a space because that is what i split by)
            # Key was originally records[i][1:]
            accpro20_dict[records[i].split(' ')[0][1:]] = [int(x) for x in records[i + 1].split(' ')]

    return accpro20_dict