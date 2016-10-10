import os
import os.path as op

import ssbio.databases.ncbi

from ssbio import utils
date = utils.Date()

class ATLAS():
    """Class to represent an ATLAS workflow to carry out multi-strain comparisons

    Main steps are:
    1. Strain-specific model construction based on orthologous genes & systems modeling
    2. Phylogenetic analysis to pick out important genes
    3. GEM-PRO of the "base strain"
    4. Structure property calculation & integrated structural systems analysis

    Each step may generate a report and also request additional files if something is missing
    """

    def __init__(self, base_strain_name, root_dir):
        self.root_dir = root_dir

        # atlas_dir - directory where all ATLAS analysis will be carried out
        # TODO: should _ATLAS be appended to the name?
        self.atlas_dir = op.join(root_dir, base_strain_name + '_ATLAS')

        # model_files - directory where base strain GEM and new models will be stored
        self.model_files = op.join(self.atlas_dir, 'model_files')

    def prep_folders(self):
        """Prepare all folders for an ATLAS project
        """
        # data - directory where all data will be stored
        self.data = op.join(self.atlas_dir, 'data')

        # notebooks - directory where ipython notebooks will be stored for manual analyses
        self.notebooks = op.join(self.atlas_dir, 'notebooks')

        # # figures - directory where all figures will be stored
        self.figures = op.join(self.atlas_dir, 'figures')

        # seq_files - sequence related files are stored here
        self.seq_files = op.join(self.atlas_dir, 'sequence_files')

        for directory in [self.atlas_dir, self.data, self.notebooks, self.figures,
                          self.model_files, self.seq_files]:
            if not op.exists(directory):
                os.makedirs(directory)
                print('[PREP] Created directory: {}'.format(directory))
            else:
                print('[PREP] Directory already exists: {}'.format(directory))

    def input_genome_ids(self, infile):
        """Loads a line-delimited file of NCBI GIs or RefSeq complete genome IDs for the ATLAS pipeline.

        Args:
            infile: path to file containing a list of genome identifiers, separated by newlines

        """
        with open(infile) as f:
            lines = f.readlines()

    def download_all_protein_sequences(genome_accession_or_id, email, outdir='', outfile=''):
        """Download the entire set of amino acid sequences from protein-encoding genes in a genome from NCBI.

        Saves a FASTA file in the optional directory specified.

        Args:
            genome_accession_or_id: RefSeq complete genome ID (e.g. "NC_000913") or NCBI GI (e.g. "556503834")
            email: mandatory email so NCBI knows who is accessing the data
            outdir: optional output directory (default is the current directory)

        Returns:
            Path to downloaded FASTA file.

        """

        # path and filename parsing
        if outfile:
            outfile = op.join(outdir, '{}.faa'.format(outfile))
        else:
            # if no outfile is specified, default is "$DATE-$GI.faa"
            outfile = op.join(outdir, '{}-{}.faa'.format(date.short_date, genome_accession_or_id))

        # return the path to the file if it was already downloaded
        if op.exists(outfile):
            print('[INFO] Genome FASTA file already exists at {}'.format(outfile))
            return outfile

        Entrez.email = email

        # convert to NCBI GI
        convert_to_gi = Entrez.read(Entrez.esearch(db="nucleotide", term=genome_accession_or_id, retmode="xml"))
        gi = convert_to_gi['IdList']

        # download only the CDS protein FASTA (https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/)
        record = Entrez.efetch(db="nucleotide", id=gi, rettype="fasta_cds_aa", retmode="fasta")

        # save the FASTA file
        with open(outfile, 'wt') as f:
            f.write(record.read())

        print('[INFO] Saved FASTA file at {}'.format(outfile))
        return outfile