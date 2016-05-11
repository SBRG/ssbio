import getpass
import os
import os.path as op

from ssbio.sequence import fasta as fasta


class ITASSERPrep():
    '''
    Prepares a sequence for ITASSER runs.
    '''

    def __init__(self, ident, seq_str,
                 root_dir, itasser_loc, itlib_loc, data_dir='', light='true', runtype='local', print_exec=False,
                 binding_site_pred=False, ec_pred=False, go_pred=False,
                 slurm_email='', slurm_username='', slurm_walltime='48:00:00',
                 slurm_queue='shared'):
        self.ident = ident
        self.seq_str = seq_str

        self.root_dir = root_dir
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        if not data_dir:
            self.data_dir = self.prep_folders(seq_str)
        elif data_dir:
            orig_data_dir = self.prep_folders(seq_str)
            self.data_dir = op.join(data_dir, op.basename(orig_data_dir))

        self.print_exec = print_exec
        self.runtype = runtype
        self.light = light

        additional_options = ''
        if binding_site_pred:
            additional_options += '-LBS true '
        if ec_pred:
            additional_options += '-EC true '
        if go_pred:
            additional_options += '-GO true '
        self.additional_options = additional_options

        if runtype == 'local' or runtype == 'torque':
            self.prep_script_local(itasser_loc=itasser_loc,
                                   itlib_loc=itlib_loc)

        if runtype == 'slurm':
            self.slurm_email = slurm_email
            self.slurm_username = slurm_username
            self.slurm_walltime = slurm_walltime
            self.slurm_queue = slurm_queue

            self.prep_script_slurm(itasser_loc=itasser_loc,
                                   itlib_loc=itlib_loc)

    def prep_folders(self, seq):
        """Take in a sequence string and prepares the folder for it to run ITASSER
        """
        itasser_dir = op.join(self.root_dir, self.ident)

        if not op.exists(itasser_dir):
            os.makedirs(itasser_dir)

        fasta.write_fasta_file(seq_str=seq,
                               ident='seq',
                               extension='fasta',
                               outpath=itasser_dir)
        return itasser_dir

    def prep_script_local(self, itasser_loc, itlib_loc, java_home='/usr/'):
        script_file = '{}.sh'.format(self.ident)
        outfile = os.path.join(self.root_dir, script_file)

        itasser = {'executable': op.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
                   'pkgdir': itasser_loc,
                   'libdir': itlib_loc,
                   'seqname': self.ident,
                   'datadir': self.data_dir,
                   'usrname': getpass.getuser(),
                   'java_home': java_home,
                   'additional_options': self.additional_options,
                   'light': self.light}

        script = open(outfile, 'w')
        script.write('#!/bin/bash -l\n')

        if self.runtype == 'torque':
            script.write('#PBS -l walltime=12:00:00\n')
            script.write('#PBS -q regular\n')
            script.write('#PBS -N {i[seqname]}\n'.format(i=itasser))
            script.write('#PBS -o {i[seqname]}.out\n'.format(i=itasser))
            script.write('#PBS -e {i[seqname]}.err\n'.format(i=itasser))

        script.write(("{i[executable]} "
                      "-pkgdir {i[pkgdir]} "
                      "-libdir {i[libdir]} "
                      "-seqname {i[seqname]} "
                      "-datadir {i[datadir]} "
                      "-usrname {i[usrname]} "
                      "-java_home {i[java_home]} "
                      "{i[additional_options]}"
                      "-light {i[light]}\n\n").format(i=itasser))
        script.close()

        os.chmod(outfile, 0o755)

        if self.print_exec and self.runtype=='local':
            print('nohup ./{} > {}.out &'.format(op.basename(outfile), os.path.join(self.root_dir, self.ident)),
                  end='\n\n')

        if self.print_exec and self.runtype == 'torque':
            print('qsub {}'.format(op.basename(outfile), os.path.join(self.root_dir, self.ident)),
                  end='; ')

        return outfile

    def prep_script_slurm(self,
                          itasser_loc,
                          itlib_loc,
                          java_home='${JAVA_HOME}'
                          ):
        script_file = '{}.slm'.format(self.ident)
        outfile = os.path.join(self.root_dir, script_file)

        itasser = {'executable': op.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
                   'pkgdir': itasser_loc,
                   'libdir': itlib_loc,
                   'seqname': self.ident,
                   'datadir': self.data_dir,
                   'usrname': self.slurm_username,
                   'java_home': java_home,
                   'light': self.light,
                   'additional_options': self.additional_options,
                   'slurm_walltime': self.slurm_walltime,
                   'slurm_queue': self.slurm_queue,
                   'slurm_email': self.slurm_email}

        slurm = open(outfile, 'w')

        slurm.write('#!/bin/bash -l\n')
        slurm.write('#SBATCH -p shared\n')
        slurm.write('#SBATCH -t {i[slurm_walltime]}\n'.format(i=itasser))
        slurm.write('#SBATCH --mem=8GB\n')
        slurm.write('#SBATCH -J {i[seqname]}\n'.format(i=itasser))
        slurm.write('#SBATCH --mail-user {i[slurm_email]}\n'.format(i=itasser))
        slurm.write('#SBATCH -o {i[seqname]}.out\n'.format(i=itasser))
        slurm.write('#SBATCH -e {i[seqname]}.err\n'.format(i=itasser))
        slurm.write('#SBATCH --mail-type=BEGIN\n')
        slurm.write('#SBATCH --mail-type=END\n')
        slurm.write(('{i[executable]} '
                     '-pkgdir {i[pkgdir]} '
                     '-libdir {i[libdir]} '
                     '-seqname {i[seqname]} '
                     '-datadir {i[datadir]} '
                     '-usrname {i[usrname]} '
                     '-java_home {i[java_home]} '
                     '{i[additional_options]}'
                     '-light {i[light]}\n\n').format(i=itasser))

        slurm.close()

        os.chmod(outfile, 0o755)

        if self.print_exec:
            print('sbatch {}'.format(op.basename(outfile)), end='; ')

        return outfile


if __name__ == '__main__':
    pass

    # TODO: make this an executable script to
    # 1) ask for global I-TASSER locations
    # 2) ask for working directory
    # 3) take in multiple inputs and prepare them for I-TASSER runs
    #     a) input types
    #         i) a single FASTA file with single or multiple sequences
    #         ii) multiple FASTA files contained in the working directory
    #         iii) a dataframe with IDs and sequences
    #         iv) a sequence string and an ID (and optional additional identifiers)
    #     b) types of runs
    #         i) NERSC slurm (sbatch) inputs
    #         ii) local torque (qsub) inputs
    #         iii) simple executable background scripts
    # 4) Output executable scripts or submit things to the queue

    # root = '/home/nathan/projects/GEM-PRO/cyano/'
    # files = glob.glob(os.path.join(root,'*.faa'))
    # for f in files:
    #     identifier = os.path.splitext(os.path.basename(f))[0]
    #     ip = ITASSERPrep(ident=identifier, root_dir='/home/nathan/projects/GEM-PRO/cyano')
    #
    #     sequence = sl.seq_loader(f, is_file=True)
    #     data_dir = ip.prep_folders(sequence)
    #     ip.prep_script_local(itasser_loc='/home/nathan/software/I-TASSER4.4',
    #                          itlib_loc='/home/nathan/software/ITLIB',
    #                          datadir=data_dir)

    # ip = ITASSERPrep(ident='W5EP13', root_dir='/home/nathan/projects/GEM-PRO/cyano/')
    #
    # sequence = sl.seq_loader('/home/nathan/Downloads/W5EP13.faa', is_file=True)
    # data_dir = ip.prep_folders(sequence)
    # ip.prep_script_local(itasser_loc='/home/nathan/software/I-TASSER4.4',
    #                      itlib_loc='/home/nathan/software/ITLIB',
    #                      datadir=data_dir)
