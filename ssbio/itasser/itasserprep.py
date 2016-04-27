import os
import getpass
from seqtools import SeqTools
from seqwriter import SeqWriter

st = SeqTools()


class ITASSERPrep():
    '''
    Prepares sequences for ITASSER runs.
    '''

    def __init__(self, ident, root_dir):
        self.ident = ident
        self.root_dir = root_dir
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)

    def prep_folders(self, seq):
        '''
        Takes in a Biopython SeqRecord object,
        prepares the folders for it to run ITASSER
        '''
        itasser_dir = os.path.join(self.root_dir, self.ident)

        if not os.path.exists(itasser_dir):
            os.makedirs(itasser_dir)

        sw.write_fasta_file(seq_records=seq,
                            ident='seq',
                            extension='fasta',
                            outpath=itasser_dir)
        return itasser_dir

    def prep_script_local(self, itasser_loc, itlib_loc, datadir,
                          java_home='/usr/', light='true', print_exec=True):
        script_file = '{}.sh'.format(self.ident)
        outfile = os.path.join(self.root_dir, script_file)

        itasser = {'executable': os.path.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
                   'pkgdir': itasser_loc,
                   'libdir': itlib_loc,
                   'seqname': self.ident,
                   'datadir': datadir,
                   'usrname': getpass.getuser(),
                   'java_home': java_home,
                   'light': light}

        script = open(outfile, 'w')
        script.write('#!/bin/bash -l\n')
        script.write("{i[executable]} -pkgdir {i[pkgdir]} -libdir {i[libdir]} -seqname {i[seqname]} -datadir {i[datadir]} -usrname {i[usrname]} -java_home {i[java_home]} -light {i[light]}\n\n".format(i=itasser))
        script.close()

        if print_exec:
            print('nohup ./{} > {}.out &'.format(script_file, self.ident))

        return outfile



if __name__ == '__main__':
    import glob
    sw = SeqWriter()

    TODO: make this an executable script to
    1) ask for global I-TASSER locations
    2) ask for working directory
    3) take in multiple inputs and prepare them for I-TASSER runs
        a) input types
            i) a single FASTA file with single or multiple sequences
            ii) multiple FASTA files contained in the working directory
            iii) a dataframe with IDs and sequences
            iv) a sequence string and an ID (and optional additional identifiers)
        b) types of runs
            i) NERSC slurm (sbatch) inputs
            ii) local torque (qsub) inputs
            iii) simple executable background scripts
    4) Output executable scripts or submit things to the queue

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

    ip = ITASSERPrep(ident='W5EP13', root_dir='/home/nathan/projects/GEM-PRO/cyano/')

    sequence = sl.seq_loader('/home/nathan/Downloads/W5EP13.faa', is_file=True)
    data_dir = ip.prep_folders(sequence)
    ip.prep_script_local(itasser_loc='/home/nathan/software/I-TASSER4.4',
                         itlib_loc='/home/nathan/software/ITLIB',
                         datadir=data_dir)
