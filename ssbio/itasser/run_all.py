


import os
import shutil
import subprocess

thedir = '.'
folders = [name for name in os.listdir(
    thedir) if os.path.isdir(os.path.join(thedir, name))]
folders = sorted(folders, reverse=True)
for_ssb3 = folders[:len(folders) / 2]

for fo in for_ssb3:
    coach = open('%s_coach.sh' % fo, 'w')

    coach.write('#!/bin/bash\n')
    coach.write('#PBS -l walltime=05:20:00\n')
    coach.write('#PBS -q regular\n')
    coach.write('#PBS -N %s\n' % fo)
    coach.write('perl ~/software/I-TASSER4.4/I-TASSERmod/runCOACH.pl -pkgdir /home/nathan/software/I-TASSER4.4 -libdir /home/nathan/software/ITLIB -protname %s -model model1.pdb -datadir /home/nathan/projects/GEM-PRO/yome/all_test/%s -GO true\n\n' % (fo, fo))

    coach.close()

    # subprocess.call('qsub %s_coach.sh;' % (fo), shell=True)
    print('qsub %s_coach.sh;' % (fo)),
