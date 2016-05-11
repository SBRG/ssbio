#!/usr/bin/env python

# this is to run dock.py on a folder of PDBs. it splits up the pdbs into multiple folders to make it run faster.

import os
import re
import glob
import shutil
import subprocess
import numpy as np
import time

######### MAKE SURE THESE ARE RIGHT!! ######################
dock_script = "/home/nathan/projects/scripts/DOCK/dock.py"
docking_folder = "/home/nathan/projects/DOCK/G6PD/2_MD/2_WT_2BH9/2_NADP/"
md_frames = docking_folder + 'main/'
ligand = "G6Q/G6Q_ideal.mol2"
cofactors = "NPD"
bindingsite = "171,205,239,258,360,365,395"
desc = 'G6PD_WT_NADP_G6Q'
split = 8
#############################################################


today = time.strftime("%y%m%d")
logs_folder = '/home/nathan/projects/DOCK/logs'
os.chdir(md_frames)
# get all files in the main trajectory folder
orig_pdbs = sorted(glob.glob("*.*.pdb"))

# divide total number of pdb frames by number we wanna split by
# save each list of files into a list (make a list of lists)
to_split = np.array_split(orig_pdbs,split)

# for each of these lists
for i,sub in enumerate(to_split):
    # add already generated files that may correspond with frames
    sub_list = list(sub)
    for f in sub_list:
        frame_num = f.split('.')[1]
        other_files = glob.glob('*')
        for other_file in other_files:
            match = re.search('([.]|[-])%s([.]|[_])' % frame_num, other_file)
            if match and other_file != f:
                sub_list.append(other_file)

    # make a new folder
    split_folder = md_frames + 'split_' + str(i)

    if not os.path.exists(split_folder):
        os.makedirs(split_folder)

    # move these files into that folder
    for ff in sub_list:
        shutil.move(ff, split_folder)

    # start running dock.py in the background!
    return_code = subprocess.call('nohup %s %s %s %s %s > %s/%s-%s &' % (dock_script,split_folder,ligand,cofactors,bindingsite,logs_folder,today,desc), shell=True)
    if return_code == 0:
        print 'Running DOCK on %s' % split_folder

print 'DOCK is now running in the background for %s.' % desc