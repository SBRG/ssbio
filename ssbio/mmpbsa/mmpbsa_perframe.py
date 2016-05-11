#! /usr/bin/env python

import os
import jinja2
import time
import yaml
import subprocess
import glob


# IMPORTANT: input the directories of your stuff here
# ROOT: location of all files
ROOT = '/home/nathan/projects'
# MMPBSA_FOLDER: this is the "main" MMPBSA folder
MMPBSA_FOLDER = os.path.join(ROOT, 'MMPBSA')
# IMPORTANT: this is the template for making an input script for mmpbsa running per frame
MMPBSA_PERFRAME_TEMPLATE = os.path.join(ROOT, 'templates', 'mmpbsa_perframe.jinja')
# IMPORTANT: this is where the template for the script calling MMPBSA.py is located
MMPBSA_RUNNER_TEMPLATE = os.path.join(ROOT, 'templates', 'mmpbsa_run.jinja')
# IMPORTANT: here are the ligand parameters
PARAMETERS = os.path.join(ROOT, 'parameters')
LIGANDS_PATH = os.path.join(PARAMETERS, 'ligands')


# input your proteins, modifications, types, ligands, and cluster choices here
# EXAMPLE:
# protein - G6PD:
#     modification - 1_WT_2BH9:
#         type - 2_NADP:
#             location of mmpbsa receptor prep files (receptor in gas) - receptor_prep: /home/nathan/projects/MMPBSA/G6PD/1_WT_2BH9/2_NADP/0_receptor_prep
#             name of receptor in gas - receptor_prmtop_file: receptor_gas.prmtop
#             the number of residues + your ligand in the complex  -- num_residues: 492
#             the pdb id of your protein - basename: 2BH9
#             indented list of ligands - LIGANDS:
#                 ligand1 - BG6:
#                     name of the folder which contains your ligand parameters in LIGANDS_PATH - param_loc: BG6_minus2
#                     indented list of all clusters and their frames from the MD trajectory - CLUSTERS:
#                         '111': '298'
#                         '110': '199'
#                         '100': '045'
#                         '001': '448'
#                         '000': '010'
#                     indented list of the trajectory you want to use for MMPBSA and the number of frames - TRAJ:
#                         '111':
#                             - '2BH9_WT_BG6_C111_0ns_11ns_weeded_step_100.traj'
#                             - 70
#                         #'110': ''
#                         '100':
#                             - 'G6PD_BG6_minus2_100_1ns_7ns_weeded_step_100.traj'
#                             - 42
#                         '001':
#                             - 'G6PD_BG6_minus2_001_1ns_10ns_weeded_step_100.traj'
#                             - 66
#                         '000':
#                             - '2BH9_WT_BG6_C000_0ns_11ns_weeded_step_100.traj'
#                             - 70
#                  ligand2 - XXX:
#                      .....

mmpbsa_prep = yaml.load("""
---
G6PD:
    1_WT_2BH9:
        2_NADP:
            receptor_prep: /home/nathan/projects/MMPBSA/G6PD/1_WT_2BH9/2_NADP/0_receptor_prep
            receptor_prmtop_file: receptor_gas.prmtop
            num_residues: 492
            basename: 2BH9
            LIGANDS:
                BG6:
                    param_loc: BG6_minus2
                    CLUSTERS:
                        '111': '298'
                        '110': '199'
                        '100': '045'
                        '001': '448'
                        '000': '010'
                    TRAJ:
                        '111':
                            - '2BH9_WT_BG6_C111_0ns_11ns_weeded_step_100.traj'
                            - 70
                        #'110': ''
                        '100':
                            - 'G6PD_BG6_minus2_100_1ns_7ns_weeded_step_100.traj'
                            - 42
                        '001':
                            - 'G6PD_BG6_minus2_001_1ns_10ns_weeded_step_100.traj'
                            - 66
                        '000':
                            - '2BH9_WT_BG6_C000_0ns_11ns_weeded_step_100.traj'
                            - 70

    2_SNP_2BH9:
        2_NADP:
            receptor_prep: /home/nathan/projects/MMPBSA/G6PD/2_SNP_2BH9/2_NADP/0_receptor_prep
            receptor_prmtop_file: receptor_gas.prmtop
            num_residues: 492
            basename: 2BH9
            LIGANDS:
                BG6:
                    param_loc: BG6_minus2
                    CLUSTERS:
                        '111': '038'
                        '011': '643'
                        '010': '100'
                        '001': '056'
                        '000': '170'
                    TRAJ:
                        '111':
                            - '2BH9_SNP_BG6_C111_0ns_11ns_weeded_step_100.traj'
                            - 70
                        '011':
                            - 'G6PD_BG6_minus2_011_1ns_10ns_weeded_step_100.traj'
                            - 66
                        '010':
                            - 'G6PD_BG6_minus2_010_1ns_10ns_weeded_step_100.traj'
                            - 66
                        #'001': ''
                        '000':
                            - '2BH9_SNP_BG6_C000_0ns_11ns_weeded_step_100.traj'
                            - 70
""")


def mmpbsa_maker(template_file, mmpbsa_dict):
    '''
    Creates MMPBSA input script (to specify gb, pb, number of frames, other options)
    Requires template file: 'mmpbsa_perframe.jinja'
    Returns filename of input script
    '''
    # allows loading of templates off filesystems
    templateLoader = jinja2.FileSystemLoader( searchpath="/" )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # load the template
    template = templateEnv.get_template( template_file )

    # process the template to produce our final text.
    outputText = template.render( mmpbsa_dict )

    file_out = mmpbsa_dict['filename']
    with open(file_out, "w") as f:
        f.write(outputText)

    return file_out


def mmpbsa_runner(template, mmpbsa_run_dict):
    '''
    Creates script to call MMPBSA.py and points to all input files (receptor, complex, ligand prmtops)
    Requires template file: 'mmpbsa_run.jinja'
    Returns filename of input script
    '''
    # allows loading of templates off filesystems
    templateLoader = jinja2.FileSystemLoader( searchpath="/" )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # load the template
    MMPBSA_TEMPLATE = template
    template = templateEnv.get_template( MMPBSA_TEMPLATE )

    # process the template to produce our final text.
    outputText = template.render( mmpbsa_run_dict )

    file_out = mmpbsa_run_dict['date'] + '-' + mmpbsa_run_dict['mmpbsa_input_params'] + '.in'
    with open(file_out, "w") as f:
        f.write(outputText)

    return file_out


# options to use for running MMPBSA perframe (leave this the same)
pb_perframe = yaml.load("""
---
folder: 'pb_perframe'
filename: mmpbsa_pb_perframe.in
gen_options:
    - keep_files=0
    - verbose=2
pb_options:
    - fillratio=4.0
    - inp=1
    - radiopt=0
""")


if __name__ == '__main__':

    TODAY = time.strftime("%y%m%d")

    # a lot of this is to just first grab all the names
    for PROTEIN,mods in mmpbsa_prep.iteritems():

        for MOD,types in mods.iteritems():
            # IMPORTANT: double check if this folder is correct
            f1 = os.path.join(MMPBSA_FOLDER, PROTEIN, MOD)

            for TYPE,info in types.iteritems():
                # IMPORTANT: double check if this folder is correct
                f2 = os.path.join(f1, TYPE)

                basename = info['basename']
                num_residues = info['num_residues']
                ligands = info['LIGANDS']
                receptor_prep = info['receptor_prep']
                receptor_prmtop_file = info['receptor_prmtop_file']

                for LIGAND,lig_info in ligands.iteritems():
                    # IMPORTANT: double check if this folder is correct
                    f3 = os.path.join(f2, LIGAND)

                    clusters = lig_info['CLUSTERS']
                    for_mmpbsa_trajs = lig_info['TRAJ']
                    param_loc = lig_info['param_loc']

                    for CLUSTER,FRAME in clusters.iteritems():
                        main_folder = os.path.join(f3, CLUSTER)

                        # IMPORTANT: double check if these folders are correct
                        complex_prep = os.path.join(main_folder, '0_complex_prep')
                        md = os.path.join(main_folder, '1_md')
                        mmpbsa = os.path.join(main_folder, '2_mmpbsa')

                        os.chdir(mmpbsa)

                        # only run MMPBSA if there is a trajectory to work with!
                        if for_mmpbsa_trajs:

                            if CLUSTER in for_mmpbsa_trajs.keys():
                                traj_to_run_mmpbsa_on = for_mmpbsa_trajs[CLUSTER][0]
                                numframes = for_mmpbsa_trajs[CLUSTER][1]

                                mmpbsa_dict = pb_perframe
                                mmpbsa_param = mmpbsa_dict['folder']

                                # run under a new folder under the main MMPBSA folder ("pb_perframe")
                                perframe_folder = os.path.join(mmpbsa, mmpbsa_param)
                                if not os.path.exists(perframe_folder):
                                    os.mkdir(perframe_folder)

                                # and for each of the frames, also make a folder for the input script and results
                                for f in range(1, numframes):
                                    each_frame_folder = os.path.join(perframe_folder, '%.3d' % f)

                                    if not os.path.exists(each_frame_folder):
                                        os.mkdir(each_frame_folder)

                                    print 'CREATING:', each_frame_folder

                                    # make the input file for that parameter
                                    # specify startframe and endframe to be 1,1, 2,2, 3,3, etc (one by one!)
                                    mmpbsa_dict['start_frame'] = f
                                    mmpbsa_dict['end_frame'] = f

                                    # the input file is made and deleted in /tmp/
                                    os.chdir('/tmp/')
                                    mmpbsa_input = mmpbsa_maker(MMPBSA_PERFRAME_TEMPLATE, mmpbsa_dict)

                                    print mmpbsa_input

                                    names = {'protein':PROTEIN,
                                             'mod':MOD,
                                             'type':TYPE,
                                             'ligand':LIGAND,
                                             'cluster':CLUSTER,
                                             'mmpbsa_param':mmpbsa_param,
                                             'receptor_prep_folder':receptor_prep,
                                             'complex_prep_folder':complex_prep,
                                             'ligands_path':LIGANDS_PATH,
                                             'md_run_folder':md,
                                             'param_loc':param_loc,
                                             'num_residues':num_residues,
                                             'for_mmpbsa':traj_to_run_mmpbsa_on}

                                    # finally make the script to call MMPBSA.py (calls it using nohup)
                                    # IMPORTANT: make sure the names point to the correct files - especially prmtop paths!
                                    mmpbsa_run_template = """
---
input_scripts_dir: /tmp
mmpbsa_input_params: mmpbsa_{mmpbsa_param}
complex_solv_prmtop_path: {complex_prep_folder}/{protein}_{ligand}_{cluster}.prmtop
complex_gas_prmtop_path: {complex_prep_folder}/{protein}_{ligand}_{cluster}_GAS.prmtop
receptor_gas_prmtop_path: {receptor_prep_folder}/receptor_gas.prmtop
ligand_gas_prmtop_path: {ligands_path}/{param_loc}/{ligand}.prmtop
trajectory_path: {md_run_folder}/{for_mmpbsa}
"""
                                    mmpbsa_run_template = mmpbsa_run_template.format(**names)
                                    mmpbsa_run_yaml = yaml.load(mmpbsa_run_template)
                                    mmpbsa_run_yaml['date'] = TODAY

                                    os.chdir(each_frame_folder)
                                    files = glob.glob('*')
                                    if 'FINAL_RESULTS_mmpbsa_%s.dat' % mmpbsa_param in files:
                                        print 'RESULTS PREVIOUSLY CALCULATED IN:', each_frame_folder
                                        continue
                                    mmpbsa_runner_file = mmpbsa_runner(MMPBSA_RUNNER_TEMPLATE, mmpbsa_run_yaml)

                                    print 'RUNNING MMPBSA IN:', each_frame_folder
                                    os.chmod(mmpbsa_runner_file, 0755)
                                    retva = subprocess.call('./%s' % mmpbsa_runner_file, shell=True)
                                    print retva

                                    # IMPORTANT: adjust the sleep time (in seconds) so you don't kill the computer
                                    # 5 to 10 recommended to be safe, pb should take around 15 seconds? to run (test it on a frame first!)
                                    time.sleep(10)