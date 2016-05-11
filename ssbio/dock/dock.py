#! /usr/bin/env python

import os
import re
import sys

from pkg_resources import resource_filename
amb = resource_filename('ssbio.etc', 'vdw_AMBER_parm99.defn')
f1 = resource_filename('ssbio.etc', 'flex.defn')
f2 = resource_filename('ssbio.etc', 'flex_drive.tbl')

def dockprep(file_name, prefix):
    outfile = '%s_prep.mol2' % prefix

    prep = open("prep.py", "w")
    prep.write('import chimera\n')
    prep.write('from DockPrep import prep\n')
    prep.write('models = chimera.openModels.list(modelTypes=[chimera.Molecule])\n')
    prep.write('prep(models)\n')
    prep.write('from WriteMol2 import writeMol2\n')
    prep.write('writeMol2(models, "%s")\n' % outfile)
    prep.close()

    cmd = 'chimera --nogui %s prep.py' % file_name
    os.system(cmd)
    os.remove('prep.py')
    os.remove('prep.pyc')

    return outfile

def protein_only_and_noH(file_name, keep_these, prefix):
    receptor_mol2 = '%s_receptor.mol2' % prefix
    receptor_noh = '%s_receptor_noH.pdb' % prefix
    prot_only = open("prly.com", "w")
    prot_only.write('open %s\n' % file_name)

    keep_str = 'delete ~protein'
    for res in keep_these:
        keep_str += ' & ~:%s ' % res
    keep_str = keep_str.strip() + '\n'
    prot_only.write(keep_str)

    prot_only.write('write format mol2 0 %s\n' % receptor_mol2)
    prot_only.write('delete element.H\n')
    prot_only.write('write format pdb 0 %s\n' % receptor_noh)
    prot_only.close()

    cmd = 'chimera --nogui prly.com'
    os.system(cmd)
    os.remove('prly.com')

    return receptor_mol2, receptor_noh

def dms_maker(file_name, prefix):
    dms = '%s_receptor.dms' % prefix

    cmd = 'dms %s -n -w 1.4 -o %s' % (file_name, dms)
    os.system(cmd)

    return dms

def sphgen(file_name, prefix):
    sph = "%s_receptor.sph" % prefix

    in_sph = open("INSPH", "w")
    in_sph.write("%s\n" % file_name)
    in_sph.write("R\n")
    in_sph.write("X\n")
    in_sph.write("0.0\n")
    in_sph.write("4.0\n")
    in_sph.write("1.4\n")
    in_sph.write("%s\n" % sph)
    in_sph.close()

    cmd = "sphgen_cpp"
    os.system(cmd)
    os.remove('INSPH')

    return sph

def binding_site_mol2(file_name, residues, prefix_x):
    '''
    This function will take in a .pdb file (preferably the _recpetor_noH.pdb file)
    and a string of residues (ex: 144,170,199) and delete all other residues in the
    .pdb file. It then saves the coordinates of the selected residues as a .mol2 file.
    This is necessary for Chimera to select spheres within the radius of the binding
    site.
    '''

    prefix = prefix_x + '_' + 'binding_residues'

    mol2_maker = open('%s_make_mol2.py' % prefix, 'w')
    mol2_maker.write('#! /usr/bin/env python\n')
    mol2_maker.write('from chimera import runCommand\n')
    mol2_maker.write('runCommand("open %s")\n' % file_name)
    mol2_maker.write('runCommand("delete ~:%s")\n' % residues)
    mol2_maker.write('runCommand("write format mol2 resnum 0 %s.mol2")\n' % prefix)
    mol2_maker.write('runCommand("close all")')
    mol2_maker.close()

    cmd = "chimera --nogui %s_make_mol2.py" % prefix
    os.system(cmd)
    os.remove("%s_make_mol2.py" % prefix)
    os.remove("%s_make_mol2.pyc" % prefix)

    return '%s.mol2' % prefix

def sphere_selector_using_residues(sph_file, res_mol2, radius, prefix):
    selsph = '%s_selected_spheres_using_binding_residues.sph' % prefix

    cmd = "sphere_selector %s %s %d" % (sph_file, res_mol2, radius)
    rename = "mv selected_spheres.sph %s" % selsph

    os.system(cmd)
    os.system(rename)

    return selsph

def split_sph(file_name, prefix):
    selsph = '%s_selected_spheres.sph' % prefix

    sph_file = open("%s" % file_name, "r")
    text = sph_file.read()
    sph_file.seek(0)
    lines = sph_file.readlines()
    sph_file.close()
    paragraphs = re.split("cluster    ...  number of spheres in cluster   ...\n", text)

    sel_sph = open(selsph, "w")
    sel_sph.write(lines[1])
    sel_sph.write(paragraphs[1])
    sel_sph.close()

    return selsph

def showbox(file_name, prefix):
    boxfile = "%s_box.pdb" % prefix

    in_box = open("%s_box.in" % prefix, "w")
    in_box.write("Y\n")
    in_box.write("0\n")
    in_box.write("%s\n" % file_name)
    in_box.write("1\n")
    in_box.write("%s" % boxfile)
    in_box.close()

    cmd = "showbox < %s_box.in" % prefix
    os.system(cmd)

    return boxfile

def grid(receptor_file, box_file, amber_params, prefix):
    in_name = "%s_grid.in" % prefix
    out_name = "%s_grid.out" % prefix

    in_grid = open("%s" % in_name, "w")
    grid_text = """compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  yes
contact_cutoff_distance        4.5

energy_score                   yes
energy_cutoff_distance         9999

atom_model                     all
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4

bump_filter                    yes
bump_overlap                   0.75

receptor_file                  %s
box_file                       %s
vdw_definition_file            %s
score_grid_prefix              %s_grid
""" % (receptor_file, box_file, amber_params, prefix)

    in_grid.write(grid_text)
    in_grid.close()

    cmd = "grid -i %s -o %s" % (in_name, out_name)
    os.system(cmd)

    return out_name

def do_dock6_rigid(ligand_path, sph_path, grid_path, amber, flex1, flex2, prefix):
    ligand_name = os.path.basename(args.ligand).split('.')[0]
    in_name = "%s_%s_dock.in" % (prefix, ligand_name)
    out_name = "%s_%s_dock.out" % (prefix, ligand_name)

    in_dock = open("%s" % in_name, "w")
    dock_text = """ligand_atom_file                                             %s
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           %s
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
flexible_ligand                                              no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       %s
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        yes
simplex_coefficient_restraint                                10.0
atom_model                                                   all
vdw_defn_file                                                %s
flex_defn_file                                               %s
flex_drive_file                                              %s
ligand_outfile_prefix                                        %s_%s_rigid
write_orientations                                           yes
num_scored_conformers                                        20
rank_ligands                                                 no
""" % (ligand_path, sph_path, grid_path.split('.')[0], amber, flex1, flex2, prefix, ligand_name)

    in_dock.write(dock_text)
    in_dock.close()

    cmd = "dock6 -i %s -o %s" % (in_name, out_name)
    os.system(cmd)

def do_dock6_flexible(ligand_path, sph_path, grid_path, amber, flex1, flex2,prefix):
    ligand_name = os.path.basename(args.ligand).split('.')[0]
    in_name = "%s_%s_flexdock.in" % (prefix, ligand_name)
    out_name = "%s_%s_flexdock.out" % (prefix, ligand_name)

    in_dock = open("%s" % in_name, "w")
    dock_text = """ligand_atom_file                                             %s
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           %s
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
flexible_ligand                                              yes
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          100
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100
use_clash_overlap                                            no
write_growth_tree                                            no
bump_filter                                                  yes
bump_grid_prefix                                             %s
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       %s
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        yes
simplex_coefficient_restraint                                10.0
atom_model                                                   all
vdw_defn_file                                                %s
flex_defn_file                                               %s
flex_drive_file                                              %s
ligand_outfile_prefix                                        %s_%s_flexible
write_orientations                                           no
num_scored_conformers                                        20
rank_ligands                                                 no
""" % (ligand_path, sph_path, grid_path.split('.')[0], grid_path.split('.')[0], amber, flex1, flex2, prefix, ligand_name)

    in_dock.write(dock_text)
    in_dock.close()

    cmd = "dock6 -i %s -o %s -v" % (in_name, out_name)
    os.system(cmd)

def do_dock6_amberscore(ligand_path, sph_path, grid_path, amber, flex1, flex2,prefix):
    ligand_name = os.path.basename(args.ligand).split('.')[0]
    in_name = "%s_%s_amberscore.in" % (prefix, ligand_name)
    out_name = "%s_%s_amberscore.out" % (prefix, ligand_name)

    in_dock = open("%s" % in_name, "w")
    dock_text = """ligand_atom_file                                             {}.amber_score.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
use_internal_energy                                          no
flexible_ligand                                              no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
descriptor_score_primary                                     no
descriptor_score_secondary                                   no
gbsa_zou_score_primary                                       no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_primary                                   no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_primary                                no
SASA_descriptor_score_secondary                              no
amber_score_primary                                          yes
amber_score_secondary                                        no
amber_score_receptor_file_prefix                             {}
amber_score_movable_region                                   ligand
amber_score_minimization_rmsgrad                             0.01
amber_score_before_md_minimization_cycles                    100
amber_score_md_steps                                         3000
amber_score_after_md_minimization_cycles                     100
amber_score_gb_model                                         5
amber_score_nonbonded_cutoff                                 18.0
amber_score_temperature                                      300.0
amber_score_abort_on_unprepped_ligand                        yes
ligand_outfile_prefix                                        output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
""".format()

    in_dock.write(dock_text)
    in_dock.close()

    cmd = "dock6 -i %s -o %s -v" % (in_name, out_name)
    os.system(cmd)


if __name__ == '__main__':

    import glob
    import argparse
    import shlex

    # load inputs from command line
    p = argparse.ArgumentParser(
        description='Run the DOCK steps on a folder of structure files. To run in the background, execute using nohup: nohup dock.py $BASENAME $NUMFRAMES /path/to/structures/ /path/to/parameters/ --ligand /path/to/ligand.mol2 --cofactors $COFACTORS --residues $RESIDUES > /path/to/logs/$LOGNAME &')
    p.add_argument('basename', help='base filename that you used to name your files')
    p.add_argument('numframes', help='total number of frames from your trajectory', type=int)
    p.add_argument('folder', help='path to folder with your structure files')
    p.add_argument('params', help='path to folder with parameter files')
    p.add_argument('--ligand', help='path to file of your ligand that you want to dock')
    p.add_argument('--cofactors', help='comma-separated string of cofactors that you want to keep while docking (e.g. SAM,DNC,WAT)')
    p.add_argument('--residues', help='comma-separated string of the binding residues')
    p.add_argument('--radius', help='radius around binding residues to dock to (default 9 A)', type=int, default=9)
    p.add_argument('--redock', help='run DOCK again for the specified ligand, even if docking files exist', default=False)
    args = p.parse_args()

    # file paths for docking parameters
    amb = os.path.join(args.params, 'vdw_AMBER_parm99.defn')
    f1 = os.path.join(args.params, 'flex.defn')
    f2 = os.path.join(args.params, 'flex_drive.tbl')

    
    print(args)
    # loading current files
    os.chdir(args.folder)
    #pdbs = glob.glob('{}-*.pdb'.format(args.basename))
    current_files = os.listdir(os.getcwd())

    # ligand name
    if args.ligand:
        ligandname = os.path.basename(args.ligand)

    # cofactors
    if args.cofactors:
        cofactors_list = shlex.split(args.cofactors)
    else:
        cofactors_list = []

    print('***************PARAMETERS***************')
    print('FULL LIST: {0}'.format(vars(args)))
    if args.ligand:
        print('LIGAND: {0}'.format(ligandname))
    if args.cofactors:
        print('COFACTORS: {0}'.format(cofactors_list))
    if args.residues:
        print('BINDING RESIDUES: {0}'.format(args.residues))
        print('RADIUS: {0}'.format(args.radius))

    counter = 1
    for frame in range(1,args.numframes+1):
        # just a simple counter
        print(str(counter) + '/' + str(args.numframes))
        counter += 1

        # file_prefix = '{0}-{1:03d}'.format(args.basename, frame)
        file_prefix = '{0}'.format(args.basename)
        print(file_prefix)

        # DOCKPREP
        # example: 3bwm-440_prep.mol2
        pdb = '{0}.pdb'.format(file_prefix)
        prepped_check = '%s_prep.mol2' % file_prefix
        if prepped_check in current_files:
            print('***DOCKPREP PREVIOUSLY RUN***')
            prepped_file = prepped_check
        else:
            print('RUNNING: DOCKPREP')
            prepped_file = dockprep(pdb, file_prefix)

        # ISOLATE RECEPTOR
        # example: 3bwm-440_receptor.mol2, 3bwm-440_receptor_noH.pdb
        receptor_check = '%s_receptor.mol2' % file_prefix
        receptor_noH_check = '%s_receptor_noH.pdb' % file_prefix
        if receptor_check in current_files and receptor_noH_check in current_files:
            print('***RECEPTOR FILES PREVIOUSLY GENERATED***')
            receptor, receptor_noH = receptor_check, receptor_noH_check
        else:
            print('RUNNING: ISOLATE RECEPTOR')
            receptor, receptor_noH = protein_only_and_noH(prepped_file, cofactors_list, file_prefix)

        # DMS
        # example: 3bwm-440_receptor.dms
        dms_check = '%s_receptor.dms' % file_prefix
        if dms_check in current_files:
            print('***SURFACE PREVIOUSLY GENERATED***')
            dms = dms_check
        else:
            print('RUNNING: DMS')
            dms = dms_maker(receptor_noH, file_prefix)

        # SPHGEN
        # example: 3bwm-440_receptor.sph
        sph_check = '%s_receptor.sph' % file_prefix
        if sph_check in current_files:
            print('***SPHERES PREVIOUSLY GENERATED***')
            sph = sph_check
        else:
            print('RUNNING: SPHGEN')
            sph = sphgen(dms, file_prefix)

        # SPHERE_SELECTOR
        # first choose binding site and save it as separate .mol2
        # example: 3BWY-418_binding_residues.mol2
        binding_site_mol2_file_check = '%s_binding_residues.mol2' % file_prefix
        if binding_site_mol2_file_check in current_files:
            print('***BINDING SITE RESIDUES ALREADY DEFINED***')
            binding_site_mol2_file = binding_site_mol2_file_check
        else:
            print('RUNNING: BINDING SITE MOL2')
            binding_site_mol2_file = binding_site_mol2(receptor_noH, args.residues, file_prefix)

        # then select the spheres based on these binding residues
        # example: 3bwm-440_selected_spheres_using_binding_residues.sph
        sel_sph_check = '%s_selected_spheres_using_binding_residues.sph' % file_prefix
        if sel_sph_check in current_files:
            print('***SPHERES ALREADY SELECTED***')
            sel_sph = sel_sph_check
        else:
            print('RUNNING: SPHERE_SELECTOR')
            sel_sph = sphere_selector_using_residues(sph, binding_site_mol2_file, args.radius, file_prefix)

        # SHOWBOX
        # example: 3bwm-440_box.pdb
        box_check = '%s_box.pdb' % file_prefix
        if box_check in current_files:
            print('***BOX PREVIOUSLY MADE***')
            box = box_check
        else:
            print('RUNNING: SHOWBOX')
            box = showbox(sel_sph, file_prefix)

        # GRID
        # example: 3bwm-440_grid.out
        gr_check = '%s_grid.out' % file_prefix
        if gr_check in current_files:
            print('***GRID PREVIOUSLY CALCULATED***')
            gr = gr_check
        else:
            print('RUNNING: GRID')
            gr = grid(receptor, box, amb, file_prefix)

        # DOCK
        if args.ligand:
            dock6_flexible_check = '%s_%s_flexible_scored.mol2' % (file_prefix, ligandname.split('.')[0])
            if dock6_flexible_check in current_files and not args.redock:
                print('***DOCK PREVIOUSLY RUN***')
            else:
                print('RUNNING: DOCK')
                do_dock6_flexible(args.ligand,sel_sph,gr,amb,f1,f2,file_prefix)

    print('***DOCKING COMPLETE***')
