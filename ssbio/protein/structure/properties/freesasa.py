import subprocess
import ssbio.utils
import os
import os.path as op
from collections import OrderedDict


def run_freesasa(infile, outfile, include_hetatms=True, outdir=None, force_rerun=False):
    """Run freesasa on a PDB file, output using the NACCESS RSA format.
    
    Args:
        infile (str): Path to PDB file (only PDB file format is accepted) 
        outfile (str): Path or filename of output file
        include_hetatms (bool): If heteroatoms should be included in the SASA calculations
        outdir (str): Path to output file if not specified in outfile
        force_rerun (bool): If freesasa should be rerun even if outfile exists

    Returns:
        str: Path to output SASA file

    """
    if not outdir:
        outdir = ''

    outfile = op.join(outdir, outfile)

    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        if op.exists(outfile):
            os.remove(outfile)
        if include_hetatms:
            shell_command = 'freesasa --format=rsa --hetatm {} -o {}'.format(infile, outfile)
        else:
            shell_command = 'freesasa --format=rsa {} -o {}'.format(infile, outfile)
        command = subprocess.Popen(shell_command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True)
        out, err = command.communicate()

    return outfile


def parse_rsa_data(rsa_outfile, ignore_hets=True):
    """Process a NACCESS or freesasa RSA output file. Adapted from Biopython NACCESS modele.
    
    Args:
        rsa_outfile (str): Path to RSA output file
        ignore_hets (bool): If HETATMs should be excluded from the final dictionary. This is extremely important
            when loading this information into a ChainProp's SeqRecord, since this will throw off the sequence matching.

    Returns:
        dict: Per-residue dictionary of RSA values

    """

    naccess_rel_dict = OrderedDict()

    with open(rsa_outfile, 'r') as f:
        for line in f:
            if line.startswith('RES'):
                res_name = line[4:7]
                chain_id = line[8]
                resseq = int(line[9:13])
                icode = line[13]
                res_id = (' ', resseq, icode)
                all_atoms_abs = line[16:22].strip()
                all_atoms_rel = line[23:28].strip()
                side_chain_abs = line[29:35].strip()
                side_chain_rel = line[36:41].strip()
                main_chain_abs = line[42:48].strip()
                main_chain_rel = line[49:54].strip()
                non_polar_abs = line[55:61].strip()
                non_polar_rel = line[62:67].strip()
                all_polar_abs = line[68:74].strip()
                all_polar_rel = line[75:80].strip()

                if all_atoms_rel =='N/A' and main_chain_rel =='N/A' and all_polar_rel =='N/A' and non_polar_rel =='N/A' and side_chain_rel =='N/A' and ignore_hets:
                    continue

                naccess_rel_dict[(chain_id, res_id)] = {
                    'res_name'      : res_name,
                    'all_atoms_abs' : ssbio.utils.conv_to_float(all_atoms_abs, inf_str='N/A'),
                    'all_atoms_rel' : ssbio.utils.conv_to_float(all_atoms_rel, inf_str='N/A'),
                    'side_chain_abs': ssbio.utils.conv_to_float(side_chain_abs, inf_str='N/A'),
                    'side_chain_rel': ssbio.utils.conv_to_float(side_chain_rel, inf_str='N/A'),
                    'main_chain_abs': ssbio.utils.conv_to_float(main_chain_abs, inf_str='N/A'),
                    'main_chain_rel': ssbio.utils.conv_to_float(main_chain_rel, inf_str='N/A'),
                    'non_polar_abs' : ssbio.utils.conv_to_float(non_polar_abs, inf_str='N/A'),
                    'non_polar_rel' : ssbio.utils.conv_to_float(non_polar_rel, inf_str='N/A'),
                    'all_polar_abs' : ssbio.utils.conv_to_float(all_polar_abs, inf_str='N/A'),
                    'all_polar_rel' : ssbio.utils.conv_to_float(all_polar_rel, inf_str='N/A')}

    return naccess_rel_dict