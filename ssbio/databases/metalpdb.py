import glob
import os.path as op
from ssbio.protein.structure.structprop import StructProp


def get_metalpdb_info(metalpdb_lig_file):
    """Parse a MetalPDB .lig file and return a tuple of the chain ID it represents, along with metal binding information.

    Args:
        metalpdb_lig_file (str): Path to .lig file

    Returns:
        tuple: (str, dict) of the chain ID and the parsed metal binding site information

    """

    pdb_metals = ['CU', 'ZN', 'MN', 'FE', 'MG', 'CO', 'SE', 'YB', 'SF4', 'FES', 'F3S', 'NI', 'FE2']

    # Information to collect
    coordination_number = 0
    endogenous_ligands = []
    exogenous_ligands = []

    # Load the structure
    ss = StructProp(ident='metalpdb', structure_path=metalpdb_lig_file, file_type='pdb')

    # This lig file should just be for one chain
    chain_id = op.basename(metalpdb_lig_file)[5]
    metal_id = (op.basename(metalpdb_lig_file).split('_')[2], op.basename(metalpdb_lig_file).split('_')[3])

    for r in ss.parse_structure().first_model.get_residues():
        return_id = (r.get_id(), r.get_resname())
        #         print(r.resname)
        # Binding partners
        ## check if residue is a normal one (not a HETATM, WAT, or the metal that is identified)
        if r.get_id()[0] != ' ':
            if not r.resname.strip() in pdb_metals and r.resname != 'HOH':
                #                 print('appended', r.resname)
                exogenous_ligands.append(return_id)
        else:
            endogenous_ligands.append(return_id)

        # Coordination number
        for a in r.get_atom():
            if not a.element in pdb_metals:
                coordination_number += 1

    infodict = {metal_id: {'endogenous_ligands' : endogenous_ligands,
                           'exogenous_ligands'  : exogenous_ligands,
                           'coordination_number': coordination_number}}

    return chain_id, infodict