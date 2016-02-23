if __name__ == '__main__':
    from Bio import PDB
    p = PDB.PDBParser()
    structure = p.get_structure("1mot", "test_structures/1mot.pdb")
    model = structure[0]
    dssp = PDB.DSSP(model, "test_structures/1mot.pdb")
    # DSSP data is accessed by a tuple (chain_id, res_id)

    a_key = list(dssp)[2]

    #print(a_key)
    dssp[a_key]
    # residue object, secondary structure, solvent accessibility,
    # relative accessiblity, phi, psi
    # print(calc_sasa('1mot.pdb'))

