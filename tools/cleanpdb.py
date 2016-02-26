def final_mutation_function(corrected_df, pdb_id, pdb_file):
    # default is to keep atoms with alternate location ID 'A'
    alt = 'A'
    # default is to not keep "alt" identifier
    keep_alt = 0

    output = pdb_id + '_modified.pdb'
#     output = pdb_id + '.pdb'

    #read in structure to parser
    parser = PDB.PDBParser()
    struct = parser.get_structure(pdb_id, pdb_file)

    # create PDBIO class for writing
    io = PDB.PDBIO()
    io.set_structure(struct)
    # save it
    io.save('/tmp/%s_no_disorder.pdb' % pdb_id, select=NotDisordered(alt))

    # choose only parts of the protein to mutate
    df = corrected_df[corrected_df.type == 'mutation']
    # load structure (no disorder) into biopython & perform modifications
    structure_new = setup_pdb_for_amber(df, pdb_id, '/tmp/%s_no_disorder.pdb' % pdb_id)
    #  write out new structure
    write_to_pdb(structure_new, output)

    return output

def setup_pdb_for_amber_monomer(df, name, file_name):

    '''takes in a PDB or homology file and
    (1) performs site-directed mutagenesis given a user-defined sequence change,
    (2) removes all hetero atoms (e.g. water and metals/cofactors) and
    (3) removes disordered atoms
    MOD: keep only 1 chain (df is already the "corrected" df with only the gene of interest anyway)
    '''

    parser = PDBParser()
    structure = parser.get_structure(name, file_name)
    model = structure[0]
#     chains = [i.id for i in model.child_list]

    # get mutations, if any
    mutation_df = df[df.type == 'mutation']

    # choosing a representative chain
    sorted_by_match_and_chain = df[df['type']=='match'].sort(['count','chain'], ascending=[False, True]).set_index(['chain','count']).index.tolist()
    chain_id = str(sorted_by_match_and_chain[0][0])

    # operate only on this chain
    chain = model[chain_id]
    residue_list = chain.child_list

    for residue in list(residue_list):
        # mutate residues according to sequence alignment data
#             print residue.id[1]
        if len(mutation_df)>0:
            if residue.id[1] in df[df.id_b.str.endswith(chain_id)].id_b_start.tolist():

#                 print name, residue.id[1], chain_id
#                 print df
                res_w = df[df.id_b.str.endswith(chain_id)][df.id_b_start == residue.id[1]].id_b_aa.values[0]
                res_mut = df[df.id_b.str.endswith(chain_id)][df.id_b_start == residue.id[1]].id_a_aa.values[0]

                if res_mut not in AAdict2.keys():
                    warnings.warn('***UNKNOWN AMINO ACID IN UNIPROT SEQUENCE. SKIPPING***') ### TO CHANGE!!!!!!!!!!!!!!!
                    continue

                print df[df.id_b.str.endswith(chain_id)][df.id_b_start == residue.id[1]].values
                print '\n'
                print residue.id[1], residue.resname.upper()," | ", res_w, "mutate to:", res_mut

                # Remove all atoms except protein backbone atoms:
                for atom in residue.child_dict.keys():
                    if atom in keep_atom_list:
                        pass
                    else:
                        residue.detach_child(atom)

                residue.resname = AAdict2[res_mut]
                # if residue is non-standard, make it a non-heteroatom
                if residue.id[0] != ' ':
                    resnum = residue.id[1]
                    residue.id = (' ',resnum, ' ')
                # check that its the residue you expect and rename the residue:
#                 if AAdict[residue.resname] == res_w:#AAdict[res_w]:
#                     print 'matches \n'
#                     residue.resname = AAdict2[res_mut]
#                     print residue.resname
#                 else:
#                     raise ValueError("Unrecognised residue %r" % residue.resname)

    # Remove all hydrogens and heteroatom residues (e.g. WAT or metal) from structure:
    for residue in list(chain):
        id = residue.id
        if id[0] != ' ':
            chain.detach_child(id)
        if len(chain) == 0:
            model.detach_child(chain.id)
        for atom in residue.get_list():
        #print residue.resname, residue.id[1], atom.element
            if atom.element == 'H':
                residue.detach_child(atom.id)


    return chain

def write_to_pdb(struct, new_file_name):
    w = PDBIO()
    w.set_structure(struct)
    w.save(new_file_name)
