def pdb_chain_stoichiometry_biomolone(pdbid):
    '''
    Takes in a PDB ID and returns the stoichiometry of the chains in biological assembly 1 as a dictionary.
    Steps taken are:
    1) download PDB and parse header, make biomolecule if provided
    2) count how many times each chain appears in biomolecule #1
    3) convert chain id to uniprot id
    4) return final dictionary
    Input: pdbid - a 4 character PDB ID
    Output: a dictionary of {(ChainID,UniProtID): # occurences}
    '''

    try:
        # change to biomol=True if you want biological assemblies
        pdb, header = pr.parsePDB(pdbid, header=True, biomol=True, secondary=False)

        # if multiple biomols choose the first one, also if NMR choose the first coord set
        if type(pdb) == list:
            pdb = pdb[0]

    except IOError:
        return None, None
    # if there are no biological assemblies
    except ValueError:
        pdb, header = pr.parsePDB(pdbid, header=True, biomol=False, secondary=False)

    # count the occurences of each chain
    chain_stoich = defaultdict(int)
    hier = pdb.getHierView()
    for chain in hier:
        chain_stoich[chain.getChid()] += 1

    # DBREF entry in PDB file sometimes contains obsolete UniProt entries
    # chain_to_uniprot = {}
    # for chain in header['polymers']:
        # for dbref in chain.dbrefs:
            # if dbref.database.lower() == 'uniprot':
                # chain_to_uniprot[chain.chid] = dbref.accession

    # convert chain IDs to uniprot IDs
    chain_to_uniprot = {}
    for chain in header['polymers']:
        try:
            chain_to_uniprot[chain.chid] = sifts_pdb_chain_to_uniprot(pdbid.lower(), chain.chid)
        except KeyError:
            chain_to_uniprot[chain.chid] = ['PDB-'+chain.chid]

    # keep both chain ID and uniprot ID (this is the final dictionary)
    combined = {}
    for k,v in chain_to_uniprot.iteritems():
        for uni in v:
            combined[(k,uni)] = chain_stoich[k]

    return combined

def pisa_complex_information(pdb_id):

    pdb_id = pdb_id.lower()
    pisa = {}

    # request the xml file for a PDB
    r = requests.post('http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimers.pisa?%s' % pdb_id)
    soup = BeautifulSoup(r.content)

    # if PISA can't calculate an assembly...
    if 'not found' in str(soup.pdb_entry.status.string) or 'No symmetry' in str(soup.pdb_entry.status.string) or 'Overlapping structures' in str(soup.pdb_entry.status.string):
        pisa[pdb_id] = 'ERROR'
        return pisa

    # if it is a monomer...
    num_complexes = int(soup.total_asm.string)
    if num_complexes == 0:
        pisa[pdb_id] = 'MONOMER'
        return pisa

    # otherwise, get the complex information!
    elif num_complexes > 0:

        # all "assembly sets" (see PISA sets for more info)
        sets = soup.findAll('asm_set')

        for s in sets:

            set_id = int(s.ser_no.string)

            # all assemblies
            complexes = s.findAll('assembly')

            for cplx in complexes:

                ############################################################################################
                # this part tells you the actual composition of the predicted complex (chains and ligands) #
                parts = cplx.findAll('molecule')

                chains = defaultdict(int)

                for part in parts:
                    part_id = part.chain_id.string
                    if part_id.startswith('['):
                        part_id = 'LIG_' + part_id.split(']')[0].strip('[')
                    chains[str(part_id)] += 1

                ligands = {}

                for key in chains.keys():
                    if key.startswith('LIG_'):
                        ligands[str(key.split('_')[1])] = chains.pop(key)

                chains_final = {}
                for k,v in chains.iteritems():
                    try:
                        chain_to_uniprot = sifts_pdb_chain_to_uniprot(pdb_id.lower(), k)
                    except KeyError:
                        chain_to_uniprot = ['PDB-' + k]
                    for m in chain_to_uniprot:
                        chains_final[(str(k), str(m))] = v
                ############################################################################################

                # this part give you something to add to a dataframe
                adder = {}

                cplx_id = int(cplx.id.string)
                cplx_composition = str(cplx.composition.string)

                d_g_diss = float(cplx.diss_energy.string)
                d_g_int = float(cplx.int_energy.string)

                pdb_biomol = int(cplx.r350.string)

                if d_g_diss >= 0:
                    stable = True
                else:
                    stable = False

                adder['cplx_composition'] = cplx_composition.strip()
                adder['cplx_chains'] = chains_final
                adder['cplx_ligands'] = ligands
                adder['stable'] = stable
                adder['d_g_diss'] = d_g_diss
                adder['d_g_int'] = d_g_int
                adder['pdb_biomol'] = pdb_biomol

                pisa[(set_id,cplx_id)] = adder

        return pisa
