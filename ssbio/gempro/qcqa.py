import requests
import cachetools

SEVEN_DAYS = 60 * 60 * 24 * 7

@cachetools.func.ttl_cache(maxsize=500, ttl=SEVEN_DAYS)
def get_best_structures(uniprot_id):
    """Use the PDBe REST service to query for the best PDB structures for a UniProt ID.

    More information found here: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html
    Link used to retrieve results: https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/:accession
    The list of PDB structures mapping to a UniProt accession sorted by coverage of the protein and, if the same, resolution.

    Args:
        uniprot_id: a valid UniProt ID

    Returns:
        A rank-ordered list of dictionaries, which contain these keys:
        pdb_id: the PDB ID which maps to the UniProt ID
        chain_id: the specific chain of the PDB which maps to the UniProt ID
        coverage: the percent coverage of the entire UniProt sequence
        resolution: the resolution of the structure
        start: the structure residue number which maps to the start of the mapped sequence
        end: the structure residue number which maps to the end of the mapped sequence
        unp_start: the sequence residue number which maps to the structure start
        unp_end: the sequence residue number which maps to the structure end
        experimental_method: type of experiment used to determine structure
        tax_id: taxonomic ID of the protein

    """
    r = requests.get('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{}'.format(uniprot_id))
    if not r.status_code == 200:
        return None

    # evaluate false, true, and null correctly (http://stackoverflow.com/questions/1083250/running-json-through-pythons-eval)
    best_structures = eval(r.text, {'false': False, 'true': True, 'null': None})

    # return a rank ordered list
    return best_structures[uniprot_id]


        # below is code from GEM-PRO stage 3, QC/QA
        # import ranking as rk
        #
        #
        # def ranker(ranks, reverse_flag):
        #     """Rank a list of tuples (id, property) by the numeric property.
        #
        #     Args:
        #         ranks: list of tuples: [(identifier (str), property (float)), ...]
        #         reverse_flag: rank in descending order if True
        #
        #     Returns:
        #
        #     """
        #     # TODO: reverse_flag is confusing
        #     ranked = {}
        #
        #     idents = [x[0] for x in ranks]
        #     prop = [x[1] for x in ranks]
        #     rank = list(rk.Ranking(prop, reverse=reverse_flag))
        #
        #     for i,v in enumerate(rank):
        #         ranked[idents[i]] = v[0]
        #
        #     return ranked
        #
        #
        # def ranker_averager(ro1, ro2):
        #     avg_rank = {}
        #     for k in ro1:
        #         avg = (float(ro1[k]) + float(ro2[k])) / 2
        #         avg_rank[k] = avg
        #
        #     order_by_avg_rank = sorted(avg_rank.iteritems(), key=operator.itemgetter(1))
        #     return order_by_avg_rank
        #
        #
        # def remap( x, oMin, oMax, nMin, nMax ):
        #     """function to map to a 0 to 1 scale
        #     http://stackoverflow.com/questions/929103/convert-a-number-range-to-another-range-maintaining-ratio
        #
        #     Args:
        #         x:
        #         oMin:
        #         oMax:
        #         nMin:
        #         nMax:
        #
        #     Returns:
        #
        #     """
        #
        #     #range check
        #     if oMin == oMax:
        #         # print "Warning: Zero input range"
        #         return None
        #
        #     if nMin == nMax:
        #         # print "Warning: Zero output range"
        #         return None
        #
        #     #check reversed input range
        #     reverseInput = False
        #     oldMin = min( oMin, oMax )
        #     oldMax = max( oMin, oMax )
        #     if not oldMin == oMin:
        #         reverseInput = True
        #
        #     #check reversed output range
        #     reverseOutput = False
        #     newMin = min( nMin, nMax )
        #     newMax = max( nMin, nMax )
        #     if not newMin == nMin :
        #         reverseOutput = True
        #
        #     portion = (x-oldMin)*(newMax-newMin)/(oldMax-oldMin)
        #     if reverseInput:
        #         portion = (oldMax-x)*(newMax-newMin)/(oldMax-oldMin)
        #
        #     result = portion + newMin
        #     if reverseOutput:
        #         result = newMax - portion
        #
        #     return result
        #
        #
        # def jaccard_index_resnums(residue_list1, residue_list2):
        #     union = sorted(list(set(residue_list1).union(residue_list2)))
        #     residue_list1_final = [0] * len(union)
        #     residue_list2_final = [0] * len(union)
        #
        #     for i, v in enumerate(union):
        #         if v in residue_list1:
        #             residue_list1_final[union.index(v)] = v
        #         if v in residue_list2:
        #             residue_list2_final[union.index(v)] = v
        #     return jaccard_similarity_score(residue_list1_final, residue_list2_final)
        #
        #
        # def rank_by_sequence(canonical_sequence, pdb_ids):
        #     for pdb_id in pdb_ids:
        #
        #         pdb_fasta = SEQ_PDB_FILES + pdb_id + '.faa'
        #         if not os.path.exists(pdb_fasta):
        #             warnings.warn('***WARNING: No PDB FASTA file found for %s. Please check if step 3.3 was run.***' % pdb_id)
        #             seq_pdb_align_errors.append((aligner_id, pdb_id))
        #             continue
        #
        #         # RANKING USING ALIGNMENT STATISTICS
        #         alignment_filename = "%s_%s_align.txt" % (aligner_id, pdb_id)
        #         if not os.path.exists(SEQ_ALIGN_FILES + alignment_filename):
        #             try:
        #                 #                 print 'running alignment'
        #                 alignment_filename = ssbio.sequence.align.run_alignment(aligner_id, seq_fasta, pdb_id, pdb_fasta)
        #             except:
        #                 seq_pdb_align_errors.append((aligner_id, pdb_id))
        #                 warnings.warn('***ERROR: with alignment of %s and %s***' % (aligner_id, pdb_id))
        #                 alignment_scores[pdb_id] = 0
        #                 alignment_coverage[pdb_id] = 0
        #                 alignment_coverage_sim[pdb_id] = 0
        #                 continue
        #         # getting the "raw" alignment dataframe (includes chains that might not match the gene uniprot)
        #         # also does not include pdb_start/pdb_stop yet
        #
        #
        #         pdb_alignment_plus_start_stop_df_path = '%s_%s_align_table.csv' % (aligner_id, pdb_id)
        #         if not os.path.exists(SEQ_ALIGN_FILES + pdb_alignment_plus_start_stop_df_path):
        #             pdb_alignment_df = ssbio.sequence.align.get_alignment_df(alignment_filename)
        #
        #             # ADDING PDB_START AND PDB_STOP ##############################
        #             # split the chain ID into another column
        #             separate_chain = pdb_alignment_df.join(pdb_alignment_df['id_b'].apply(lambda x: pd.Series(x.split('.')[1])))
        #             separate_chain = separate_chain.rename(columns={0: 'chain'})
        #             # adding pdb_start and pdb_stop as empty columns
        #             separate_chain['id_b_start'] = np.nan
        #             separate_chain['id_b_stop'] = np.nan
        #             # getting the pdb start residues for all chains using the previously generated dataframe
        #             pdb_res_starts_slice = DF_03A_PDB_STARTS[DF_03A_PDB_STARTS.pdb_id == pdb_id]
        #             chains_and_starts = dict(pdb_res_starts_slice.set_index(['chain_id', 'start_res']).index.tolist())
        #             # now adding in pdb_start and pdb_stop
        #             for chain in separate_chain.chain.unique():
        #                 #                 print 'adding pdb starts', pdb_id, chain
        #                 subset = separate_chain[separate_chain.chain == chain]
        #                 if chain == 'NA':
        #                     chain = nan
        #                 pdb_start = chains_and_starts[chain]
        #                 for idx, row in subset.iterrows():
        #                     if row['type'] == 'insertion':
        #                         adder = row['count'] - 1
        #                         pdb_stop = adder + pdb_start
        #                         separate_chain.loc[idx, 'id_b_start'] = pdb_start
        #                         separate_chain.loc[idx, 'id_b_stop'] = pdb_stop
        #                         pdb_start = pdb_stop + 1
        #                     elif row['type'] != 'deletion':
        #                         adder = row['id_a_stop'] - row['id_a_start']
        #                         pdb_stop = adder + pdb_start
        #                         separate_chain.loc[idx, 'id_b_start'] = pdb_start
        #                         separate_chain.loc[idx, 'id_b_stop'] = pdb_stop
        #                         pdb_start = pdb_stop + 1
        #             # separate_chain should now have pdb_start and pdb_stop
        #             # save this...
        #             separate_chain.to_csv(SEQ_ALIGN_FILES + pdb_alignment_plus_start_stop_df_path)
        #             # END ADDING PDB_START AND PDB_STOP ##############################
        #         else:
        #             separate_chain = pd.read_csv(pdb_alignment_plus_start_stop_df_path, index_col=0)
