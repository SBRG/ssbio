import warnings
from Bio import SeqIO


import pandas as pd
def gene_to_uniprot(gene, isoforms, mapping_df, isoform_df, metadata_df):
    '''
    Input: a gene, its isoforms - 314, [1,2]
    Output: the mapped UniProt IDs to each of the isoforms as a dictionary {'314.1': 'O75106-1', '314.2': 'O75106-2'}

    Checks if there are UniProt ID mappings for a gene and its isoforms. A couple cases to consider:
    1. (Perfect) There are the same number of isoforms in UniProt as the entry query. The correct ones map to the correct numbers (ie .1 -> -1 and .2 -> -2)
    2. (Imperfect) The isoform IDs in UniProt are odd (-A and -B, or -3 is the main isoform)
        ## MAYBE NOT TODO: how to check this one? is there an isoform comment? probably have to manually read that. but there don't seem to be any!! yay!
        Output: {'314.1': 'O75106-A', '314.2': 'O75106-B'}
    3. (Missing) There are a lesser number UniProt isoforms for the gene.
        3a. Because the other isoform is in an unreviewed entry
            Output: {'314.1': 'O75106-1', '314.2': 'A0SDF98JH-1'}
        3b. Because UniProt doesn't have it
            Output: {'314.1': 'O75106-1', '314.2': None}
            # TODO: is this the best way to do it?
    4. (Missing) There is no UniProt entry
        Output: None
    5. (Other) There are multiple reviewed entries. These must be manually handled.

    ##TODO: please draw a diagram or something to explain what's going on here
    '''
    g_to_u_df = mapping_df[mapping_df.m_gene == gene]
    g_to_u_all = g_to_u_df.u_uniprot_acc.unique().tolist()
    g_to_u_reviewed = g_to_u_df[g_to_u_df.u_reviewed == True].u_uniprot_acc.unique().tolist()
    g_to_u_unreviewed = g_to_u_df[g_to_u_df.u_reviewed == False].u_uniprot_acc.unique().tolist()

    # if we cannot map to reviewed or unreviewed, we must give up
    if len(g_to_u_reviewed) == 0 and len(g_to_u_unreviewed) == 0:
#         print 'ah'
        warnings.warn('No reviewed or unreviewed UniProt entries for gene %s and its isoforms %s' % (gene, isoforms), UserWarning)
        return None

    # now check the number of isoforms
    if len(g_to_u_reviewed) == 1:
        g_to_u_reviewed_match = g_to_u_reviewed[0]

        # load the isoforms
        mapped_isoforms_df = isoform_df[isoform_df.u_uniprot_acc == g_to_u_reviewed_match]
        mapped_uni_isoforms = mapped_isoforms_df.u_isoform_id.tolist()
        mapping_dict = {}

        # check if the number of isoforms available is greater than or equal to the number of expected isoforms
        if len(mapped_uni_isoforms) >= len(isoforms):
            # for each of the recon 2 isoform ids...
            for isoform_id in isoforms:
                # create a hypothetical uniprot isoform
                ssb_uni_isoform = g_to_u_reviewed_match + '-' + isoform_id
                # check to see if it exists
                if ssb_uni_isoform in mapped_uni_isoforms:
                    mapping_dict[gene + '.' + isoform_id] = ssb_uni_isoform
                # if it doesn't, and there is only one isoform, that means there is an issue with the metadata retrieval. Take the original uniprot ID as the mapping
                elif len(isoforms) == 1:
                    mapping_dict[gene + '.' + isoform_id] = g_to_u_reviewed_match + '-' + isoform_id
                else:
                    warnings.warn('TODO: unknown error', UserWarning)
                    return None
            return mapping_dict

        # if the number of isoforms available is less than the number of expected isoforms..
        # we first check the unreviewed uniprots
        elif len(g_to_u_unreviewed) >= 1:
            # get the sequence of the reviewed uniprot (from metadata DF just to make sure..)
            rev_uni_seqs = isoform_df[isoform_df.u_uniprot_acc == g_to_u_reviewed_match].u_seq.unique().tolist()
            # for the unreviewed entries, check if they are different sequences!
            mini_len_dict = {}
            # make sure to make a copy of the list since we are removing things..
            for unrev_uni in g_to_u_unreviewed[:]:
                unrev_uni_seq = metadata_df.loc[unrev_uni].u_seq
#                 print unrev_uni, unrev_uni_seq
                if unrev_uni_seq in rev_uni_seqs:
                    # remove the same sequences
#                     print 'REMOVE', unrev_uni
                    g_to_u_unreviewed.remove(unrev_uni)
                else:
                    # record the length for later use
#                     print 'KEEP', unrev_uni, len(unrev_uni_seq)
                    mini_len_dict[unrev_uni] = len(unrev_uni_seq)

            # now if our collection of unreviewed and reviewed entries is greater than the number of isoforms...
            if len(g_to_u_unreviewed) + len(mapped_uni_isoforms) >= len(isoforms):
                # match the reviewed isoforms to the genes like above..
                # and if the hypothetical uniprot isoform does not exist, match it to the longest unreviewed isoform

                # first sort the list of unreviewed isoforms..
                sorted_unrev_uni = [y[0] for y in sorted(mini_len_dict.items(), key=lambda x: x[1])]
                mapping_dict = {}
#                 print gene, isoforms, mapped_uni_isoforms,mini_len_dict,sorted_unrev_uni
                for isoform_id in isoforms:
                    # first map to reviewed entries
                    ssb_uni_isoform = g_to_u_reviewed_match + '-' + isoform_id
                    if ssb_uni_isoform in mapped_uni_isoforms:
                        mapping_dict[gene + '.' + isoform_id] = ssb_uni_isoform

                    # then map to longest unreviewed
                    else:
                        if len(sorted_unrev_uni) == 0:
                            break
                        # print sorted_unrev_uni
                        mapping_dict[gene + '.' + isoform_id] = sorted_unrev_uni[-1]
                        sorted_unrev_uni.pop(-1)
#                 print mapping_dict
                return mapping_dict

            # if our collection of unreviewed and reviewed entries is still less than the number of isoforms...
            # just map as many as we can
            else:
                for isoform_id in isoforms:
                    # create a hypothetical uniprot isoform
                    ssb_uni_isoform = g_to_u_reviewed_match + '-' + isoform_id
                    # check to see if it exists
                    if ssb_uni_isoform in mapped_uni_isoforms:
                        mapping_dict[gene + '.' + isoform_id] = ssb_uni_isoform
                    # if it doesn't, and there is only one isoform, that means there is an issue with the metadata retrieval. Take the original uniprot ID as the mapping
                    elif len(isoforms) == 1:
                        mapping_dict[gene + '.' + isoform_id] = g_to_u_reviewed_match
                    else:
                        mapping_dict[gene + '.' + isoform_id] = None
                return mapping_dict
        # if we just have a reviewed entry, and it has not enough isoforms
        # just map as many as we can
        else:
            for isoform_id in isoforms:
                # create a hypothetical uniprot isoform
                ssb_uni_isoform = g_to_u_reviewed_match + '-' + isoform_id
                # check to see if it exists
                if ssb_uni_isoform in mapped_uni_isoforms:
                    mapping_dict[gene + '.' + isoform_id] = ssb_uni_isoform
                # if it doesn't, and there is only one isoform, that means there is an issue with the metadata retrieval. Take the original uniprot ID as the mapping
                elif len(isoforms) == 1:
                    mapping_dict[gene + '.' + isoform_id] = g_to_u_reviewed_match
                else:
                    mapping_dict[gene + '.' + isoform_id] = None
            return mapping_dict
    elif len(g_to_u_reviewed) == 0 and len(g_to_u_unreviewed) >= len(isoforms):
        mapping_dict = {}
        mini_len_dict = {}
        # get the length of the unreviewed sequences
        for unrev_uni in g_to_u_unreviewed[:]:
            unrev_uni_seq = metadata_df.loc[unrev_uni].u_seq
            mini_len_dict[unrev_uni] = len(unrev_uni_seq)
            
        sorted_unrev_uni = [y[0] for y in sorted(mini_len_dict.items(), key=lambda x: x[1])]
        for isoform_id in isoforms:
            mapping_dict[gene + '.' + isoform_id] = sorted_unrev_uni[-1]
            sorted_unrev_uni.pop(-1)
        return mapping_dict
        
        

import glob
SEQ_SUNPRO_FILES = '/home/nathan/projects/GEM-PRO/SUNPRO/refseq_sequences/'
available_refseq_files = glob.glob(SEQ_SUNPRO_FILES + '*.faa')
available_refseq = {}
for f in available_refseq_files:
    refseq_id = f.split('.faa')[0].split('/')[-1].split('.')[0]
    available_refseq[refseq_id] = f.split('/')[-1]

def gene_to_refseq(gene, isoforms, mapping_df):#, metadata_df):
    '''
    Input: a gene, its isoforms (314, [1,2]), mapping_df (the gene to refseq df with no isoforms), metadata_df (the refseq df with sequence and whatnot)
    Output: the mapped RefSeq IDs to each of the isoforms as a dictionary {'314.1': 'NP_033720', '314.2': 'NP_001149'}

    Checks if there are RefSeq ID mappings for a gene and its isoforms. A couple cases to consider:
    1. (Perfect) There are the same (or greater) number of isoforms in RefSeq as the entry query. The correct ones map to the correct transcript names (ie .1 -> -001 and .2 -> -002)
    2. (Imperfect) The isoform IDs in RefSeq are odd (-201 and -202, or -5 is the main isoform)
        ## NOTTODO: how to check this one? is there an isoform comment? seems like we have to do this manually...
        ## ANSWER: just map by length of sequence! https://www.biostars.org/p/19029/ - longest isoform is "canonical"
        Output: {'314.1': 'NP_033720', '314.2': 'NP_001149'}
    3. (Missing) There is only one RefSeq isoform for the mapped entry when there should be more.
        Output: None
    4. (Missing) There is no RefSeq entry
        Output: None
    5. (Other) There are multiple entries. These must be manually handled.
    '''

    isoforms = sorted(isoforms)
    g_to_many_df = mapping_df[mapping_df.m_gene_noiso == gene]
    g_to_r = g_to_many_df[pd.notnull(g_to_many_df.bm_refseq)].bm_refseq.unique().tolist()

    if len(g_to_r) == 0:
#         warnings.warn('WARNING: no mapped RefSeq IDs for gene %s and isoforms %s' % (gene, isoforms))
        return None

    # load the sequences of the isoforms, record their length. assign them to the isoforms based on length
    refseq_lengths = {}
    refseq_seqs = []
    for mapped_r in g_to_r:
        if mapped_r in available_refseq.keys():
            mapped_r_file = available_refseq[mapped_r]
            handle = open(SEQ_SUNPRO_FILES + mapped_r_file, "rU")
            for record in SeqIO.parse(handle, "fasta") :
                if record.seq not in refseq_seqs:
                    refseq_lengths[mapped_r] = len(record.seq)
                    refseq_seqs.append(record.seq)
            handle.close()
    sorted_refseqs = [y[0] for y in sorted(refseq_lengths.items(), key=lambda x: x[1])]
    if len(refseq_lengths.keys()) >= len(isoforms):
        mapping_dict = {}
        for iso in isoforms:
            mapping_dict[gene + '.' + iso] = sorted_refseqs.pop(-1)
        return mapping_dict
    
    # if the number of refseq ids is less than the number of isoforms
    # we should not remove the sequences if they are the same - they are just probably protein sequences with different UTRs or something
    elif len(refseq_lengths.keys()) < len(isoforms):
        
        # redo the above code but add even if the sequence is in there
        refseq_lengths = {}
        transcript_names = {}
        refseq_seqs = []
        for mapped_r in g_to_r:
            if mapped_r in available_refseq.keys():
                mapped_r_file = available_refseq[mapped_r]
                handle = open(SEQ_SUNPRO_FILES + mapped_r_file, "rU")
                for record in SeqIO.parse(handle, "fasta") :
                    refseq_lengths[mapped_r] = len(record.seq)
                    transcript_names[mapped_r] = g_to_many_df[g_to_many_df.bm_refseq == mapped_r].bm_name.unique()[0]
                    refseq_seqs.append(record.seq)
                handle.close()
#         sorted_refseqs = [y[0] for y in sorted(refseq_lengths.items(), key=lambda x: x[1])]
        # assign based on transcript name?
        sorted_refseqs = [y[0] for y in sorted(transcript_names.items(), key=lambda x: x[1], reverse=True)]
        
        mapping_dict = {}
        for iso in isoforms:
            if len(sorted_refseqs) == 0:
                mapping_dict[gene + '.' + iso] = None
            else:
                mapping_dict[gene + '.' + iso] = sorted_refseqs.pop(-1)
        return mapping_dict
