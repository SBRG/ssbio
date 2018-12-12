import json
import os.path as op
import pandas as pd
import ssbio.utils
import pickle
from collections import defaultdict
from tqdm import tqdm_notebook as tqdm

filter_proteome = [False, ]

filter_observations = ['A22.CCCP', 'A22.PYOCYANIN', 'CCCP.2', 'CHIR090.CCCP', 'DOXYCYCLINE.PARAQUAT', 'DOXYCYCLINE.PMS',
                       'DOXYCYCLINE.PYOCYANIN', 'FUMARATE.40C', 'FUMARATE.40MM', 'FUMARATE.A22', 'FUMARATE.CEFSULODIN',
                       'FUMARATE.PARAQUAT', 'FUMARATE.TOBRAMYCIN', 'NACL.PARAQUAT', 'NOSALT.FUMARATE', 'NOSALT.PARAQUAT',
                       'PARAQUAT.10UM', 'PH5.FUMARATE', 'PMS.PROCAINE', 'PYOCYANIN.0P2UM', 'PYOCYANIN.10', 'PYOCYANIN.1UM',
                       'TRIMETHOPRIM.PYOCYANIN', 'PYOCYANIN.10UM', 'pathotype_simple',
                       'ros_simulated', 'pathotype_simple2', 'pathotype_simple3'] #'isolation_source', 'pathotype', ]

filter_ismem = [('vizrecon','membrane'),
                ('vizrecon','inner_membrane'),
                ('vizrecon','outer_membrane'),]
                #('tmhmm','membrane')]

filter_notmem = [#('vizrecon','non_membrane'),
                 ('vizrecon','periplasm'),
                 ('vizrecon','cytosol'), ]
                 #('vizrecon','extracellular'),]
                 #('tmhmm','non_membrane')]

filter_ismem_subseq = [['all', 'acc_3D', 'metal_2_5D', 'metal_3D', 'tm_2D', 'tm_3D',
                        'csa_2_5D', 'sites_2_5D'], ]  # need the comma to treat it like one set
filter_notmem_subseq = [['all', 'disorder_2D', 'ss_disorder_2D', 'disorder_3D', 'ss_disorder_3D', 'acc_2D',
                         'acc_3D', 'surface_3D', 'metal_2_5D', 'metal_3D', 'dna_2_5D', 'csa_2_5D', 'sites_2_5D'], ]
filter_all_subseq = [['all', 'disorder_2D', 'ss_disorder_2D', 'disorder_3D', 'ss_disorder_3D', 'acc_2D',
                      'tm_2D', 'tm_3D',
                      'acc_3D', 'surface_3D', 'metal_2_5D', 'metal_3D', 'dna_2_5D', 'csa_2_5D', 'sites_2_5D'], ]

filter_subseq_suffixes = {'all': ['aa_%_ord', 'aa_%_dis', 'aa_%_M', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg', 'aa_%_Y'],
                          'disorder_2D': ['aa_%_ord', 'aa_%_dis'],
                          'disorder_3D': ['aa_%_ord', 'aa_%_dis'],
                          'ss_disorder_2D': ['aa_%_ord', 'aa_%_dis'],
                          'ss_disorder_3D': ['aa_%_ord', 'aa_%_dis'],
                          'acc_2D': ['aa_%_M', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg', 'aa_%_Y'],
                          'acc_3D': ['aa_%_M', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg', 'aa_%_Y'],
                          'surface_3D': ['aa_%_M', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg', 'aa_%_Y'],
                          'metal_2_5D': ['aa_%_M', 'aa_%_bulk', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg'],
                          'metal_3D': ['aa_%_M', 'aa_%_bulk', 'aa_%_C', 'aa_%_chrg', 'aa_%_carb', 'aa_%_poschrg', 'aa_%_negchrg'],
                          'tm_2D': ['aa_%_M', 'aa_%_tmstab', 'aa_%_tmunstab'],
                          'tm_3D': ['aa_%_M', 'aa_%_tmstab', 'aa_%_tmunstab'],
                          'dna_2_5D': ['aa_%_dis', 'aa_%_ord'],
                          'csa_2_5D': ['aa_%_chrg', 'aa_%_M', 'aa_%_C'],
                          'sites_2_5D': ['aa_%_chrg', 'aa_%_M', 'aa_%_C']}

# filter_subseq_3D_suffixes = {'disorder_3D': ['aa_count_ord', 'aa_count_dis'],
#                           'ss_disorder_3D': ['aa_count_ord', 'aa_count_dis'],
#                           'acc_3D': ['aa_count_M', 'aa_count_C', 'aa_count_chrg', 'aa_count_carb', 'aa_count_poschrg', 'aa_count_negchrg'],
#                           'surface_3D': ['aa_count_M', 'aa_count_C', 'aa_count_chrg', 'aa_count_carb', 'aa_count_poschrg', 'aa_count_negchrg'],
#                           'metal_3D': ['aa_count_M', 'aa_count_bulk', 'aa_count_C', 'aa_count_chrg', 'aa_count_carb', 'aa_count_poschrg', 'aa_count_negchrg'],
#                           'tm_3D': ['aa_count_M', 'aa_count_tmstab', 'aa_count_tmunstab'],
#                           'dna_2_5D': ['aa_count_dis', 'aa_count_ord'],
#                           'csa_2_5D': ['aa_count_chrg', 'aa_count_M', 'aa_count_C'],
#                           'sites_2_5D': ['aa_count_chrg', 'aa_count_M', 'aa_count_C']}
#
# filter_subseq_2D_suffixes = {'disorder_2D': ['aa_count_ord', 'aa_count_dis'],
#                           'ss_disorder_2D': ['aa_count_ord', 'aa_count_dis'],
#                           'acc_2D': ['aa_count_M', 'aa_count_C', 'aa_count_chrg', 'aa_count_carb', 'aa_count_poschrg', 'aa_count_negchrg'],
#                           'metal_2_5D': ['aa_count_M', 'aa_count_bulk', 'aa_count_C', 'aa_count_chrg', 'aa_count_carb', 'aa_count_poschrg', 'aa_count_negchrg'],
#                           'tm_2D': ['aa_count_M', 'aa_count_tmstab', 'aa_count_tmunstab'],
#                           'dna_2_5D': ['aa_count_dis', 'aa_count_ord'],
#                           'csa_2_5D': ['aa_count_chrg', 'aa_count_M', 'aa_count_C'],
#                           'sites_2_5D': ['aa_count_chrg', 'aa_count_M', 'aa_count_C']}


def get_protein_feather_paths(protgroup, memornot, protgroup_dict, protein_feathers_dir, core_only_genes=None):
    """
    protgroup example: ('subsystem', 'cog_primary', 'H')
    memornot example: ('vizrecon', 'membrane')
    protgroup_dict example: {'databases': {'redoxdb': {'experimental_sensitive_cys': ['b2518','b3352','b2195','b4016'], ...}}}
    """
    prots_memornot = protgroup_dict['localization'][memornot[0]][memornot[1]]

    if protgroup[0] == 'localization':
        if protgroup[2] != 'all':
            if memornot[1] in ['membrane', 'inner_membrane', 'outer_membrane'] and protgroup[2] not in ['membrane', 'inner_membrane', 'outer_membrane']:
                return []
            if memornot[1] not in ['membrane', 'inner_membrane', 'outer_membrane'] and protgroup[2] in ['membrane', 'inner_membrane', 'outer_membrane']:
                return []

    prots_group = protgroup_dict[protgroup[0]][protgroup[1]][protgroup[2]]
    prots_filtered = list(set(prots_group).intersection(prots_memornot))
    if core_only_genes:
        prots_filtered = list(set(prots_filtered).intersection(core_only_genes))
    return [op.join(protein_feathers_dir, '{}_protein_strain_properties.fthr'.format(x)) for x in prots_filtered if op.exists(op.join(protein_feathers_dir, '{}_protein_strain_properties.fthr'.format(x)))]


def get_observed_strains_and_df(observation, observation_dict):
    """
    observation example: 'ros_simulated'
    observation_dict example: {'ros_simulated': [['NT12204_755', 'wt'], ['NT12120_270', 'wt'], ...] ...}
    """
    observed_df = pd.DataFrame.from_records(observation_dict[observation], columns=['strain','phenotype']).set_index('strain')
    observed_strains = observed_df.index.tolist()
    return observed_strains, observed_df


def get_interested_subsequences(subsequences):#, subseqdim):
    """
    subsequences example: ['disorder_2D', 'ss_disorder_2D', ...]
    filter_subseq_suffixes example: {'disorder_2D': ['aa_count_ord', 'aa_count_dis'], ... } -- defined above
    """
    keep_indices = []
    # if subseqdim == 'ALLD':
    #     filter_subseq_suffixes = filter_subseq_suffixes
    # elif subseqdim == '3D':
    #     filter_subseq_suffixes = filter_subseq_3D_suffixes
    # elif subseqdim == '2D':
    #     filter_subseq_suffixes = filter_subseq_2D_suffixes
    # else:
    #     raise ValueError('ALLD, 3D, or 2D only!')

    for subseq in subsequences:
        if subseq == 'all':
            keep_indices.extend([x for x in filter_subseq_suffixes[subseq]])
        else:
            keep_indices.extend([subseq + '_' + x for x in filter_subseq_suffixes[subseq]])
    return keep_indices


def load_feather(protein_feather, length_filter_pid=None, copynum_scale=False, copynum_df=None):
    """Load a feather of amino acid counts for a protein.

    Args:
        protein_feather (str): path to feather file
        copynum_scale (bool): if counts should be multiplied by protein copy number
        copynum_df (DataFrame): DataFrame of copy numbers

    Returns:
        DataFrame: of counts with some aggregated together

    """
    protein_df = pd.read_feather(protein_feather).set_index('index')

    # Combine counts for residue groups
    from ssbio.protein.sequence.properties.residues import _aa_property_dict_one, EXTENDED_AA_PROPERTY_DICT_ONE
    aggregators = {
        'aa_count_bulk'    : {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['Bulky'],
                              'subseqs' : ['metal_2_5D', 'metal_3D']},
        'aa_count_carb'    : {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['Carbonylation susceptible'],
                              'subseqs' : ['metal_2_5D', 'metal_3D', 'acc_2D', 'acc_3D', 'surface_3D']},
        'aa_count_chrg'    : {'residues': _aa_property_dict_one['Charged'],
                              'subseqs' : ['metal_2_5D', 'metal_3D', 'csa_2_5D', 'sites_2_5D', 'acc_2D', 'acc_3D',
                                           'surface_3D']},
        'aa_count_poschrg' : {'residues': _aa_property_dict_one['Basic'],
                              'subseqs' : ['metal_2_5D', 'metal_3D', 'acc_2D', 'acc_3D', 'surface_3D']},
        'aa_count_negchrg' : {'residues': _aa_property_dict_one['Acidic'],
                              'subseqs' : ['metal_2_5D', 'metal_3D', 'acc_2D', 'acc_3D', 'surface_3D']},
        'aa_count_tmstab'  : {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['TM stabilizing'],
                              'subseqs' : ['tm_2D', 'tm_3D']},
        'aa_count_tmunstab': {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['TM to Thr stabilizing'],
                              'subseqs' : ['tm_2D', 'tm_3D']},
        'aa_count_dis'     : {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['Disorder promoting'],
                              'subseqs' : ['disorder_2D', 'ss_disorder_2D', 'disorder_3D', 'ss_disorder_3D',
                                           'dna_2_5D']},
        'aa_count_ord'     : {'residues': EXTENDED_AA_PROPERTY_DICT_ONE['Order promoting'],
                              'subseqs' : ['disorder_2D', 'ss_disorder_2D', 'disorder_3D', 'ss_disorder_3D',
                                           'dna_2_5D']}}

    # Do combination counts for all types of subsequences
    for suffix, info in aggregators.items():
        agg_residues = info['residues']
        for prefix in info['subseqs']:
            to_add_idxes = []
            for agg_res in agg_residues:
                to_add_idx = prefix + '_aa_count_' + agg_res
                if to_add_idx in protein_df.index:
                    to_add_idxes.append(to_add_idx)
            subseq_agged_col = protein_df.loc[to_add_idxes, :].sum()  # Add each residue series
            protein_df.loc[prefix + '_' + suffix] = subseq_agged_col  # Append to df

    ## REMOVE OTHER STRAINS WITH DELETIONS (use float -- length_filter_pid=0.8 to get only strains with >80% length
    ## alternative to atlas2.calculate_residue_counts_perstrain wt_pid_cutoff param -- works a little differently just considering length
    if length_filter_pid:
        keep_cols = protein_df.loc['aa_count_total'][protein_df.loc['aa_count_total'] > protein_df.at['aa_count_total', 'K12'] * length_filter_pid].index
        protein_df = protein_df[keep_cols]

    # Multiply by proteomics copy number?
    if copynum_scale:
        if not isinstance(copynum_df, pd.DataFrame):
            raise ValueError('Please supply copy numbers')
        protein_id = op.basename(protein_feather).split('_protein')[0]
        if protein_id in copynum_df.index:
            copynum = copynum_df.at[protein_id, 'copynum']
            if copynum > 0:  # TODO: currently keeping one copy of proteins with 0, is that ok?
                protein_df = protein_df * copynum

    return protein_df


def get_proteome_counts_simple(prots_filtered_feathers, outpath, length_filter_pid=None,
                               copynum_scale=False, copynum_df=None,
                               force_rerun=False):
    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
        big_strain_counts_df = pd.DataFrame()
        first = True
        for feather in prots_filtered_feathers:
            loaded = load_feather(protein_feather=feather, length_filter_pid=length_filter_pid,
                                  copynum_scale=copynum_scale,
                                  copynum_df=copynum_df)
            if first:
                big_strain_counts_df = pd.DataFrame(index=loaded.index, columns=loaded.columns)
                first = False
            big_strain_counts_df = big_strain_counts_df.add(loaded, fill_value=0)
        if len(big_strain_counts_df) > 0:
            big_strain_counts_df.astype(float).reset_index().to_feather(outpath)
        return big_strain_counts_df
    else:
        return pd.read_feather(outpath).set_index('index')


# def get_proteome_counts_simple_sc(sc, prots_filtered_feathers, outpath, length_filter_pid=None,
#                                     copynum_scale=False, copynum_df=None,
#                                     force_rerun=False):
#     import ssbio.utils
#     if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
#         protein_feathers_final_rdd = sc.parallelize(prots_filtered_feathers)
#         mapper = protein_feathers_final_rdd.map(lambda x: load_feather(protein_feather=x, length_filter_pid=None,
#                                                                        copynum_scale=copynum_scale,
#                                                                        copynum_df=copynum_df))
#         big_strain_counts_df = mapper.reduce(lambda df1, df2: df1.add(df2, fill_value=0))
#         big_strain_counts_df.astype(float).reset_index().to_feather(outpath)
#         return big_strain_counts_df
#     else:
#         return pd.read_feather(outpath).set_index('index')


def get_proteome_counts_impute_missing(prots_filtered_feathers, outpath, length_filter_pid=None,
                                       copynum_scale=False, copynum_df=None,
                                       force_rerun=False):
    """Get counts, uses the mean feature vector to fill in missing proteins for a strain"""
    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
        big_strain_counts_df = pd.DataFrame()
        first = True
        for feather in prots_filtered_feathers:
            loaded = load_feather(protein_feather=feather, length_filter_pid=length_filter_pid,
                                  copynum_scale=copynum_scale,
                                  copynum_df=copynum_df)
            if first:
                big_strain_counts_df = pd.DataFrame(index=_all_counts, columns=loaded.columns)
                first = False

            new_columns = list(set(loaded.columns.tolist()).difference(big_strain_counts_df.columns))
            if new_columns:
                for col in new_columns:
                    big_strain_counts_df[col] = big_strain_counts_df.mean(axis=1)

            not_in_loaded = list(set(big_strain_counts_df.columns).difference(loaded.columns.tolist()))
            if not_in_loaded:
                for col in not_in_loaded:
                    big_strain_counts_df[col] = big_strain_counts_df[col] + loaded.mean(axis=1)

            big_strain_counts_df = big_strain_counts_df.add(loaded, fill_value=0)

        if len(big_strain_counts_df) > 0:
            big_strain_counts_df.astype(float).reset_index().to_feather(outpath)
        return big_strain_counts_df
    else:
        return pd.read_feather(outpath).set_index('index')


def get_proteome_percentages(counts_df, outpath, force_rerun=False):
    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
        big_strain_percents_df = pd.DataFrame(columns=counts_df.columns)
        for strain in counts_df.columns:
            totals = list(filter(lambda x: x.endswith('total'), counts_df[strain].index))
            for t in totals:
                counts = t.rsplit('_', 1)[0]
                aa_counts = list(filter(lambda x: (x.startswith(counts) and x not in totals), counts_df[strain].index))
                for aa_count in aa_counts:
                    big_strain_percents_df.at[aa_count.replace('count', '%'), strain] = counts_df[strain][aa_count]/counts_df[strain][t]

        big_strain_percents_df.astype(float).reset_index().to_feather(outpath)
    else:
        big_strain_percents_df = pd.read_feather(outpath).set_index('index')

    big_strain_percents_df.index.name = None
    return big_strain_percents_df


def get_proteome_correct_percentages(prots_filtered_feathers, outpath, length_filter_pid=None,
                                     copynum_scale=False, copynum_df=None,
                                     force_rerun=False):
    """Get counts and normalize by number of proteins, providing percentages"""
    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
        prot_tracker = defaultdict(int)
        big_strain_counts_df = pd.DataFrame()
        first = True
        for feather in prots_filtered_feathers:
            loaded = load_feather(protein_feather=feather, length_filter_pid=length_filter_pid,
                                  copynum_scale=copynum_scale,
                                  copynum_df=copynum_df)

            if first:
                big_strain_counts_df = pd.DataFrame(columns=loaded.columns)
                first = False
            tmp_df = pd.DataFrame(columns=loaded.columns)
            for strain in loaded.columns:
                prot_tracker[strain] += 1
                totals = list(filter(lambda x: x.endswith('total'), loaded[strain].index))
                for t in totals:
                    counts = t.rsplit('_', 1)[0]
                    aa_counts = list(
                        filter(lambda x: (x.startswith(counts) and x not in totals), loaded[strain].index))
                    for aa_count in aa_counts:
                        tmp_df.at[aa_count.replace('count', '%'), strain] = loaded[strain][aa_count] / \
                                                                            loaded[strain][t]
            big_strain_counts_df = big_strain_counts_df.add(tmp_df, fill_value=0)

        for c, total in prot_tracker.items():
            big_strain_counts_df.loc[:, c] /= total

        if len(big_strain_counts_df) > 0:
            big_strain_counts_df.astype(float).reset_index().to_feather(outpath)
        return big_strain_counts_df
    else:
        return pd.read_feather(outpath).set_index('index')


def remove_correlated_feats(df):
    tmp = df.T
    # Remove columns with no variation
    nunique = tmp.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    tmp.drop(cols_to_drop, axis=1, inplace=True)

    perc_spearman = scipy.stats.spearmanr(tmp)
    abs_corr = np.subtract(np.ones(shape=perc_spearman.correlation.shape),
                           np.absolute(perc_spearman.correlation))
    np.fill_diagonal(abs_corr, 0)
    abs_corr_clean = np.maximum(abs_corr,
                                abs_corr.transpose())  # some floating point mismatches, just make symmetric
    clustering = linkage(squareform(abs_corr_clean), method='average')
    clusters = fcluster(clustering, .1, criterion='distance')
    names = tmp.columns.tolist()
    names_to_cluster = list(zip(names, clusters))
    indices_to_keep = []
    ### Extract models closest to cluster centroids
    for x in range(1, len(set(clusters)) + 1):
        # Create mask from the list of assignments for extracting submatrix of the cluster
        mask = np.array([1 if i == x else 0 for i in clusters], dtype=bool)

        # Take the index of the column with the smallest sum of distances from the submatrix
        idx = np.argmin(sum(abs_corr_clean[:, mask][mask, :]))

        # Extract names of cluster elements from names_to_cluster
        sublist = [name for (name, cluster) in names_to_cluster if cluster == x]

        # Element closest to centroid
        centroid = sublist[idx]
        indices_to_keep.append(centroid)

    return df.loc[df.index.isin(indices_to_keep)]


def get_simple_sigdict(prots_filtered_feathers, subsequences, observation, observation_dict,
                       remove_corr=True, force_rerun=False):
    sigdict = {'less': defaultdict(list),
               'greater': defaultdict(list)}

    for p in prots_filtered_feathers:

        p_id = op.basename(p).split('_')[0]
        outpath = op.join(op.dirname(p), '{}_protein_strain_properties_percfilt.fthr'.format(p_id))

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=outpath):
            p_df = load_feather(protein_feather=p, length_filter_pid=0.8)
            p_perc_df = get_proteome_percentages(counts_df=p_df,
                                                 outpath=outpath,
                                                 force_rerun=force_rerun)
        else:
            p_perc_df = pd.read_feather(outpath).set_index('index')

        # Clean data first
        keep_subsequences = get_interested_subsequences(subsequences=subsequences)
        p_perc_df2 = p_perc_df.loc[p_perc_df.index.isin(keep_subsequences)]
        p_perc_df2 = p_perc_df2.astype(float).fillna(0)
        p_perc_df2 = p_perc_df2.loc[(p_perc_df2 > 0).any(axis=1)]

        if remove_corr:
            try:
                p_perc_df2 = remove_correlated_feats(p_perc_df2)
            except:
                continue
        all_features = p_perc_df2.index.tolist()

        # Add observations
        keep_strains, observations_df = get_observed_strains_and_df(observation=observation,
                                                                    observation_dict=observation_dict)
        feat_obs_df = p_perc_df2.T.join(observations_df, how='inner')

        # Split into 2 groups
        if observation == 'isolation_source':
            continue
        elif observation == 'ros_simulated':
            wt_df = feat_obs_df[feat_obs_df.phenotype == 'wt']
            mut_df = feat_obs_df[feat_obs_df.phenotype != 'wt']
        elif 'pathotype' in observation:
            wt_df = feat_obs_df[feat_obs_df.phenotype == 'Other']
            mut_df = feat_obs_df[feat_obs_df.phenotype != 'Other']
        else:
            wt_df = feat_obs_df[feat_obs_df.phenotype == 'No growth']
            mut_df = feat_obs_df[feat_obs_df.phenotype == 'Growth']

        if len(wt_df) == 0 or len(mut_df) == 0:
            continue

        # Mann-whitney test of each feature
        for alt in ['less', 'greater']:
            for feat in all_features:
                try:
                    z_stat, p_val = mannwhitneyu(wt_df[feat], mut_df[feat], alternative=alt)
                except ValueError:
                    continue
                if p_val < 0.01:
                    # Store differentiating features for this protein in a dictionary
                    sigdict[alt][p_id].append((feat, p_val))

    return sigdict


from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np
import scipy.stats
from scipy.spatial.distance import squareform
from sklearn import preprocessing
from sklearn.decomposition import PCA
from scipy.stats import ks_2samp
import itertools
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
from itertools import combinations
from sklearn.metrics.pairwise import euclidean_distances  # http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.euclidean_distances.html
import scipy.stats.stats as st
from scipy.stats import mannwhitneyu
from numpy.random import choice

plt.ioff()  # Turn interactive plotting off
sns.set(rc={'figure.figsize': (16, 8.5)})
sns.set_context('talk')
sns.set_style('ticks')
sns.set_palette('Set2')


class PCAMultiROS():
    def __init__(self, features_df, observations_df, plot_title, observation_colname='phenotype'):
        self.features_df = features_df
        self.observations_df = observations_df
        self.observation_colname = observation_colname
        self.plot_title = plot_title

        self.num_components = None
        self.pca = None
        self.pc_names_list = None
        self.pc_names_dict = None
        self.principal_df = None
        self.principal_observations_df = None
        self.markers = None

    def clean_data(self, keep_features=None, remove_correlated_feats=True):
        self.features_df = self.features_df.astype(float).fillna(0)
        self.features_df = self.features_df.loc[(self.features_df > 0).any(axis=1)]

        if keep_features:
            self.features_df = self.features_df.loc[self.features_df.index.isin(keep_features)]

        if remove_correlated_feats:
            tmp = self.features_df.T

            # Remove columns with no variation
            nunique = tmp.apply(pd.Series.nunique)
            cols_to_drop = nunique[nunique == 1].index
            tmp.drop(cols_to_drop, axis=1, inplace=True)

            perc_spearman = scipy.stats.spearmanr(tmp)
            abs_corr = np.subtract(np.ones(shape=perc_spearman.correlation.shape),
                                   np.absolute(perc_spearman.correlation))
            np.fill_diagonal(abs_corr, 0)
            abs_corr_clean = np.maximum(abs_corr,
                                        abs_corr.transpose())  # some floating point mismatches, just make symmetric
            clustering = linkage(squareform(abs_corr_clean), method='average')
            clusters = fcluster(clustering, .1, criterion='distance')
            names = tmp.columns.tolist()
            names_to_cluster = list(zip(names, clusters))
            indices_to_keep = []
            ### Extract models closest to cluster centroids
            for x in range(1, len(set(clusters)) + 1):
                # Create mask from the list of assignments for extracting submatrix of the cluster
                mask = np.array([1 if i == x else 0 for i in clusters], dtype=bool)

                # Take the index of the column with the smallest sum of distances from the submatrix
                idx = np.argmin(sum(abs_corr_clean[:, mask][mask, :]))

                # Extract names of cluster elements from names_to_cluster
                sublist = [name for (name, cluster) in names_to_cluster if cluster == x]

                # Element closest to centroid
                centroid = sublist[idx]
                indices_to_keep.append(centroid)

            self.features_df = self.features_df.loc[self.features_df.index.isin(indices_to_keep)]

    def run_pca(self, whiten=True):
        # Normalize
        for_pca_df = self.features_df.T
        for_pca_df_scaled = pd.DataFrame(preprocessing.scale(for_pca_df), columns=for_pca_df.columns)

        # Run PCA
        self.num_components = min(len(for_pca_df.T.columns), len(for_pca_df.T.index))
        pca = PCA(n_components=self.num_components, whiten=whiten)
        pca_fit = pca.fit_transform(for_pca_df_scaled)
        self.pc_names_list = ['PC{} ({:.0%})'.format(x + 1, pca.explained_variance_ratio_[x]) for x in
                                  range(self.num_components)]
        self.pc_names_dict = {k.split(' ')[0]: k for k in self.pc_names_list}
        principal_df = pd.DataFrame(data=pca_fit, columns=self.pc_names_list, index=for_pca_df.index)
        principal_df.index.name = 'strain'

        self.principal_df = principal_df
        self.pca = pca
        # self.principal_observations_df = self.principal_df.join(self.observations_df, how='inner')
        #
        # # Make iterable list of markers
        # mks = itertools.cycle(["<", "+", "o", 'D', 'x', '^', '*', '8', 's', 'p', 'v', 'X', '_', 'h'])
        # self.markers = [next(mks) for i in range(len(self.principal_observations_df[self.observation_colname].unique()))]

    def make_biplot(self, pc_x=1, pc_y=2, outpath=None, dpi=150, custom_markers=None, custom_order=None):
        if not custom_order:
            custom_order = sorted(self.observations_df[self.observation_colname].unique().tolist())
        if not custom_markers:
            custom_markers = self.markers
        plot = sns.lmplot(data=self.principal_observations_df,
                              x=self.principal_observations_df.columns[pc_x - 1],
                              y=self.principal_observations_df.columns[pc_y - 1],
                              hue=self.observation_colname,
                              hue_order=custom_order,
                              fit_reg=False,
                              size=6,
                              markers=custom_markers,
                            scatter_kws={'alpha': 0.5})
        plot = (plot.set(title='PC{} vs. PC{}'.format(pc_x, pc_y)))
        if outpath:
            plot.savefig(outpath, dpi=dpi)
        else:
            plt.show()
        plt.close()

    def make_pairplot(self, num_components_to_plot=4, outpath=None, dpi=150):
        # Get columns
        components_to_plot = [self.principal_observations_df.columns[x] for x in range(num_components_to_plot)]

        # Plot
        plot = sns.pairplot(data=self.principal_observations_df, hue=self.observation_colname,
                                vars=components_to_plot, markers=self.markers, size=4)
        plt.subplots_adjust(top=.95)
        plt.suptitle(self.plot_title)

        if outpath:
            plot.fig.savefig(outpath, dpi=dpi)
        else:
            plt.show()
        plt.close()

    def make_3Dplot(self, outpath=None, dpi=150):
        import mpl_toolkits.mplot3d
        fig = plt.figure(1, figsize=(8, 6))
        ax = mpl_toolkits.mplot3d.Axes3D(fig)
        for key, group in self.principal_observations_df.groupby(self.observation_colname):
            ax.plot(group[self.principal_observations_df.columns[0]],
                    group[self.principal_observations_df.columns[1]],
                    group[self.principal_observations_df.columns[2]],
                    'o', alpha=0.5, label=key)

        # Make simple, bare axis lines through space:
        xAxisLine = ((min(self.principal_observations_df[self.principal_observations_df.columns[0]]),
                      max(self.principal_observations_df[self.principal_observations_df.columns[0]])), (0, 0), (0, 0))
        ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r', alpha=0.4)
        yAxisLine = ((0, 0), (min(self.principal_observations_df[self.principal_observations_df.columns[1]]),
                              max(self.principal_observations_df[self.principal_observations_df.columns[1]])), (0, 0))
        ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r', alpha=0.4)
        zAxisLine = ((0, 0), (0, 0), (min(self.principal_observations_df[self.principal_observations_df.columns[2]]),
                                      max(self.principal_observations_df[self.principal_observations_df.columns[2]])))
        ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r', alpha=0.4)

        ax.set_title("PC1 vs. PC2 vs. PC3")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_zlabel("PC3")
        ax.legend(loc='upper left', fontsize=15)

        if outpath:
            fig.savefig(outpath, dpi=dpi)
        else:
            plt.show()
        plt.close()

    def make_contribplot(self, pc_to_look_at=1, sigadder=0.01, outpath=None, dpi=150, return_top_contribs=False):
        """Make a plot showing contributions of properties to a PC"""
        cont = pd.DataFrame(self.pca.components_, columns=self.features_df.index, index=self.pc_names_list)
        tmp_df = pd.DataFrame(cont.iloc[pc_to_look_at - 1]).reset_index().rename(columns={'index': 'Property'})
        tmp_df['Contribution'] = tmp_df.iloc[:, 1] ** 2
        tmp_df = tmp_df[tmp_df['Contribution'] > 1 / len(
                cont.iloc[0]) + sigadder]  # Alter sigadder to just plot more/less significant contributors
        tmp_df['Sign'] = np.where(tmp_df.iloc[:, 1] >= 0, 'Positive', 'Negative')
        tmp_df = tmp_df.sort_values(by='Contribution', ascending=False)

        fig, ax = plt.subplots(figsize=(30, 10))
        sns.barplot(data=tmp_df, y='Property', x='Contribution', hue='Sign', dodge=False, ax=ax, hue_order=['Positive', 'Negative'],
                        palette=sns.color_palette("coolwarm", 2))

        # Random formatting crap
        self._change_height(ax, .6)  # Make bars thinner
        ax.set_title('{} contributors'.format(self.pc_names_list[pc_to_look_at - 1]))
        legend = plt.legend(loc=8, bbox_to_anchor=(1.2, .8), ncol=1, title='Sign', fontsize=10)
        plt.setp(legend.get_title(), fontsize=12)
        plt.gcf().subplots_adjust(left=.5, right=.65)
        if outpath:
            fig.savefig(outpath, dpi=dpi)
        else:
            plt.show()
        plt.close()

        if return_top_contribs:
            return tmp_df.Property.values.tolist()

    def _change_height(self, ax, new_value):
        """Make bars in horizontal bar chart thinner"""
        for patch in ax.patches:
            current_height = patch.get_height()
            diff = current_height - new_value

            # we change the bar height
            patch.set_height(new_value)

            # we recenter the bar
            patch.set_y(patch.get_y() + diff * .5)

    def get_pca_ks_stats(self, maxrange=5):
        """Get a dictionary of PC#: K-S test stat for each """
        pc_to_phenotype_pairs = {}
        num_components = self.principal_observations_df.shape[1]
        if num_components < maxrange:
            maxrange = num_components

        phenotypes = self.principal_observations_df.phenotype.unique().tolist()
        for i in range(0, maxrange):
            phenotype_pair_to_ks = {}
            for p1, p2 in combinations(phenotypes, 2):
                p1_pc = self.principal_observations_df[self.principal_observations_df.phenotype == p1].iloc[:,i].as_matrix()
                p2_pc = self.principal_observations_df[self.principal_observations_df.phenotype == p2].iloc[:,i].as_matrix()
                phenotype_pair_to_ks[(p1, p2)] = ks_2samp(p1_pc, p2_pc)
            pc_to_phenotype_pairs[i + 1] = phenotype_pair_to_ks

        return pc_to_phenotype_pairs


def get_intra_inter_distances(feat_df, obs_df, normalize=False, plot=False):
    # Drop zero rows, transpose, set type to float, and fill any missing values
    feat_df = feat_df.loc[(feat_df != 0).any(axis=1)].T.astype(float).fillna(0)
    if normalize:
        feat_df = pd.DataFrame(preprocessing.scale(feat_df), columns=feat_df.columns, index=feat_df.index)
    feat_obs_df = feat_df.join(obs_df, how='inner')

    obs_intra_distances = {}
    obs_inter_distances = {}

    obs_feat_vectors = {}

    for phenotype in feat_obs_df.phenotype.unique():
        feat_vectors = feat_obs_df[feat_obs_df.phenotype == phenotype].drop(columns=feat_obs_df.columns[-1]).as_matrix()
        obs_feat_vectors[phenotype] = feat_vectors

        # Intra-distances
        intra_distances_calc = pdist(feat_vectors)
        obs_intra_distances[phenotype] = intra_distances_calc

    # Inter-distances
    # Randomly sample from 1...?  uncomment next 2 lines and use obs_feat_vector instead of obs_feat_vector_choice
    # if you want to try that
    #     obs_shortest = min([x.shape[0] for x in obs_feat_vectors.values()])
    #     obs_feat_vector_choice = {k:(v[choice(v.shape[0], obs_shortest, replace=False), :] if v.shape[0]!=obs_shortest else v) for k,v in obs_feat_vectors.items()}
    for inter1, inter2 in combinations(obs_feat_vectors, 2):
        obs_inter_distances[(inter1, inter2)] = euclidean_distances(obs_feat_vectors[inter1],
                                                                    obs_feat_vectors[inter2])

    if plot:
        df = pd.DataFrame()
        for k, v in obs_intra_distances.items():
            ser = pd.Series(v)
            df[k] = ser

        for k, v in obs_inter_distances.items():
            ser = pd.Series(v.flatten())
            df[str(k)] = ser

        plotter = pd.melt(df, value_vars=df.columns.tolist())
        plotter = plotter[pd.notnull(plotter.value)]

        sns.violinplot(data=plotter, x='variable', y='value')

    return obs_intra_distances, obs_inter_distances


def compute_skew_stats(intra, inter):
    """Returns two dictionaries reporting (skew, skew_pval) for all groups"""
    # Intra (within a group) stats
    intra_skew = {}
    for k, v in intra.items():
        skew = st.skew(v)
        try:
            skew_zstat, skew_pval = st.skewtest(v)
        except ValueError:  # if sample size too small
            skew_zstat, skew_pval = (0, 1)
        intra_skew[k] = (skew, skew_zstat, skew_pval)

    # Inter (between groups) stats
    inter_skew = {}
    for k, v in inter.items():
        # Inter skew stats
        skew_sep = st.skew(v.flatten())
        try:
            skew_sep_zstat, skew_sep_pval = st.skewtest(v.flatten())
        except ValueError:
            skew_sep_zstat, skew_sep_pval = (0, 1)
        inter_skew['-'.join(k)] = (skew_sep, skew_sep_zstat, skew_sep_pval)

        # Significance of difference between intra and inter distributions
        for intra_key in k:
            try:
                separation_zstat, separation_pval = mannwhitneyu(intra[intra_key],
                                                                 v.flatten(),
                                                                 alternative='less')
            except ValueError:  # All numbers are identical in mannwhitneyu
                separation_zstat, separation_pval = (0, 1)
            inter_skew['{}<{}'.format(intra_key, '-'.join(k))] = (separation_zstat, separation_pval)

    return intra_skew, inter_skew


def run_all(protgroup, memornot, subsequences, observation, proteomescale, base_outdir,
            protgroup_dict, protein_feathers_dir, observation_dict, copynum_df, date,
            errfile, subseqdim, statfile_duo, statfile_trio, cutoff_num_proteins=0, core_only_genes=None,
            impute_counts=True, length_filter_pid=.8,
            force_rerun_counts=False,
            sc=None, run_percentages=True, force_rerun_percentages=False,
            run_pca_and_stats=True, remove_correlated_feats=True, save_plots=False, force_rerun_pca=False):
    import ssbio.utils

    # Need to set multiprocessing limit for scipy/numpy stuff if parallelizing anything
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    # First, filter down the protein group to the membrane/nonmembrane definition
    prots_filtered_feathers = get_protein_feather_paths(protgroup=protgroup, memornot=memornot,
                                                        protgroup_dict=protgroup_dict,
                                                        protein_feathers_dir=protein_feathers_dir,
                                                        core_only_genes=core_only_genes)
    num_proteins = len(prots_filtered_feathers)
    if num_proteins <= cutoff_num_proteins:
        return

    # Make output directories
    if proteomescale:
        protscale = 'proteome_scaled'
    else:
        protscale = 'proteome_unscaled'
    outdir_d0 = ssbio.utils.make_dir(op.join(base_outdir, protscale))
    outdir_d1 = ssbio.utils.make_dir(op.join(outdir_d0, '-'.join(memornot)))
    outdir_final = ssbio.utils.make_dir(op.join(outdir_d1, '-'.join(protgroup)))
    outdir_observ = ssbio.utils.make_dir(op.join(outdir_final, observation))
    outdir_observ_subseqdim = ssbio.utils.make_dir(op.join(outdir_observ, subseqdim))

    # Then load the protein feathers and add them all together to represent a "proteome"
    # if sc:
    #     big_strain_counts_df = get_proteome_counts_sc(sc=sc, prots_filtered_feathers=prots_filtered_feathers,
    #                                                   outpath=op.join(outdir_final,
    #                                                                   '{}-subsequence_proteome.fthr'.format(date)),
    #                                                   copynum_scale=proteomescale, copynum_df=copynum_df,
    #                                                   force_rerun=force_rerun_counts)
    # else:
    #     big_strain_counts_df = get_proteome_counts(prots_filtered_feathers=prots_filtered_feathers,
    #                                                outpath=op.join(outdir_final,
    #                                                                '{}-subsequence_proteome.fthr'.format(date)),
    #                                                copynum_scale=proteomescale, copynum_df=copynum_df,
    #                                                force_rerun=force_rerun_counts)
    # if len(big_strain_counts_df) == 0:
    #     with open(errfile, "a") as myfile:
    #         myfile.write('COUNT ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
    #     return

    # Impute averages for missing counts
    if impute_counts:
        big_strain_counts_df = get_proteome_counts_impute_missing(prots_filtered_feathers=prots_filtered_feathers,
                                                                  outpath=op.join(outdir_final,
                                                                                  '{}-subsequence_proteome_IMP.fthr'.format(
                                                                                      date)),
                                                                  length_filter_pid=length_filter_pid,
                                                                  copynum_scale=proteomescale, copynum_df=copynum_df,
                                                                  force_rerun=force_rerun_counts)

        big_strain_percents_df = get_proteome_percentages(counts_df=big_strain_counts_df,
                                                          outpath=op.join(outdir_final,
                                                                          '{}-subsequence_proteome_perc_IMP.fthr'.format(
                                                                                  date)),
                                                          force_rerun=force_rerun_percentages)

    # Divide by totals to get percentages in a new dataframe
    else:
        try:
            big_strain_percents_df = get_proteome_correct_percentages(prots_filtered_feathers=prots_filtered_feathers,
                                                                      outpath=op.join(outdir_final,
                                                                                      '{}-subsequence_proteome_perc.fthr'.format(
                                                                                          date)),
                                                                      length_filter_pid=length_filter_pid,
                                                                      force_rerun=force_rerun_percentages)
        except:
            with open(errfile, "a") as myfile:
                myfile.write('PERCENTAGES ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(
                    observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            return

    # Stop here if only percentages desired
    if not run_pca_and_stats:
        return

    pca_pickle = op.join(outdir_observ_subseqdim, '{}-subsequence_pca.pckl'.format(date))
    if ssbio.utils.force_rerun(flag=force_rerun_pca, outfile=pca_pickle):

        # Then, get filters for the columns to available strains in an observation for percentage df
        keep_strains, observations_df = get_observed_strains_and_df(observation=observation,
                                                                    observation_dict=observation_dict)
        if keep_strains:
            big_strain_percents_df = big_strain_percents_df[big_strain_percents_df.columns[big_strain_percents_df.columns.isin(keep_strains)]]

        # Then, get filters for rows of the loaded feathers for interested subsequences
        keep_subsequences = get_interested_subsequences(subsequences=subsequences)

        # Some numbers: number of observations
        num_obs = observations_df.phenotype.value_counts().to_dict()
        if len(num_obs) < 2:  # If only one observation, what are we trying to compare?? nothing really
            return
        observations_string = ';'.join('{}:{}'.format(key, val) for key, val in num_obs.items())

        # Some numbers: number of features
        num_feats = len(big_strain_percents_df)

        # Make an unwieldy title
        big_title = 'LOC={0}; PROTGROUP={1}; PHENOTYPE={2}; PROTSCALE={3}; SUBSEQDIM={4};\n' \
                    'NUMPROTS={5}; NUMFEATS={6}; NUMSTRAINS={7}'.format('-'.join(memornot),
                                                                        '-'.join(protgroup),
                                                                        str(observation),
                                                                        str(proteomescale),
                                                                        subseqdim,
                                                                        num_proteins,
                                                                        num_feats,
                                                                        observations_string)

        # Run PCA and make plots
        runner = PCAMultiROS(features_df=big_strain_percents_df, observations_df=observations_df, plot_title=big_title)
        try:
            runner.clean_data(keep_features=keep_subsequences, remove_correlated_feats=remove_correlated_feats)
        except:
            with open(errfile, "a") as myfile:
                myfile.write(
                    'CLEAN ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            return
        try:
            runner.run_pca()
        except:
            with open(errfile, "a") as myfile:
                myfile.write(
                    'PCA ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            return

        if save_plots:
            try:
                runner.make_biplot(pc_x=1, pc_y=2,
                                   outpath=op.join(outdir_observ_subseqdim, '{}-subsequence_biplot_1_2.png'.format(date)))
            except:
                with open(errfile, "a") as myfile:
                    myfile.write(
                        'PCA BIPLOT ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            try:
                runner.make_pairplot(num_components_to_plot=4,
                                     outpath=op.join(outdir_observ_subseqdim, '{}-subsequence_pairplot.png'.format(date)))
            except:
                with open(errfile, "a") as myfile:
                    myfile.write(
                        'PCA PAIRPLOT ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            try:
                runner.make_3Dplot(outpath=op.join(outdir_observ_subseqdim, '{}-subsequence_3Dplot.png'.format(date)))
            except:
                with open(errfile, "a") as myfile:
                    myfile.write(
                        'PCA 3DPLOT ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
            try:
                runner.make_contribplot(pc_to_look_at=1, sigadder=0.01,
                                        outpath=op.join(outdir_observ_subseqdim, '{}-subsequence_contribplot.png'.format(date)))
            except:
                with open(errfile, "a") as myfile:
                    myfile.write(
                            'PCA CONTRIBPLOT ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(
                                observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
        with open(pca_pickle, 'wb') as f:
            pickle.dump(runner, f)
    else:
        with open(pca_pickle, 'rb') as f:
            runner = pickle.load(f)

    # Get stats
    try:
        ks_stats = runner.get_pca_ks_stats()
    except:
        with open(errfile, "a") as myfile:
            myfile.write(
                    'STAT K-S ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(
                            observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
        return

    try:
        # perc_intra, perc_inter = get_intra_inter_distances(feat_df=runner.features_df,
        #                                                    obs_df=runner.observations_df,
        #                                                    normalize=True,
        #                                                    plot=False)
        # perc_stats = compute_skew_stats(intra=perc_intra, inter=perc_inter)

        pca_intra, pca_inter = get_intra_inter_distances(feat_df=runner.principal_df.T,
                                                         obs_df=runner.observations_df,
                                                         normalize=False,
                                                         plot=False)
        skew_stats = compute_skew_stats(intra=pca_intra, inter=pca_inter)

    except:
        with open(errfile, "a") as myfile:
            myfile.write('STAT DIST ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(
                observation) + '\t' + str(proteomescale) + '\t' + subseqdim + "\n")
        return

    # with open(statfile, 'a') as statf:
    #     statf.write('-'.join(memornot) + '\t' + '-'.join(protgroup) + '\t' + str(observation) + '\t' + str(
    #         proteomescale) + '\t' + subseqdim + '\t' + str(num_proteins) + '\t' + str(num_feats) + '\t' + str(
    #         num_obs) + '\t' + str(perc_stats) + '\t' + str(skew_stats) + '\t' + str(ks_stats) + '\n')

    skew_pval_significant = True
    skew_worst_pval = -1
    ks_pval_significant = True
    ks_worst_pval = -1

    # Search for clusters not close to WT
    if observation == 'ros_simulated':
        for k, v in skew_stats[1].items():
            if 'wt' in k and '<' in k:
                pval = v[1]
                if pval > skew_worst_pval:
                    skew_worst_pval = pval
                if pval > 0.05:
                    skew_pval_significant = False
        for k, v in ks_stats[1].items():
            if 'wt' in k:
                if v.pvalue > ks_worst_pval:
                    ks_worst_pval = v.pvalue
                if v.pvalue > 0.05:
                    ks_pval_significant = False
        # Combination stats
        if ks_pval_significant and skew_pval_significant:
            if ks_worst_pval != -1 and skew_worst_pval != -1:
                with open(statfile_trio, 'a') as f:
                    f.write(str((protgroup, memornot, subsequences, observation, proteomescale)) + '\t' + str(ks_worst_pval) + '\t' + str(skew_worst_pval) + '\n')

    # Search for G/NG clusters and pathotype_simple ones
    else:
        for k, v in skew_stats[1].items():
            if '<' in k:
                pval = v[1]
                if pval > skew_worst_pval:
                    skew_worst_pval = pval
                if pval > 0.05:
                    skew_pval_significant = False
        for k, v in ks_stats[1].items():
            if v.pvalue > ks_worst_pval:
                ks_worst_pval = v.pvalue
            if v.pvalue > 0.05:
                ks_pval_significant = False

        # Combination stats
        if ks_pval_significant and skew_pval_significant:
            if ks_worst_pval != -1 and skew_worst_pval != -1:
                with open(statfile_duo, 'a') as f:
                    f.write(str((protgroup, memornot, subsequences, observation, proteomescale)) + '\t' + str(ks_worst_pval) + '\t' + str(skew_worst_pval) + '\n')


def run_all2(protgroup, memornot, subsequences, base_outdir,
            protgroup_dict, protein_feathers_dir, date, errfile, impute_counts=True,
            cutoff_num_proteins=0, core_only_genes=None,
            length_filter_pid=.8, remove_correlated_feats=True,
            force_rerun_counts=False, force_rerun_percentages=False, force_rerun_pca=False):
    """run_all but ignoring observations before pca"""
    import ssbio.utils

    # Need to set multiprocessing limit for scipy/numpy stuff if parallelizing anything
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    # First, filter down the protein group to the membrane/nonmembrane definition
    prots_filtered_feathers = get_protein_feather_paths(protgroup=protgroup, memornot=memornot,
                                                        protgroup_dict=protgroup_dict,
                                                        protein_feathers_dir=protein_feathers_dir,
                                                        core_only_genes=core_only_genes)
    num_proteins = len(prots_filtered_feathers)
    if num_proteins <= cutoff_num_proteins:
        return

    # Make output directories
    protscale = 'proteome_unscaled'
    outdir_d0 = ssbio.utils.make_dir(op.join(base_outdir, protscale))
    outdir_d1 = ssbio.utils.make_dir(op.join(outdir_d0, '-'.join(memornot)))
    outdir_final = ssbio.utils.make_dir(op.join(outdir_d1, '-'.join(protgroup)))

    if impute_counts:
        big_strain_counts_df = get_proteome_counts_impute_missing(prots_filtered_feathers=prots_filtered_feathers,
                                                                  outpath=op.join(outdir_final,
                                                                                  '{}-subsequence_proteome_IMP.fthr'.format(
                                                                                      date)),
                                                                  length_filter_pid=length_filter_pid,
                                                                  force_rerun=force_rerun_counts)

        big_strain_percents_df = get_proteome_percentages(counts_df=big_strain_counts_df,
                                                          outpath=op.join(outdir_final,
                                                                          '{}-subsequence_proteome_perc_IMP.fthr'.format(
                                                                                  date)),
                                                          force_rerun=force_rerun_percentages)
        pca_pickle = op.join(outdir_final, '{}-subsequence_pca.pckl'.format(date))
    # Divide by totals to get percentages in a new dataframe
    else:
        try:
            big_strain_percents_df = get_proteome_correct_percentages(prots_filtered_feathers=prots_filtered_feathers,
                                                                      outpath=op.join(outdir_final,
                                                                                      '{}-subsequence_proteome_perc_AVG.fthr'.format(
                                                                                              date)),
                                                                      length_filter_pid=length_filter_pid,
                                                                      force_rerun=force_rerun_percentages)
            pca_pickle = op.join(outdir_final, '{}-subsequence_pca_AVG.pckl'.format(date))
        except:
            with open(errfile, "a") as myfile:
                myfile.write('PERCENTAGES ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + "\n")
            return



    if ssbio.utils.force_rerun(flag=force_rerun_pca, outfile=pca_pickle):

        # Then, get filters for rows of the loaded feathers for interested subsequences
        keep_subsequences = get_interested_subsequences(subsequences=subsequences)

        # Some numbers: number of features
        num_feats = len(big_strain_percents_df)

        # Make an unwieldy title
        big_title = 'LOC={0}; PROTGROUP={1};\n' \
                    'NUMPROTS={2}; NUMFEATS={3}'.format('-'.join(memornot),
                                                                        '-'.join(protgroup),
                                                                        num_proteins,
                                                                        num_feats)

        # Run PCA and make plots
        runner = PCAMultiROS(features_df=big_strain_percents_df, observations_df=pd.DataFrame(), plot_title=big_title)
        try:
            runner.clean_data(keep_features=keep_subsequences, remove_correlated_feats=remove_correlated_feats)
        except:
            with open(errfile, "a") as myfile:
                myfile.write(
                    'CLEAN ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + "\n")
            return
        # try:
        runner.run_pca()
        # except:
        #     with open(errfile, "a") as myfile:
        #         myfile.write(
        #             'PCA ERR: ' + '-'.join(memornot) + '\t' + '-'.join(protgroup) + "\n")
        #     return
        with open(pca_pickle, 'wb') as f:
            pickle.dump(runner, f)
    else:
        with open(pca_pickle, 'rb') as f:
            runner = pickle.load(f)


def run_all_simple(protgroup, memornot, subsequences, observation, protgroup_dict, protein_feathers_dir, observation_dict, statfile, cutoff_num_proteins=0):

    # Need to set multiprocessing limit for scipy/numpy stuff if parallelizing anything
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    # First, filter down the protein group to the membrane/nonmembrane definition
    prots_filtered_feathers = get_protein_feather_paths(protgroup=protgroup, memornot=memornot,
                                                        protgroup_dict=protgroup_dict,
                                                        protein_feathers_dir=protein_feathers_dir)
    num_proteins = len(prots_filtered_feathers)
    if num_proteins <= cutoff_num_proteins:
        return

    # WT HAVE MORE OF THESE...
    sigdict = get_simple_sigdict(prots_filtered_feathers=prots_filtered_feathers,
                                        subsequences=subsequences,
                                        observation=observation,
                                        observation_dict=observation_dict)

    sig_predf = []
    for protein, feat_and_pvals in sigdict['greater'].items():
        for feat, pval in feat_and_pvals:
            sig_predf.append((protein, feat, pval))

    sig_df1 = pd.DataFrame(sig_predf, columns=['protein', 'feature', 'pval']).sort_values(by='pval', ascending=True)

    if len(sig_df1) > 0:
        num_prots1 = len(sig_df1.protein.unique())
        num_feats1 = len(sig_df1.feature.unique())

        # Very signficant features
        sig_df2 = sig_df1[sig_df1.pval <= sig_df1.pval.quantile(.25)]
        num_prots2 = len(sig_df2.protein.unique())
        num_feats2 = len(sig_df2.feature.unique())

        # Very significant features that are explainable
        sig_df3 = sig_df2[sig_df2.feature.isin(expected_in_wt)]
        num_prots3 = len(sig_df3.protein.unique())
        num_feats3 = len(sig_df3.feature.unique())

        wt_perc_sig_feats = num_feats3 / num_feats2
        wt_perc_sig_prots = num_prots3 / num_proteins
    else:
        wt_perc_sig_feats = 0
        wt_perc_sig_prots = 0


    # MUT HAVE MORE OF THESE...
    sig_predf = []
    for protein, feat_and_pvals in sigdict['less'].items():
        for feat, pval in feat_and_pvals:
            sig_predf.append((protein, feat, pval))

    sig_df1 = pd.DataFrame(sig_predf, columns=['protein', 'feature', 'pval']).sort_values(by='pval', ascending=True)

    if len(sig_df1) > 0:
        num_prots1 = len(sig_df1.protein.unique())
        num_feats1 = len(sig_df1.feature.unique())

        # Very signficant features
        sig_df2 = sig_df1[sig_df1.pval <= sig_df1.pval.quantile(.25)]
        num_prots2 = len(sig_df2.protein.unique())
        num_feats2 = len(sig_df2.feature.unique())

        # Very significant features that are explainable
        sig_df3 = sig_df2[sig_df2.feature.isin(expected_in_mut)]
        num_prots3 = len(sig_df3.protein.unique())
        num_feats3 = len(sig_df3.feature.unique())

        mut_perc_sig_feats = num_feats3 / num_feats2
        mut_perc_sig_prots = num_prots3 / num_proteins
    else:
        mut_perc_sig_feats = 0
        mut_perc_sig_prots = 0

    with open(statfile, 'a') as f:
        f.write(str((protgroup, memornot, observation)) + '\t' + str(wt_perc_sig_feats) + '\t' + str(wt_perc_sig_prots) + '\t' + str(mut_perc_sig_feats) + '\t' + str(mut_perc_sig_prots) + '\t' + str(num_proteins) + '\n')


expected_in_wt = ['aa_%_dis', 'disorder_2D_aa_%_dis', 'ss_disorder_2D_aa_%_dis', 'disorder_3D_aa_%_dis', 'ss_disorder_3D_aa_%_dis', 'dna_2_5D_aa_%_ord', 'aa_%_C', 'acc_2D_aa_%_C', 'acc_3D_aa_%_C', 'surface_3D_aa_%_C', 'metal_2_5D_aa_%_C', 'metal_3D_aa_%_C', 'csa_2_5D_aa_%_C', 'sites_2_5D_aa_%_C', 'aa_%_carb', 'acc_2D_aa_%_carb', 'acc_3D_aa_%_carb', 'surface_3D_aa_%_carb', 'metal_2_5D_aa_%_carb', 'metal_3D_aa_%_carb', 'aa_%_chrg', 'acc_2D_aa_%_chrg', 'acc_3D_aa_%_chrg', 'surface_3D_aa_%_chrg', 'metal_2_5D_aa_%_chrg', 'metal_3D_aa_%_chrg', 'csa_2_5D_aa_%_chrg', 'sites_2_5D_aa_%_chrg', 'acc_2D_aa_%_poschrg', 'acc_3D_aa_%_poschrg', 'surface_3D_aa_%_poschrg', 'metal_2_5D_aa_%_poschrg', 'metal_3D_aa_%_poschrg', 'acc_2D_aa_%_negchrg', 'acc_3D_aa_%_negchrg', 'surface_3D_aa_%_negchrg', 'tm_2D_aa_%_tmunstab', 'tm_3D_aa_%_tmunstab', 'surface_3D_aa_%_Y','aa_%_Y','acc_3D_aa_%_Y','acc_2D_aa_%_Y']
expected_in_mut = ['aa_%_ord', 'disorder_2D_aa_%_ord', 'ss_disorder_2D_aa_%_ord', 'disorder_3D_aa_%_ord', 'ss_disorder_3D_aa_%_ord', 'dna_2_5D_aa_%_dis', 'aa_%_M', 'acc_2D_aa_%_M', 'acc_3D_aa_%_M', 'surface_3D_aa_%_M', 'metal_2_5D_aa_%_M', 'metal_3D_aa_%_M', 'csa_2_5D_aa_%_M', 'sites_2_5D_aa_%_M', 'metal_2_5D_aa_%_negchrg', 'metal_3D_aa_%_negchrg', 'metal_2_5D_aa_%_bulk', 'metal_3D_aa_%_bulk', 'tm_2D_aa_%_M', 'tm_3D_aa_%_M', 'tm_2D_aa_%_tmstab', 'tm_3D_aa_%_tmstab']


_all_counts = ['aa_count_A', 'aa_count_C', 'aa_count_D', 'aa_count_E', 'aa_count_F', 'aa_count_G', 'aa_count_H',
               'aa_count_I', 'aa_count_K', 'aa_count_L', 'aa_count_M', 'aa_count_N', 'aa_count_P', 'aa_count_Q',
               'aa_count_R', 'aa_count_S', 'aa_count_T', 'aa_count_V', 'aa_count_W', 'aa_count_Y', 'aa_count_total',
               'acc_2D_aa_count_A', 'acc_2D_aa_count_C', 'acc_2D_aa_count_D', 'acc_2D_aa_count_E', 'acc_2D_aa_count_F',
               'acc_2D_aa_count_G', 'acc_2D_aa_count_H', 'acc_2D_aa_count_I', 'acc_2D_aa_count_K', 'acc_2D_aa_count_L',
               'acc_2D_aa_count_M', 'acc_2D_aa_count_N', 'acc_2D_aa_count_P', 'acc_2D_aa_count_Q', 'acc_2D_aa_count_R',
               'acc_2D_aa_count_S', 'acc_2D_aa_count_T', 'acc_2D_aa_count_V', 'acc_2D_aa_count_W', 'acc_2D_aa_count_Y',
               'acc_2D_aa_count_carb', 'acc_2D_aa_count_chrg', 'acc_2D_aa_count_negchrg', 'acc_2D_aa_count_poschrg',
               'acc_2D_aa_count_total', 'acc_3D_aa_count_A', 'acc_3D_aa_count_C', 'acc_3D_aa_count_D',
               'acc_3D_aa_count_E', 'acc_3D_aa_count_F', 'acc_3D_aa_count_G', 'acc_3D_aa_count_H', 'acc_3D_aa_count_I',
               'acc_3D_aa_count_K', 'acc_3D_aa_count_L', 'acc_3D_aa_count_M', 'acc_3D_aa_count_N', 'acc_3D_aa_count_P',
               'acc_3D_aa_count_Q', 'acc_3D_aa_count_R', 'acc_3D_aa_count_S', 'acc_3D_aa_count_T', 'acc_3D_aa_count_V',
               'acc_3D_aa_count_W', 'acc_3D_aa_count_Y', 'acc_3D_aa_count_carb', 'acc_3D_aa_count_chrg',
               'acc_3D_aa_count_negchrg', 'acc_3D_aa_count_poschrg', 'acc_3D_aa_count_total', 'csa_2_5D_aa_count_A',
               'csa_2_5D_aa_count_C', 'csa_2_5D_aa_count_D', 'csa_2_5D_aa_count_E', 'csa_2_5D_aa_count_F',
               'csa_2_5D_aa_count_G', 'csa_2_5D_aa_count_H', 'csa_2_5D_aa_count_I', 'csa_2_5D_aa_count_K',
               'csa_2_5D_aa_count_L', 'csa_2_5D_aa_count_M', 'csa_2_5D_aa_count_N', 'csa_2_5D_aa_count_P',
               'csa_2_5D_aa_count_Q', 'csa_2_5D_aa_count_R', 'csa_2_5D_aa_count_S', 'csa_2_5D_aa_count_T',
               'csa_2_5D_aa_count_V', 'csa_2_5D_aa_count_W', 'csa_2_5D_aa_count_Y', 'csa_2_5D_aa_count_chrg',
               'csa_2_5D_aa_count_total', 'disorder_2D_aa_count_A', 'disorder_2D_aa_count_C', 'disorder_2D_aa_count_D',
               'disorder_2D_aa_count_E', 'disorder_2D_aa_count_F', 'disorder_2D_aa_count_G', 'disorder_2D_aa_count_H',
               'disorder_2D_aa_count_I', 'disorder_2D_aa_count_K', 'disorder_2D_aa_count_L', 'disorder_2D_aa_count_M',
               'disorder_2D_aa_count_N', 'disorder_2D_aa_count_P', 'disorder_2D_aa_count_Q', 'disorder_2D_aa_count_R',
               'disorder_2D_aa_count_S', 'disorder_2D_aa_count_T', 'disorder_2D_aa_count_V', 'disorder_2D_aa_count_W',
               'disorder_2D_aa_count_Y', 'disorder_2D_aa_count_dis', 'disorder_2D_aa_count_ord',
               'disorder_2D_aa_count_total', 'disorder_3D_aa_count_A', 'disorder_3D_aa_count_C',
               'disorder_3D_aa_count_D', 'disorder_3D_aa_count_E', 'disorder_3D_aa_count_F', 'disorder_3D_aa_count_G',
               'disorder_3D_aa_count_H', 'disorder_3D_aa_count_I', 'disorder_3D_aa_count_K', 'disorder_3D_aa_count_L',
               'disorder_3D_aa_count_M', 'disorder_3D_aa_count_N', 'disorder_3D_aa_count_P', 'disorder_3D_aa_count_Q',
               'disorder_3D_aa_count_R', 'disorder_3D_aa_count_S', 'disorder_3D_aa_count_T', 'disorder_3D_aa_count_V',
               'disorder_3D_aa_count_W', 'disorder_3D_aa_count_Y', 'disorder_3D_aa_count_dis',
               'disorder_3D_aa_count_ord', 'disorder_3D_aa_count_total', 'dna_2_5D_aa_count_A', 'dna_2_5D_aa_count_C',
               'dna_2_5D_aa_count_D', 'dna_2_5D_aa_count_E', 'dna_2_5D_aa_count_F', 'dna_2_5D_aa_count_G',
               'dna_2_5D_aa_count_H', 'dna_2_5D_aa_count_I', 'dna_2_5D_aa_count_K', 'dna_2_5D_aa_count_L',
               'dna_2_5D_aa_count_M', 'dna_2_5D_aa_count_N', 'dna_2_5D_aa_count_P', 'dna_2_5D_aa_count_Q',
               'dna_2_5D_aa_count_R', 'dna_2_5D_aa_count_S', 'dna_2_5D_aa_count_T', 'dna_2_5D_aa_count_V',
               'dna_2_5D_aa_count_W', 'dna_2_5D_aa_count_Y', 'dna_2_5D_aa_count_dis', 'dna_2_5D_aa_count_ord',
               'dna_2_5D_aa_count_total', 'metal_2_5D_aa_count_A', 'metal_2_5D_aa_count_C', 'metal_2_5D_aa_count_D',
               'metal_2_5D_aa_count_E', 'metal_2_5D_aa_count_F', 'metal_2_5D_aa_count_G', 'metal_2_5D_aa_count_H',
               'metal_2_5D_aa_count_I', 'metal_2_5D_aa_count_K', 'metal_2_5D_aa_count_L', 'metal_2_5D_aa_count_M',
               'metal_2_5D_aa_count_N', 'metal_2_5D_aa_count_P', 'metal_2_5D_aa_count_Q', 'metal_2_5D_aa_count_R',
               'metal_2_5D_aa_count_S', 'metal_2_5D_aa_count_T', 'metal_2_5D_aa_count_V', 'metal_2_5D_aa_count_W',
               'metal_2_5D_aa_count_Y', 'metal_2_5D_aa_count_bulk', 'metal_2_5D_aa_count_carb',
               'metal_2_5D_aa_count_chrg', 'metal_2_5D_aa_count_negchrg', 'metal_2_5D_aa_count_poschrg',
               'metal_2_5D_aa_count_total', 'metal_3D_aa_count_bulk', 'metal_3D_aa_count_carb',
               'metal_3D_aa_count_chrg', 'metal_3D_aa_count_negchrg', 'metal_3D_aa_count_poschrg',
               'notdeep_3D_aa_count_A', 'notdeep_3D_aa_count_C', 'notdeep_3D_aa_count_D', 'notdeep_3D_aa_count_E',
               'notdeep_3D_aa_count_F', 'notdeep_3D_aa_count_G', 'notdeep_3D_aa_count_H', 'notdeep_3D_aa_count_I',
               'notdeep_3D_aa_count_K', 'notdeep_3D_aa_count_L', 'notdeep_3D_aa_count_M', 'notdeep_3D_aa_count_N',
               'notdeep_3D_aa_count_P', 'notdeep_3D_aa_count_Q', 'notdeep_3D_aa_count_R', 'notdeep_3D_aa_count_S',
               'notdeep_3D_aa_count_T', 'notdeep_3D_aa_count_V', 'notdeep_3D_aa_count_W', 'notdeep_3D_aa_count_Y',
               'notdeep_3D_aa_count_total', 'sites_2_5D_aa_count_A', 'sites_2_5D_aa_count_C', 'sites_2_5D_aa_count_D',
               'sites_2_5D_aa_count_E', 'sites_2_5D_aa_count_F', 'sites_2_5D_aa_count_G', 'sites_2_5D_aa_count_H',
               'sites_2_5D_aa_count_I', 'sites_2_5D_aa_count_K', 'sites_2_5D_aa_count_L', 'sites_2_5D_aa_count_M',
               'sites_2_5D_aa_count_N', 'sites_2_5D_aa_count_P', 'sites_2_5D_aa_count_Q', 'sites_2_5D_aa_count_R',
               'sites_2_5D_aa_count_S', 'sites_2_5D_aa_count_T', 'sites_2_5D_aa_count_V', 'sites_2_5D_aa_count_W',
               'sites_2_5D_aa_count_Y', 'sites_2_5D_aa_count_chrg', 'sites_2_5D_aa_count_total',
               'ss_disorder_2D_aa_count_A', 'ss_disorder_2D_aa_count_C', 'ss_disorder_2D_aa_count_D',
               'ss_disorder_2D_aa_count_E', 'ss_disorder_2D_aa_count_F', 'ss_disorder_2D_aa_count_G',
               'ss_disorder_2D_aa_count_H', 'ss_disorder_2D_aa_count_I', 'ss_disorder_2D_aa_count_K',
               'ss_disorder_2D_aa_count_L', 'ss_disorder_2D_aa_count_M', 'ss_disorder_2D_aa_count_N',
               'ss_disorder_2D_aa_count_P', 'ss_disorder_2D_aa_count_Q', 'ss_disorder_2D_aa_count_R',
               'ss_disorder_2D_aa_count_S', 'ss_disorder_2D_aa_count_T', 'ss_disorder_2D_aa_count_V',
               'ss_disorder_2D_aa_count_W', 'ss_disorder_2D_aa_count_Y', 'ss_disorder_2D_aa_count_dis',
               'ss_disorder_2D_aa_count_ord', 'ss_disorder_2D_aa_count_total', 'ss_disorder_3D_aa_count_A',
               'ss_disorder_3D_aa_count_C', 'ss_disorder_3D_aa_count_D', 'ss_disorder_3D_aa_count_E',
               'ss_disorder_3D_aa_count_F', 'ss_disorder_3D_aa_count_G', 'ss_disorder_3D_aa_count_H',
               'ss_disorder_3D_aa_count_I', 'ss_disorder_3D_aa_count_K', 'ss_disorder_3D_aa_count_L',
               'ss_disorder_3D_aa_count_M', 'ss_disorder_3D_aa_count_N', 'ss_disorder_3D_aa_count_P',
               'ss_disorder_3D_aa_count_Q', 'ss_disorder_3D_aa_count_R', 'ss_disorder_3D_aa_count_S',
               'ss_disorder_3D_aa_count_T', 'ss_disorder_3D_aa_count_V', 'ss_disorder_3D_aa_count_W',
               'ss_disorder_3D_aa_count_Y', 'ss_disorder_3D_aa_count_dis', 'ss_disorder_3D_aa_count_ord',
               'ss_disorder_3D_aa_count_total', 'surface_3D_aa_count_A', 'surface_3D_aa_count_C',
               'surface_3D_aa_count_D', 'surface_3D_aa_count_E', 'surface_3D_aa_count_F', 'surface_3D_aa_count_G',
               'surface_3D_aa_count_H', 'surface_3D_aa_count_I', 'surface_3D_aa_count_K', 'surface_3D_aa_count_L',
               'surface_3D_aa_count_M', 'surface_3D_aa_count_N', 'surface_3D_aa_count_P', 'surface_3D_aa_count_Q',
               'surface_3D_aa_count_R', 'surface_3D_aa_count_S', 'surface_3D_aa_count_T', 'surface_3D_aa_count_V',
               'surface_3D_aa_count_W', 'surface_3D_aa_count_Y', 'surface_3D_aa_count_carb', 'surface_3D_aa_count_chrg',
               'surface_3D_aa_count_negchrg', 'surface_3D_aa_count_poschrg', 'surface_3D_aa_count_total',
               'tm_2D_aa_count_tmstab', 'tm_2D_aa_count_tmunstab', 'tm_3D_aa_count_tmstab', 'tm_3D_aa_count_tmunstab']
