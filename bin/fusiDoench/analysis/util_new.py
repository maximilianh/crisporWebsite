"""
New modules with utils from the earlier util module which are needed for crispor.
"""

import numpy as np
import pandas
import scipy.stats as st
import scipy as sp

def spearmanr_nonan(x,y):
    '''
    same as scipy.stats.spearmanr, but if all values are unique, returns 0 instead of nan
    (Output: rho, pval)
    '''
    r, p = st.spearmanr(x, y)
    if np.isnan(p):
        if len(np.unique(x))==1 or len(np.unique(y))==1:
            print "WARNING: spearmanr is nan due to unique values, setting to 0"
            p = 0.0
            r = 0.0
        else:
            raise Exception("found nan spearman")
    assert not np.isnan(r)
    return r, p

def concatenate_feature_sets(feature_sets):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    '''
    assert feature_sets != {}, "no feature sets present"
    F = feature_sets[feature_sets.keys()[0]].shape[0]
    for set in feature_sets.keys():
        F2 = feature_sets[set].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (feature_sets.keys()[0], set)

    N = feature_sets[feature_sets.keys()[0]].shape[0]
    inputs = np.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for set in feature_sets.keys():
        inputs_set = feature_sets[set].values
        dim[set] = inputs_set.shape[1]
        dimsum += dim[set]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[set].columns.tolist())

    return inputs, dim, dimsum, feature_names

def get_data(data, y_names, organism="human", target_gene=None):
    outputs = pandas.DataFrame()
    '''
    this is called once for each gene (aggregating across cell types)
    y_names are cell types
    e.g. call: X_CD13, Y_CD13 = get_data(cd13, y_names=['NB4 CD13', 'TF1 CD13'])
    '''

    #generate ranks for each cell type before aggregating to match what is in Doench et al
    thresh = 0.8
    for y_name in y_names: # for each cell type
        y = pandas.DataFrame(data[y_name])
        # these thresholds/quantils are not used:
        y_rank, y_rank_raw, y_threshold, y_quantiles = get_ranks(y, thresh=thresh, flip=False, col_name=y_name)
        y_rank.columns = [y_name + " rank"]
        y_rank_raw.columns = [y_name + " rank raw"]
        y_threshold.columns = [y_name + " threshold"]

        outputs = pandas.concat([outputs, y, y_rank, y_threshold, y_rank_raw], axis=1)


    #aggregated rank across cell types
    average_activity = pandas.DataFrame(outputs[[y_name for y_name in y_names]].mean(1))
    average_activity.columns = ['average activity']

    average_rank_from_avg_activity = get_ranks(average_activity,  thresh=thresh, flip=False, col_name='average activity')[0]
    average_rank_from_avg_activity.columns = ['average_rank_from_avg_activity']
    average_threshold_from_avg_activity = (average_rank_from_avg_activity > thresh)*1
    average_threshold_from_avg_activity.columns = ['average_threshold_from_avg_activity']

    average_rank = pandas.DataFrame(outputs[[y_name + ' rank' for y_name in y_names]].mean(1))
    average_rank.columns = ['average rank']
    # higher ranks are better (when flip=False as it should be)
    average_threshold = (average_rank > thresh)*1
    average_threshold.columns = ['average threshold']

    # undo the log2 trafo on the reads per million, apply rank trafo right away
    average_rank_raw = pandas.DataFrame(outputs[[y_name+' rank raw' for y_name in y_names]].mean(1))
    average_rank_raw.columns = ['average rank raw']
    outputs = pandas.concat([outputs, average_rank, average_threshold, average_activity, average_rank_raw, average_rank_from_avg_activity, average_threshold_from_avg_activity], axis=1)

    # import ipdb; ipdb.set_trace()

    #sequence-specific computations
    #features = featurize_data(data)
    #strip out featurization to later
    features = pandas.DataFrame(data['30mer'])

    if organism is "human":
        target_gene = y_names[0].split(' ')[1]

    outputs['Target gene'] = target_gene
    outputs['Organism'] = organism

    features['Target gene'] = target_gene
    features['Organism'] = organism
    features['Strand'] = pandas.DataFrame(data['Strand'])
    
    return features, outputs

def impute_gene_position(gene_position):
    '''
    Some amino acid cut position and percent peptide are blank because of stop codons, but
    we still want a number for these, so just set them to 101 as a proxy
    '''

    gene_position['Percent Peptide'] = gene_position['Percent Peptide'].fillna(101.00)

    if 'Amino Acid Cut position' in gene_position.columns:
        gene_position['Amino Acid Cut position'] = gene_position['Amino Acid Cut position'].fillna(gene_position['Amino Acid Cut position'].mean())

    return gene_position

def get_ranks(y, thresh=0.8, prefix="", flip=False, col_name='score'):
    """
    y should be a DataFrame with one column
    thresh is the threshold at which to call it a knock-down or not
    col_name = 'score' is only for V2 data
    flip should be FALSE for both V1 and V2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """

    if prefix is not None:
        prefix = prefix + "_"

    #y_rank = y.apply(ranktrafo)
    y_rank = y.apply(sp.stats.mstats.rankdata)
    y_rank /= y_rank.max()

    if flip:
        y_rank = 1.0 - y_rank # before this line, 1-labels where associated with low ranks, this flips it around (hence the y_rank > thresh below)
        # we should NOT flip (V2), see README.txt in ./data

    y_rank.columns = [prefix + "rank"]
    y_threshold = (y_rank > thresh)*1

    y_threshold.columns = [prefix + "threshold"]

    # JL: undo the log2 transform (not sure this matters?)
    y_rank_raw = (2**y).apply(sp.stats.mstats.rankdata)
    y_rank_raw /= y_rank_raw.max()
    if flip:
        y_rank_raw = 1.0 - y_rank_raw
    y_rank_raw.columns = [prefix + "rank raw"]
    assert ~np.any(np.isnan(y_rank)), "found NaN ranks"

    # divides into quantiles, but not used:
    y_quantized = pandas.DataFrame(pandas.qcut(y[col_name], 5, labels=np.arange(5.0)).astype(float)) # quantized vector
    y_quantized.columns = [prefix + "quantized"]

    return y_rank, y_rank_raw, y_threshold, y_quantized
