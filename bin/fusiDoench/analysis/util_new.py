"""
New modules with utils from the earlier util module which are needed for crispor.
"""

import numpy as np
import scipy.stats as st

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
