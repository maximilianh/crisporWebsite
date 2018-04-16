import os

import numpy as np
import pandas as pd
import scipy.stats

from azimuth.model_comparison import predict


def predict_moreno(seq):
    return predict(seq, aa_cut=None, percent_peptide=None, model_file="moreno_model.pkl")


def test_moreno_predict(dataset_file):
    """
    use azimuth-moreno model on mouse/zebra-fish data
    """

    # test on alena data
    data = pd.read_csv(dataset_file, delimiter="\t", index_col=[1])
    data["30mer"] = [np.nan] * data.shape[0]
    for i in range(data.shape[0]):
        row = data.iloc[i]
        seq = row['seq']

        if len(seq) == 20:
            twentymer = seq
        elif len(seq) == 23:
            # the rest is PAM
            twentymer = seq[:20]
        else:
            raise Exception("%s: %s" % (i, twentymer))

        hundredmer = row['longSeq100Bp']
        pos = hundredmer.find(twentymer)

        if pos == -1:
            raise Exception("%s not found in 100mer" % twentymer)

        assert twentymer == hundredmer[pos:pos + 20], "%s: %s != %s" % (i, twentymer, hundredmer[pos:pos + 20])
        thirtymer = hundredmer[pos - 4:pos + 26]
        data.loc[data.iloc[i].name, "30mer"] = thirtymer

    sequences = data['30mer'].values

    y_pred = predict_moreno(sequences)
    y_true = data["modFreq"].values

    corr, _ = scipy.stats.spearmanr(y_pred, y_true)
    print "moreno spearmanr on %s = %s" % (dataset_file, corr)
    return corr


if __name__ == "__main__":
    assert np.allclose(0.312935544002, test_moreno_predict("alenaAll.scores.tab.txt"))
    assert np.allclose(0.583916446247, test_moreno_predict("teboulVivo_mm9.scores.tab.txt"))
