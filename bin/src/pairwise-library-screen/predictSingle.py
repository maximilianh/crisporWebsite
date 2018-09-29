"""
based on predict_activity_single.py
by Max Haeussler

rewritten to be faster when calculating scores for a single guide len and for
many off-targets

"""

base_pairings = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
mismatch_ids = ['dT:rC', 'dT:rG', 'dT:rU', 'dG:rA', 'dG:rG', 'dG:rU', 'dC:rA', 'dC:rC', 'dC:rU', 'dA:rA', 'dA:rC', 'dA:rG']

mm_dict = {_: __ for _, __ in zip(base_pairings, mismatch_ids)}

from os.path import dirname, join
import pandas as pd
from collections import defaultdict

class SaCas9Scorer:

    def __init__(self, guide_len):
        self.guide_len =guide_len
        model_dir = join(dirname(__file__), 'models/')
        num_mismatches = range(20, 25)
        if guide_len not in num_mismatches:
            raise Exception('Model not found for this guide length.')

        guideLenDict = dict()
        model_file = model_dir + 'rstan_mult_%s_coeffs.txt' % guide_len
        model_df = pd.read_table(model_file)

        for i in range(1, guide_len + 1):
            guideLenDict[i] = defaultdict(dict)

        for row in model_df.iterrows():
            row_values = row[1]
            guideLenDict[guide_len - row_values['Position'] + 1][row_values['Mismatch']] = row_values['Median']
        self.guideLenDict = guideLenDict

    def calcScore(self, guide_seq, target_seq):
        assert(len(guide_seq)==self.guide_len)
        assert(len(target_seq)==self.guide_len)
        if '-' in guide_seq or '-' in target_seq:
            raise Exception('This model does not support gaps/bulges.')

        # Add up penalties
        pred_activity = 1
        guideLenDict = self.guideLenDict

        for i, [dna, rna] in enumerate(zip(target_seq, guide_seq)):
            mm_pos = i + 1

            if mm_pos == 1:  # Mismatches at at the 5' position are not truly determined by the model
                continue

            if dna != rna:
                if dna + rna in mm_dict:
                    mm_coeff = guideLenDict[mm_pos][mm_dict[dna + rna]]
                    pred_activity *= mm_coeff

        return pred_activity

if __name__=="__main__":
    scorer = SaCas9Scorer(21)
    assert(scorer.calcScore("GCAACCACAAACCCACGAGGG", "GCAACCACAAACCCACGAGGG")==1.0)
    got = scorer.calcScore("GCAACCACAAACCCACGAGGG", "ACAAACACATACCCACAAGGA")
    exp = 0.041293256
    print got, exp

    got = scorer.calcScore("GGGTGAGTGAGTGTGTGCGTG", "GAGCGAGTGGGTGTGTGCGTG")
    exp = 0.134419832
    print got, exp

    scorer = SaCas9Scorer(20)
    got = scorer.calcScore("TGTGGGTGAGTGTGTGCGTG", "TGTGAGTGAGTGTGTGCGTG")
    exp = 0.314432049304524
    print got, exp
