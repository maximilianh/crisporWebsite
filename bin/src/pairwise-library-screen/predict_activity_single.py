#!/usr/bin/env python3

"""
Return the SaCas9 Specificity score for a guide sequence -- target site pair

./predict_activity_single.py guide_seq target_seq

Output (tab-delimited):

guide_seq target_seq pred_activity_score num_mismatches

"""


from __future__ import print_function
from sys import argv
import pandas as pd
from collections import defaultdict

model_dir = 'models/'

num_mismatches = range(20, 25)


# Check input

guide_seq, target_seq = argv[1], argv[2]

guide_seq = guide_seq.replace('U', 'T')

if len(target_seq) > len(guide_seq):
    target_seq = target_seq[0:len(guide_seq)]  # This deals with target sequences that include the PAM

pair_key = (guide_seq, target_seq)

guide_len = len(guide_seq)
if guide_len not in num_mismatches:
    raise Exception('Model not found for this guide length.')

if '-' in guide_seq or '-' in target_seq:
    raise Exception('This model does not support gaps/bulges.')


# Define mismatches

base_pairings = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
mismatch_ids = ['dT:rC', 'dT:rG', 'dT:rU', 'dG:rA', 'dG:rG', 'dG:rU', 'dC:rA', 'dC:rC', 'dC:rU', 'dA:rA', 'dA:rC', 'dA:rG']

mm_dict = {_: __ for _, __ in zip(base_pairings, mismatch_ids)}

# Read models

par_dict = defaultdict(dict)
for guide_len in num_mismatches:

    model_file = model_dir + 'rstan_mult_%s_coeffs.txt' % guide_len
    print (model_file)
    model_df = pd.read_table(model_file)

    par_dict[guide_len] = {}

    for i in range(1, guide_len + 1):
        par_dict[guide_len][i] = defaultdict(dict)

    for row in model_df.iterrows():
        row_values = row[1]
        par_dict[guide_len][guide_len - row_values['Position'] + 1][row_values['Mismatch']] = row_values['Median']


# Add up penalties

print("len", guide_len)
print(par_dict[20])
pred_activity = 1
num_mm = 0


for i, [dna, rna] in enumerate(zip(target_seq, guide_seq)):
    mm_pos = i + 1

    if mm_pos == 1:  # Mismatches at at the 5' position are not truly determined by the model
        continue
    if dna != rna:
        if dna + rna in mm_dict:
            print(dna+rna)
            mm_coeff = par_dict[guide_len][mm_pos][mm_dict[dna + rna]]
            num_mm += 1
            print(mm_coeff)
            pred_activity *= mm_coeff

print(guide_seq, target_seq, pred_activity, num_mm, sep='\t')

