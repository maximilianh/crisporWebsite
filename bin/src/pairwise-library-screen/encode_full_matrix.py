#!/usr/bin/env python 
# Converts a data frame processed by ./filter_dataset.py into a format that can be used by train_rstan_model.R
# Usage: ./encode_full_matrix.py input_df.csv out_matrix.txt

from __future__ import print_function
import pandas as pd
import numpy as np
from sys import argv

in_file = argv[1]
out_file = argv[2]

df = pd.read_csv(in_file)

out_handle = open(out_file, 'w')

max_len = 24

mismatches = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']

output_size = 24 * len(mismatches)

for row in df.iterrows():
    out_vector = np.zeros(output_size)

    row_values = row[1]

    guide_seq = row_values['guide']

    target_seq = row_values['target']
    target_seq = target_seq[len(target_seq) - len(guide_seq):]

    to_print = []
    

    num_mm = 0
    for i, [dna, rna] in enumerate(zip(target_seq, guide_seq)):

        if dna != rna:

            mm_index = len(mismatches) * i + mismatches.index(dna + rna)
            out_vector[mm_index] = 1
            num_mm += 1
    
    to_print = [str(row_values['Avg Weighted Indel Rate']), str(row_values['Off:On Target Ratio']),
                str(row_values['GuideGroupID']), str(len(guide_seq)), guide_seq, target_seq, str(num_mm)]

    to_print += list(out_vector)

    print('\t'.join([str(_) for _ in to_print]), file=out_handle)


out_handle.close()
