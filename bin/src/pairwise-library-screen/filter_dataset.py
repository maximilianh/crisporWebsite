#!/usr/bin/env python

# Filters the Pairwise Library dataset before fitting the model 
# This script was originally tested only on "Day 3" data

# Usage: ./filter_dataset.py input_file.csv min_rbcs out_file.csv

from __future__ import print_function
from sys import argv
import pandas as pd

input_file = argv[1]
min_rbc_threshold = int(argv[2])
out_file = argv[3]

df = pd.read_csv(input_file, dtype = {'Negative Control Target': int})

# Filter out negative controls

df = df[df['Negative Control Target'] == 0]

# Filter by number of barcodes
df = df[df['rBCs BR1'] >= min_rbc_threshold]
df = df[df['rBCs BR2'] >= min_rbc_threshold]

# Exclude guides with length 25 (which were not included in any analyses)

df = df[df['guidelen'] < 25]

df.to_csv(out_file)
