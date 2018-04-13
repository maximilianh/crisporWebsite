#Calculates the Rule set 2 score for the given 30-mer
#Input: 1. 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
#       2. Amino acid cut position, for full model prediction only
#       3. Percent peptide, for full model prediction only
#Output: Rule set 2 score

import pandas as pd
import csv, argparse, sys
import pickle
import model_comparison

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq',
        type=str,
        help='30-mer')
    parser.add_argument('--aa-cut',
        type=int,
        default=None,
        help='Amino acid cut position of sgRNA')
    parser.add_argument('--per-peptide',
        type=float,
        default=None,
        help='Percentage of protein cut by sgRNA')
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    seq = args.seq.upper()
    if len(seq)!=30: 
        print "Please enter a 30mer sequence."
        sys.exit(1)
    aa_cut = args.aa_cut
    per_peptide = args.per_peptide
    model_file_1 = '../saved_models/V3_model_nopos.pickle'
    model_file_2 = '../saved_models/V3_model_full.pickle'
    if (aa_cut == None) or (per_peptide == None):
        model_file = model_file_1
    else:
        model_file = model_file_2
    try:
        with open(model_file, 'rb') as f:
            model= pickle.load(f)    
    except:
        raise Exception("could not find model stored to file %s" % model_file)
    if seq[25:27] == 'GG':
        score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)
        print 'Rule set 2 score: %.4f'% (score)
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'