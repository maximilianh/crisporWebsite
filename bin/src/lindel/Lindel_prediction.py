#!/usr/bin/env python 
# Author: Will Chen
''' 
1. All functions are tested under python3.5 and python 3.6
2. Add Lindel folder to your python path.
3. y_hat is the prediction of all ~450 classes of indels <30bp.
4. fs is the frameshift ratio for this sequence.
5. Input should be 65bp (30 bp upstream and 35 bp downstream of the cleavage site)
usage: pyton Lindel_predction.py your_sequence_here your_file_name_here(can be gene name or guide name you designed)
'''
import Lindel, os, sys
from Lindel.Predictor import *
import pickle as pkl

weights = pkl.load(open(os.path.join(Lindel.__path__[0], "Model_weights.pkl"),'rb'))
prerequesites = pkl.load(open(os.path.join(Lindel.__path__[0],'model_prereq.pkl'),'rb'))
seq = sys.argv[1].upper() #input your sequence here
filename = sys.argv[2]
try:
    y_hat, fs = gen_prediction(seq,weights,prerequesites)
    filename += '_fs=' + str(round(fs,3))+'.txt'
    rev_index = prerequesites[1]
    pred_freq = {}
    for i in range(len(y_hat)):
        if y_hat[i]!=0:
            pred_freq[rev_index[i]] = y_hat[i]
    pred_sorted = sorted(pred_freq.items(), key=lambda kv: kv[1],reverse=True)
    write_file(seq,pred_sorted,pred_freq,filename)
except ValueError:
    print ('Error: No PAM sequence is identified.Please check your sequence and try again')
