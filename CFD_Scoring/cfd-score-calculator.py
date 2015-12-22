#Calculates the Cutting Frequency Determination score
#Requirements: 1. Pickle file with mismatch scores in working directory
#              2. Pickle file containing PAM scores in working directory 
#Input: 1. 23mer WT sgRNA sequence
#       2. 23mer Off-target sgRNA sequence
#Output: CFD score
import pickle
import argparse
import re
import numpy as np

def get_parser():
    parser = argparse.ArgumentParser(description='Calculates CFD score')
    parser.add_argument('--wt',
        type=str,
        help='WT 23mer sgRNA sequence')
    parser.add_argument('--off',
        type=str,
        help='Off-target 23mer sgRNA sequence')
    return parser

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Unpickle mismatch scores and PAM scores
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open('pam_scores.pkl','rb'))
        return (mm_scores,pam_scores)
    except: 
        raise Exception("Could not find file with mismatch scores or PAM scores")

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

if __name__ == '__main__':
    args = get_parser().parse_args()
    mm_scores,pam_scores = get_mm_pam_scores()
    wt = args.wt
    off = args.off
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:-3]
        cfd_score = calc_cfd(wt,sg,pam)
        print "CFD score: "+str(cfd_score)



