 # wrapper script to score predicted sites
 # SP: PAM NGG
 # ST1: PAM NNAGAAW
 # Takes an input file of FASTA sequences and generates scored gRNA sequences
 # Need biopython, numpy and scipy installed

from __future__ import division

import sys
import optparse
import os
import operator
import re
import numpy as np
import scipy
import platform
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

import scipy.stats as ss
import time

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 4:
    raise Exception, "Usage: python scoreMySites.py [inputFASTAFile] [organism(Hg/Mm)] [SP/ST1] [outputHeader]"

# handles
infile = args[0]
organism = args[1]
species = args[2]
header = args[3]

# know which files to use
if species=='SP' and organism=='Hg':
	model = '293T.HiSeq.SP.Nuclease.100.SVM.Model.txt'
	dist = 'Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt'
	pam = 'NGG'
elif species=='SP' and organism=='Mm':
	model = '293T.HiSeq.SP.Nuclease.100.SVM.Model.txt'
	dist = 'Mm10.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt'
	pam = 'NGG'
elif species=='ST1' and organism=='Hg':
	model = '293T.HiSeq.ST1.Nuclease.100.V2.SVM.Model.txt'
	dist = 'Hg19.RefFlat.Genes.75bp.NoUTRs.ST1Sites.SVMOutput.txt'
	pam = 'NNAGAAW'

elif species=='ST1' and organism=='Mm':
	model = '293T.HiSeq.ST1.Nuclease.100.V2.SVM.Model.txt'
	dist = 'Mm10.RefFlat.Genes.75bp.NoUTRs.ST1Sites.SVMOutput.txt'
	pam = 'NNAGAAW'	
else:
	raise Exception, "Invalid selection! Choose Hg/Mm and SP/ST1"


# file names
gRNAFile = header + '.putativeGRNASequences.fasta'
svmInputFile = header + '.SVMInput.txt'
svmOutputFile = header + '.SVMOutput.txt'
finalOutput = header + '.FinalOutput.txt'

# first generate the sites
print 'Time: ' + str(time.ctime())
print 'Generating putative gRNA sites...'
runID = 'python identifyPutativegRNASites.py ' + infile + ' ' + pam + ' ' + gRNAFile
p = subprocess.Popen(runID,shell=True)
p.communicate()

# next generate the SVM input file
print 'Time: ' + str(time.ctime())
print 'Generating SVM input file from gRNA sequences...'
runSVMGen = 'python generateSVMFile.FASTA.py ' + gRNAFile + ' ' + svmInputFile
p = subprocess.Popen(runSVMGen,shell=True)
p.communicate()

# run the SVM
print 'Time: ' + str(time.ctime())
print 'Running classification using SVM-Light'

# SVM Classify
currPlatform = platform.system()
if currPlatform == 'Windows':
	runSVMClassify = 'svm_classify -v 0 ' + svmInputFile + ' ' + model + ' ' + svmOutputFile
else:
	runSVMClassify = './svm_classify -v 0 ' + svmInputFile + ' ' + model + ' ' + svmOutputFile
p = subprocess.Popen(runSVMClassify,shell=True)
p.communicate()

# write the final outputHeader
print 'Time: ' + str(time.ctime())
print 'Converting scores to ranks based on global ' + species + ' score distribution'
runMakeTable = 'python makeFinalTable.py ' + gRNAFile + ' ' + svmOutputFile + ' ' + dist + ' ' + finalOutput
print runMakeTable
p = subprocess.Popen(runMakeTable,shell=True)
p.communicate()

