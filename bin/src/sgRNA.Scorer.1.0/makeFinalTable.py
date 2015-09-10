from __future__ import division

import sys
import optparse
import os
import operator
import re
import numpy as np
import scipy

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

import matplotlib.pyplot as pyplot
import scipy.stats as ss
import time

# fasta entr

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 4:
    raise Exception, "Usage: python makeFinalTable.py [gRNAFile] [SVMScoresForThisFile] [AllPossibleScores] [outputFile]"

# assign file handles
gRNAFile = open(args[0])
svmThis = open(args[1],'r')
svmAll = open(args[2],'r')
outfile = open(args[3],'w')

# first through go all scores and get the max and min
allData = []
for line in svmAll:
	line = line.rstrip('\r\n')
	allData.append(float(line))
svmAll.close()

# now get the scores of the sequences
scoreArray = []
for line in svmThis:
	line = line.rstrip('\r\n')
	scoreArray.append(float(line))
svmThis.close()

allData = np.array(allData)

# go through the fasta file
outfile.write('SeqID\tSequence\tScore\tRank\n')

index = 0
for record in SeqIO.parse(gRNAFile,'fasta'):
	# get the percentile
	#oldPerc = ss.percentileofscore(allData,scoreArray[index])
	percentile = 100.0*(allData[allData < scoreArray[index]].size / float(allData.size))
        #print "old, new", oldPerc, percentile
	# let everyone know where we are in the process
	if index % 100==0:
		print 'Finished ' + str(index) + ' sequences'
		print 'Time: ' + str(time.ctime())
	outfile.write(str(record.id) + '\t' + str(record.seq) + '\t' + str(scoreArray[index])+ '\t' + str(percentile) + '\n')
	index += 1
outfile.close()
