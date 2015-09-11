from __future__ import division

import sys
import optparse
import os
import operator
import re

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 2:
    raise Exception, "Usage: python generateSVMFile.FASTA.py [gRNASequences.Good] [svmInputFile]\n"

# file handles
goodgRNAs = open(args[0])
svmFile = open(args[1],'w')

# binary encoding
encoding = defaultdict(str)
encoding['A'] = '0001'
encoding['C'] = '0010'
encoding['T'] = '0100'
encoding['G'] = '1000'

# store the guideRNAs
for record in SeqIO.parse(goodgRNAs,'fasta'):
	sequence = str(record.seq).upper()
	x = 0
	tw = '-1'
	# end index
	if len(sequence)==27:
		endIndex = 22
	else:
		endIndex = 21
	while x < endIndex:
		y = 0
		while y < 4:
			tw = tw + ' ' + str(x+1) + str(y+1) + ':' + encoding[sequence[x]][y]
			y += 1
		x += 1
	svmFile.write(tw + '\n')
svmFile.close()
