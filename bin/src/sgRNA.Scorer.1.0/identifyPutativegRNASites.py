# script to identify gRNA sites in a list of FASTA files
# different recognition sites:
# ST1: NNNNNNNNNNNNNNNNNNNN NNAGAAW -> W is A or T
# NM1: NNNNNNNNNNNNNNNNNNNN NNNNGATT
# SP: NNNNNNNNNNNNNNNNNNNN NGG
# Edited Dec 10, 2013 to correct the lower case to upper case. Also to add teh 


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
if len(args) != 3:
	raise Exception, "Usage: python identifyPutativegRNASites.py [fastaFile] [RecognitionSite] [outputFile]"

# open the list of files
summaryfile = open(args[2],'w')

# recognition site
recognitionSite = args[1]

# go through fasta file
fastaReads = open(args[0])
count = 1
totalgRNAsFound = 0
gRNASeen = defaultdict(str)

# process recognition site
spacer = '\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S\S'
sites = []
sitesRC = []
for char in recognitionSite:
	if len(sites)==0:
		site = '\S'
		sites.append(site)
		sitesRC.append(site)
	else:
		myList = []
		myListRC = []
		if char=='N':
			for site in sites:
				site1 = site + '\S'
				myList.append(site1)
			for site in sitesRC:
				site1 = site + '\S'
				myListRC.append(site1)

		elif char=='W':
			for site in sites:
				site1 = site + 'A'
				site2 = site + 'T'
				myList.append(site1)
				myList.append(site2)
			for site in sitesRC:
				site1 = 'A' + site
				site2 = 'T' + site
				myListRC.append(site1)
				myListRC.append(site2)
		else:
			charS = Seq(char)
			charRC = str(charS.reverse_complement())
			for site in sites:
				site1 = site + char
				myList.append(site1)
			for site in sitesRC:
				site1 = charRC + site
				myListRC.append(site1)
		sites = myList
		sitesRC = myListRC

# target site length
spacerSize = len(spacer) + len(sites[0])

for record in SeqIO.parse(fastaReads,'fasta'):
	seqID = record.id
	strSeq = str(record.seq)
	end = len(strSeq)-spacerSize
	index = 0
	while index < end:
		for toFind in sites:
			toFind = spacer + toFind
			myIter = re.finditer(toFind,strSeq[index:])
			for gRNA in myIter:
				if not 'N' in str(gRNA.group(0)) and str(gRNA.group(0)) not in gRNASeen:
					summaryfile.write('>' + str(seqID) + '_Plus_' + str(gRNA.start()) + '\n' + str(gRNA.group(0)) + '\n')
					count = count + 1
					totalgRNAsFound = totalgRNAsFound + 1
					gRNASeen[str(gRNA.group(0))] = 'Y'

		for toFindRC in sitesRC:
			toFindRC = toFindRC + spacer
			myIterRC = re.finditer(toFindRC,strSeq[index:])
			for gRNA in myIterRC:
				if not 'N' in str(gRNA.group(0)) and str(gRNA.group(0)) not in gRNASeen:
					seq = Seq(gRNA.group(0))
					seqRC = seq.reverse_complement()
					summaryfile.write('>' + str(seqID) + '_Minus_' + str(gRNA.start()) + '\n' + str(seqRC) + '\n')
					count = count + 1
					totalgRNAsFound = totalgRNAsFound + 1
					gRNASeen[str(gRNA.group(0))] = 'Y'
		index += 1
print 'Total number of gRNAs: ' + str(totalgRNAsFound)    
summaryfile.close()
