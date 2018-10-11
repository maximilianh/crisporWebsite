'''
Inputs: library_features.csv that describes the pairwise library
Output: library_features_annotated.csv the includes additional annotations
	annotate the genomic region type, the mismatch position and type(s), and other design features

Josh Tycko
'''

import pandas as pd
import os

PAIRINGS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',
    'Y': 'R',
    'R': 'Y'
}

def rev_comp(bases):
    return ''.join([PAIRINGS[c] for c in bases[::-1]])

libdf = pd.read_csv('library_features.csv')
libdf.set_index('LibMemberNumber')

#annotate region type: intron vs exon
def regiontype(df):
	if 'intron' in str(df['comment']):
		return 'Intron'
	elif 'exon' in str(df['comment']):
		return 'Exon'
	elif 'bless' in str(df['comment']):
		return 'Exon'
libdf['Genomic Region'] = libdf.apply(regiontype,axis = 1)

#annotate mismatch position and type
def mismatch_pos_type(df):
	if df['Negative Control Target'] == 1 or df['On Target'] == 1:
		return None,None,None,None,None,None
	guide = str(df['guide']).upper()
	target = str(df['target']).upper()

	if df['Bulge Target'] == 1:
		return None,None,None,None,None,None

	for pos in range(len(guide)):
		pos += 1
		if guide[-pos] == target[-pos]:
			continue
		else:
			pos1 = pos
			RNApos1 = guide[-pos]
			DNApos1 = target[-pos]
			if RNApos1 == 'T':
				RNApos1 = 'U'
			DNApos1 = rev_comp(DNApos1)
			if '2MM' in str(df['comment']) or '2TV' in str(df['comment']):
				for pos in range(pos1+1,len(guide)):
					if guide[-pos] == target[-pos]:
						continue
					else:
						pos2 = pos
						RNApos2 = guide[-pos]
						DNApos2 = target[-pos]
						if RNApos2 == 'T':
							RNApos2 = 'U'
						DNApos2 = rev_comp(DNApos2)	
						return pos1, RNApos1, DNApos1, pos2, RNApos2, DNApos2
				return pos1, RNApos1, DNApos1, None,None,None
			else:
				return pos1, RNApos1, DNApos1, None,None,None
	return None, None, None, None, None, None

def bulge_pos_type(df):
	guide = str(df['guide']).upper()
	target = str(df['target']).upper()
	if df['Negative Control Target'] == 1 or df['On Target'] == 1:
		return None,None,None,None
	if df['Bulge Target'] == 1:
		if len(target) == 24: #i.e. deletion
			for pos in range(len(guide)):
				pos += 1
				if guide[-pos] == target[-pos]:
					continue
				else:
					delpos = pos
					delRNA = guide[-pos]
					if delRNA == 'T':
						delRNA = 'U'
					return delpos, delRNA, None, None

		if len(target) == 26: #i.e. insertion
			for pos in range(len(guide)):
				pos += 1
				if guide[-pos] == target[-pos]:
					continue
				else:
					inspos = pos
					insDNA = rev_comp(target[-pos])
					return None, None, inspos, insDNA

	else:
		return None, None, None, None


libdf['mismatch_tuple'] = libdf.apply(mismatch_pos_type, axis = 1)

mismatch_colList = ['MM1 Position', 'MM1 RNA', 'MM1 DNA', 'MM2 Position', 'MM2 RNA', 'MM2 DNA']
for n, col in enumerate(mismatch_colList):
	libdf[col] = libdf['mismatch_tuple'].apply(lambda mismatch_tuple: mismatch_tuple[n])

libdf['bulge_tuple'] = libdf.apply(bulge_pos_type, axis = 1)

bulge_colList = ['RNA Bulge Position', 'RNA Bulge', 'DNA Bulge Position', 'DNA Bulge']
for n, col in enumerate(bulge_colList):
	libdf[col] = libdf['bulge_tuple'].apply(lambda bulge_tuple: bulge_tuple[n])

libdf.to_csv('library_features_annotated.csv')
