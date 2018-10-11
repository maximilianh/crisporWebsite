'''
SaCas9 Specificity high throughput pairwise library screen analysis

Inputs: library design CSV, stitched Illumina reads FASTQ (post-assembly by PEAR and split into chunks for multiprocessing)
Outputs: Hamming Code representation, rBCs representation, Recombination rates, indel rate overall

Josh Tycko, Patrick Hsu
'''
import pandas as pd
from Bio import SeqIO
import collections
from collections import defaultdict
from collections import Counter
import hamstringPDH
import csv
import multiprocessing
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import timeit
import operator
from itertools import imap
import numpy as np
import os
import re
import difflib
import datetime
import math
import copy
import pickle
import gc

start_time = timeit.default_timer()

libdf = pd.read_csv('library_output.csv')
path = os.sys.argv[1]
os.chdir(path)
filename = os.sys.argv[2]

HC_obj = libdf.Hammingbarcode.to_dict()
HC_dict = {v:k for k,v in HC_obj.items()}
HC_dict = {k:0 for k,v in HC_dict.items()}
guide_obj = libdf.guide.to_dict() 
guide_dict = {v:0 for k,v in guide_obj.items()} 
target_obj = libdf.target.to_dict()
target_dict = {v:0 for k,v in target_obj.items()} 
targuidedict = libdf.set_index(['target']).guide.to_dict()

HC_guide_hash = libdf.set_index(['Hammingbarcode']).guide.to_dict() #HC as key, guide as value
HC_target_hash = libdf.set_index(['Hammingbarcode']).target.to_dict() #HC as key, target as value
HC_synth_hash = libdf.set_index(['Hammingbarcode']).libraryseq.to_dict() #HC as key, synthesis seq as value

PAMbuffer = 'CAGGGTGTAACAT'
tracr = 'GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGATTTTTT'
right_gibson = 'AAGCTTGGCGTAACTAGATC'
left_gibson = 'TGTGGAAAGGACGAAACACC'
HC_length = 15
#target_len = 25
Day = filename[4:6]
print 'Day ', Day

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def linch_dict_divider(raw_dict, num):
    list_result = []
    len_raw_dict = len(raw_dict)
    if len_raw_dict > num:
        base_num = len_raw_dict / num
        addr_num = len_raw_dict % num
        for i in range(num):
            this_dict = dict()
            keys = list()
            if addr_num > 0:
                keys = raw_dict.keys()[:base_num + 1]
                addr_num -= 1
            else:
                keys = raw_dict.keys()[:base_num]
            for key in keys:
                this_dict[key] = raw_dict[key]
                del raw_dict[key]
            list_result.append(this_dict)

    else:
        for d in raw_dict:
            this_dict = dict()
            this_dict[d] = raw_dict[d]
            list_result.append(this_dict)
    return list_result

def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return 1*np.array(list(imap(ne, str1, str2)))

def op_ver(opcodes):
    ops = [x[0][0] for x in opcodes]
    if len(ops) % 2:
        if not (ops[0] == 'i' and ops[-1] == 'i' and set(ops[1::2]) == set(['e'])):
            return False
        else:
            proc_ops = [(x[0][0], x[3], x[4]-x[3], x[1], x[1]-x[2]) for x in opcodes[2:-2:2]]
            #type of mutation, position in ref seq, length of mutation in ref, position in read, length of mutation in read
            return proc_ops
    else:
        return False
        
def fastq_process(fastq_filename):
	handle = open(fastq_filename,"rU")
	records = SeqIO.parse(handle,"fastq")
	chunk_output = map(seq_process,records)
	seq_output, good_read_HCrBCindel = zip(*chunk_output)

	assembled_counters = [sum(l) for l in zip(*seq_output)]
	good_read_HClistfinal = []
	rBC_list = []
	indel_list = []
	for k in range(len(good_read_HCrBCindel)):
		if good_read_HCrBCindel[k][0] is not None:
			good_read_HClistfinal.append(good_read_HCrBCindel[k][0])
			rBC_list.append(good_read_HCrBCindel[k][1])
			indel_list.append(good_read_HCrBCindel[k][2])
	del(seq_output, good_read_HCrBCindel, chunk_output)
	print 'finished ' + fastq_filename, (timeit.default_timer() - start_time), 'seconds run time'

	return assembled_counters, good_read_HClistfinal, rBC_list, indel_list

def seq_process(record):
	#initialize goodreadHC_counter, all counter variables
	readnum, HC_matched, HC_unmatched, HC_corrected, HC_uncorrected, PAM_match, tracr_match, L_Gibson_arm, R_Gibson_arm, rBC_count = (0,) * 10
	guides_correct, guides_matched, guides_unmatched, guides_recombined, target_correct, target_matched, target_unmatched,target_recombined, G_T_HC_correct, G_Tr_rBC_P_HC_correct, GT_no_mismatch, GT_no_mismatch_cut,GT_no_mismatch_not_cut = (0,) * 13
	indel, err, far_from_cut, full_err, miscall, strange_miscall, high_difflib= (0,) * 7
	good_read_counter = None
	rBC = None
	readnum += 1
	seq = str(record.seq)

	#check for PAM
	PAMlocation = seq.find(PAMbuffer[0:10])
	if PAMlocation > 0:		
		PAM_match += 1

	#find right arm
	loc = seq.find(right_gibson[0:10]) 
	if loc > 0: 
		R_Gibson_arm += 1
	HC = seq[loc-HC_length:loc]
	HCraw = HC

	if HC_dict.get(HC) is not None:
		HC_matched += 1
	else:
		HC_unmatched +=1

		#try to error-correct the Hamming code
		decode = hamstringPDH.decodeHamming(HC)
		if decode['chksum'] != 'ok' and HC_dict.get(decode['nucleotide']) is not None:
			HC = decode['nucleotide']
			HC_corrected += 1
		else:
			HC_uncorrected += 1

	if HC_matched == 1 or HC_corrected == 1:
		#look for rBC
		#TTNNNNNNNNCTCC = tracr-rBC-spacer
		rBCanchor = re.search('TT[AGTC]{7,9}CTCC',seq) #flexibility in length of rBC boosts depth
		if rBCanchor:
			start,end = rBCanchor.span()
			rBC_count += 1
			rBC = seq[start+2:end-4] #6-8bp rBC, 4bp spacer

		#look for guide
		guide_designed = HC_guide_hash.get(HC)
		guide_len = len(guide_designed)

		location_GibL = seq.find(left_gibson[-10:]) + 10
		if location_GibL > 10:
			L_Gibson_arm += 1

		guide_putative = seq[location_GibL:location_GibL+guide_len]

		if guide_dict.get(guide_putative) is not None:
			guides_matched += 1
			if guide_putative == guide_designed:
				guides_correct += 1
			else:
				guides_recombined += 1 #i.e. when there is a perfect guide but with another guide's HC
		else:
			guides_unmatched += 1

		#look for target
		target_designed = HC_target_hash.get(HC)
		target_len = len(target_designed)
		if PAM_match == 1:
			target_putative = seq[PAMlocation - target_len: PAMlocation]
			if target_dict.get(target_putative) is not None:
				target_matched += 1
				if target_putative == target_designed:
					target_correct += 1 #when you have a perfect target but with another target's HC
				else:
					target_recombined += 1
			else:
				target_unmatched += 1

		#check for a perfectly good read
		if guides_correct == 1 and target_correct == 1:
			G_T_HC_correct += 1

		#check for reads passing filters (any day)
		if guides_correct == 1 and PAM_match == 1 and rBC != None:
			
			#look for tracr
			if tracr in seq or difflib.SequenceMatcher(None, tracr, seq[location_GibL+guide_len:location_GibL+guide_len+82]).real_quick_ratio() > .95:
				tracr_match += 1 #allowing for 1-2 mismatch in tracr
				G_Tr_rBC_P_HC_correct += 1 #note HC was already checked
				good_read_counter = HC

				if target_correct != 1:
					start,end = [PAMlocation - 30, PAMlocation + 15]
					target_read = seq[start:end]
					cutsite = PAMlocation-3-start
					target_ref = rBC+'CTCC'+target_designed+PAMbuffer+HCraw+right_gibson[0:9]
					s = difflib.SequenceMatcher(None, target_read, target_ref, autojunk=False)
					a = s.get_opcodes()

					if len(a) < 15: #filter out any reads with more than 6 indels + mismatches
					    indel_list = op_ver(a)
					    # indel list = [type of mutation, position in ref seq, length of mutation in ref, position in read, length of mutation in read] (see 'proc_ops' for details)

					    if not indel_list:
					        err += 1

					    elif sum(abs(x[3]-cutsite) for x in indel_list)/len(indel_list) > 22:
					    	far_from_cut += 1
					    elif set.union(set(x[2] for x in indel_list), set(x[4] for x in indel_list), set(x[0] for x in indel_list)) == set(['r', 1, -1]):
					        miscall += 1

					    elif set(x[0] for x in indel_list) == set('r'):
					        strange_miscall += 1

					    else:
					        indel += 1

					else:
					    high_difflib += 1

		#For looking at on-target sites
		if guide_designed in target_designed and G_Tr_rBC_P_HC_correct == 1:
			GT_no_mismatch += 1
			if indel == 1:
				GT_no_mismatch_cut += 1
			elif target_correct == 1:
				GT_no_mismatch_not_cut += 1

	#return all counters
	return (readnum, HC_matched, HC_unmatched, HC_corrected, HC_uncorrected, PAM_match, tracr_match, R_Gibson_arm, L_Gibson_arm, guides_correct, guides_matched, guides_unmatched, guides_recombined, target_correct, target_matched, target_unmatched, target_recombined, G_T_HC_correct, G_Tr_rBC_P_HC_correct, rBC_count, indel, err, full_err, miscall, strange_miscall, high_difflib, far_from_cut, GT_no_mismatch, GT_no_mismatch_cut,GT_no_mismatch_not_cut), (good_read_counter, rBC, indel)

def indelratecompute(HC_rBC_dict):
	HC_rBC_indelrate_dict = dict()
	for k in HC_rBC_dict:
		#HC, Read Num, Indel Read Num, Unique rBC pre-val, Unique rBC post-val, Unweighted Indel Rate, UIR stdev, Weighted Indel Rate, WIR stdev
		rBCnum_preval = len(HC_rBC_dict[k])
		invalid_rBC = []
		for r in HC_rBC_dict[k]:
			if r not in HC_validated_rBCs[k]:
				if float(HC_rBC_dict[k][r][0])/HC_rBC_dict[k][r][1] <= .6: #validate this new rBC if >40% perfect reads
					continue
				else:
					invalid_rBC.append(r)
		for r in invalid_rBC:
			del HC_rBC_dict[k][r] 
		R = sum(HC_rBC_dict[k][r][1] for r in HC_rBC_dict[k])
		IR = float(sum(HC_rBC_dict[k][r][0] for r in HC_rBC_dict[k])) 
		rBCnum = len(HC_rBC_dict[k])
		if rBCnum == 0:
			continue
		UIR = IR/R
		UIR_SD = math.sqrt(UIR*(1-UIR)) #stdev of Bernoulli distribution
		WIRlist = []
		for r in HC_rBC_dict[k]:
			WIRlist.append(float(HC_rBC_dict[k][r][0])/HC_rBC_dict[k][r][1])
		WIR = float(sum(WIRlist))/rBCnum
		WIR_SD = np.std(WIRlist)
		HC_rBC_indelrate_dict[k] = R, IR, rBCnum_preval, rBCnum, UIR, UIR_SD, WIR, WIR_SD
		
	return HC_rBC_indelrate_dict

def print_dict(d, file_name, header = None):
	with open(file_name, 'w') as f_name:
		if header is not None:
			for x in header:
				f_name.write(str(x) + ',')
			f_name.write('\n')
		for k,v in d.items():
			if isinstance(v,tuple):
				v = str(v).replace("(","")
				v = str(v).replace(")","")
			f_name.write('{},{}\n'.format(k,v))
		# w = csv.DictWriter(f_name, d.keys())
		# w.writerow(d)

def multiprocessor(fastq_list):
	#initiate the parallel processing, assemble the results, print outputs
	os.chdir('SplitFastq')

	#start n worker processes
	p = multiprocessing.Pool(processes=36) 
	total_output = p.map(fastq_process, fastq_list)
	p.close()
	gc.collect()
	print 'p.map Complete', (timeit.default_timer() - start_time), 'seconds run time'
	os.chdir('..')

	seq_output, good_read_HClist, rBC_list, indel_list = zip(*total_output)
	good_read_HClist = [item for sublist in good_read_HClist for item in sublist] #append the chunk-lists one after another
	rBC_list = [item for sublist in rBC_list for item in sublist]
	indel_list = [item for sublist in indel_list for item in sublist]
	assembled_counters = [sum(l) for l in zip(*seq_output)]
	assembled_counters = [float(i) for i in assembled_counters]
	del(seq_output)

	readnum = assembled_counters[0]
	good_HC_count = assembled_counters[1]+assembled_counters[3]

	HC_rBC_dict = defaultdict(dict)
	for k in range(len(good_read_HClist)):
		if good_read_HClist[k] in HC_rBC_dict and rBC_list[k] in HC_rBC_dict[good_read_HClist[k]]: #if HC and rBC already in dict
			HC_rBC_dict[good_read_HClist[k]][rBC_list[k]] = HC_rBC_dict[good_read_HClist[k]][rBC_list[k]][0] + indel_list[k], HC_rBC_dict[good_read_HClist[k]][rBC_list[k]][1] + 1
		else: #if HC-rBC combo isn't in dict
			HC_rBC_dict[good_read_HClist[k]][rBC_list[k]] = indel_list[k], 1 #new rBC: indel, first read
	del(good_read_HClist)
	HC_rBC_dict.pop(None,None)
	numHC = len(HC_rBC_dict)

	print_dict(HC_rBC_dict,filename+'_'+str(datetime.date.today())+'HC_rBC_alldata.csv',['HC','rBC','Indel Reads','Reads', 'etc.'])

	if Day == '00':
		print 'Checking for rBCs with intact target sites at Day 0', (timeit.default_timer() - start_time), '(s)'
		Day0num_rBClist = []
		for HC in HC_rBC_dict.keys():
			for r in HC_rBC_dict[HC].keys():
				#delete the rBCs that have >10% indel reads
				if float(HC_rBC_dict[HC][r][0])/HC_rBC_dict[HC][r][1] > 0.1: 
					HC_rBC_dict[HC].pop(r,None)

			#delete that HC that only had rBCs with indel reads		
			if len(HC_rBC_dict[HC]) == 0:
				HC_rBC_dict.pop(HC,None) 

			#write list of validated rBCs	
			HC_rBC_dict[HC] = HC_rBC_dict[HC].keys() 
			Day0num_rBClist.append(len(HC_rBC_dict[HC]))
		print_dict(HC_rBC_dict,filename+'_'+str(datetime.date.today())+'HC_validated_rBCs.csv',['HC','Validated rBCs'])
		save_obj(HC_rBC_dict,filename+'_'+str(datetime.date.today())+'HC_validated_rBCs')


	else:
		print 'indel rate compute start', (timeit.default_timer() - start_time), 'seconds run time'

		HC_rBC_dict2list = linch_dict_divider(HC_rBC_dict,36)

		# start n worker processes
		p = multiprocessing.Pool(processes=36) 
		HC_rBC_indelrate_dict = p.map(indelratecompute,HC_rBC_dict2list)
		p.close()
		HC_rBC_indelrate_dict = dict((k,v) for d in HC_rBC_indelrate_dict for (k,v) in d.items())
		print 'indel rate compute finish', (timeit.default_timer() - start_time), 'seconds run time'
		print_dict(HC_rBC_indelrate_dict,filename+'_'+str(datetime.date.today())+'HC_rBC_indelRate_data.csv',['HC','Good Reads','Indel Reads', 'rBCs pre-validation','rBCs','Unweighted Indel Rate','UIR StDev', 'Weighted Indel Rate', 'WIR StDev'])
		numHCpostval = len(HC_rBC_indelrate_dict)

		#make pandas data frame of HC dict
		indel_df = pd.DataFrame.from_dict(HC_rBC_indelrate_dict, 'index')
		indel_df.index.name = 'Hammingbarcode'
		indel_df.columns = ['Good Reads','Indel Reads', 'rBCs pre-validation', 'rBCs','Unweighted Indel Rate','UIR StDev', 'Weighted Indel Rate', 'WIR StDev']
		indel_df.to_pickle(filename+'_'+str(datetime.date.today())+'_indel_df')
		print indel_df['Weighted Indel Rate'].describe()
		print '\n'

	sum_dict = collections.OrderedDict()
	sum_dict['readnum']= readnum
	sum_dict['HC_matched']= 100*assembled_counters[1]/readnum
	sum_dict['HC_corrected']= round(100*assembled_counters[3]/assembled_counters[2],2)
	sum_dict['HC_uncorrected']= 100*assembled_counters[4]/readnum
	sum_dict['PAM_match']= 100*assembled_counters[5]/readnum
	sum_dict['R Gibson_arm match']= 100*assembled_counters[7]/readnum

	sum_dict['# of found HCs ']= numHC
	sum_dict['# of validated HCs ']= numHCpostval

	sum_dict['guide correct / good HC reads']= round(100*assembled_counters[9]/good_HC_count,2) #divide by number of matching HCs + corrected HCs
	sum_dict['guide match']= round(100*assembled_counters[10]/good_HC_count,2)
	sum_dict['guide unmatch']= round(100*assembled_counters[11]/good_HC_count,2)
	sum_dict['guide recombined']= round(100*assembled_counters[12]/good_HC_count,2)
	sum_dict['target correct']= round(100*assembled_counters[13]/good_HC_count,2) #divide by number of matching HCs + corrected HCs
	sum_dict['target match']= round(100*assembled_counters[14]/good_HC_count,2)
	sum_dict['target unmatch']= round(100*assembled_counters[15]/good_HC_count,2)
	sum_dict['target recombined']= round(100*assembled_counters[16]/good_HC_count,2)
	sum_dict['HC reads with rBC']= round(100*assembled_counters[19]/good_HC_count,2)
	sum_dict['Guide-Target-HC correct / good HC reads']= round(100*assembled_counters[17]/good_HC_count,2)
	sum_dict['Guide-tracr-rBC-PAM-HC correct / good HC reads']= round(100*assembled_counters[18]/good_HC_count,2)
	sum_dict['Guide-tracr-rBC-PAM-HC correct / total reads ']= round(100*assembled_counters[18]/readnum,2)
	sum_dict['G-Tr-rBC-PAM-HC reads with indel']= round(100*assembled_counters[20]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with err']= round(100*assembled_counters[21]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with far_from_cut']= round(100*assembled_counters[26]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with full err']= round(100*assembled_counters[22]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with miscall']= round(100*assembled_counters[23]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with strange_miscall']= round(100*assembled_counters[24]/assembled_counters[18],2)
	sum_dict['G-Tr-rBC-PAM-HC reads with high_difflib']= round(100*assembled_counters[25]/assembled_counters[18],2)
	sum_dict['GT pair with no mismatch reads']= assembled_counters[27]
	sum_dict['GT_no_mismatch with indel rate'] = round(100*assembled_counters[28]/assembled_counters[27],2)
	sum_dict['GT_no_mismatch with uncut target'] = round(100*assembled_counters[29]/assembled_counters[27],2)

	if Day == '00':
		sum_dict['# of HCs after removing rBCs with indel reads,'] = len(HC_rBC_dict) 
		sum_dict['Max # indel-free rBC per HC,']= max(Day0num_rBClist)
		sum_dict['Average # indel-free rBCs,'] = round(float(sum(Day0num_rBClist))/len(HC_rBC_dict),2)
		sum_dict['Median # indel-free rBC per HC,']= np.median(Day0num_rBClist)
		sum_dict['StDev # indel-free rBC per HC,']= round(np.std(Day0num_rBClist),2)
	else:
		sum_dict['Max # rBC per HC,']= indel_df['rBCs'].max()
		sum_dict['Average # rBC per HC,']= round(indel_df['rBCs'].mean(),2)
		sum_dict['Median # rBC per HC,']= indel_df['rBCs'].median()
		sum_dict['StDev # rBC per HC,']= round(indel_df['rBCs'].std(),2)
		sum_dict['Median Weighted Indel,']= indel_df['Weighted Indel Rate'].median()
		sum_dict['StDev Weighted Indel,']= round(indel_df['Weighted Indel Rate'].std(),2)

	sum_dict['Run Time (s)'] = (timeit.default_timer() - start_time)
	print_dict(sum_dict,filename+'_'+str(datetime.date.today())+'SummaryAnalysis.csv')
	for x in sum_dict:
		print x, sum_dict[x]

def fastqlist():
	#go to subfolder where we've split fastq into smaller chunks, output list of file names
	fastq_list = os.listdir('SplitFastq')
	return fastq_list

def mergerBC_dicts(dict_a, dict_b):
    merged_dict = dict_a
    for key, value in dict_b.iteritems():
        try:
            merged_dict[key] = list(set(merged_dict[key] + value))
        except (KeyError, AttributeError, TypeError):
            merged_dict[key] = value
    return merged_dict

#check we have validated rBC file from baseline Day 0 timepoint    
if Day != '00': 
	validated_rBC_filenameBR3 = '../Sample_BR3D00/BR3D00_HC_validated_rBCs' #delete timepoint string from the file you want to use
	HC_validated_rBCsBR3 = load_obj(validated_rBC_filenameBR3)
	validated_rBC_filenameBR4 = '../Sample_BR4D00/BR4D00_HC_validated_rBCs' #delete timepoint string from the file you want to use
	HC_validated_rBCsBR4 = load_obj(validated_rBC_filenameBR4)
	HC_validated_rBCs = mergerBC_dicts(HC_validated_rBCsBR3,HC_validated_rBCsBR4) #join both bioreps' validated rBCs
	del(HC_validated_rBCsBR3, HC_validated_rBCsBR4)

multiprocessor(fastqlist())
print(timeit.default_timer() - start_time), 'seconds run time'