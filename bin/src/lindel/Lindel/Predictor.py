from Bio import pairwise2
import numpy as np
import scipy.sparse as sparse
import re
import json

def gen_indel(sequence,cut_site):
    nt = ['A','T','C','G']
    up = sequence[0:cut_site]
    down = sequence[cut_site:]
    dmax = min(len(up),len(down))
    uniqe_seq ={}
    #dlen,dstart,dstop = 
    for dstart in range(1,cut_site+3):
        for dlen in range(1,dmax):
            if len(sequence) > dlen+dstart > cut_site-2:
                seq = sequence[0:dstart]+sequence[dstart+dlen:]
                try: uniqe_seq[seq] += 1
                except KeyError: uniqe_seq[seq] = 1 
    for base in nt:
        seq = sequence[0:cut_site]+base+sequence[cut_site:]
        try: uniqe_seq[seq] += 1
        except KeyError: uniqe_seq[seq] = 1 
        for base2 in nt:
            seq = sequence[0:cut_site] + base + base2 + sequence[cut_site:]
            try: uniqe_seq[seq] += 1
            except KeyError: uniqe_seq[seq] = 1 
    uniq_align = []
    for key in uniqe_seq:
        ali = pairwise2.align.globalms(key, sequence,5,-4,-13,-0.3)
        array = [None]*15
        if len(sequence)-len(key)==1:
            array[0:2] = ali[0][0], ali[0][1]
        else:
            array[0:2] = ali[-1][0], ali[-1][1]
        array[8] = cut_site-17
        uniq_align.append(array)
    uniq_align = assign_indel(np.array(uniq_align))
    uniq_align = label_mh(uniq_align,4)
    for read in uniq_align:
        if read[-2]=='mh':
            merged=[]
            for i in range(0,read[2]+1):
                merged.append((read[10]-i,read[11]))
            read[-1] = merged
    return uniq_align

def label_mh(sample,mh_len):
    for k in range(len(sample)):
        read = sample[k]
        sample[k][-2] = None
        if read[9] == 'del':
            idx = read[8] + int(read[10]) +17
            idx2 = idx + read[11]
            x = mh_len if read[11] > mh_len else read[11]
            for i in range(x,0,-1):
                if read[1][idx-i:idx] == read[1][idx2-i:idx2] and i <= read[11]:
                    sample[k][-2] = 'mh'
                    sample[k][2] = i
                    break
            if sample[k][-2]!='mh':
                sample[k][2]=0
    return sample

def create_feature_array(ft,uniq_indels):
    #require the features and label 
    ft_array = np.zeros(len(ft))
    for read in uniq_indels:
        if read[-2] == 'mh':
            mh = str(read[10]) + '+' + str(read[11]) + '+' + str(read[2])
            try:
                ft_array[ft[mh]] = 1
            except KeyError:
                pass
        else:
            pt = str(read[10]) + '+' + str(read[11]) + '+' + str(0)
            try:
                ft_array[ft[pt]]=1
            except KeyError:
                pass
    return ft_array
#convert to single and di-nucleotide hotencode
def onehotencoder(seq):
    nt= ['A','T','C','G']
    head = []
    l = len(seq)
    for k in range(l):
        for i in range(4):
            head.append(nt[i]+str(k))

    for k in range(l-1):
        for i in range(4):
            for j in range(4):
                head.append(nt[i]+nt[j]+str(k))
    head_idx = {}
    for idx,key in enumerate(head):
        head_idx[key] = idx
    encode = np.zeros(len(head_idx))
    for j in range(l):
        encode[head_idx[seq[j]+str(j)]] =1.
    for k in range(l-1):
        encode[head_idx[seq[k:k+2]+str(k)]] =1.
    return encode

def create_label_array(lb,ep_freq,seq):
    lb_array = np.zeros(len(lb))
    for pt in ep_freq[seq]['del']:
        lb_array[lb[pt]] = ep_freq[seq]['del'][pt]
    for pt in ep_freq[seq]['ins']:
        lb_array[lb[pt]] = ep_freq[seq]['ins'][pt]
    return lb_array


def assign_indel(array):
    for i in range(len(array)):
        read  = array[i,0]
        ref   = array[i,1]
        ss,cs = array[i,8], array[i,8] +17 # Define start site of the target
        s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
        if '-' in read:
            pattern = s2[1]
            loc = list(re.search(r"\b"+pattern, read).span())
            y1 = set(range(cs-15,cs+20))
            x=set(range(loc[0],loc[1]))
            if x.intersection(y1):
                array[i,9:12] = 'del',loc[0]-cs,loc[1]-loc[0]
    # Label insertion
        elif '-' in ref:
            pattern = s1[1]
            loc = list(re.search(r"\b"+pattern+r"\b", ref).span())
            y1 = set(range(cs-15,cs+15))
            x=set(range(loc[0],loc[1]))
            if x.intersection(y1):
                array[i,9:13] = 'ins',loc[0]-cs,loc[1]-loc[0],read[loc[0]:loc[1]]
    # Step 6 For 1bp deletion, alignment index changed
    for i in range(len(array)):
        if array[i,9]=='del' and array[i,11]==1:
            read  = array[i,0]
            ref   = array[i,1]
            s1, s2 = read.replace('-',''),ref.replace('-','')
            ali = pairwise2.align.globalms(s1, s2,5,-4,-13,-0.3)
            if len(ali)>1:
                array[i,0:2] = ali[0][0],ali[0][1]
                read  = array[i,0]
                ref   = array[i,1]
                ss = array[i,8]
                cs = ss + 17 # Define start site of the target
                s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
                if '-' in read:
                    pattern = s2[1]
                    loc = list(re.search(r"\b"+pattern, read).span())
                    y1 = set(range(cs-15,cs+15))
                    x=set(range(loc[0],loc[1]))
                    if x.intersection(y1):
                        array[i,9:12] = 'del',loc[0]-cs,loc[1]-loc[0]

    for i in range(len(array)):
        if array[i,10]>0:
            read  = array[i,0]
            ref   = array[i,1]
            s1, s2 = read.replace('-',''),ref.replace('-','')
            ali = pairwise2.align.globalms(s1, s2,5,-4,-13,-0.3)
            if len(ali)>1:
                shift = array[i,10]
                if len(ali)>shift:
                    if array[i,11]>1 or array[i,9]=='ins':
                        array[i,0:2] = ali[-1-shift][0],ali[-1-shift][1]
                    else:
                        array[i,0:2] = ali[0+shift][0],ali[0+shift][1]
                elif len(ali)<=shift:
                    if array[i,11]>1:
                        array[i,0:2] = ali[0][0],ali[0][1]
                    else:
                        array[i,0:2] = ali[-1][0],ali[-1][1]
                read  = array[i,0]
                ref   = array[i,1]
                ss = array[i,8]
                cs = ss +17 # Define start site of the target
                s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
                if '-' in read:
                    pattern = s2[1]
                    loc = list(re.search(r"\b"+pattern, read).span())
                    y1 = set(range(cs-5,cs+5))
                    x=set(range(loc[0],loc[1]))
                    if x.intersection(y1):
                        array[i,9:12] = 'del',loc[0]-cs,loc[1]-loc[0]
    # Label insertion
                elif '-' in ref:
                    pattern = s1[1]
                    loc = list(re.search(r"\b"+pattern+r"\b", ref).span())
                    y1 = set(range(cs-5,cs+5))
                    x=set(range(loc[0],loc[1]))
                    if x.intersection(y1):
                        array[i,9:13] = 'ins',loc[0]-cs,loc[1]-loc[0],read[loc[0]:loc[1]]
    return array

def gen_prediction(seq,wb,prereq):
    pam = {'AGG':0,'TGG':0,'CGG':0,'GGG':0}
    guide = seq[13:33]
    if seq[33:36] not in pam:
        raise Exception('Error: No NGG at position 33 (0-based). Guide: %s' % guide)
    w1,b1,w2,b2,w3,b3 = wb
    label,rev_index,features,frame_shift = prereq
    indels = gen_indel(seq,30)
    input_indel = onehotencoder(guide)
    input_ins   = onehotencoder(guide[-6:])
    input_del   = np.concatenate((create_feature_array(features,indels),input_indel),axis=None)
    cmax = gen_cmatrix(indels,label) # combine redundant classes
    dratio, insratio = softmax(np.dot(input_indel,w1)+b1)
    ds  = softmax(np.dot(input_del,w2)+b2)
    ins = softmax(np.dot(input_ins,w3)+b3)
    y_hat = np.concatenate((ds*dratio,ins*insratio),axis=None) * cmax
    return (y_hat,np.dot(y_hat,frame_shift))

def softmax(weights):
    return (np.exp(weights)/sum(np.exp(weights)))

def gen_cmatrix(indels,label):
    combine = []
    for s in indels:
        if s[-2] == 'mh':
            tmp = []
            for k in s[-1]:
                try:
                    tmp.append(label['+'.join(list(map(str,k)))])
                except KeyError:
                    pass
            if len(tmp)>1:
                combine.append(tmp)
    temp = np.diag(np.ones(557), 0)
    for key in combine:
        for i in key[1:]:
            temp[i,key[0]] = 1
            temp[i,i]=0    
    return (sparse.csr_matrix(temp))

def write_json(seq,array,freq):
    sequences,frequency,indels = [],[],[]
    ss = 13
    sequences.append(seq[0:30] + ' | '+ seq[30:60])
    frequency.append(0)
    indels.append('')
    for i in range(len(array)):
        pt = array[i][0]
        try:
            idx1,dl = map(int,pt.split('+'))
            idx1 += ss+17
            idx2 = idx1 + dl
            cs = ss+17
            if idx1 < cs:
                if idx2>=cs:
                    s = seq[0:idx1]+'-'*(cs-idx1) + ' ' + '|' + ' ' + '-'*(idx2-cs)+seq[idx2:]
                else:
                    s = seq[0:idx1]+'-'*(idx2-idx1) + seq[idx2:cs]  + ' ' + '|' + ' ' + seq[cs:]
            elif idx1 > cs:
                s = seq[0:cs]+' ' + '|' + ' '+ seq[cs:idx1]+'-'*int(dl)+seq[idx2:]
            else:
                s = seq[0:idx1]+ ' ' + '|' + ' ' +'-'*int(dl)+seq[idx2:]
            indels.append('D' + str(dl) + '  ' +str(idx1-30))            
        except ValueError:
            idx1 = int(pt.split('+')[0])
            if pt!='3':
                bp = pt.split('+')[1]
                il = str(idx1)
                indels.append('I' +il +'+' + bp)
            else:
                bp ='X' # label any insertion >= 3bp as X
                il = '>=3'
                indels.append('I3' + '+' + bp)
            s = seq[0:ss+17]+' '+bp+' '*(2-len(bp))+seq[ss+17:]
        sequences.append(s)
        frequency.append("{0:.2f}".format(freq[pt]*100))
    output = [{"Sequence": s, "Frequency": f, "Indels": i} for s,f,i in zip(sequences,frequency,indels)]
    return (json.dumps(output, indent=1))


def write_file(seq,array,freq,fname):
    sequences,frequency,indels = [],[],[]
    ss = 13
    sequences.append(seq[0:30] + ' | '+ seq[30:60])
    frequency.append('0')
    indels.append('')
    for i in range(len(array)):
        pt = array[i][0]
        try:
            idx1,dl = map(int,pt.split('+'))
            idx1 += ss+17
            idx2 = idx1 + dl
            cs = ss+17
            if idx1 < cs:
                if idx2>=cs:
                    s = seq[0:idx1]+'-'*(cs-idx1) + ' ' + '|' + ' ' + '-'*(idx2-cs)+seq[idx2:]
                else:
                    s = seq[0:idx1]+'-'*(idx2-idx1) + seq[idx2:cs]  + ' ' + '|' + ' ' + seq[cs:]
            elif idx1 > cs:
                s = seq[0:cs]+' ' + '|' + ' '+ seq[cs:idx1]+'-'*int(dl)+seq[idx2:]
            else:
                s = seq[0:idx1]+ ' ' + '|' + ' ' +'-'*int(dl)+seq[idx2:]
            indels.append('D' + str(dl) + '  ' +str(idx1-30))            
        except ValueError:
            idx1 = int(pt.split('+')[0])
            if pt!='3':
                bp = pt.split('+')[1]
                il = str(idx1)
                indels.append('I' +il +'+' + bp)
            else:
                bp ='X' # label any insertion >= 3bp as X
                il = '>=3'
                indels.append('I3' + '+' + bp)
            s = seq[0:ss+17]+' '+bp+' '*(2-len(bp))+seq[ss+17:]
        sequences.append(s)
        frequency.append("{0:.8f}".format(freq[pt]*100))
    f0 = open(fname,'w')
    for s,f,i in zip(sequences,frequency,indels):
        f0.write(s+'\t'+f + '\t'+i +'\n')
    f0.close()

def iter_results(seq,array,freq):
    " added by Max, based on write_file "
    sequences,frequency,indels = [],[],[]
    ss = 13
    sequences.append(seq[0:30] + ' | '+ seq[30:60])
    #frequency.append('0')
    frequency.append(0)
    indels.append('')
    for i in range(len(array)):
        pt = array[i][0]
        try:
            idx1,dl = map(int,pt.split('+'))
            idx1 += ss+17
            idx2 = idx1 + dl
            cs = ss+17
            if idx1 < cs:
                if idx2>=cs:
                    s = seq[0:idx1]+'-'*(cs-idx1) + ' ' + '|' + ' ' + '-'*(idx2-cs)+seq[idx2:]
                else:
                    s = seq[0:idx1]+'-'*(idx2-idx1) + seq[idx2:cs]  + ' ' + '|' + ' ' + seq[cs:]
            elif idx1 > cs:
                s = seq[0:cs]+' ' + '|' + ' '+ seq[cs:idx1]+'-'*int(dl)+seq[idx2:]
            else:
                s = seq[0:idx1]+ ' ' + '|' + ' ' +'-'*int(dl)+seq[idx2:]
            indels.append('D' + str(dl) + '  ' +str(idx1-30))
        except ValueError:
            idx1 = int(pt.split('+')[0])
            if pt!='3':
                bp = pt.split('+')[1]
                il = str(idx1)
                indels.append('I' +il +'+' + bp)
            else:
                bp ='X' # label any insertion >= 3bp as X
                il = '>=3'
                indels.append('I3' + '+' + bp)
            s = seq[0:ss+17]+' '+bp+' '*(2-len(bp))+seq[ss+17:]
        sequences.append(s)
        #frequency.append("{0:.8f}".format(freq[pt]*100))
        frequency.append(freq[pt]*100)

    for s,f,i in zip(sequences,frequency,indels):
        yield (f, s, i)
