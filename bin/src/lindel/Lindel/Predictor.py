import numpy as np
import scipy.sparse as sparse
import re
import json

def gen_indel(sequence,cut_site):
    '''This is the function that used to generate all possible unique indels and 
    list the redundant classes which will be combined after'''
    nt = ['A','T','C','G']
    up = sequence[0:cut_site]
    down = sequence[cut_site:]
    dmax = min(len(up),len(down))
    uniqe_seq ={}
    for dstart in range(1,cut_site+3):
        for dlen in range(1,dmax):
            if len(sequence) > dlen+dstart > cut_site-2:
                seq = sequence[0:dstart]+sequence[dstart+dlen:]
                indel = sequence[0:dstart] + '-'*dlen + sequence[dstart+dlen:]
                array = [indel,sequence,13,'del',dstart-30,dlen,None,None,None]
                try: 
                    uniqe_seq[seq]
                    if dstart-30 <1:
                        uniqe_seq[seq] = array
                except KeyError: uniqe_seq[seq] = array
    for base in nt:
        seq = sequence[0:cut_site]+base+sequence[cut_site:]
        indel = sequence[0:cut_site]+'-'+sequence[cut_site:]
        array = [sequence,indel,13,'ins',0,1,base,None,None]
        try: uniqe_seq[seq] = array
        except KeyError: uniqe_seq[seq] = array
        for base2 in nt:
            seq = sequence[0:cut_site] + base + base2 + sequence[cut_site:]
            indel = sequence[0:cut_site]+'--'+sequence[cut_site:]
            array = [sequence,indel,13,'ins',0,2,base+base2,None,None]
            try: uniqe_seq[seq] = array
            except KeyError:uniqe_seq[seq] = array
    uniq_align = label_mh(list(uniqe_seq.values()),4)
    for read in uniq_align:
        if read[-2]=='mh':
            merged=[]
            for i in range(0,read[-1]+1):
                merged.append((read[4]-i,read[5]))
            read[-3] = merged
    return uniq_align

def label_mh(sample,mh_len):
    '''Function to label microhomology in deletion events'''
    for k in range(len(sample)):
        read = sample[k]
        if read[3] == 'del':
            idx = read[2] + read[4] +17
            idx2 = idx + read[5]
            x = mh_len if read[5] > mh_len else read[5]
            for i in range(x,0,-1):
                if read[1][idx-i:idx] == read[1][idx2-i:idx2] and i <= read[5]:
                    sample[k][-2] = 'mh'
                    sample[k][-1] = i
                    break
            if sample[k][-2]!='mh':
                sample[k][-1]=0
    return sample


def create_feature_array(ft,uniq_indels):
    '''Used to create microhomology feature array 
       require the features and label 
    '''
    ft_array = np.zeros(len(ft))
    for read in uniq_indels:
        if read[-2] == 'mh':
            mh = str(read[4]) + '+' + str(read[5]) + '+' + str(read[-1])
            try:
                ft_array[ft[mh]] = 1
            except KeyError:
                pass
        else:
            pt = str(read[4]) + '+' + str(read[5]) + '+' + str(0)
            try:
                ft_array[ft[pt]]=1
            except KeyError:
                pass
    return ft_array


def onehotencoder(seq):
    '''convert to single and di-nucleotide hotencode'''
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


def gen_prediction(seq,wb,prereq):
    '''generate the prediction for all classes, redundant classes will be combined'''
    pam = {'AGG':0,'TGG':0,'CGG':0,'GGG':0}
    guide = seq[13:33]
    if seq[33:36] not in pam:
        #return ('Error: No PAM sequence is identified.', None)
        raise Exception('Error for guide %s: the PAM %s (pos 33-36 in input seq) is not supported by Lindel.' %
            (guide, seq[33:36]))

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
    ''' Combine redundant classes based on microhomology, matrix operation'''
    combine = []
    for s in indels:
        if s[-2] == 'mh':
            tmp = []
            for k in s[-3]:
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
    """ added by Max, based on write_file, yields rows instead of writing to file,
    order is: frequency, sequence, indel-descr. """
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
        #frequency.append("{0:.8f}".format(freq[pt]*100))
        frequency.append(freq[pt]*100)

    for s,f,i in zip(sequences,frequency,indels):
        yield (f, s, i)
