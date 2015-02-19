#!/usr/bin/python2.7 
from math import exp
from re import findall

seq='GGAGGAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAACCGGTGGCG' # sequence is recommend within 60~80 bases.
length_weight=20.0
left=30 # Insert the position expected to be broken. 
right=len(seq)-int(left)
print "position of breakpoint", left
print 'length of seq = '+str(len(seq))
file_temp=open("1.before removing duplication.txt", "w") 
for k in range(2,left)[::-1]:
    for j in range(left,left+right-k+1): 
        for i in range(0,left-k+1):
            if seq[i:i+k]==seq[j:j+k]:
                length=j-i
                file_temp.write(seq[i:i+k]+'\t'+str(i)+'\t'+str(i+k)+'\t'+str(j)+'\t'+str(j+k)+'\t'+str(length)+'\n')
file_temp.close()

### After searching out all microhomology patterns, duplication should be removed!! 
f1=open("1.before removing duplication.txt", "r")
s1=f1.read()
f2=open("2.all microhomology patterns.txt", "w") #After removing duplication 
f2.write(seq+'\t'+'microhomology\t'+'deletion length\t'+'score of a pattern\n')
if s1!="":
    list_f1=s1.strip().split('\n')
    sum_score_3=0
    sum_score_not_3=0

    for i in range(len(list_f1)):
        n=0
        score_3=0
        score_not_3=0
        line=list_f1[i].split('\t')
        scrap=line[0]
        left_start=int(line[1])
        left_end=int(line[2])
        right_start=int(line[3])
        right_end=int(line[4])
        length=int(line[5])

        for j in range(i):
            line_ref=list_f1[j].split('\t')
            left_start_ref=int(line_ref[1])
            left_end_ref=int(line_ref[2])
            right_start_ref=int(line_ref[3])
            right_end_ref=int(line_ref[4])

            if (left_start >= left_start_ref) and (left_end <= left_end_ref) and (right_start >= right_start_ref) and (right_end <= right_end_ref):
                if (left_start - left_start_ref)==(right_start - right_start_ref) and (left_end - left_end_ref)==(right_end - right_end_ref):
                    n+=1
            else: pass

        if n == 0:
            if (length % 3)==0:
                length_factor = round(1/exp((length)/(length_weight)),3)
                num_GC=len(findall('G',scrap))+len(findall('C',scrap))
                score_3=100*length_factor*((len(scrap)-num_GC)+(num_GC*2))
            elif (length % 3)!=0:
                length_factor = round(1/exp((length)/(length_weight)),3)
                num_GC=len(findall('G',scrap))+len(findall('C',scrap))
                score_not_3=100*length_factor*((len(scrap)-num_GC)+(num_GC*2))

            f2.write(seq[0:left_end]+'- '*length+seq[right_end:]+'\t'+scrap+'\t'+str(length)+'\t'+str(100*length_factor*((len(scrap)-num_GC)+(num_GC*2)))+'\n')
        sum_score_3+=score_3
        sum_score_not_3+=score_not_3
    print sum_score_3, sum_score_not_3
    print 'Microhomology score = ' + str(sum_score_3+sum_score_not_3)
    print 'Out-of-frame score = ' + str((sum_score_not_3)*100/(sum_score_3+sum_score_not_3))
f1.close()
f2.close()
