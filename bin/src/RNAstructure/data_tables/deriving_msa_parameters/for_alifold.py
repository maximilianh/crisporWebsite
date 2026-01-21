#--------------------------------------------------------------------------

# From rna.tstack.dg
#Data Arrangement:
#For the 4x4 tables, each starting pair needs its own table.
#This is the AA table:
#                       5' --> 3'
#                           AX
#                           UY
#                       3' <-- 5'
#                           Y->
#               ----------------------------
#       (X)     A       C       G       U
#        |      ----------------------------
#        V
#       (A)     -0.8    -1.0    -0.8    -1.0
#       (C)     -0.6    -0.7    -0.6    -0.7
#       (G)     -0.8    -1.0    -0.8    -1.0
#       (U)     -0.6    -0.8    -0.6    -0.8


# From rna.int11.dg
#Data Arrangement:
#For the 4x4 tables, each starting pair needs its own table.
#This is the AA table:
#                       5' --> 3'
#                           X
#                          A A
#                          U U
#                           Y
#                       3' <-- 5'
#                           Y->
#               ----------------------------
#       (X)     A       C       G       U
#        |      ----------------------------
#        V
#       (A)     1.9     1.9     1.9     1.9
#       (C)     1.9     1.9     1.9     1.9
#       (G)     1.9     1.9     -0.7    1.9
#       (U)     1.9     1.9     1.9     1.5

#An energy of "." is "infinity," i.e. not allowed.

#--------------------------------------------------------------------------

import numpy as np

#Watson-Crick base pairs which constitute secondary structure
wc_bp = set(['AU','UA', 'GC', 'CG', 'GU', 'UG'])

#Non-canonical base pairs which constitute secondary structure
nc_bp = set(['AA','AG', 'AC', 'UU', 'UC', 'GA', 'GG', 'CA', 'CU', 'CC'])


#--------------------------------------------------------------------------
'''
Stacks 

| bp, ... stacking

     5' --> 3'
	A G
	| |	
	U C
     3' <-- 5'
	AUGC

     5' --> 3'
	C U
	| |
	G A
     3' <-- 5'
	CGUA

 AUGC is identical to CGAU

'''

stacks = []
for i in wc_bp:
	for j in nc_bp:
		if i+j not in stacks and (i+j)[::-1] not in stacks:
			stacks.append(i+j)

stacks24 = []
for i in stacks:
	if i[2]==i[3]:
		stacks24.append(i)

'''
stacks

['GUAA', 'GUAC', 'GUAG', 'GUUU', 'GUCA', 'GUGG', 'GUCC', 'GUCU', 'GUGA', 'GUUC', 'CGAA', 'CGAC', 'CGAG', 'CGUU', 'CGCA', 'CGGG', 'CGCC', 'CGCU', 'CGGA', 'CGUC', 'GCAA', 'GCAC', 'GCAG', 'GCUU', 'GCCA', 'GCGG', 'GCCC', 'GCCU', 'GCGA', 'GCUC', 'AUAA', 'AUAC', 'AUAG', 'AUUU', 'AUCA', 'AUGG', 'AUCC', 'AUCU', 'AUGA', 'AUUC', 'UGAA', 'UGAC', 'UGAG', 'UGUU', 'UGCA', 'UGGG', 'UGCC', 'UGCU', 'UGGA', 'UGUC', 'UAAA', 'UAAC', 'UAAG', 'UAUU', 'UACA', 'UAGG', 'UACC', 'UACU', 'UAGA', 'UAUC']
'''


'''
                       5' --> 3'
                           X
                          A G
                          U C
                           Y
                       3' <-- 5'

		AUXYGC = AUXY + CGYX 

split_into_stacks returns a 1D numnpy array where two stacks 
which are present in six_char strings are set 1 and rest are set 0

'''

def split_into_stacks(six_char, stks):
	p = np.zeros(len(stks), dtype = float)
	p[stks.index(six_char[:4])] += 1
	p[stks.index(six_char[2:8][::-1])] += 1
	return p

def into_stacks(stk_char, stks):
	p = np.zeros(len(stks), dtype = float)
	p[stks.index(stk_char)] += 1
	return p
	

#--------------------------------------------------------------------------
# XY contains possible combinations of nucleotides as they appear in tables in rna.int11.dg. 
# See the sample table on the top of the file.

XY = ['AA', 'AC', 'AG', 'AU', 
      'CA', 'CC', 'CG', 'CU', 
      'GA', 'GC', 'GG', 'GU', 
      'UA', 'UC', 'UG', 'UU']

double_stacks = []
g = open('rna.int11.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 2:
		for j in XY:
			double_stacks.append(f[i].rstrip().split()[0:][0] + f[i+1].rstrip().split()[0:][0] + j + f[i].rstrip().split()[0:][1] + f[i+1].rstrip().split()[0:][1])
		i = i + 2
	else:
		i = i + 1	
g.close()


tstacks = []
g = open('rna.tstack.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 1:
#		print f[i], f[i+1]
		for j in XY:
			tstacks.append(f[i].strip()[0] + f[i+1].strip()[0] + j)
		i = i + 2
	else:
		i = i + 1	
g.close()

#--------------------------------------------------------------------------
# Open rna.int11.dg
# Extract all dg values and populate double_stack_dg in same order as in double_stacks
# Generate dictionary deouble_stack_dg

double_stack_dg = []
g = open('rna.int11.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 5:
		array = []	
		array.extend(map(float, f[i].split()[1:]))
		array.extend(map(float, f[i+1].split()[1:]))
		array.extend(map(float, f[i+2].split()[1:]))
		array.extend(map(float, f[i+3].split()[1:]))
		double_stack_dg.extend(array)
		i = i + 4
	i = i + 1
g.close()


double_stack_dict = {}
for i in range(len(double_stacks)):
	double_stack_dict[double_stacks[i]] = double_stack_dg[i]
	
# delete those keys which have wc XY
keys_with_wc_XY = []
for i  in double_stack_dict:
	if i[2:4] in {'AU','UA','GC','CG','GU','UG'}:
		keys_with_wc_XY.append(i)
for i in keys_with_wc_XY:
	del double_stack_dict[i]

not_repeated_keys = []
for i in double_stack_dict:
	if (i[::-1] not in not_repeated_keys) and (i not in not_repeated_keys):
		not_repeated_keys.append(i)
unique_dict = {}
for i in not_repeated_keys:
	unique_dict[i]= double_stack_dict[i]

#--------------------------------------------------------------------------
# Open rna.tstack.dg
# Extract all dg values and populate tstack_dg in same order as in tstacks
# Generate dictionary tstack_dg

tstack_dg = []
g = open('rna.tstack.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 5:
		array = []	
		array.extend(map(float, f[i].split()[1:]))
		array.extend(map(float, f[i+1].split()[1:]))
		array.extend(map(float, f[i+2].split()[1:]))
		array.extend(map(float, f[i+3].split()[1:]))
		tstack_dg.extend(array)
		i = i + 4
	i = i + 1
g.close()

tstack_dict = {}
for i in range(len(tstacks)):
	tstack_dict[tstacks[i]] = tstack_dg[i]

unique_tstack_dict = {}
for i in tstack_dict:
	if (i not in unique_tstack_dict) and (i[::-1] not in unique_tstack_dict) and (i[2:4] not in {'AU','UA','GC','CG','GU','UG'}):
		unique_tstack_dict[i] = tstack_dict[i] 

#--------------------------------------------------------------------------

# create not_repeated_keys and create unique_dict
# AUXYGC and CGYXUA are same!

#print unique_dict

#for i in unique_dict:
#	if (i[2] == i[3]) and i != i[::-1]:
#		print i[0:4]+i[0:2][::-1], i[2:][::-1]+i[4:]
#		print double_stack_dict[i[0:4]+i[0:2][::-1]], double_stack_dict[i[2:][::-1]+i[4:]]
#		print (double_stack_dict[i[0:4]+i[0:2][::-1]]+ double_stack_dict[i[2:][::-1]+i[4:]])/2
#		print i, double_stack_dict[i], (double_stack_dict[i[0:4]+i[0:2][::-1]]+ double_stack_dict[i[2:][::-1]+i[4:]])/2.0
#		print i, unique_dict[i]
#		print i,(unique_dict[i[0:4]+i[0:2][::-1]]+ unique_dict[i[2:][::-1]+i[4:]])/2.0

#Generate numpy 2D array A 
#Columns are various stacks
#Rows are various 1x1 internal loops approximated as combination of 2 stacks
#stacks which are there in the internal loop will 1 rest of the stacks witll be 0

A = [np.array(np.zeros(len(stacks)))]
for i in unique_dict:
		A = np.append(A, [np.array(split_into_stacks(i,stacks))], axis = 0)
for i in unique_tstack_dict:
		A = np.append(A, [np.array(into_stacks(i,stacks))], axis = 0)
A = np.delete(A, 0, 0)

#Ax = b
# Array b contains the total energy of loops/double stacks
b = []
for i in unique_dict:
		b.append(unique_dict[i])
for i in unique_tstack_dict:
		b.append(unique_tstack_dict[i])


x, residuals, rank, s = np.linalg.lstsq(A,b)

k = 0
for i in stacks:
	print i, round(x[k],1)
	k = k + 1


'''
k = 0
for i in stacks:
 	if i[2]==i[3]:	
		print i, round(x[k],2)
	k = k + 1

print residuals
'''

'''
for i in unique_dict:
	if i[2] != i[3] and i[:2] == i[4:][::-1]:
		print i, round(x[stacks.index(i[0:4])]+x[stacks.index(i[2:][::-1])],2)
'''
'''
for i in unique_dict:
	if i[2:4] == 'CA': 
		print i, unique_dict[i]
	elif i[2:4] == 'AC':
		print i[::-1], unique_dict[i]
'''
'''
for i in unique_dict:
	if (i[2] == i[3]) and i != i[::-1]:
		print i, round(x[stacks.index(i[0:4])]+x[stacks.index(i[2:][::-1])],2)
'''
'''
for i in unique_dict:
	if (i[2] != i[3]):
		print i,unique_dict[i] ,round(x[stacks.index(i[0:4])]+x[stacks.index(i[2:][::-1])],2)

'''
