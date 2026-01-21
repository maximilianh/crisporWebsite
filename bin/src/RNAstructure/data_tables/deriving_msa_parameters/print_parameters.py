'''
#		5' --> 3'
		   AX
		   UY
#		3' <-- 5'
#		    Y 
# 	---------------------------- 
	A	C	G	U    
# 	----------------------------  
A	.	.	.	-0.9
C	.	.	-2.2	.
G	.	-2.1	.	-0.6
U	-1.1	.	-1.4	.
'''

nuc_list = ['A', 'C', 'G', 'U', '-']

wc_bp	 = ['AU','UA','GC','CG','UG','GU']
nc_bp	 = ['AA','AG', 'AC', 'UU', 'UC', 'GA', 'GG', 'CA', 'CU', 'CC']
lgap_bp	 = ['-A','-U','-G','-C']
rgap_bp	 = ['A-','U-','G-','C-']
gap_bp 	 = lgap_bp + rgap_bp

dgap	 = '--'

comb = wc_bp + nc_bp + gap_bp + ['--']
temp = []


def print_stack_dg(bp, big_dic):
	print "#\t\t5' --> 3'"
	print "\t\t   " + bp[0] + "X"
	print "\t\t   " + bp[1] + "Y"
	print "#\t\t3' <-- 5'"
	print "#  \t---------------------------------"
	print "\tA       C       G       U       -"
	print "#  \t---------------------------------"

	for i in nuc_list:
		print i+"\t",
		for j in nuc_list:
			print str(big_dic[bp+i+j]) + "\t",
		print ""
	print "\n"	
		
big_dic = {}

for i in comb:
	for j in comb:
		if i in nc_bp and j in nc_bp:
			# nc:nc
			big_dic.update({i+j:.2})
		if i in nc_bp and j in gap_bp:
			# nc:-N
			big_dic.update({i+j:.6})
		if i in nc_bp and j == dgap:
			# nc:--
			big_dic.update({i+j:0.0})

		if i in wc_bp and j == dgap:
			# wc:--
			big_dic.update({i+j:0.0})
		if i == dgap and j == dgap:
			# --:--
			big_dic.update({i+j:0.0})
		if i in rgap_bp and j in rgap_bp:
			# -N:-N
			big_dic.update({i+j:1.2})
		if i in lgap_bp and j in lgap_bp:
			# N-:N-
			big_dic.update({i+j:1.2})
		if i in lgap_bp and j in rgap_bp:
			# -N:N-
			big_dic.update({i+j:0.0})	
		if i in rgap_bp and j in lgap_bp:
			# N-:-N
			big_dic.update({i+j:0.0})	
		
		if i in gap_bp and j == dgap:
			# N-:--
			big_dic.update({i+j:.6})
		if i in gap_bp and j in wc_bp:
			# N-:wc
			big_dic.update({i+j:.6})

f = open('wc_nc_stack_dg','r')
g = f.readlines()
for i in g:
	big_dic[i.split()[0]] = float(i.split()[1])
f.close()

f = open('wc_wc_stack_dg','r')
g = f.readlines()
for i in g:
	big_dic[i.split()[0]] = float(i.split()[1])
f.close()

temp_dic = {}
for i in big_dic:
	if i[::-1] not in big_dic:
			temp_dic.update({i[::-1]:big_dic[i]})
big_dic.update(temp_dic)
	
#print big_dic
for i in comb:
	print_stack_dg(i, big_dic)
