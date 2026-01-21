#--------------------------------------------------------------------------
#Watson-Crick base pairs which constitute secondary structure
wc_bp = set(['AU','UA', 'GC', 'CG', 'GU', 'UG'])

wc_stacks = []

for i in wc_bp:
	for j in wc_bp:
		wc_stacks.extend([i+j])
#--------------------------------------------------------------------------
# XY contains possible combinations of nucleotides as they appear in tables in rna.int11.dg. 
# See the sample table on the top of the file.

XY = ['AA', 'AC', 'AG', 'AU', 
      'CA', 'CC', 'CG', 'CU', 
      'GA', 'GC', 'GG', 'GU', 
      'UA', 'UC', 'UG', 'UU']

stacks = []
g = open('rna.stack.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 1:
		for j in XY:
			stacks.append(f[i].strip()[0] + f[i+1].strip()[0] + j)
		i = i + 2
	else:
		i = i + 1	
g.close()
#--------------------------------------------------------------------------
# Extract all dg values and populate double_stack_dg in same order as in double_stacks
# Generate dictionary deouble_stack_dg

stack_dg = []
g = open('rna.stack.dg')
f = g.readlines()
i = 0
while (i < len(f)):
	if f[i][0] != '#' and len(f[i].split()) == 5:
		array = []	
		array.extend(f[i].split()[1:])
		array.extend(f[i+1].split()[1:])
		array.extend(f[i+2].split()[1:])
		array.extend(f[i+3].split()[1:])
		stack_dg.extend(array)
		i = i + 4
	i = i + 1
g.close()

for i in range(len(stacks)):
	if stacks[i] in wc_stacks:
		print stacks[i], float(stack_dg[i])
#--------------------------------------------------------------------------
