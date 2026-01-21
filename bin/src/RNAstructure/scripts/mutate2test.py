#imports
import os
from argparse import ArgumentParser
# import timeit
import blosum as bl
# import numpy as np
import math
# import Bio
from Bio.Seq import Seq
import subprocess

import tqdm



#global variables
codons = {} #This initilizes the codon dictionary containing all codons one-letter symbols as keys, their 3-nt triplets as values.


#the function reads in the codons from the codon table file. It stores the codon info in the global variable codons.
#file: name of codonTable file (with full path)
def getCodonTable(file):
	#open codonTable.txt file
	with open(file) as fp: 
		for line in fp:
			temp = line.split() #spliting line because each entry is tab-delimeted
			if(len(temp) != 0): #removing empty lines
				for i in range(0, len(temp), 3): #each codon has 3 entries: 3-nt triplet, the codon one letter key, the codon 3 letter key 
					if(temp[i+1] not in codons.keys()):  #if the codon has been added before, you append the 3-nt triplet to the exlready existing value list corresponding to this key
						value = [temp[i]]
						codons[temp[i+1]] = value
					else: #if codon has not been added to the dictionary yet, make a new entry
						newValue = list(codons.get(temp[i+1])) + [temp[i]]
						codons[temp[i+1]] = newValue
	fp.close()
	# for key, value in codons.items():
	# 	print(key, value)


#this function reads in a nucleotide triplet and return the codon it represents. 
#codon: 3 nucleotide triplet
#return: codon key represting triplet
def getCodon(codon):
	codon = codon.replace("U", "T").upper() #make sure all 3-nt triplets have Ts since codon Table has Ts
	finalKey = ""
	for key, value in codons.items():
		if codon in value: #check if 3-nt triplet is found in a codon corresponding 3-nt list
			finalKey = key
			break
	return finalKey #return the codon key that corresponds to the 3-nt triplet

#this function takes in a codonKey and returns the conservative codon keys that have a score higher than a threshold
#codonSymbol: codon key that I want to get conservative codon changes for.
#blosumMatrix: the type of blosum matrix to use (default is 62).  You can choose from BLOSUM 45, 50, 62, 80 and 90.
#thresholdScore: the homology score threshold used; only output codons with a score higher than this value.
#return: list of codons shown to be coservative changes to codonSymbol
def getAltAAs(codonSymbol, blosumMatrix, thresholdScore=0):
	altList = []
	alts = blosumMatrix[codonSymbol] #retreive matrix information using Biopython
	for key, value in alts.items():
		if(value > thresholdScore):
			#Removing Redundant Special Chars: B, J, Z, X, *
			if(key not in ["X", "B", "Z", "J", "*"]): #X is any amino acid, B is D or E, and Z is N or Q (they are keys representing multiple other keys) #the original AAs will already be respresented 
				altList.append(key)

	# print(altList)
	return altList #returning a list of codon keys that are considered acceptable conservative amino acid changes (by passing the score threshold provided)


#this function counts the number of nt changes between 2 sequences
#this function is no longer used in the final code but can be utilitzed in troubleshooting, testing, or analysis
#original: sequence 1; most likley WT sequence
#mutated: sequence 2; most likely the mutated sequence
#return: number of nt changes
def countMutations(original, mutated):
	x = list(original)
	y = list(mutated)

	count = 0
	for i in range(0, len(x)): #iterate over each nucleotide in original and mutated
		if(x[i] != y[i]): #Check if each original nuleotide matches its correponding mutated nucleotide
			count += 1 #if they don't match, add to the count tally

	return count #return count tally


#this function reads 2 sequences (most likley wt and a mutant), translates them and compares them and outputs location of the AA change found! 
#this is used in final output when using blosum to record the amino acid changes following this format: original codon key + location + mutated codon key
#seq1: first sequence (most likley WT)
#seq2: second sequence (most likely mutant)
#orfPosition: location of where first full codon starts in respect with the full ORF and the used sequence
#return: original codon + location + mutated codon
def annotateAltAAs(seq1, seq2, orfPosition):
	changes = ""

	#to remove any nts at the end that don't belong to a codon
	seq1 = seq1[orfPosition:]
	seq2 = seq2[orfPosition:]
	# print(orfPosition, len(seq1))

	outOfContext = len(seq1)%3
	if(outOfContext > 0):
		seq1=seq1[:-outOfContext]
		seq2=seq2[:-outOfContext]
	# print(outOfContext,len(seq1))


	seq1 = Seq(seq1) #since I use Biopything to translate, I convert sequence string to a Seq object
	seq2 = Seq(seq2)

	aa1 = seq1.translate() #this calls biopython's translate function to translate seq1 to its amino acid sequence
	aa2 = seq2.translate()

	for i in range(0, len(aa1)): #iterate over each amino acid correponsind to both seq1 and seq2
		if(aa1[i] != aa2[i]): #check if any amino acids in aa1 don't match to their corrensoding partner in aa2
			ntLocation = (i*3)+len(seq1[:orfPosition]) + 1 #annotate nucleotide location corresponding to that amino acid
			# print((i*3)+len(seq1[orfPosition:]) + 1, (i*3),len(seq1[orfPosition:]), 1)
			changes+= aa1[i]+str(ntLocation)+aa2[i]+";" #annotate names

	return changes[:-1] #return amino acid chanes


#this function calculates the probability of structure formation of model/target structure in a sequence
#mut_seq: the sequence for which the probability is calculated (can be WT or mutant sequence)
#modelStructure: the target structure
#disruptions: a description of the sequence
#disrupting: True means that any non-accepted pairing found will be removed before probability calculation (most likley used in disrupting mutation)
#disrupting: False means that the exact structure is used for probability calculation even if there are mismatch pairs introduced (most likley used in restoring mutation)
#returns: probability of structure formation as a float
def calcPstructure(mut_seq, modelStructure, disruptions, disrupting=False):
	mut_seq = mut_seq.replace("U", "T") #make sure sequence has T if using Us

	#make a temporary ct/dot file with the input mut_seq and modelStructure
	ctfile = "tempct_Pstructure.ct"
	dotfile = open("tempdot_Pstructure.dot", "w")
	dotfile.write(">"+ disruptions + "\n")
	dotfile.write(mut_seq+"\n")
	
	#if disrupting = True; all base-pairs that in the model structure that cannot form given the inputted sequence are removed to make sure probability will not be 0
	if(disrupting):
		stack = []
		pair1 = []
		pair2 = []

		mutStruct = list(modelStructure)

		#here I figure out the locations of the pairing nucleotides in pair1 and their partners in pair2
		for i in range(0, len(modelStructure)):
			if(modelStructure[i] == "("):
				stack.append(i)
			elif(modelStructure[i] == ")"):
				pair1.append(stack[-1])
				stack.pop()
				pair2.append(i)
		#here I check whether the pairing nucleotides results in an acceptable/canonical base-pairing or not
		for i in range(0, len(pair1)):
			acceptedPairing = ["AT", "TA", "GC", "CG", "GT", "TG"]
			mutPairing = mut_seq[pair1[i]] + mut_seq[pair2[i]]			
			if(mutPairing not in acceptedPairing): #if pairing is not acceptable, I make sure these 2 nucleotides are now unpaired in the structure
				mutStruct[pair1[i]] = "."
				mutStruct[pair2[i]] = "."
		dotfile.write("".join(mutStruct)+"\n") #add edited structure to the temporary dot file
	else:
		dotfile.write(modelStructure+"\n") #add original structure to the temporary dot file. Any any not acceptable base-pairs will result in a probability of 0.
	dotfile.close()

	subprocess.run(['dot2ct', 'tempdot_Pstructure.dot', ctfile], capture_output = True) #convert temporary dot file to a ct file needed to run Efn2.

	efn2 = 0
	ensembleE = 0
	#gasConstant
	RT = 310.15 * 0.0019872 # R = 1.9872 * 10^-3 Kcal/mol #RT = 0.61633008

	subprocess.run(['efn2', '--simple', ctfile, 'mut_efn2.txt'], capture_output = True)

	#extract the free energy calculated using efn2 from the output file. Information is found in first line.
	with open("mut_efn2.txt") as fp:
		efn2 = float(fp.readline().split()[4])

	#create a temporary fasta file needed to run EnsembleEnergy.
	seqfile = "temp_seqfile.fasta"
	seqfile_output = open(seqfile, "w")
	seqfile_output.write(">" + disruptions + "\n")
	seqfile_output.write(mut_seq+ "\n")
	seqfile_output.close()

	result = subprocess.run(['EnsembleEnergy', '--sequence', seqfile], capture_output = True, text = True)
	#extract ensemble energy outputted. The information is found in the line starting with Ensemble energy for
	for line in result.stdout.split("\n"):
		if(line.startswith("Ensemble energy for")):
			ensembleE = float(line.split()[4]) #kcal/mol

	#Pstructure = exp(-(DGstructure + (Ensemble folding free energy))/RT)
	Pstructure = math.exp((efn2 - ensembleE)/-RT)

	# if(Pstructure > 10):
	# 	print( Pstructure ,mut_seq, modelStructure, disruptions, disrupting)
	# return float("%.3g" %Pstructure)
	return Pstructure

#this function removes the parts of the structure that do not belong in the full open reading frame (the out of context chunks in the begining) in a dot file
#this function is not used in the final mutate2test run. It can be used in custom analysis.
#dotfile: the dotfile with the path that has the structure and where the truncated structure will be stored
#ntNumOutOfFrame: the number of nts at the begining that don't belong to the first full codon.
def frameDotFile(dotfile, ntNumOutOfFrame):
	header = ""
	seq = ""
	struct = ""
	#read dot file and extract sequence and structure
	with open(dotfile) as fp:
		header = fp.readline()
		seq = fp.readline()
		struct = fp.readline()
	fp.close()

	#create a new dot file that contains the trimmed sequence and structure
	output = open(dotfile[:-3]+"-framed.dot", "w")
	output.write(header)
	output.write(seq[ntNumOutOfFrame:])
	output.write(struct[ntNumOutOfFrame:])
	# print(ntNumOutOfFrame)
	# print(seq[ntNumOutOfFrame:])
	# print(struct[ntNumOutOfFrame:])
	output.close()
	return dotfile[:-3]+"-framed.dot"


#this function (NED version) converts indices from 0-indexed to 1-indexed and it adds the length of the begining leftout segment to the indices!!
#Given that mutatet2test referes to all mutations in 0-indexing, the final outputted mutations use 1-indexing (like a biologist). This function is used for the conversion.
#RNAsegment: the WT sequence
#mutationName: 111-CTT:108-GCC is an example or 111-CTT,115-CCC:108-GCC,99-GGG as a multi mutation example
#leftOutSegment: the length of the nts at the begining that did not belong to a full codon, they are added to correct and reflect for the full sequence inputted
#nonORFrun: only used to signfiy that the nonORF function is used here and only 1 nt at a time is used and stored mutation
#return: the new name with the corrected indices
def modifyIndices_NED(RNAsegment, mutationName, leftOutSegment, nonORFrun=False):

	#Using the mutationName, I split the disruptions and the restorations
	disruptions = mutationName.split(":")[0].split(",")
	restorations = mutationName.split(":")[1].split(",")

	newName = ""
	for i in range(0, len(disruptions)): #iterate over all disruptions, converts to 1-indexing, append the WT triplet to the final name
		disruption = disruptions[i].split("-")
		newLoc = str(int(disruption[0])+len(leftOutSegment)+1)

		wtCodon = ""
		if(nonORFrun): #if nonORF run, I use single nucleotides not 3-nucleotide triplets
			wtCodon = RNAsegment[int(disruption[0])]
		else:
			wtCodon = RNAsegment[int(disruption[0]):int(disruption[0])+3]

		newName += wtCodon +newLoc+disruption[1]+","

	newName= newName[:-1]+":"
	for i in range(0, len(restorations)):  #iterate over all restorations, converts to 1-indexing, append the WT triplet to the final name
		restoration = restorations[i].split("-")
		newLoc = str(int(restoration[0])+len(leftOutSegment)+1)

		wtCodon = ""
		if(nonORFrun): #if nonORF run, I use single nucleotides not 3-nucleotide triplets
			wtCodon = RNAsegment[int(restoration[0])]
		else:
			wtCodon = RNAsegment[int(restoration[0]):int(restoration[0])+3]

		newName += wtCodon+ newLoc+restoration[1]+","

	#111-CTT:108-GCC
	# print(mutationName,newName[:-1])
	return newName[:-1]

#this function (Pstructure version) converts indices from 0-indexed to 1-indexed and it adds the length of the begining leftout segment to the indices!! 
#the mutation naming is a bit different when using NED and Pstructure - that's why a different function is needed
#Given that mutatet2test referes to all mutations in 0-indexing, the final outputted mutations use 1-indexing (like a biologist). This function is used for the conversion.
#RNAsegment: the WT sequence
#mutation: is the full row list representing the mutation where the disrupting and restoring mutations fall in first and third location in list 
#leftOutSegment: the length of the nts at the begining that did not belong to a full codon, they are added to correct and reflect for the full sequence inputted
#return: the new mutation row valye with the corrected indices for disrupting and restoring mutations
def modifyIndices_Pstructure(RNAsegment, mutation, leftOutSegment, nonORFrun=False):
	#Using the mutation entry, I split the disruptions and the restorations
	disruptions = mutation[0].split(",")
	restorations = mutation[2].split(",")

	newMutation = list(mutation)

	newDisruption = ""
	for i in range(0, len(disruptions)):#iterate over all disruptions, converts to 1-indexing, append the WT triplet to the final name
		disruption = disruptions[i].split("-")
		newLoc = str(int(disruption[0])+len(leftOutSegment)+1)

		wtCodon = ""
		if(nonORFrun):#if nonORF run, I use single nucleotides not 3-nucleotide triplets
			wtCodon = RNAsegment[int(disruption[0])]
		else:
			wtCodon = RNAsegment[int(disruption[0]):int(disruption[0])+3]

		newDisruption += wtCodon+ newLoc+disruption[1]+","

	newMutation[0] = newDisruption[:-1]

	newRestoration= ""
	for i in range(0, len(restorations)): #iterate over all restorations, converts to 1-indexing, append the WT triplet to the final name
		restoration = restorations[i].split("-")
		newLoc = str(int(restoration[0])+len(leftOutSegment)+1)
		wtCodon = ""
		if(nonORFrun):#if nonORF run, I use single nucleotides not 3-nucleotide triplets
			wtCodon = RNAsegment[int(restoration[0])]
		else:
			wtCodon = RNAsegment[int(restoration[0]):int(restoration[0])+3]

		newRestoration += wtCodon+ newLoc+restoration[1]+","

	newMutation[2] = newRestoration[:-1]
	#111-CTT:108-GCC

	# print("mutation ",mutation, "newMutations", newMutation)
	# print(RNAsegment, mutation, newMutation)
	return newMutation #return new mutation entry with the final mutation names/indices

#this function outputs 2 drawings for each sequence:
#1- uses the sequence, calculates its parition function and uses base-pairing probabilites to color-annotate model structure (to see effect of mutation on bp probabilities in target)
#2- takes in the sequence, predicts its structure using MaxEpect and color-annotates structure using partition function base pair probabilities
#seq: is the sequence for which the prediction is done
#mutationName: represents the codon indices and what they are mutated to
#outputPSFileName_path: the path and the name of the final ps file outputted
#modelDotFile: the path of the model/target dot file where it will be used to draw drawing 1
def visualizeMutation(seq, mutationName, outputPSFileName_path, modelDotFile):
	#create a temporary fasta file that will be used in partition function calculations
	fastaFile = open("disuptionTemp.fasta", "w")
	fastaFile.write(">mutation is " + mutationName+"\n")
	fastaFile.write(seq+"\n")
	fastaFile.close()

	#run programs needed to predict and draw structure.
	subprocess.run(['MaxExpect', '--sequence', 'disuptionTemp.fasta', 'disuptionTemp.ct'], capture_output = True)
	subprocess.run(['partition', 'disuptionTemp.fasta', 'disuptionTemp.psf'], capture_output = True)
	subprocess.run(['draw', 'disuptionTemp.ct', '-p','disuptionTemp.psf', outputPSFileName_path[:-3]+'_MEA.ps'], capture_output = True)
	subprocess.run(['draw', modelDotFile, '-p','disuptionTemp.psf', outputPSFileName_path[:-3]+'_bpProbEffectOnModelStruct.ps'], capture_output = True)


#this function iterates over all mutation pairs found and stored in mutMatrix and reformats them if needed, visualizes top mutation, and outputs mutations to a txt file 
#mutMatrix: a 2D array where each row represents a mutation pair
#leftOutSegment: nts found in the begining that do not belong to a full codon
#RNAsegment: WT sequence without leftOutSegment
#modelDotFile: the path for the dotfile which contains the WT sequence and structure and will be used for color-annotating mutant base-pairs on mutations
#outputFile:the name and path of txt file that will have the final mutations outputted to it
#orfPosition: location of where first full codon starts in respect with the full ORF and the used sequence
#checkAltAAs: T/F if the user chose to check blosum alternative AAs
#maxMutations: this is the number of mutations the user wants to output from each mutational iteration(single, double, triple mutations.. etc). Defaults to outputting all.
#nonORFrun: only used to signfiy that the nonORF function is used here and only 1 nt at a time is used and stored mutation
def storeRankedMutations_NED(mutMatrix, leftOutSegment, RNAsegment, modelDotFile, modelED, outputFile, orfPosition=0, checkAltAAs=False, maxMutations = None, nonORFrun=False):
	output = open(outputFile, "w") #create final output file handler
	#outputWithSeqs = open(outputFile[:-4]+"withSeq.txt", "w")

	output.write("Model NED = " + str(modelED) + "\n") #append model ED to the output file.

	#append column headers to the output file.
	if(not checkAltAAs):
		output.write("Index \t Disruptions \t DisruptedNED \t\t Restorations \t RestoringNED \t DisruptedNTs# \t\t DisruptedSeq \t RestoredSeq \n")
	else: #if altAA is checked, I add another colum to annotate the alternative codon restoring mutations. I append WT amino acid + location + mutated amino acid.
		output.write("Index \t Disruptions \t DisruptedNED \t\t Restorations \t RestoringNED \t DisruptedNTs# \t\t DisruptedSeq \t RestoredSeq \t AlternativeAA \n")

	#create folders/subfolders where visualized mutations will be stored. These will created in the same working directory as the input ct file
	output_path = os.path.dirname(outputFile)+"/"
	if not os.path.exists(output_path+"topHitDrawings/"):
		os.makedirs(output_path+"topHitDrawings/")
	if not os.path.exists(output_path+"topHitDrawings/totalNEDrank/"):
		os.makedirs(output_path+"topHitDrawings/totalNEDrank/")
	if not os.path.exists(output_path+"topHitDrawings/nt-NEDrank/"):
		os.makedirs(output_path+"topHitDrawings/nt-NEDrank/")

	disruptedNTcount = 0 #count the number of nts that are disrupted (aka have an ED higher than the threshold)
	for k in range(0, len(mutMatrix)):

		modifiedName = modifyIndices_NED(RNAsegment, mutMatrix[k][0],leftOutSegment, nonORFrun=nonORFrun) #this modifies the final mutation names to be 1-indexed and include the WT codon, location of change, mutated codon
		#print(modifiedName,mutMatrix[k][0])

		####this adds the number of disrupted NT to the end of the row for each mutations
		########this is needed for storing and ranking the mutations according to the number of nts disrupted
		for p in range(3, len(mutMatrix[k])):
			if(mutMatrix[k][p] == 1):
				disruptedNTcount += 1
		mutMatrix[k].append(disruptedNTcount)

		#here we're creating a full disrupting and restoring sequences that will be used in final output file.
		newRNAsegment  = RNAsegment
		disruptions = modifiedName.split(":")[0].split(",") #iterate over disruptions
		for i in range(0, len(disruptions)):
		
			disuptionLoc = 0
			disruptionCodon = ""
			if(nonORFrun): #if nonORF run, I use single nucleotides not 3-nucleotide triplets
				disuptionLoc = int(disruptions[i][1:-1]) - (1 + len(leftOutSegment)) #this is to make sure it's back to zero-indexed and without the leftoutsegmeent
				disruptionCodon = disruptions[i][-1]

				newRNAsegment = newRNAsegment[:disuptionLoc] + disruptionCodon + newRNAsegment[disuptionLoc+1:]

			else:
				disuptionLoc = int(disruptions[i][3:-3]) - (1 + len(leftOutSegment)) #this is to make sure it's back to zero-indexed and without leftoutsegmeent
				disruptionCodon = disruptions[i][-3:]

				newRNAsegment = newRNAsegment[:disuptionLoc] + disruptionCodon + newRNAsegment[disuptionLoc+3:]



		disruptedSeq =leftOutSegment +newRNAsegment
		####this adds the number of disrupted NT to the end of the row for each mutations
		mutMatrix[k].append(disruptedSeq)		

		#each line represnts a mutation pair: disrupting and restoring.
		disruptNED_rounded = round(float(mutMatrix[k][1]), 2)
		line = str(k+1) + "\t" + modifiedName.split(":")[0] + "\t" + str(disruptNED_rounded)  + "\t\t" #I append the disrupting sequence

		restorations = modifiedName.split(":")[1].split(",") #create the full restoring sequence
		for i in range(0, len(restorations)): #iterating over all restorations
			restorationLoc = 0
			restorationCodon = ""

			if(nonORFrun): #if nonORF run, I use single nucleotides not 3-nucleotide triplets
				disuptionLoc = int(disruptions[i][1:-1]) - (1 + len(leftOutSegment)) #this is to make sure it's back to zero-indexed and without the leftoutsegmeent
				restorationLoc = int(restorations[i][1:-1]) - (1 + len(leftOutSegment))
				restorationCodon = restorations[i][-1]

				newRNAsegment = newRNAsegment[:restorationLoc] + restorationCodon + newRNAsegment[restorationLoc+1:]

			else:
				restorationLoc = int(restorations[i][3:-3]) - (1 + len(leftOutSegment))
				restorationCodon = restorations[i][-3:]

				newRNAsegment = newRNAsegment[:restorationLoc] + restorationCodon + newRNAsegment[restorationLoc+3:]

		restoredSeq = leftOutSegment +newRNAsegment
		# print(disruptedSeq, restoredSeq)

		####this adds the number of disrupted NT to the end of the row for each mutations
		mutMatrix[k].append(restoredSeq)		

		#append restorting mutation to final line representing that mutation pair
		restoreNED_rounded = round(float(mutMatrix[k][2]), 2)
		line += modifiedName.split(":")[1] + "\t" + str(restoreNED_rounded) + "\t" + str(disruptedNTcount) +"\t\t"+disruptedSeq +"\t" + restoredSeq
		
		#if the user specified the use of alternative restoring amino acids, I annotate the changes and add to the final output file.
		if(checkAltAAs):
			line+= "\t" + annotateAltAAs(leftOutSegment+RNAsegment, restoredSeq, orfPosition) #this is because the leftOutSegment is excluded from both WT & restored

		#lineWithoutSeqs = str(k) + "\t" + mutMatrix[k][0].split(":")[0] + "\t" + str(mutMatrix[k][1]) + "\t\t" + mutMatrix[k][0].split(":")[1] + "\t"+ str(mutMatrix[k][2]) +"\n"
		# outputWithSeqs.write(line)
		output.write(line+"\n")
		line =""
		disruptedNTcount = 0

		#if the user specified a maxMutation number, I only output that number of mutations in the final output file.
		if(k == 0):
			#here I only draw mutations of the top ranking mutations
			mut = outputFile.split("/")[-1].split("-MutationsRanked")[0]
			visualizeMutation(disruptedSeq, modifiedName, output_path+"topHitDrawings/totalNEDrank/"+mut+"topDisruption.ps", modelDotFile)
			visualizeMutation(restoredSeq, modifiedName,  output_path+"topHitDrawings/totalNEDrank/"+mut+"topRestoration.ps", modelDotFile)

		if(maxMutations == k):
			break

		
	output.close()

	####this is to sort according to the number of dirsupted nts, and then store
	output = open(outputFile[:-4]+"-disruptedNTsNum.txt", "w")
	#outputWithSeqs = open(outputFile[:-4]+"withSeq.txt", "w")

	output.write("Model NED = " + str(modelED) + "\n")
	output.write("Index \t Disruptions \t DisruptedNED \t\t Restorations \t RestoringNED \t DisruptedNTs# \t\t DisruptedSeq \t RestoredSeq \n")

	# I sort mutations in an descending order using the number of disupted nts, then in a descending order using the disruption NED, then in a ascending order using the restoration NED 
	mutMatrixSortedAgain = sorted(mutMatrix, key=lambda x: (-x[-3], -x[1], x[2]))
	for k in range(0, len(mutMatrixSortedAgain)):
		modifiedName = modifyIndices_NED(RNAsegment, mutMatrixSortedAgain[k][0],leftOutSegment, nonORFrun)
		
		disruptNED_rounded = round(float(mutMatrixSortedAgain[k][1]), 2)
		restoreNED_rounded = round(float(mutMatrixSortedAgain[k][2]), 2)

		line = str(k+1) + "\t" + modifiedName.split(":")[0] + "\t" + str(disruptNED_rounded)  + "\t\t" + modifiedName.split(":")[1] + "\t" + str(restoreNED_rounded) + "\t" + str(mutMatrixSortedAgain[k][-3]) +"\t\t"+ mutMatrixSortedAgain[k][-2] +"\t" + mutMatrixSortedAgain[k][-1] 
		#if the user specified the use of alternative restoring amino acids, I annotate the changes and add to the final output file.
		if(checkAltAAs):
			line+= "\t" + annotateAltAAs(leftOutSegment+RNAsegment, restoredSeq, orfPosition) #this is because the leftOutSegment is excluded from both WT & restored

		output.write(line+"\n")
		line =""
		
		#if the user specified a maxMutation number, I only output that number of mutations in the final output file.
		if(k == 0):
			# visualizeMutation(mutMatrixSortedAgain[k][-2] , modifiedName, outputFile[:-4]+"topNTnumDisruption.ps", modelDotFile)
			# visualizeMutation(mutMatrixSortedAgain[k][-1], modifiedName, outputFile[:-4]+"topNTnumRestoration.ps", modelDotFile)
			mut = outputFile.split("/")[-1].split("-MutationsRanked")[0]

			#here I only draw mutations of the top ranking mutations
			visualizeMutation(mutMatrixSortedAgain[k][-2] , modifiedName, output_path+"topHitDrawings/nt-NEDrank/"+mut+"topDisruption.ps", modelDotFile)
			visualizeMutation(mutMatrixSortedAgain[k][-1], modifiedName,  output_path+"topHitDrawings/nt-NEDrank/"+mut+"topRestoration.ps", modelDotFile)


		if(maxMutations == k):
			break

	output.close()



#this function iterates over all mutation pairs found and stored in mutMatrix and reformats them if needed, visualizes top mutation, and outputs mutations to a txt file 
#mutations: a 2D array where each row represents a mutation pair
#leftOutSegment: nts found in the begining that do not belong to a full codon
#modelDotFile: the path for the dotfile which contains the WT sequence and structure and will be used for color-annotating mutant base-pairs on mutations
#modelPstuct: the probability of WT sequence forming model/target structure
#outputFile:the name and path of txt file that will have the final mutations outputted to it
#RNAsegment: WT sequence without leftOutSegment
#orfPosition: location of where first full codon starts in respect with the full ORF and the used sequence
#checkAltAAs: T/F if the user chose to check blosum alternative AAs
#maxMutations: this is the number of mutations the user wants to output from each mutational iteration(single, double, triple mutations.. etc). Defaults to outputting all.
#return: sorted mutations based on Pstructure
def storeRankedMutation_Pstructure(mutations, leftOutSegment, modelDotFile, modelPstruct, outputFile, RNAsegment, orfPosition=0, checkAltAAs=False, maxMutations = None, nonORFrun=False):
	#print(mutations)
	output = open(outputFile, "w")  #create final output file handler
	#outputWithSeqs = open(outputFile[:-4]+"withSeq.txt", "w")

	output.write("Model Pstructure = " + str(modelPstruct) + "\n")

	#append column headers to the output file.
	if(not checkAltAAs):
		output.write("Index \t Disruptions \t DisruptedPstructure \t Restorations \t RestoringPstructure \t DisruptedSeq \t RestoredSeq \n")
	else:
		output.write("Index \t Disruptions \t DisruptedPstructure \t Restorations \t RestoringPstructure \t DisruptedSeq \t RestoredSeq \t AlternativeAA \n")


	# I sort mutations in an ascending order using the disruption probability, then in a descending order using the restoration probability 
	mutations_sorted = sorted(mutations, key=lambda x: (-x[1], x[3]))
	# print(mutations_sorted)

	#create folders/subfolders where visualized mutations will be stored. These will created in the same working directory as the input ct file
	output_path = os.path.dirname(outputFile)+"/"
	if not os.path.exists(output_path+"topHitDrawings/"):
		os.makedirs(output_path+"topHitDrawings/")

	for k in range(0, len(mutations)): #iterate over all mutations

		modifiedMutation = modifyIndices_Pstructure(RNAsegment, mutations[k],leftOutSegment, nonORFrun) #convert mutations to 1-indexing and append WT codon changes to the final name.

		line = str(k+1) #each line represnts a mutation pair: disrupting and restoring.
		for i in range(0, len(modifiedMutation)):
			if(i == 1):
				disruptP = float("%.3g" %modifiedMutation[i])
				line = line+ "\t" + str(disruptP)
				pass
			elif(i == 3):
				restoreP = float("%.3g" %modifiedMutation[i])				
				line = line+ "\t" + str(restoreP)				
				pass
			else:
				line = line+ "\t" + str(modifiedMutation[i])

		#if the user specified the use of alternative restoring amino acids, I annotate the changes and add to the final output file.
		if(checkAltAAs):
			line+= "\t" + annotateAltAAs(leftOutSegment+RNAsegment, modifiedMutation[-1], orfPosition)

		line += "\n"

		output.write(line) #write final mutation line to the output file.
		line =""

		#if the user specified a maxMutation number, I only output that number of mutations in the final output file.
		if(k == 0):
			#visualizeMutation(seq, mutationName, outputPSFileName_path, modelDotFile):
			mut = outputFile.split("/")[-1].split("-MutationsRanked")[0]
			visualizeMutation(modifiedMutation[-2], modifiedMutation[0], output_path+"topHitDrawings/"+mut+"topDisruption.ps", modelDotFile)
			visualizeMutation(modifiedMutation[-1], modifiedMutation[2],  output_path+"topHitDrawings/"+mut+"topRestoration.ps", modelDotFile)
		
		if(maxMutations == k):
			break
		
	output.close()

	# print(mutations)
	# print(mutations_sorted)
	return mutations_sorted



#this function is run for each NED mutation where I iterate over each nt with its individual ED and determine if it is affected or not (1 or 0)
#mutationMatrix: the matrix for all successfull mutational pairs that were found using NED. It's 2D list: a list of lists.
#oneEDthreshold: is the individual ED score based on which I determine if the individual nt is disrupted or not
#return: the new matrix with zeros and ones represening individual nts disruption for each mutation
def convertEDMatrixToZeroAndOnes(mutMatrix, oneEDthreshold):

	WT = mutMatrix[0] #the first line is always added as is because it contains the WT sequence information 
	newMatrix = []
	newMatrix.append(WT)
	#print(WT)

	for i in range(1, len(mutMatrix)): #iterate over all mutations (not including WT) lists
		row = []
		#add the name, disruption NED, and restoring NED in the begining of the list as is.
		row.append(mutMatrix[i][0])
		row.append(mutMatrix[i][1])
		row.append(mutMatrix[i][2])

		#iterate over every nucleotide ED in each mutation. The nucleotide NED information start from the 4th location of the list.
		for j in range(3, len(mutMatrix[i])):
			#print(mutMatrix[i][j], WT[j])
			EDdifference = mutMatrix[i][j] - WT[j] #Used to check if the ED difference is considered disruption or not.
			if(EDdifference >= oneEDthreshold):
				row.append(1) #1 means disrupted nucleotide
			else:
				row.append(0) # 0 means not disrupted

		newMatrix.append(row) #append 0/1 row entry to new matrix

	return newMatrix #return the ne wmatrix with zeros and ones.


#this function checks if there is any overlap in nt diruption for 2 mutations
#this function if no longer used in the final mutate2test. It can be used later in combination with convertEDMatrixToZeroAndOnes
#if we'd like to limit the number of combined mutations tried. such that no 2 mutations with exact nt disruption overlap are combined.
#row1: mutation #1
#row2: mutation #2
#return: True or False if it finds overlap or not
def checkOverlap(row1, row2):
	# row1Unique = False
	# row2Unique = False

	row1= row1[3:]
	row2=row2[3:]
	if(row1 == row2):
		return True
	else:
		return False

#this function calls EDcalculator and calculates NED using start and end locations if provided, and EDs for each nucleotide
#seq: the sequence to calculate ED for if adopting structure
#structure: the structure file that EDcalculator will calculate ED for.
#desiredStart: start location of local region to place mutations at (default is none, looking everywhere at structure)
#desiredEnd: end location of local region to place mutations at (default is none, looking everywhere)
#return: return NED and a list of EDs for each nuleotide!
def calculateNED(seq, structure, desiredStart, desiredEnd, smp):
	mutFile = open("tempMut.dot", "w")
	mutFile.write("> mutation \n")
	mutFile.write(seq + "\n")
	mutFile.write(structure + "\n")
	mutFile.close()

	result = ""
	EDcalculator_version = ""
	if(smp):
		EDcalculator_version = "EDcalculator-smp"
	else:
		EDcalculator_version = "EDcalculator"

	if((desiredStart != None) and (desiredEnd != None)): #this is run when desired start and end are identified 
		
		#since EDcalculator appends to file if present, I need to remove old file if present to avoid confusion
		if os.path.exists("nucfileOutputTemp.txt"):
			os.remove("nucfileOutputTemp.txt")
		result = subprocess.run([EDcalculator_version, "tempMut.dot", '-s', str(desiredStart), '-e', str(desiredEnd), '--nucfile','nucfileOutputTemp.txt'], capture_output = True, text = True)
		#print(result.stdout)
	elif(desiredStart != None):
		#since EDcalculator appends to file if present, I need to remove old file if present to avoid confusion
		if os.path.exists("nucfileOutputTemp.txt"):
			os.remove("nucfileOutputTemp.txt")
		result = subprocess.run([EDcalculator_version, "tempMut.dot", '-s', str(desiredStart),'--nucfile','nucfileOutputTemp.txt'], capture_output = True, text = True)
	elif(desiredEnd != None):
		if os.path.exists("nucfileOutputTemp.txt"): #since EDcalculator appends to file if present, I need to remove old file if present to avoid confusion
			os.remove("nucfileOutputTemp.txt")
		result = subprocess.run([EDcalculator_version, "tempMut.dot", '-e', str(desiredEnd), '--nucfile','nucfileOutputTemp.txt'], capture_output = True, text = True)
	else:
		if os.path.exists("nucfileOutputTemp.txt"): #since EDcalculator appends to file if present, I need to remove old file if present to avoid confusion
			os.remove("nucfileOutputTemp.txt")
		result = subprocess.run([EDcalculator_version, "tempMut.dot", '--nucfile','nucfileOutputTemp.txt'], capture_output = True, text = True)


	#here I extract the NED calculated for the file. That info. is found in the line that starts with Structure 
	ED = 0
	for line in result.stdout.split("\n"):
		if(line.startswith("Structure")):
			# print(l)
			temp = line.split()
			# ED = round(float(temp[7]), 2)
			ED = float(temp[7])
	
	#here I extract the ED for each nucleotide and append to a new list where first ED corresponds to the ED for the first nucleotide.. etc.
	newRow = []
	with open("nucfileOutputTemp.txt") as fpNT:
		#remove first 2 lines that do not contain information for the nuceotide EDs
		fpNT.readline()
		fpNT.readline()

		#iterate over each line where each line represents a nucleotide, extract the ED, append to list
		for lineNT in fpNT:
			if(not lineNT.isspace()):
				tempNT = lineNT.split()
				newRow.append(float(tempNT[1]))
	fpNT.close()

	return ED, newRow

#this function combines 2 mutations (using NED), and then evaluates if it's with acceptable NED ranges or not. 
#row1: first mutation info
#row2: second mutation info
#leftOutSegment: nts found in the begining that do not belong to a full codon
#RNAsegment: WT sequence without leftOutSegment
#modelED: the ED calculated for the WT seq and WT target structure which will be used for evaluating success of combined mutation
#desiredStart: start location of local region to place mutations at (default is none, looking everywhere at structure)
#desiredEnd: end location of local region to place mutations at (default is none, looking everywhere)
#acceptedEDrestoringDiff: the max ED increase beyond WT ED that is acceptable
#return: T/F flag of whether combined mutation is accepted && combined mutation info.
def combineAndCalc_NED(row1, row2, leftOutSegment, RNAsegment, WTstructure, modelED, desiredStart = None , desiredEnd = None, acceptedEDrestoringDiff = 0, nonORFrun=False, smp=False, verbose=False):
	#I extract the mutation names that contain mutated codon location and mutated 3-nt triplet for both mutations to be combined
	row1Name = row1[0]
	row2Name = row2[0]

	#the default flag is False indicating that the default for 2 mutations to be combined is that they don't fall within acceptable ranges unless otherwise indicated later. 
	validMutationFlag = False

	#this has the new row that represents the combined mutation
	newRow = []


	#here I list all combined disruptions and restorations.
	combinedDisruptions = ','.join( sorted(row1Name.split(":")[0].split(",") + row2Name.split(":")[0].split(",")))
	combinedRestorations = ','.join( sorted(row1Name.split(":")[1].split(",") + row2Name.split(":")[1].split(",")))

	#this is the new name for the combined mutations where I merge the combined mutation disruptions and restorations together.
	newName = combinedDisruptions + ":" + combinedRestorations

	newRow.append(newName) #first item in the muattion list is the new merged name

	newRNAsegment = RNAsegment #the newRNAsegement will contain the changed mutations in the sequence (change are done using the following lines)
	disruptions = row1Name.split(":")[0].split(",") + row2Name.split(":")[0].split(",") #the comma will be used to annotate mutliple disruptions or restorations
	for i in range(0, len(disruptions)):
		disuptionLoc = int(disruptions[i].split("-")[0])
		disruptionCodon = disruptions[i].split("-")[1]

		if(nonORFrun): #f nonORF theme is used, each disruption is only 1 nt change.
			newRNAsegment = newRNAsegment[:disuptionLoc] + disruptionCodon + newRNAsegment[disuptionLoc+1:]
		else:
			newRNAsegment = newRNAsegment[:disuptionLoc] + disruptionCodon + newRNAsegment[disuptionLoc+3:]
		#print(newRNAsegment, disruptionCodon)

	
	ED, nt_NED_disrupt = calculateNED(leftOutSegment+newRNAsegment, WTstructure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption

	restorations = ""
	restoringED = 0.0

	highestED = max([float(row1[1]), float(row2[1])]) #I figure out the highest disruption NED between the 2 combined mutations
	if(ED > highestED): #only accept combined mutations where they have a disruption NED higher than the individual mutations
		
		# add to the disrupting sequence the restoring changes resulting in the restoring sequence!!!!
		restorations = row1Name.split(":")[1].split(",") + row2Name.split(":")[1].split(",")
		for i in range(0, len(restorations)):
			restorationLoc = int(restorations[i].split("-")[0])
			restorationCodon = restorations[i].split("-")[1]
			
			# print(restorationCodon)
			if(nonORFrun): #if nonORF theme is used, each restoration is only 1 nt change.
				newRNAsegment = newRNAsegment[:restorationLoc] + restorationCodon + newRNAsegment[restorationLoc+1:]
			else:
				newRNAsegment = newRNAsegment[:restorationLoc] + restorationCodon + newRNAsegment[restorationLoc+3:]

		restoringED,_ = calculateNED(leftOutSegment+newRNAsegment, WTstructure, desiredStart, desiredEnd, smp)
		
		restoringED_diff = restoringED - modelED
		# print(disruptions, restorations,restoringED, leftOutSegment+newRNAsegment)
		if(restoringED_diff <= acceptedEDrestoringDiff):

			# print(row1[1], row2[1], ED, row1[0], row2[0])
			# #print("restoringED" , restoringED, "acceptedEDrestoringDiff", acceptedEDrestoringDiff, restoringED_diff, modelED)
			# print("restoringED" , restoringED)
			
			validMutationFlag = True
			newRow.append(ED)
			newRow.append(restoringED)
			newRow = newRow + nt_NED_disrupt
	
	if(verbose): #print out mutation as it is found if verbose is selected
		print(">Found disruption at ", disruptions, "with NED =" , ED, "and restoration at ", restorations, "with NED =",restoringED, ". Restoring seq =", newRNAsegment)

	return validMutationFlag, newRow

#this function combines 2 mutations (NED), and then evaluates if it's with acceptable ranges or not
#row1: first mutation info
#row2: second mutation info
#leftOutSegment: nts found in the begining that do not belong to a full codon
#RNAsegment: WT sequence without leftOutSegment
#modelPstuct: the probability of WT sequence forming model/target structure
#disruptionFactor: factor used to determine if disruptin is successful (when disruptionP >= wtP * disruptionFactor)
#restorationFactor: factor used to determine if restoration is successful (when restorationP <= wtP * restorationFactor)
#return: T/F flag of whether combined mutation is accepted && combined mutation info.
def combineAndCalc_Pstruct(row1, row2, leftOutSegment, RNAsegment, WTstructure, modelPstruct, disruptionFactor = 50, restorationFactor = 3, nonORFrun=False, verbose=False):
	#I extract the mutation names that contain mutated codon location and mutated 3-nt triplet for both mutations to be combined
	# row1Name = row1[0]
	# row2Name = row2[0]

	#the default flag is False indicating that the default for 2 mutations to be combined is that they don't fall within acceptable ranges unless otherwise indicated later. 
	validMutationFlag = False

	#this has the new row that represents the combined mutation
	newRow = []

	#here I list all combined disruptions and restorations.
	# print(combinedDisruptions, combinedRestorations)
	combinedDisruptions = ','.join( sorted(row1[0].split(",") + row2[0].split(",")))
	combinedRestorations = ','.join( sorted(row1[2].split(",") + row2[2].split(",")))
	
	#these are the inital WT sequence to be changed later to reflect disruptions or restorations
	disruptiveSegment = RNAsegment
	restoredSegment = RNAsegment
	
	#combined mutation name that represents the combined/merged mutation
	#print(combinedDisruptions, combinedRestorations)
	tempMutName = combinedDisruptions + ":" + combinedRestorations
	#newRow.append(newName)

	#the comma will be used to annotate mutliple disruptions or restorations
	for i in range(0, len(combinedDisruptions.split(","))):
		#print(combinedDisruptions[i])
		disuptionLoc = int(combinedDisruptions.split(",")[i].split("-")[0])
		disruptionCodon = combinedDisruptions.split(",")[i].split("-")[1]

		if(nonORFrun):#if nonORF theme is used, each restoration is only 1 nt change.
			disruptiveSegment = disruptiveSegment[:disuptionLoc] + disruptionCodon + disruptiveSegment[disuptionLoc+1:]
		else:
			disruptiveSegment = disruptiveSegment[:disuptionLoc] + disruptionCodon + disruptiveSegment[disuptionLoc+3:]
		#print(newRNAsegment, disruptionCodon)
	
	#calculate probability of structure for the disruption
	disruptivePstruct = calcPstructure(leftOutSegment + disruptiveSegment, WTstructure, "disruptingSeq", disrupting=True)

	#here I calculate the disruption factor increase to help determine if it falls within acceptable ranges later.
	# disruptiveSecCodonFactorIncrease = 0
	if(disruptivePstruct != 0):
		disruptiveFactorIncrease =  modelPstruct * 1.0 / disruptivePstruct

	lowestP = min([float(row1[1]), float(row2[1])]) #check lowest disruption probability
	#here I selectd combined mutations that have a disruption factor higher that the threshold and has lower probability lower than individual disruption probabilities
	if((disruptiveFactorIncrease >= disruptionFactor) and (disruptivePstruct < lowestP)):

		#I add restoring mutations to the disruptive sequence
		restoredSegment = disruptiveSegment
		for i in range(0, len(combinedRestorations.split(","))):
			restorationLoc = int(combinedRestorations.split(",")[i].split("-")[0]) # I extract the location of the mutations
			restorationCodon = combinedRestorations.split(",")[i].split("-")[1] # I extract the mutated change (3nt triplet in ORF theme, and 1nt in nonORF theme)

			if(nonORFrun):#if nonORF theme is used, each restoration is only 1 nt change.
				restoredSegment = restoredSegment[:restorationLoc] + restorationCodon + restoredSegment[restorationLoc+1:]
			else:
				# print(restorationLoc, restorationCodon)
				restoredSegment = restoredSegment[:restorationLoc] + restorationCodon + restoredSegment[restorationLoc+3:]

		#calculate the probability of the restoring probability 
		restoredPstruct = calcPstructure( leftOutSegment+restoredSegment, WTstructure, "restoringSeq", disrupting=False)
		#check if the disruption falls within acceptable ranges		
		if( restoredPstruct >= (modelPstruct / (1.0 * restorationFactor) )):
			newRow = [combinedDisruptions, disruptivePstruct, combinedRestorations, restoredPstruct, leftOutSegment+disruptiveSegment, leftOutSegment+restoredSegment]
			# newRow = [combinedDisruptions, disruptivePstruct, combinedRestorations, restoredPstruct, disruptiveSegment, restoredSegment]
			validMutationFlag = True
	if(verbose): #print out mutation as it is found if verbose is selected
		print(">Found disruption at ", combinedDisruptions, "with probability =" , disruptivePstruct, "and restoration at ", combinedRestorations, "with probability =",restoredPstruct, ". Restoring seq =", restoredSegment)

	# print(newRow)
	return validMutationFlag, newRow


#this function tries all mutation combinations of inputted mutations when using NED as the objective function
#this function combines 2 mutational matrices; could be singleMatrices vs multiple or single vs single.. etc
#mutMatrix1 and mutMatrix2 can be of single mutations or of mutliple mutations
#leftOutSegment: nts found in the begining that do not belong to a full codon
#zeroOneMatrix1: the first mutational matrix that will be combined -  it's in zero/one mode
#zeroOneMatrix2: the second mutational matrix that will be combined -  it's in zero/one mode
#RNAsegment: WT sequence without leftOutSegment
#WTstructure: the model/target structure
#desiredStart: start location of local region to place mutations at (default is none, looking everywhere at structure)
#desiredEnd: end location of local region to place mutations at (default is none, looking everywhere)
#acceptedEDrestoringDiff: the max ED increase beyond WT ED that is acceptable
#disruptiveEDLowestDiff: the lowest threshold of ED increase from WT to be considered disrupting
def combineMutations_NED(zeronOneMatrix1, zeronOneMatrix2, leftOutSegment, RNAsegment, WTstructure, desiredStart = None, desiredEnd = None, acceptedEDrestoringDiff = 0, disruptiveEDLowestDiff = 0.3, nonORFrun=False, smp=False, verbose=False):
	modelED = zeronOneMatrix1[0][1] #here I extract the model ED. The model/wt info is always the first sublist in the mtarix

	#here I initilize a new matrix used to store the combinedMatrix information
	combinedMutMatrix = []
	combinedMutMatrix.append(zeronOneMatrix1[0]) #here I add the model/wt info as is to the combined matrix.
	triedMutations = [] #this list stored all tried mutations to help avoid redundancy
	
	####this is for combining every 2 mutations together
	###I start from 1 since the first row is the WT 
	for i in range(1, len(zeronOneMatrix1)): #I iterate over all the mutations in zeronOneMatrix1
		row1 = zeronOneMatrix1[i] #extract mutation info.

		for j in range(1, len(zeronOneMatrix2)): #I iterate over all the mutations in zeronOneMatrix2
			#111-ACG:39-GGG
			row2=zeronOneMatrix2[j] #extract mutation info.

			#the next few lines make sure that the 2 mutations change totally different codons. 

			#here I add all disrupting and restoring changes into one list
			mutationsWithCodons = row2[0].split(":")[0].split(",") + row1[0].split(":")[0].split(",") + row2[0].split(":")[1].split(",") + row1[0].split(":")[1].split(",")
			mutationList = []
			for p in range(0, len(mutationsWithCodons)):
				mutationList.append(mutationsWithCodons[p].split("-")[0])

			#here I create a set copy of the mutation list. Any redundant elements will be present only once.
			mutationSet = set(mutationList)
			# print(row1, row2)
			# print(mutationSet, mutationList)

			#here I create the new merged name that will be added to the triedMutation list.
			row1Name = row1[0]
			row2Name = row2[0]
			combinedDisruptions = ','.join( sorted(row1Name.split(":")[0].split(",") + row2Name.split(":")[0].split(",")))
			combinedRestorations = ','.join( sorted(row1Name.split(":")[1].split(",") + row2Name.split(":")[1].split(",")))
			tempMutName = combinedDisruptions + ":" + combinedRestorations
			#print(tempMutName)
			#print(triedMutations, tempMutName)
			
			#here I make sure the mutation combination has not be tried before
			if(tempMutName not in triedMutations):

				#this means that both mutations target different codons changes - if the length of the mutation list and set is the same
				if(len(mutationList) == len(mutationSet)):
					##### checkOverlap is not used. Can be used later to narrow down the number of mutations to be tried! 
					# if(checkOverlap(row1, row2)):

					# combine and calculate effects
					combinationOutput = combineAndCalc_NED(row1, row2, leftOutSegment, RNAsegment, WTstructure, modelED, desiredStart, desiredEnd, acceptedEDrestoringDiff, nonORFrun, smp=smp, verbose=verbose)
					if(combinationOutput[0]):
						combinedMutMatrix.append(combinationOutput[1])
						triedMutations.append(tempMutName)

	# this checks if no combined mutations were found
	if(len(combinedMutMatrix) == 1):
		print("No combined mutations meeting the required thresholds are found!")
		return []
	else:
		#return the converted 0/1 matrix of the combined mutation
		combinedMutMatrixBinary = convertEDMatrixToZeroAndOnes(combinedMutMatrix, disruptiveEDLowestDiff)
		return combinedMutMatrixBinary

#this function tries all mutation combinations of inputted mutations when using probability as the objective function
#this function combines 2 mutational matrices; could be singleMatrices vs multiple or single vs single.. etc
#mutMatrix1 and mutMatrix2 can be of single mutations or of mutliple mutations
#mutations1: the first mutational matrix that will be combined -  it's in zero/one mode
#mutations2: the second mutational matrix that will be combined -  it's in zero/one mode
#leftOutSegment: nts found in the begining that do not belong to a full codon
#RNAsegment: WT sequence without leftOutSegment
#modelPstuct: the probability of WT sequence forming model/target structure
#WTstructure: the model/target structure
#disruptionFactor: factor used to determine if disruptin is successful (when disruptionP >= wtP * disruptionFactor)
#restorationFactor: factor used to determine if restoration is successful (when restorationP <= wtP * restorationFactor)
def combineMutations_Pstructure(mutations1, mutations2, leftOutSegment, RNAsegment, modelPstruct, WTstructure, disruptionFactor = 50, restorationFactor = 3, nonORFrun=False, verbose=False):
	#here I initilize a new matrix used to store the combinedMatrix information
	combinedMutMatrix = []	
	triedMutations = [] #this list stored all tried mutations to help avoid redundancy

	####this is for combining every 2 mutations together
	###I start from 1 since the first row is the WT 
	for i in range(1, len(mutations1)): #I iterate over all the mutations in mutations1
		row1 = mutations1[i] #extract mutation info.
		for j in range(1, len(mutations2)): #I iterate over all the mutations in mutations2

			row2=mutations2[j] #extract mutation info.

			#the next few lines make sure that the 2 mutations change totally different codons. 

			#here I add all disrupting and restoring changes into one list
			mutationsWithCodons = row1[0].split(",") + row1[2].split(",") + row2[0].split(",") + row2[2].split(",") 
			mutationList = []
			for p in range(0, len(mutationsWithCodons)):
				mutationList.append(mutationsWithCodons[p].split("-")[0])
			
			#here I create a set copy of the mutation list. Any redundant elements will be present only once.
			mutationSet = set(mutationList)
			# print(mutationList, mutationSet)
			# print(row1, row2)
			# print(mutationSet, mutationList)

			#here I create the new merged name that will be added to the triedMutation list.
			combinedDisruptions = ','.join( sorted(row1[0].split(",") + row2[0].split(",")))
			combinedRestorations = ','.join( sorted(row1[2].split(",") + row2[2].split(",")))
			tempMutName = combinedDisruptions + ":" + combinedRestorations
			#print(tempMutName)
			#print(triedMutations, tempMutName)

			#here I make sure the mutation combination has not be tried before
			if(tempMutName not in triedMutations):
				
				#this means that both mutations target different codons changes - if the length of the mutation list and set is the same
				if(len(mutationList) == len(mutationSet)):
					#combine and calculate effects 
					combinationOutput = combineAndCalc_Pstruct(row1, row2, leftOutSegment, RNAsegment, WTstructure, modelPstruct, disruptionFactor, restorationFactor, nonORFrun, verbose)
					if(combinationOutput[0]):
						combinedMutMatrix.append(combinationOutput[1])
						triedMutations.append(tempMutName)
						# print(row1, row2, leftOutSegment, RNAsegment, WTstructure, modelPstruct, disruptionFactor, restorationFactor)
	
	# this checks if no combined mutations were found
	if(len(combinedMutMatrix) == 1):
		print("No combined mutations meeting the required thresholds are found!")
		return []
	else:
		# print(combinedMutMatrix)
		return combinedMutMatrix


#this function starts the mutational analysis using NED objective function with the option of checking alternative amino acids in restoring mutations
#this function is considered the "main" function for running ORF sequences using the NED as the objective function
#modelStructureFile: the path to the target structure file in dot-bracket format or ct format
#fullORF_file: a text file that only contains the full ORF sequence and nothing else
#desiredStart: start location of local region to place mutations at (default is none, looking everywhere at structure)
#desiredEnd: end location of local region to place mutations at (default is none, looking everywhere)
#numberOfMutations: this is the number of max mutations the user wants to explore, default is 3
#acceptedEDrestoringDiff: the max ED increase beyond WT ED that is acceptable
#disruptiveEDLowestDiff: the lowest threshold of ED increase from WT to be considered disrupting
#blosumNum: blosum matrix used to get AA alignment scores, defaults to 62
#checkAltAAs: T/F if the user chose to check blosum alternative AAs
#verbose: T/F is user wants to print mutations as they are found
def getMutations_single_multiple_NED_altAAs(path, modelStructFile, fullORF_file, desiredStart = None , desiredEnd = None, numberOfMutations =3, acceptedEDrestoringDiff = 0, disruptiveEDLowestDiff = 0.3, blosumNum=62, checkAltAAs=False, verbose=False, smp=False, quiet=False):
	
	#here I create a new temp folder where I store background files needed for calculations but not used for final output
	if not os.path.exists(path+'temp/'):
		os.makedirs(path+'temp/')
	os.chdir(path+"temp/") #this makes it the working directory

	# print(path)
	matrix = bl.BLOSUM(blosumNum, default=0) #here I extract the blosum matrix using the blosum library

	#extract full ORF sequence
	fullORF = ""
	with open(path+fullORF_file) as fp:
		fp.readline()
		for line in fp:
			fullORF+= line.strip()
	fp.close()

	#if user inputs a ct file, convert it to a dot file
	# modelStruct = ""
	if (modelStructFile.endswith(".ct")):
		result = subprocess.run(['ct2dot', path + modelStructFile, '1', path + modelStructFile[:-3]+".dot"], capture_output = True)
	modelFile = modelStructFile.split(".")[0] + ".dot"


	#extract WT sequence and model structure
	RNAsegment = ""
	Structure = ""
	with open(path+modelFile) as fp:
		header = fp.readline().strip()
		RNAsegment = fp.readline().strip()
		Structure = fp.readline().strip()
	fp.close()


	#find location of where first full codon is found
	orfFound = fullORF.find(RNAsegment) 
	orfPosition = orfFound % 3

	#this is to make sure framing is right 
	if(orfFound == -1):
		print("Your target sequence cannot be found in the provided full ORF! Please check sequences and re-run again.")
		quit()
	elif(orfPosition == 1):
		orfPosition = 2
	elif(orfPosition == 2):
		orfPosition = 1


	leftOutSegment = RNAsegment[:orfPosition] #left outsegment represents the nts at the begining of the sequence that don't form a full codon
	RNAsegment = RNAsegment[orfPosition:] #this removes any out of frame nts at the begining of the sequence

	# print(fullORF.find(RNAsegment), orfPosition,leftOutSegment)

	nutEDMatrix = [] #initilize the matrix where mutations info will be stored

	#calculate NED for model structure
	modelED, nt_NED_model= calculateNED(leftOutSegment+RNAsegment, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption

	#create the WT mutation list used for comparison later
	row = []
	row.append("WT")
	row.append(modelED)
	row.append(modelED)
	row = row + nt_NED_model


	nutEDMatrix.append(row) #append wt info to mutation matrix
	row = []
	#print(nutEDMatrix)
	if(verbose): #print out model ED if verbose is selected
		print("Model NED = " + str(modelED))

	# firstInd = 0
	# secondInd = 0

	# firstCodon = ""
	# secondCodon = ""

	# highestEDdisruption = 0
	# disruptiveSeq = ""
	# restoratingSeq = ""

	pbar = ""
	if(not quiet) and (not verbose):
		pbar = tqdm.tqdm(total=round(len(RNAsegment)/3))

	allRestoringSeqs = [] #this is used to avoid trying the same mutation combination again
	for i in range(0, len(RNAsegment), 3): #this is the first codon iteration (to check for disrupting mutations)

		if(i+3 <= len(RNAsegment)): #this is to avoid any out of context errors that I may run into at the end codon triplet if it's not fully present
			#print(i)
			codon = RNAsegment[i:i+3] #here I extract the WT codon 3-nt triplet
			codonKey = getCodon(codon) # extract codon 1 letter symbol
			# print(codon, codonKey)
			codonAlts = codons[codonKey] # here I get all codon 3-nt triplets that correspond to that codon key
			# codonAlts.remove(codon) # here I remove the wt codon 

			# print(codon, codonKey)
			for c in codonAlts: #here I iterate over all codon triplet alternatives
				firstCodonSeq = leftOutSegment + RNAsegment[:i] + c + RNAsegment[i+3:]
				# print(firstCodonSeq)

				ED, nt_NED_firstCodon = calculateNED(firstCodonSeq, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption
				
				######check if seq with c has disruptive powers on its own? If yes, continue on with the rest of the sequence!!!
				# calculate disruption difference between model ED and disrupting ED (used later to check if it's within acceptable ranges)
				disruptiveDiff = ED - modelED
				#here I only proceed if the disruption is within acceptable ranges and that the codon I tried is not 
				if((c != codon) and (disruptiveDiff >= disruptiveEDLowestDiff)):
				# if(disruptiveDiff >= disruptiveEDLowestDiff):
					for j in range(0, len(RNAsegment), 3): #this is the second codon iteration (to check for restoring mutations)
						if(j+3 <= len(RNAsegment)):
							nextCodon = RNAsegment[j:j+3] #extract 2nd WT codon triplet	to be tried for restoration					
							nextCodonKey = getCodon(nextCodon) #get the 2nd codon symbol
							# print(nextCodon, nextCodonKey)
							#here w extract the alternative codon triplets to be tried
							nextCodonAlts = []
							if(not checkAltAAs): #limit to triplets that belong to the same amino acid
								# print(nextCodonKey)
								nextCodonAlts = codons[nextCodonKey] 
							else: #include triplets that belong to the same amino acid and the conservative amino acids extracted from the blosum selection
								altAAsymbols = getAltAAs(nextCodonKey, matrix, thresholdScore=0)
								# print(nextCodonKey, altAAsymbols)
								for aa in altAAsymbols:
									nextCodonAlts += codons[aa]
							#print(nextCodon)
							#print(nextCodon, nextCodonKey, nextCodonAlts)
							for nextC in nextCodonAlts:

								if((nextC != nextCodon) and (i != j)):
									#here I create the restoring sequence 
									twoCodonSeq = ""
									if(i < j): #this makes sure that the changes are made in the correct order 
										twoCodonSeq = leftOutSegment + RNAsegment[:i] + c + RNAsegment[i+3:j] + nextC + RNAsegment[j+3:]
									else:
										twoCodonSeq = leftOutSegment + RNAsegment[:j] + nextC + RNAsegment[j+3:i] + c + RNAsegment[i+3:]

									#print(i, j, nextCodonKey, nextCodonAlts,twoCodonSeq)
									#print(twoCodonSeq)

									#this checks that the new combination has not been tried yet
									if(twoCodonSeq not in allRestoringSeqs):
										allRestoringSeqs.append(twoCodonSeq) #append the mutation to the allRestoringSeqs so that it's not tried again

										#restoration ED 
										ED0, nt_NED_twoCodonSeq = calculateNED(twoCodonSeq,Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption
											
										#get restoration difference			
										tempED_diff0 = ED0 - modelED
				 						#check if restoation is within acceptable ranges
										if(tempED_diff0 <= acceptedEDrestoringDiff):
											
											#here I check if I use the second codon to disrupt and first codon to restore (reversing mutation codons) will also produce an accpetable mutation pair
											#I need to check that it is the same amino acid becuase alternative amino acdis can only be used in restoring mutations
											#this is to make sure I can try to change codons only if the same AAs are used!
											if(getCodon(nextC) == nextCodonKey): 

												#here I create the new disrupting sequence where the second codon triplet is now tested for disruption and is considered the disrupting mutation
												secondCodonSeq = leftOutSegment + RNAsegment[:j] + nextC + RNAsegment[j+3:]
												
												ED1, nt_NED_secondCodon = calculateNED(secondCodonSeq, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption

												disruptiveDiff1 = ED1 - modelED
												if((disruptiveDiff1 >= disruptiveEDLowestDiff)):
													
													restoratingSeq = twoCodonSeq
													if(verbose): #print out mutation as it is found if verbose is selected
														# print("Switching codos for the above mutation also produced an acceptable mutation.")
														print(">Found disruption at nucleotide", j+1, "with NED =", ED1, "and restoration at nucleotide", i+1, "with NED =", ED0, ". Restoring sequence =", twoCodonSeq)
													# firstCodon = nextC
													# secondCodon = c
													# highestEDdisruption = ED1

													# disruptiveSeq = secondCodonSeq

													row.append(str(j)+"-"+ nextC+ ":"+ str(i) + "-" + c)
													row.append(ED1)
													row.append(ED0)
													row = row + nt_NED_secondCodon
													nutEDMatrix.append(row)
													row = []
													#print(disruptiveSeq, restoratingSeq)
											if(disruptiveDiff >= disruptiveEDLowestDiff):
												
												#print("switched codon ED = ", disruptiveDiff1,"disrupted ED = ", ED1 , "regular ED", disruptiveDiff)
												# restoratingSeq = twoCodonSeq

												if(verbose): #print out mutation as it is found if verbose is selected
													print(">Found disruption at nucleotide", i+1, "with NED =",ED, "and restoration at nucleotide", j+1, "with NED =",ED0, ". Restoring seq =", twoCodonSeq)
												
												row.append(str(i)+"-"+ c+ ":"+str(j) + "-"+nextC)
												row.append(ED)
												row.append(ED0)
												row = row + nt_NED_firstCodon
												
												nutEDMatrix.append(row)
												row = []


												# firstInd = i
												# secondInd = j

												# firstCodon = c
												# secondCodon = nextC

												# highestEDdisruption = ED
												# disruptiveSeq = firstCodonSeq
		if(not quiet) and (not verbose):
			pbar.update(1)

	# print(nutEDMatrix) for the flu
	# print("Number of codon changes tried = " + str(codonChangesCount))

	#here I convert the mutation ED nucleotide information to zeros and ones where 1 is disrupted nt and 0 is not 
	singleMutationMatrix = convertEDMatrixToZeroAndOnes(nutEDMatrix, disruptiveEDLowestDiff)
	#most recent matrix is another copy of the single mutation list
	mostRecentMatrix = convertEDMatrixToZeroAndOnes(nutEDMatrix, disruptiveEDLowestDiff)

	#I rank the mutations in a descending order using diruption NED, then in an asceding order using restoring NED
	signleMutSorted = sorted(singleMutationMatrix[1:], key=lambda x: (-x[1], x[2]))

	#output the single codon mutations into an output tab delimeted .txt file where each line represents a mutation
	storeRankedMutations_NED(signleMutSorted, leftOutSegment, RNAsegment, path +modelFile ,modelED, path+"1-MutationsRanked-NED.txt", orfPosition, checkAltAAs=checkAltAAs)

	#here I combine single mutations into multi-mutations if the user chooses to. The default it to try up to 3.
	for i in range(1, numberOfMutations):
		#combine single mutation matrix and most recent matrix (will be a copy of single mutations if first iteration)
		mostRecentMatrix = combineMutations_NED(mostRecentMatrix, singleMutationMatrix, leftOutSegment, RNAsegment, Structure, desiredStart, desiredEnd, acceptedEDrestoringDiff, disruptiveEDLowestDiff, smp=smp, verbose=verbose)
	
		#only combine further if the most recent matrix is not empty
		if(len(mostRecentMatrix) > 1):
			# print(str(i+1) + " codon mutations are being calculated!")
			print(str(i+1) + " codon mutations are processed!")
			multiMutSorted = sorted(mostRecentMatrix[1:], key=lambda x: (-x[1], x[2]))
			#If there was combined mutations found, ranke, and output them into a txt file
			if(multiMutSorted != 0):
				storeRankedMutations_NED(multiMutSorted, leftOutSegment, RNAsegment, path +modelFile, modelED, path+str(1+i)+"-MutationsRanked-NED.txt", orfPosition, checkAltAAs=checkAltAAs)
		else:
			# print("No more additional mutations possible!")
			print("Max mutations is " + str(i)+ ". No more additional multi-mutations possible!")

			break

	# stop = timeit.default_timer()
	# print('Time for mutliple mutations: ', stop - start)
	subprocess.run(["rm", "-rf", path+"temp/"])
	pbar.close()								


#this function starts the mutational analysis using NED objective function with the option of checking alternative amino acids in restoring mutations
#modelStructureFile: the path to the target structure file in dot-bracket format
#fullORF_file: a text file that only contains the full ORF sequence and nothing else
#numberOfMutations: this is the number of max mutations the user wants to explore, default is 3
#disruptionFactor: factor used to determine if disruptin is successful (when disruptionP >= wtP * disruptionFactor)
#restorationFactor: factor used to determine if restoration is successful (when restorationP <= wtP * restorationFactor)
#blosumNum: blosum matrix used to get AA alignment scores, defaults to 62
#checkAltAAs: T/F if the user chose to check blosum alternative AAs
#verbose: T/F is user wants to print mutations as they are found
def getMutations_single_multiple_Pstructure_altAAs(path, modelStructFile, fullORF_file, numberOfMutations =3, disruptionFactor = 50, restorationFactor = 3, blosumNum=62, checkAltAAs=False, verbose=False, quiet=False):
	# start = timeit.default_timer() ################################
	
	#here I create a new temp folder where I store background files needed for calculations but not used for final output
	if not os.path.exists(path+'temp/'):
		os.makedirs(path+'temp/')
	os.chdir(path+"temp/") #this makes it the working directory

	#here I extract the blosum matrix using the blosum library
	matrix = bl.BLOSUM(blosumNum, default=0)

	#extract orf sequence
	fullORF = ""
	with open(path+fullORF_file) as fp:
		fp.readline()
		for line in fp:
			fullORF+= line.strip()
	fp.close()

	#convert the ct file (if provided) to a dot file 
	modelStruct = ""
	if (modelStructFile.endswith(".ct")):
		subprocess.run(['ct2dot', path + modelStructFile, '1', path + modelStructFile[:-3]+".dot"], capture_output = True)
	modelFile = modelStructFile.split(".")[0] + ".dot"


	#extract the WT sequence and model structure
	RNAsegment = ""
	Structure = ""
	with open(path+modelFile) as fp:
		header = fp.readline().strip()
		RNAsegment = fp.readline().strip()
		Structure = fp.readline().strip()
	fp.close()

	#find location of where first full codon is found
	orfFound = fullORF.find(RNAsegment) 
	orfPosition = orfFound % 3

	#this is to make sure framing is right 
	if(orfFound == -1):
		print("Your target sequence cannot be found in the provided full ORF! Please check sequences and re-run again.")
		quit()
	elif(orfPosition == 1):
		orfPosition = 2
	elif(orfPosition == 2):
		orfPosition = 1

	#this contains only full codons in the sequence!
	wt = RNAsegment
	leftOutSegment = RNAsegment[:orfPosition]
	RNAsegment = RNAsegment[orfPosition:]

	# frameDotFile(path+modelFile, len(leftOutSegment))
	modelPstruct = calcPstructure(wt, Structure, "WT", disrupting=False)
	if(verbose): #print out model probability if verbose is selected
		print("Model Pstructure = " + str(modelPstruct))

	

	# codonChangesCount = 0
	# firstInd = 0
	# secondInd = 0

	# firstCodon = ""
	# secondCodon = ""

	# highestEDdisruption = 0
	# disruptiveSeq = ""
	# restoratingSeq = ""


	mutations = [] #here I store mutation pairs found to be within acceptable ranges!

	pbar = ""
	if(not quiet) and (not verbose):
		pbar = tqdm.tqdm(total=round(len(RNAsegment)/3))
	# allRestoringSeqs = []  #this is used to avoid trying the same mutation combination again
	for i in range(0, len(RNAsegment), 3): #iterate over all 3-nt triplets in sequence

		if(i+3 <= len(RNAsegment)): #this is to avoid any out of context errors that I may run into at the end codon triplet if it's not fully present
			#print(i)

			codon = RNAsegment[i:i+3] #here I extract the WT codon 3-nt triplet
			codonKey = getCodon(codon) # extract the 1 letter codon symbol representing that triplet
			#print(codon, codonKey)
			codonAlts = codons[codonKey] #get all 3-nt triplet alternatives correspoding to that key

			#codonAlts = codonAlts.pop(codonAlts.index(codon))
			#print(codon, codonKey, codonAlts)

			for c in codonAlts: #first iteration; I iterate over all alternative 3-nt triplets to check for disruption effect
				
				######check if seq with c has disruptive powers on its own? If yes, continue on with the rest of the sequence!!!
				firstCodonSeq = leftOutSegment + RNAsegment[:i] + c + RNAsegment[i+3:]
				#print(firstCodonSeq)

				#calculate disruption probability of mutant sequence forming the model structure
				disruptivePstruct = calcPstructure(firstCodonSeq, Structure, "WT", disrupting=True)

				#calculate the disruption factor increase 
				disruptiveFactorIncrease = 0
				if(disruptivePstruct != 0):
					disruptiveFactorIncrease =  modelPstruct * 1.0 / disruptivePstruct

				#make sure I don't go over the original wt 3-nt triplet 
				# I alo check that the disruption factor increase is within accpetable ranges
				if((c != codon) and (disruptiveFactorIncrease >= disruptionFactor)):

					for j in range(0, len(RNAsegment), 3): #iterate over all 3-nt triplets in sequence
						if(j+3 <= len(RNAsegment)):  #this is to avoid any out of context errors that I may run into at the end codon triplet if it's not fully present
							nextCodon = RNAsegment[j:j+3] #extract second codon 3-nt triplet				
							nextCodonKey = getCodon(nextCodon) #get its codon symbol
							#here w extract the alternative codon triplets to be tried
							nextCodonAlts = []
							if(not checkAltAAs): #limit to triplets that belong to the same amino acid
								nextCodonAlts = codons[nextCodonKey]
							else: #include triplets that belong to the same amino acid and the conservative amino acids extracted from the blosum selection
								altAAsymbols = getAltAAs(nextCodonKey, matrix, thresholdScore=0)
								# print(nextCodonKey, altAAsymbols)
								for aa in altAAsymbols:
									nextCodonAlts += codons[aa]
							# print(nextCodon, nextCodonKey, altAAsymbols,nextCodonAlts)
							#print(nextCodon)
							#print(nextCodon, nextCodonKey, nextCodonAlts)
							for nextC in nextCodonAlts: #second iteration; I iterate over all 3-nt triplets for the second codon tried
								if((nextC != nextCodon) and (i != j)):
									#here I create the restoring sequence 
									secondCodonSeq =  leftOutSegment + RNAsegment[:j] + nextC + RNAsegment[j+3:]

									disruptivePstruct_switch = 0
									disruptiveFactorIncrease_switch = 0
									#here I check if I use the second codon to disrupt and first codon to restore (reversing mutation codons) will also produce an accpetable mutation pair
									#I need to check that it is the same amino acid becuase alternative amino acdis can only be used in restoring mutations
									#this is to make sure I can try to change codons only if the same AAs are used!
									if(getCodon(nextC) == nextCodonKey):
										disruptivePstruct_switch = calcPstructure(secondCodonSeq, Structure, "WT", disrupting=True)

									if(disruptivePstruct_switch != 0):
										disruptiveFactorIncrease_switch =  modelPstruct * 1.0 / disruptivePstruct_switch

									#here I make the restoring sequence (making sure the mutated changed are in the correct order)
									twoCodonSeq = ""
									if(i < j):
										twoCodonSeq = leftOutSegment + RNAsegment[:i] + c + RNAsegment[i+3:j] + nextC + RNAsegment[j+3:]
									else:
										twoCodonSeq = leftOutSegment + RNAsegment[:j] + nextC + RNAsegment[j+3:i] + c + RNAsegment[i+3:]
									#print(i, j, nextCodonKey, nextCodonAlts,twoCodonSeq)

									#calculate restoring probability 
									restoredPstruct = calcPstructure(twoCodonSeq, Structure, "WT", disrupting=False)
									# print(modelPstruct, disruptivePstruct, restoredPstruct)

									#check if restoring probability is within acceptable ranges
									if( restoredPstruct >= (modelPstruct / (1.0 * restorationFactor) )):
											
										# restoratingSeq = twoCodonSeq
										# firstInd = i
										# secondInd = j

										# firstCodon = c
										# secondCodon = nextC

										# disruptiveSeq = firstCodonSeq

										#this mutation pair is within acceptable ranges so it's appended to final matrix
										if(verbose): #print out mutation as it is found if verbose is selected
											print(">Found disruption at nucleotide", i+1, "with probability =" , disruptivePstruct, "and restoration at nucleotide", j+1, "with probability =",restoredPstruct, ". Restoring seq =", twoCodonSeq)

										mutation = [ str(i)+"-"+c , disruptivePstruct, str(j)+"-"+nextC, restoredPstruct, firstCodonSeq, twoCodonSeq]
										mutations.append(mutation)

										if((disruptiveFactorIncrease_switch >= disruptionFactor) and (disruptiveFactorIncrease_switch != 0)):

											if(verbose): #print out mutation as it is found if verbose is selected
												print(">Found disruption at nucleotide", j+1, " with probability = " , disruptivePstruct_switch, "and restoring at nucleotide", i+1, "with probability =" , restoredPstruct, "restoring seq =", twoCodonSeq)
										
											mutation = [ str(j)+"-"+nextC , disruptivePstruct_switch, str(i)+"-"+c, restoredPstruct, secondCodonSeq, twoCodonSeq]
											mutations.append(mutation)

	#store single codon mutations 
	singleMutationMatrix = storeRankedMutation_Pstructure(mutations, leftOutSegment, path +modelFile, modelPstruct, path+"1-MutationsRanked-Pstruc.txt", RNAsegment, orfPosition, checkAltAAs=checkAltAAs)
	mostRecentMatrix = singleMutationMatrix
	
	# print(mutations) #################
	# stop = timeit.default_timer() ################
	# print('Pstructure_altAA Execution time in flu - 1 mutations: ' + str(stop - start)+ "\n") #################

	#combine mutations into multimutations if the user chooses to.
	for i in range(1, numberOfMutations):
		#combine single mutation matrix and most recent matrix (will be a copy of single mutations if first iteration)
		mostRecentMatrix = combineMutations_Pstructure(mostRecentMatrix, singleMutationMatrix, leftOutSegment, RNAsegment, modelPstruct, Structure, disruptionFactor, restorationFactor, verbose)
		# print(mostRecentMatrix, len(mostRecentMatrix))
		#only combine further if the most recent matrix is not empty
		if(len(mostRecentMatrix) > 1): #here I check if most recent matrix is not empty (given that it will always have the WT entry thus len > 1)
			#multiMutSorted = sorted(mostRecentMatrix[1:], key=lambda x: x[1], reverse=True)
			# print(str(i+1) + " codon mutations are being calculated!")
			print(str(i+1) + " codon mutations are processed!")
			if(len(mostRecentMatrix) != 0):
				storeRankedMutation_Pstructure(mostRecentMatrix, leftOutSegment, path +modelFile, modelPstruct, path+str(1+i)+"-MutationsRanked-Pstruct.txt", RNAsegment, orfPosition, checkAltAAs=checkAltAAs)
		else:
			print("Max mutations is " + str(i)+ ". No more additional multi-mutations possible!")
			break
	# stop = timeit.default_timer() ###############################
	# print('Pstructure_altAA Execution time in flu - 3 mutations: ' + str(stop - start)+ "\n") ######################
	subprocess.run(["rm", "-rf", path+"temp/"])


#this function starts the mutational analysis using NED objective function in nonORF sequences.
#modelStructureFile: the path to the target structure file in dot-bracket format or ct format
#desiredStart: start location of local region to place mutations at (default is none, looking everywhere at structure)
#desiredEnd: end location of local region to place mutations at (default is none, looking everywhere)
#numberOfMutations: this is the number of max mutations the user wants to explore, default is 3
#acceptedEDrestoringDiff: the max ED increase beyond WT ED that is acceptable
#disruptiveEDLowestDiff: the lowest threshold of ED increase from WT to be considered disrupting
#verbose: T/F is user wants to print mutations as they are found
def getMutations_single_multiple_nonORF_basepairs_NED(path, modelStructFile, desiredStart = None , desiredEnd = None, numberOfMutations =3, acceptedEDrestoringDiff = 0, disruptiveEDLowestDiff = 0.3, verbose=False, smp=False, quiet=False):
	#here I create a new temp folder where I store background files needed for calculations but not used for final output
	if not os.path.exists(path+'temp/'):
		os.makedirs(path+'temp/')
	os.chdir(path+"temp/") #this makes it the working directory

	#if user provides ct file, I convert to dot
	modelStruct = ""
	if (modelStructFile.endswith(".ct")):
		subprocess.run(['ct2dot', path + modelStructFile, '1', path + modelStructFile[:-3]+".dot"], capture_output = True)
	modelFile = modelStructFile.split(".")[0] + ".dot"


	#extract wt sequence and structure from inputted structure file
	RNAsegment = ""
	Structure = ""
	with open(path+modelFile) as fp:
		header = fp.readline().strip()
		RNAsegment = fp.readline().strip().replace("T", "U")
		Structure = fp.readline().strip()
	fp.close()


	#here I find the pairing nucleotides and their partnters
	#this is used to make sure disrupting and restoring nts are involved in base-pairing
	stack = []
	pair1 = []
	pair2 = []
	# mutStruct = list(Structure) #I convert structure string to a list of characters
	for i in range(0, len(Structure)):
		if(Structure[i] == "("):
			stack.append(i)
		elif(Structure[i] == ")"):
			pair1.append(stack[-1])
			stack.pop()
			pair2.append(i)

	nutEDMatrix = [] #this is the matrix where all the mutation pairs within acceptable ranges are stored
	modelED, nt_NED_model= calculateNED(RNAsegment, Structure, desiredStart, desiredEnd, smp) #calculate NED for model structure
	row = []
	row.append("WT")
	row.append(modelED)
	row.append(modelED)
	row = row + nt_NED_model
	nutEDMatrix.append(row) #here I store the WT/model structure information list as the first element in the mutation matrix 
	row = []
	#print(nutEDMatrix)

	if (verbose):
		print("Model NED = " + str(modelED))

	
	# firstInd = 0
	# secondInd = 0
	# firstCodon = ""
	# secondCodon = ""
	# highestEDdisruption = 0
	# disruptiveSeq = ""
	# restoratingSeq = ""

	pbar = ""
	if(not quiet) and (not verbose):
		pbar = tqdm.tqdm(total=round(len(RNAsegment)))

	allRestoringSeqs = []  #this is used to avoid trying the same mutation combination again
	#iterate over all nucleotides in the sequence
	for i in range(0, len(RNAsegment) ):
		# if(Structure[i] in ["(", ")"]):
		if(Structure[i] in ["("]): #only mutate nucleotides that fall in a base pair. Since I switch nts below. I need to only check the first base-partner
			nts = ["A", "C", "G", "U"]
			nts.remove(RNAsegment[i]) #this makes sure that I try every other possible nt other than original
			# print(nts)
			#For C, it is probably not a good idea to try U.  Because C-G to U-G will not disrupt pairing.  
			if(RNAsegment[i] == "C"):
				nts.remove("U")

			#Likewise, I would not mutate G->A if the G was in a G-U pair.  This is something to revise if you have time.
			if(RNAsegment[i] == "G"):
				nts.remove("A")

			#iterate over all possible nucleotide changes	
			for nt in nts:

				#make the first disrupting mutant 
				firstCodonSeq = RNAsegment[:i] + nt + RNAsegment[i+1:]

				#calculate the ED for that disrupting mutant
				ED, nt_NED_firstCodon= calculateNED(firstCodonSeq, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption

				#check if disruption falls within acceptable ranges
				disruptiveDiff = ED - modelED
				if(disruptiveDiff >= disruptiveEDLowestDiff):

					#here I want to explore nts to restore that are base-pairing partners with dirupting nt!!!
					pairingLoc = 0
					if( i in pair1):
						pairingLoc = pair2[pair1.index(i)]
					else:
						pairingLoc = pair1[pair2.index(i)]

					# print(pair1, pair2, i, pairingLoc)
					j = pairingLoc 
					pairingNT = RNAsegment[pairingLoc]

					nts1 = ["A", "C", "G", "U"]
					nts1.remove(pairingNT)
					# print(nts1)

					#here I iterate over the restoring nt alternatives
					for nt1 in nts1:
						
						#make the restoring sequence
						if(i < j):
							twoCodonSeq = RNAsegment[:i] + nt + RNAsegment[i+1:j] + nt1 + RNAsegment[j+1:]
						else:
							twoCodonSeq = RNAsegment[:j] + nt1 + RNAsegment[j+1:i] + nt + RNAsegment[i+1:]


						#print(twoCodonSeq)
						if(twoCodonSeq not in allRestoringSeqs): #make sure this mutation pair has not been tried before
							allRestoringSeqs.append(twoCodonSeq) # append this mutation to the allRestoringSeqs so I don't try it again

							ED0, _= calculateNED(twoCodonSeq, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption

							#restoration ED
							tempED_diff0 = ED0 - modelED
							# print(ED0, tempED_diff0)
	 						
	 						#check if restoring probability falls within acceptable ranges
							if(tempED_diff0 <= acceptedEDrestoringDiff):
								#this is the another disrupting alternative if I use the restoring nt change as disrupting
								secondCodonSeq = RNAsegment[:j] + nt1 + RNAsegment[j+1:]
								#calcuate the disrupting probability for the switched disrupted mutation
								ED1, nt_NED_secondCodon= calculateNED(secondCodonSeq, Structure, desiredStart, desiredEnd, smp) #calculate NED for merged disruption
								#check if it falls within acceptable ranges
								disruptiveDiff1 = ED1 - modelED
								if((disruptiveDiff1 >= disruptiveEDLowestDiff)):

									# restoratingSeq = twoCodonSeq
									if(verbose): #print out mutation as it is found if verbose is selected
										print(">Found disruption at nucleotide", j+1, "with NED" ,ED1, "and restoration at nucleotide", i+1, "with NED" ,ED0, ". Restoring seq =", twoCodonSeq)

									# highestEDdisruption = ED1
									# disruptiveSeq = secondCodonSeq #indicate disruption is using the second codon from the second iteration
									#append nt ED information for the disrupting mutation
									row.append(str(j)+"-"+ nt1+ ":"+ str(i) + "-" + nt)
									row.append(ED1)
									row.append(ED0)
									row = row + nt_NED_secondCodon
									
									#append the disrupting mutation found to the final mutation matrix
									nutEDMatrix.append(row)
									row = []

								if(disruptiveDiff >= disruptiveEDLowestDiff):
									
									# restoratingSeq = twoCodonSeq

									if(verbose): #print out mutation as it is found if verbose is selected
										print(">Found disruption at nucleotide", i+1, "with NED" , ED, "and restoration at nucleotide", j+1, "with NED" , ED0, ". Restoring seq =", twoCodonSeq)
									#append nt ED information for the disrupting mutation
									row.append(str(i)+"-"+ nt+ ":"+str(j) + "-"+nt1)
									row.append(ED)
									row.append(ED0)
									row = row + nt_NED_firstCodon
									#append the disrupting mutation found to the final mutation matrix
									nutEDMatrix.append(row)
									# print(row)
									row = [] #empty row
		if(not quiet) and (not verbose):
			pbar.update(1)
	#here I convert the mutation ED nucleotide information to zeros and ones where 1 is disrupted nt and 0 is not 
	singleMutationMatrix = convertEDMatrixToZeroAndOnes(nutEDMatrix, disruptiveEDLowestDiff)
	#most recent matrix is another copy of the single mutation list	
	mostRecentMatrix = convertEDMatrixToZeroAndOnes(nutEDMatrix, disruptiveEDLowestDiff)

	#I rank the mutations in a descending order using diruption NED, then in an asceding order using restoring NED
	signleMutSorted = sorted(singleMutationMatrix[1:], key=lambda x: (-x[1], x[2]))
	#output the single codon mutations into an output tab delimeted .txt file where each line represents a mutation
	storeRankedMutations_NED(signleMutSorted, "", RNAsegment, path +modelFile ,modelED, path+"1-MutationsRanked-nonORF.txt", nonORFrun=True)

	#here I combine single mutations into multi-mutations if the user chooses to. The default it to try up to 3.
	for i in range(1, numberOfMutations):
		#combine single mutation matrix and most recent matrix (will be a copy of single mutations if first iteration)
		mostRecentMatrix = combineMutations_NED(mostRecentMatrix, singleMutationMatrix, "", RNAsegment, Structure, desiredStart, desiredEnd, acceptedEDrestoringDiff, disruptiveEDLowestDiff, nonORFrun=True, smp=smp, verbose=verbose)
		#only combine further if the most recent matrix is not empty
		if(len(mostRecentMatrix) > 1):
			# print(str(i+1) + " codon mutations are being calculated!")
			print(str(i+1) + " codon mutations are processed!")
			multiMutSorted = sorted(mostRecentMatrix[1:], key=lambda x: (-x[1], x[2]))
			#print(doubleMutSorted)
			#If there was combined mutations found, rank, and output them into a txt file			
			if(multiMutSorted != 0):
				storeRankedMutations_NED(multiMutSorted, "", RNAsegment, path +modelFile, modelED, path+str(1+i)+"-MutationsRanked-nonORF.txt", nonORFrun=True)
		else:
			# print("No more additional mutations possible!")
			print("Max mutations is " + str(i)+ ". No more additional multi-mutations possible!")			
			break
	subprocess.run(["rm", "-rf", path+"temp/"])

#this function starts the mutational analysis using NED objective function in nonORF sequences.
#modelStructureFile: the path to the target structure file in dot-bracket format or ct format
#disruptionFactor: factor used to determine if disruptin is successful (when disruptionP >= wtP * disruptionFactor)
#restorationFactor: factor used to determine if restoration is successful (when restorationP <= wtP * restorationFactor)
#numberOfMutations: this is the number of max mutations the user wants to explore, default is 3
#verbose: T/F is user wants to print mutations as they are found
def getMutations_single_multiple_nonORF_basepairs_Pstructure(path, modelStructFile, disruptionFactor = 50, restorationFactor = 3,  numberOfMutations =3, verbose=False, quiet=False):

	# start = timeit.default_timer() ###################	

	#here I create a new temp folder where I store background files needed for calculations but not used for final output
	if not os.path.exists(path+'temp/'):
		os.makedirs(path+'temp/')
	os.chdir(path+"temp/") #this makes it the working directory
	
	#convert the ct file (if provided) to a dot file 
	# modelStruct = ""
	if (modelStructFile.endswith(".ct")):
		subprocess.run(['ct2dot', path + modelStructFile, '1', path + modelStructFile[:-3]+".dot"], capture_output = True)
	modelFile = modelStructFile.split(".")[0] + ".dot"

	#extract wt sequence and structure from inputted structure file
	RNAsegment = ""
	Structure = ""
	with open(path+modelFile) as fp:
		header = fp.readline().strip()
		RNAsegment = fp.readline().strip().replace("T", "U")
		Structure = fp.readline().strip()
	fp.close()

	#here I find the pairing nucleotides and their partnters
	#this is used to make sure disrupting and restoring nts are involved in base-pairing
	stack = []
	pair1 = []
	pair2 = []

	# mutStruct = list(Structure) #I convert structure string to a list of characters
	for i in range(0, len(Structure)): 
		if(Structure[i] == "("): 
			stack.append(i)
		elif(Structure[i] == ")"):
			pair1.append(stack[-1])
			stack.pop()
			pair2.append(i)
	
	#calculatr model probability 
	modelPstruct = calcPstructure(RNAsegment, Structure, "WT", disrupting=False)
	if(verbose): #print out model probability if verbose is selected
		print("Model Pstructure = " + str(modelPstruct))

	
	# firstInd = 0
	# secondInd = 0
	# firstCodon = ""
	# secondCodon = ""
	# disruptiveSeq = ""
	# restoratingSeq = ""


	mutations = [] #this is the matrix where all the mutation pairs within acceptable ranges are stored

	pbar = ""
	if(not quiet) and (not verbose):
		pbar = tqdm.tqdm(total=round(len(RNAsegment)/3))

	allRestoringSeqs = []  #this is used to avoid trying the same mutation combination again
	#iterate over all nucleotides in the sequence
	for i in range(0, len(RNAsegment) ):
		# if(Structure[i] in ["(", ")"]):
		if(Structure[i] in ["("]): #only mutate nucleotides that fall in a base pair. Since I switch nts below. I need to only check the first base-partner
			nts = ["A", "C", "G", "U"]
			nts.remove(RNAsegment[i]) #this makes sure that I try every other possible nt other than original
			# print(nts)

			#For C, it is probably not a good idea to try U.  Because C-G to U-G will not disrupt pairing.  
			if(RNAsegment[i] == "C"):
				nts.remove("U")

			#Likewise, I would not mutate G->A if the G was in a G-U pair.  This is something to revise if you have time.
			if(RNAsegment[i] == "G"):
				nts.remove("A")

			#iterate over all possible nucleotide changes	
			for nt in nts:
				#make the first disrupting mutant 				
				firstCodonSeq = RNAsegment[:i] + nt + RNAsegment[i+1:]
				#calculate probability for the disrupting mutation using the first codon disruption
				disruptivePstruct = calcPstructure(firstCodonSeq, Structure, "WT", disrupting=True)
				#extract probability factor increase
				disruptiveFactorIncrease = 0
				if(disruptivePstruct != 0):
					disruptiveFactorIncrease =  modelPstruct * 1.0 / disruptivePstruct
				#check if disrupting mutation falls within acceptable ranges
				if(disruptiveFactorIncrease >= disruptionFactor):

					###Here I want to explore nts to restore that are base-pairing partners with dirupting nt!!!
					pairingLoc = 0
					if( i in pair1):
						pairingLoc = pair2[pair1.index(i)]
					else:
						pairingLoc = pair1[pair2.index(i)]

					# print(pair1, pair2, i, pairingLoc)
					j = pairingLoc 
					pairingNT = RNAsegment[pairingLoc]

					nts1 = ["A", "C", "G", "U"]
					nts1.remove(pairingNT)
					# print(nts1)

					#here I iterate over the restoring nt alternatives					
					for nt1 in nts1:
						#make the restoring sequence						
						secondCodonSeq =  RNAsegment[:j] + nt1 + RNAsegment[j+1:]
						#calculate probability of disrupting mutation if I switch codons
						disruptivePstruct_switch = calcPstructure(secondCodonSeq, Structure, "WT", disrupting=True)
						#calculate swithced disrupting factor increatse
						disruptiveFactorIncrease_switch = 0
						if(disruptivePstruct_switch != 0):
							disruptiveFactorIncrease_switch =  modelPstruct * 1.0 / disruptivePstruct_switch


						#create the restoring sequence
						twoCodonSeq = ""
						if(i < j):
							twoCodonSeq = RNAsegment[:i] + nt + RNAsegment[i+1:j] + nt1 + RNAsegment[j+1:]
						else:
							twoCodonSeq = RNAsegment[:j] + nt1 + RNAsegment[j+1:i] + nt + RNAsegment[i+1:]
						#print(twoCodonSeq)

						#check that the restoring sequence has not been tried before 
						if(twoCodonSeq not in allRestoringSeqs):
							allRestoringSeqs.append(twoCodonSeq) #append new restoring sequence to the allRestoringSeqs so it's not tried again
							#calculate probability of restoring mutation
							restoredPstruct = calcPstructure(twoCodonSeq, Structure, "WT", disrupting=False)

							# print(modelPstruct, disruptivePstruct, restoredPstruct)
							#check if restoring mutation falls within acceptable ranges
							if( restoredPstruct >= (modelPstruct / (1.0 * restorationFactor) )):									
								# restoratingSeq = twoCodonSeq
								# firstInd = i
								# secondInd = j

								# disruptiveSeq = firstCodonSeq

								if(verbose): #print out mutation as it is found if verbose is selected
									print("> Found disruption at nucleotide", i, "with probability",disruptivePstruct, "restoration at nucleotide", j, "with probability",restoredPstruct, ". Restoring seq =", twoCodonSeq)
								
								#add mutation info. to the mutation list
								mutation = [ str(i)+"-"+nt , disruptivePstruct, str(j)+"-"+nt1, restoredPstruct, firstCodonSeq, twoCodonSeq]
								mutations.append(mutation)
								#add mutation info. to the mutation list for the switched disrupting mutation if it falls within acceptable ranges
								if(disruptiveFactorIncrease_switch >= disruptionFactor):
									if(verbose): #print out mutation as it is found if verbose is selected
										print(">Found disruption at nucleotide", j+1, "with probability",secondCodonSeq, "restoration at nucleotide", i+1, "with probability",restoredPstruct, ". Restoring seq =", twoCodonSeq)
									mutation = [ str(j)+"-"+nt1 , disruptivePstruct_switch, str(i)+"-"+nt, restoredPstruct, secondCodonSeq, twoCodonSeq]
									mutations.append(mutation)
		if(not quiet) and (not verbose):
			pbar.update(1)
	#store single codon mutations 
	singleMutationMatrix = storeRankedMutation_Pstructure(mutations, "", path+modelFile, modelPstruct, path+"1-MutationsRanked-nonORF-Pstruc.txt", RNAsegment, nonORFrun=True)
	mostRecentMatrix = singleMutationMatrix

	# stop = timeit.default_timer()
	# print('nonORF_Pstructure Execution time in VS ribozyme - 1 mutations: ' + str(stop - start)+ "\n")
	
	#here I combine single mutations into multi-mutations if the user chooses to. The default it to try up to 3.
	for i in range(1, numberOfMutations):
		#combine single mutation matrix and most recent matrix (will be a copy of single mutations if first iteration)
		mostRecentMatrix = combineMutations_Pstructure(mostRecentMatrix, singleMutationMatrix, "", RNAsegment, modelPstruct, Structure, disruptionFactor, restorationFactor, nonORFrun=True, verbose=verbose)
		#only combine further if the most recent matrix is not empty
		if(len(mostRecentMatrix) > 1):
			# print(str(i+1) + " codon mutations are being calculated!")
			print(str(i+1) + " codon mutations are processed!")
			if(len(mostRecentMatrix) != 0):
				storeRankedMutation_Pstructure(mostRecentMatrix, "", path+modelFile, modelPstruct, path+str(1+i)+"-MutationsRanked-nonORF-Pstruc.txt", RNAsegment, nonORFrun=True)
		else:
			# print("No more additional mutations possible!")
			print("Max mutations is " + str(i)+ ". No more additional multi-mutations possible!")

			break
	# stop = timeit.default_timer() ################
	# print('nonORF_Pstructure Execution time in VS ribozyme - 1 mutations: ' + str(stop - start)+ "\n") ################
	subprocess.run(["rm", "-rf", path+"temp/"])

#this function evaluates user-input mutations by calculating its NED and probability and drawing its structure 
def evaluateMutation(mutantDotFile, outputFolder, disrupting):
	#convert ct file to dot file if needed
	if (mutantDotFile.endswith(".ct")):
		result = subprocess.run(['ct2dot', outputFolder + mutantDotFile, '1', outputFolder + mutantDotFile[:-3]+".dot"], capture_output = True)
	mutantDotFile = outputFolder + mutantDotFile[:-3]+".dot"

	#extract seq and structure info from the model structure file
	with open(mutantDotFile) as fp:
		header= fp.readline().strip()
		seq = fp.readline().strip()
		structure = fp.readline()

	#create a temporary fasta file neede for partition fynction calculation	
	fastaFile = open(outputFolder + "disuptionTemp.fasta", "w")
	fastaFile.write(">mutation is " + header+"\n")
	fastaFile.write(seq+"\n")
	fastaFile.close()

	#eub partition and draw commands
	subprocess.run(['partition', outputFolder + 'disuptionTemp.fasta', outputFolder +'disuptionTemp.psf'], capture_output = True)
	subprocess.run(['draw', outputFolder +mutantDotFile[:-3]+".dot", '-p' ,outputFolder +'disuptionTemp.psf', outputFolder +mutantDotFile[:-4]+"_bpProbEffectOnModelStruct.ps"], capture_output = True)

	#calculate muatation probability
	Pstruct = calcPstructure(seq, structure, header, disrupting)
	
	#calculate mutation NED
	NED = 0
	result = subprocess.run(['EDcalculator', str(mutantDotFile)], capture_output = True, text = True)
	#NED value can be found in the line starting Structure
	for line in result.stdout.split("\n"):
		if(line.startswith("Structure")):
			temp = line.split()
			NED = round(float(temp[7]), 2)

	print("Mutation Inputted - Probability to form target structure = " + str(Pstruct) + " and NED = " + str(NED))

#this function checks if there is any common mutations between any 2 mutational output files (NED and Probability but it can also be used for any other files following our output format)
#this function is not used in the final mutate2test script. However, it can be very usefull in mutation output analysis or comparison between NED and probability
#NEDfile: path to first mutational output file to be searched
#Pstructfile: path to first mutational output file to be searched
def overlapNEDandPstrucMutations(NEDfile, Pstructfile):
	#if restored seq is the same in both NED and Pstructure, then that means all mutations are the same (even if disruptive/restorative codons are opposite)
	NEDmutations = []
	PstructMutations = []

	#here I assume restoring mutations are found in the final column
	#open NED output file and collect all restoring mutations into NEDmutationList
	with open(NEDfile) as fp:
		fp.readline()
		fp.readline()
		for line in fp:
			mut = line.split()[-1]
			NEDmutations.append(mut)
	fp.close()
	#open probability output file and collect all restoring mutations into PstructMutations
	with open(Pstructfile) as fp:
		fp.readline()
		fp.readline()
		for line in fp:
			mut = line.split()[-1]
			PstructMutations.append(mut)

	# print(NEDmutations)
	# print(PstructMutations)

	#check if there is any overlap between the 2 lists!
	mutationMutuals = [x for x in NEDmutations if x in PstructMutations] #all the x values that are in A, if the X value is in B
	#here I print the mutation index number of mutations mutual between both files/
	for i in mutationMutuals:
		print("Mutual mutations at file1 index = " + str(NEDmutations.index(i) + 1) + " and at file2 index = " + str(PstructMutations.index(i) + 1), i)

	if(len(mutationMutuals) == 0):
		print("No mututal mutations found!")

#this function checks mutational output files for any mutations found in loops
#this function is not used in the final mutate2test script. It can be used to check mutations outputted.
#outputfile: path to mutational output file to search
#modelDotFile: the path for the dotfile which contains the WT sequence and structure and will be used for color-annotating mutant base-pairs on mutations
def checkLoopMutations(outputfile, modelDotFile):
	#extract structure and sequence from model dot file
	structure = ""
	seq = ""
	with open(modelDotFile) as fp:
		header = fp.readline().strip()
		seq = fp.readline().strip()
		structure = fp.readline().strip()
	fp.close()
	#read output file to extract mutation info.
	with open(outputfile) as fp:
		fp.readline()
		fp.readline()
		for line in fp:
			temp = line.split()
			# disrupting = temp[-2]
			restoring = temp[-1]

			flag = False
			for i in range(0, len(seq)): #iterate over all nts in a sequence
				if((seq[i] != restoring[i]) and (structure[i] == ".")): #check if mutated nt is in a loop region
					flag = True
					print("Mutation found in loop is at " + str(i), line)
					break
	fp.close()

#this function calculates Probability of target structure formation for NED mutants and outputs a differnt file with a _Pstruc end with the calculated values at the last column
#this function is not used in the final mutate2test script. It can be used to check mutations outputted.
#mutFile: path to mutational output file
def calcPstructForNEDmutants(mutFile, structure):
	#open up a new output file where I can append probability for each mutation
	output = open(mutFile[:-4] + "_Pstruc.txt", "w")
	with open(mutFile) as fp:
		output.write(fp.readline())
		output.write(fp.readline().strip() + "\t disruptedPstructure \t restorativePstructure \n")
		#each line represents a mutation
		for line in fp:
			#calculate probability for both disrupting and restoring mutations
			restorativePstructure = calcPstructure(line.split()[-1], structure, "restored", disrupting=False) 
			disruptedPstructure = calcPstructure(line.split()[-2], structure, "disrupted", disrupting=True) 
			#append new probabilites to file
			output.write(line.strip()+ "\t" + str(disruptedPstructure) + "\t" + str(restorativePstructure) + "\n")
	fp.close()
	output.close()


#the function is used to run time benchmarks for the manuscript. The test cases are flu conserved model and ribozyme
#can be used for reruns
def final_benchmarkRuns():
	output = open("/Users/saraali/Desktop/LabWork/mutate2test/0-final-scripts/time.txt", "a+")

	path = "/Users/saraali/Desktop/LabWork/mutate2test/1-benchmark-flu/"
	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt",  numberOfMutations =1, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in flu without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in flu without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")


	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt",  numberOfMutations =3, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in flu without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in flu without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")


	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt",  numberOfMutations =1, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# print('NED_altAA Execution time in flu with AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# output.write('NED_altAA Execution time in flu with AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")


	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt",  numberOfMutations =3, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in flu with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in flu with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")

	
	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt", numberOfMutations = 1, disruptionFactor = 50, restorationFactor = 3, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in flu without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in flu without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt", numberOfMutations = 3, disruptionFactor = 100, restorationFactor = 3, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in flu without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in flu without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt", numberOfMutations = 1, disruptionFactor = 50, restorationFactor = 3, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in flu with AlternativeAAs - 1: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in flu with AlternativeAAs - 1: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "flu-conservedModel-original.dot", "fullORF.txt", numberOfMutations = 3, disruptionFactor = 100, restorationFactor = 3, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in flu with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in flu with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	

	path = "/Users/saraali/Desktop/LabWork/mutate2test/1-benchmark-ribozyme-1/"
	# start = timeit.default_timer()
	getMutations_single_multiple_nonORF_basepairs_NED(path, "ribozyme_seq.dot", numberOfMutations = 1, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, verbose = True)
	# stop = timeit.default_timer()
	# output.write('nonORF_NED Execution time in VS ribozyme - 1 mutation: ' + str(stop - start)+ "\n")
	# print('nonORF_NED Execution time in VS ribozyme - 1 mutation: ' + str(stop - start)+ "\n")

	# path = "/Users/saraali/Desktop/LabWork/mutate2test/1-benchmark-ribozyme/"
	# start = timeit.default_timer()
	getMutations_single_multiple_nonORF_basepairs_NED(path, "ribozyme_seq.dot", numberOfMutations = 3, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15)
	# stop = timeit.default_timer()
	# output.write('nonORF_NED Execution time in VS ribozyme - 3 mutations: ' + str(stop - start)+ "\n")
	# print('nonORF_NED Execution time in VS ribozyme - 3 mutations: ' + str(stop - start)+ "\n")
	


	# start = timeit.default_timer()
	getMutations_single_multiple_nonORF_basepairs_Pstructure(path, "ribozyme_seq.dot", disruptionFactor = 100, restorationFactor = 3, numberOfMutations = 1)
	# stop = timeit.default_timer()
	# output.write('nonORF_Pstructure Execution time in VS ribozyme - 1 mutation: ' + str(stop - start)+ "\n")
	# print('nonORF_Pstructure Execution time in VS ribozyme - 1 mutation: ' + str(stop - start)+ "\n")

	# path = "/Users/saraali/Desktop/LabWork/mutate2test/1-benchmark-ribozyme/"
	# start = timeit.default_timer()
	getMutations_single_multiple_nonORF_basepairs_Pstructure(path, "ribozyme_seq.dot", disruptionFactor = 100, restorationFactor = 3, numberOfMutations = 3)
	# stop = timeit.default_timer()
	# output.write('nonORF_Pstructure Execution time in VS ribozyme - 3 mutations: ' + str(stop - start)+ "\n")
	# print('nonORF_Pstructure Execution time in VS ribozyme - 3 mutations: ' + str(stop - start)+ "\n")
	# output.close()


	# calcPstructForNEDmutants(path+"Jan2023-NEDresults/1-MutationsRanked.txt", ".(((((.((((.(((((.......)))))....(((......)))......((((((((..(((.(((..(((....))).))).)))))))))))............)))))))))...")
	# calcPstructForNEDmutants("/Users/saraali/Desktop/LabWork/mutate2test/NEDresults/1-MutationsRanked.txt", "...((((((((.(((((.......(((.......)))....(((((((..........(((((((..(((((.(((((((.....))))))).))))).))))).)))))))))......(((.((((((.((.(((((......))))).)))))))).)))..)).))).))))))))..")
	# calcPstructForNEDmutants("/Users/saraali/Desktop/LabWork/mutate2test/NEDresults/2-MutationsRanked.txt", "...((((((((.(((((.......(((.......)))....(((((((..........(((((((..(((((.(((((((.....))))))).))))).))))).)))))))))......(((.((((((.((.(((((......))))).)))))))).)))..)).))).))))))))..")

	# overlapNEDandPstrucMutations("/Users/saraali/Desktop/LabWork/mutate2test/NEDresults/1-MutationsRanked.txt", path+"PstructResults/1-MutationsRanked.txt")
	#checkLoopMutations(path+"not-smp-najma-output.txt", path+"A2-localModel.dot")
	# checkLoopMutations("/Users/saraali/Desktop/LabWork/mutate2test/WalterMoss/Jan2023-NEDresults/1-MutationsRanked.txt",path+"conservedModel.dot")


	# stop = timeit.default_timer()
	# print('Execution time: ', stop - start)

	pass
#the function is used to design mutations for RSV conserved model. The test cases are flu conserved model and ribozyme
#can be used for reruns
def testRSV():
	output = open("/Users/saraali/Desktop/LabWork/mutate2test/0-final-scripts/time.txt", "a+")

	path = "/Users/saraali/Desktop/LabWork/mutate2test/2-test-RSV/"

	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt",  numberOfMutations =1, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in RSV without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in RSV without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")

	# structure = "...((((((((.(((((.......(((.......)))....(((((((..........(((((((..(((((.(((((((.....))))))).))))).))))).)))))))))......(((.((((((.((.(((((......))))).)))))))).)))..)).))).)))))))).."
	# calcPstructForNEDmutants(path+"1-MutationsRanked-NED.txt", structure)
	# calcPstructForNEDmutants(path+"2-MutationsRanked-NED.txt", structure)


	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt",  numberOfMutations =3, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in RSV without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in RSV without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt",  numberOfMutations =1, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# print('NED_altAA Execution time in RSV with AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# output.write('NED_altAA Execution time in RSV with AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")


	# start = timeit.default_timer()
	getMutations_single_multiple_NED_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt",  numberOfMutations =3, acceptedEDrestoringDiff=0.05, disruptiveEDLowestDiff = 0.15, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('NED_altAA Execution time in RSV with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('NED_altAA Execution time in RSV with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")


	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt", numberOfMutations = 1, disruptionFactor = 25, restorationFactor = 3, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in RSV without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in RSV without AlternativeAAs - 1 mutation: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt", numberOfMutations = 3, disruptionFactor = 25, restorationFactor = 3, blosumNum=62, checkAltAAs=False)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in RSV without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in RSV without AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt", numberOfMutations = 1, disruptionFactor = 25, restorationFactor = 3, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in RSV with AlternativeAAs - 1: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in RSV with AlternativeAAs - 1: ' + str(stop - start)+ "\n")

	# start = timeit.default_timer()
	getMutations_single_multiple_Pstructure_altAAs(path, "A2-localModel-original.dot", "RSV-L_codingRegion.txt", numberOfMutations = 3, disruptionFactor = 25, restorationFactor = 3, blosumNum=62, checkAltAAs=True)
	# stop = timeit.default_timer()
	# output.write('Pstructure_altAA Execution time in RSV with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	# print('Pstructure_altAA Execution time in RSV with AlternativeAAs - 3 mutations: ' + str(stop - start)+ "\n")
	
#this function runs an commandline interative questions to get input needed to run mutate2test by the user.
#path and model structure are inputted from the runWithFlags function.
def interactive(path, modelStructFile):
	###############Main Input####################
	print("Mutatet2test will run in a moment. These are the available options: ")
	print("(1)Output ORF structure mutations using NED as objective function.")
	print("(2)Output ORF structure mutations using Probability of structure formation as objective function.")
	print("(3)Output nonORF/non coding structure mutations using NED as objective function.")
	print("(4)Output nonORF/non coding structure mutations using Probability as objective.")
	# print("5)Evaluate a pre-exsiting mutation.")

	val = input('Please select one of the options mentioned above by entering 1, 2, 3, 4, or 5: ')
	if(val.strip() == "1"):
		#def getMutations_single_multiple_NED_altAAs(path, modelStructFile, fullORF_file, desiredStart = None , desiredEnd = None, numberOfMutations =3, acceptedEDrestoringDiff = 0, disruptiveEDLowestDiff = 0.3, blosumNum=62, checkAltAAs=False, verbose=False):
		# print("Please make sure to include the model structure dot file and the full ORF sequence file in the directory path provided below. ")
		# path = input("Enter the full path of where input files are and where output files will be stored: ")
		# modelStructFile = input("Enter the name of the dot-bracket file containing the target structure: ")
		fullORF_file = input("Enter the name of the fasta file containing the full ORF: ")
		numberOfMutations = int(input("Enter the number of codon changes in each disrupting and restoring mutation (e.g. 1, 2, 3, .. etc): "))
		
		acceptedEDrestoringDiff = float(input("Enter the restoring ED threshold (the max ED increase beyond WT ED that is acceptable, recommendation 0.05): "))
		disruptiveEDLowestDiff = float(input("Enter the disrupting ED threshold (the lowest threshold of ED increase from WT to be considered disrupting, recommendation 0.15): "))

		checkAltAAs = input("Enter Y if you would like to check alternative amino acids for restoring mutations: ")
		if(checkAltAAs == "Y"):
			checkAltAAs = True
			blosumNum = int(input("Enter number of blosum matrix to use (recommendation 62): "))
		else:
			checkAltAAs = False

		local = input("Would you like to focus mutations on a specific motif in the traget structure? Y or N? ")
		desiredStart = None
		desiredEnd = None
		if(local.upper() == "Y"):
			desiredStart = input("Enter first position index of where motif of interest starts: ")
			desiredEnd = input("Enter last position index of where motif of interest ends: ")

			if(not desiredStart.isnumeric()) and (not desiredEnd.isnumeric()):
				print("Indices provided need to be integers. Please rerun program with correct indices!")
				return
			desiredStart = int(desiredStart)-1
			desiredEnd = int(desiredEnd) - 1

		verbose = input("Enter Y if you would like the verbose option: ")
		if(verbose.upper().strip() == "Y"):
			verbose = True
		else:
			verbose = False

		getMutations_single_multiple_NED_altAAs(path, modelStructFile, fullORF_file, desiredStart = desiredStart , desiredEnd = desiredEnd, numberOfMutations =numberOfMutations, acceptedEDrestoringDiff = acceptedEDrestoringDiff, disruptiveEDLowestDiff = disruptiveEDLowestDiff, blosumNum=blosumNum, checkAltAAs=checkAltAAs, verbose=verbose)
		pass

	elif(val.strip() == "2"):
		# print("Please make sure to include the model structure dot file and the full ORF sequence file in the directory path provided below. ")
		# path = input("Enter the full path of where input files are and where output files will be stored: ")
		# modelStructFile = input("Enter the name of the dot-bracket file containing the target structure: ")
		fullORF_file = input("Enter the name of the fasta file containing the full ORF: ")
		numberOfMutations = int(input("Enter the number of codon changes in each disrupting and restoring mutation (e.g. 1, 2, 3, .. etc): "))
		disruptionFactor = float(input("Enter the multiplication factor threshold to determine if disruptin is successful (when disruption P >= wt P * disruption Factor, recommendation 100: "))
		restorationFactor = float(input("Enter the multiplication factor threshold to determine if restoration is successful (when restoration P <= wt P * restoration Factor, recommendation 3: "))

		checkAltAAs = input("Enter Y if you would like to check alternative amino acids for restoring mutations: ")
		if(checkAltAAs == "Y"):
			checkAltAAs = True
			blosumNum = int(input("Enter number of blosum matrix to use (recommendation 62): "))
		else:
			checkAltAAs = False

		verbose = input("Enter Y if you would like the verbose option: ")
		if(verbose.upper().strip() == "Y"):
			verbose = True
		else:
			verbose = False
		getMutations_single_multiple_Pstructure_altAAs(path, modelStructFile, fullORF_file, numberOfMutations =numberOfMutations, disruptionFactor = disruptionFactor, restorationFactor = restorationFactor, blosumNum=blosumNum, checkAltAAs=checkAltAAs, verbose=verbose)
		pass
	elif(val.strip() == "3"):
		# print("Please make sure to include the model structure dot file in the directory path provided below. ")
		# path = input("Enter the full path of where input files are and where output files will be stored: ")
		# modelStructFile = input("Enter the name of the dot-bracket file containing the target structure: ")
		numberOfMutations = int(input("Enter the number of codon changes in each disrupting and restoring mutation (e.g. 1, 2, 3, .. etc): "))
		acceptedEDrestoringDiff = float(input("Enter the restoring ED threshold (the max ED increase beyond WT ED that is acceptable, recommendation 0.05): "))
		disruptiveEDLowestDiff = float(input("Enter the disrupting ED threshold (the lowest threshold of ED increase from WT to be considered disrupting, recommendation 0.15): "))
	
		local = input("Would you like to focus mutations on a specific motif in the traget structure? Y or N? ")
		desiredStart = None
		desiredEnd = None
		if(local.upper() == "Y"):
			desiredStart = input("Enter first position index of where motif of interest starts: ")
			desiredEnd = input("Enter last position index of where motif of interest ends: ")

			if(not desiredStart.isnumeric()) and (not desiredEnd.isnumeric()):
				print("Indices provided need to be integers. Please rerun program with correct indices!")
				return
			desiredStart = int(desiredStart)-1
			desiredEnd = int(desiredEnd) - 1
		verbose = input("Enter Y if you would like the verbose option: ")
		if(verbose.upper().strip() == "Y"):
			verbose = True
		else:
			verbose = False
		getMutations_single_multiple_nonORF_basepairs_NED(path, modelStructFile, desiredStart = desiredStart , desiredEnd = desiredEnd, numberOfMutations =numberOfMutations, acceptedEDrestoringDiff = acceptedEDrestoringDiff, disruptiveEDLowestDiff = disrupting, verbose=verbose)
		pass
	elif(val.strip() == "4"):
		# print("Please make sure to include the model structure dot file in the directory path provided below. ")
		# path = input("Enter the full path of where input files are and where output files will be stored: ")
		# modelStructFile = input("Enter the name of the dot-bracket file containing the target structure: ")
		numberOfMutations = int(input("Enter the number of codon changes in each disrupting and restoring mutation (e.g. 1, 2, 3, .. etc): "))
		disruptionFactor = float(input("Enter the multiplication factor threshold to determine if disruptin is successful (when disruption P >= wt P * disruption Factor, recommendation 50: "))
		restorationFactor = float(input("Enter the multiplication factor threshold to determine if restoration is successful (when restoration P <= wt P * restoration Factor, recommendation 3: "))

		verbose = input("Enter Y if you would like the verbose option: ")
		if(verbose.upper().strip() == "Y"):
			verbose = True
		else:
			verbose = False
		
		getMutations_single_multiple_nonORF_basepairs_Pstructure(path, modelStructFile, disruptionFactor = disruptionFactor, restorationFactor = restorationFactor,  numberOfMutations =numberOfMutations, verbose=verbose)
		pass
	elif(val.strip() == "5"):
		# modelStructFile = input("Enter the name of the dot-bracket file with the full path containing a header, the mutant sequence, and the target structure: ")
		# path = input("Enter the full path of where input files are and where output files will be stored: ")
		disrupting = input("Is the mutatant sequence disrupting? Y or N? ")
		if(disrupting == "Y"):
			disrupting = True
		else:
			disrupting = False
		evaluateMutation(modelStructFile, path, disrupting)
		pass
	else:
		print("Invalid option selected. Please re-run and choose one of the above choices by entering: 1, 2, 3, 4, or 5!")

#this function runs an commandline flags to get input needed to run mutate2test by the user.
def runWithFlags():
	#create a new ArgumentParser object
	parser = ArgumentParser(prog='mutate2test', description='mutate2test takes an input model structure and an optional open reading frame and then outputs a set of disrupting and restoring mutations that can be used to test the model.')
	#add CTfile required flag
	parser.add_argument('CTfile', help='The full path of a CT file containing the model structure and its sequence.\nNote that a dot-bracket file can also be provided and mutate2test distinguishes between the two automatically. Make sure each file has the correct file extension, for example, ct files should have a .ct file extension.')
	# parser.add_argument('Outputfolder', help='The name of an existing directory to which output will be written.')

	#add optional flags that don't require value
	parser.add_argument("-p", "-P", "--Prob", action='store_true',help="Specify that the probability of the structure should be used as the objective function to evaluate designs.")
	parser.add_argument("-v", "-V", "--verbose", action='store_true',help="Use extra verbosity. Mutations are provided as standard out as they are found.")
	parser.add_argument("-q", "-Q", "--quiet", action='store_true',help="Suppress progress bar output.")

	parser.add_argument("-c", "-C", "--Calc", action='store_true',help="Calculate the NED and the probability for an input structure. This does not generate mutations, but evaluates a model and sequence. The output is a predicted structure and details provided to standard out. This can be helpful to evaluate the mutations provided by mutate2test or by mutations generated manually. (See notes.)")
	parser.add_argument("-f", "-F", "--flexible", action='store_true',help="When using -c, -C, --Calc, allow flexibility in base pairing. The default is to require all pairs in the model when calculating the structure probability. In other words, for the default, if a model includes a pair that cannot base pair because the nucleotides do not allow pairing, the structure probability will be 0. If -f, -F, or --flexible is specified, the probability calculation will allow pairs to be open. This is important in evaluating disrupting mutants. (See notes.)")
	parser.add_argument("-i", "-I", "--Interact", action='store_true',help="Use an interactive mode that prompts the user for required input. The default is to control the calculation using these parameters.")
	parser.add_argument("-a", "-A", "--alt", action='store_true',help="Allow restoring mutants to make conservative amino acid changes.")
	parser.add_argument("--smp", action='store_true',help="Run the parallel version of EDcalculator. EDcalculator-smp is a parallel processing version available for use on Linux and windows. It cannot be used for MacOS with an M chip.")
	#add optional flags that require added values
	parser.add_argument("-b", "-B", "--blosum", type=int, default=62,help="Specify the blosum matrix used to indentify alternative amino acids in restoring mutations if -a, -A, or --alt is used.  You can choose from BLOSUM 45, 50, 62, 80 and 90. Default is 62.")
	parser.add_argument("-o", "-O", "--ORF", default=None,help="Specify the sequence of the open reading frame (ORF) in fasta format. If the ORF is provided, mutate2test will preserve amino acid identities. The default is to assume there is no ORF. The ORF should include the full sequence of the model structure.")
	parser.add_argument("-n", "-N", "--Number", type=int,default=3,help="Specify the number of mutations that are desired for disrupting mutations (and the same number of additional mutations in the restoring structure). When an ORF is used (see -o,-O,--ORF), then this indicates the number of codon mutations (a codon mutation can contain more than one nucleotide mutation). Default is 3.")

	parser.add_argument("-d", "-D", "--Disrupt", type=float, default=None,help="Specify the threshold for disrupting mutations. The disruptions must exceed this threshold. (See notes below.). By default, the NED threshold is 0.15 and the probability (see -p, -P, --Prob) default is a factor of 100.")
	parser.add_argument("-r", "-R", "--Restore", type=float, default=None,help="Specify the threshold for restoring mutations. The mutants must be better than the threshold. (See notes below.) By default, the NED threshold is 0.05 and the probability (see -p, -P, --Prob) default is a factor of 3.")

	parser.add_argument("-s", "-S", "--Start", type=int, default=None,help="Specify the start nucleotide of the local structure to be targeted for mutations when using the NED as the objective function. The default is the 5' end, i.e. nucleotide 1")
	parser.add_argument("-e", "-E", "--End", type=int, default=None,help="Specify the end nucleotide of the local structure to be targeted for mutations when using the NED as the objective function. The default is the 3' end, i.e. nucleotide N for an N-nucleotide sequence.")
	args = parser.parse_args()

	#check if file exsists; if not, produce an error and quit
	if not os.path.exists(args.CTfile):
		parser.error("The file %s does not exist! Please make sure the correct full path to the input file in provided." % args.CTfile)
		quit()
	
	#split full path provided to path and ctfile name
	outputFolder = os.path.dirname(args.CTfile) + "/"
	inputFile = args.CTfile.split("/")[-1]
	
	#initialize disruption and restoration thresholds based on objective function
	disruptingThreshold = args.Disrupt
	if(args.Disrupt == None):
		if(args.Prob == True):
			disruptingThreshold = 100
		else:
			disruptingThreshold = 0.15

	restoringThreshold = args.Restore
	if(args.Restore == None):
		if(args.Prob == True):
			restoringThreshold = 3
		else:
			restoringThreshold = 0.05		

	#if interact is chosen, call interactive function to call commandline questions
	if(args.Interact == True):
		# interactive(args.Outputfolder, args.CTfile)
		interactive(outputFolder, inputFile)

	elif(args.Calc == True):
		evaluateMutation(inputFile, outputFolder,False)
		# evaluateMutation(args.CTfile, args.Outputfolder,False)
	elif(args.flexible == True):
		# evaluateMutation(args.CTfile, args.Outputfolder,True)
		evaluateMutation(inputFile, outputFolder,True)

	elif(args.ORF == None): #this is the the nonORF version
		if(args.Prob == False):
			# getMutations_single_multiple_nonORF_basepairs_NED(args.Outputfolder, args.CTfile, desiredStart = args.Start , desiredEnd = args.End, numberOfMutations =args.Number, acceptedEDrestoringDiff = restoringThreshold, disruptiveEDLowestDiff = disruptingThreshold, verbose=args.verbose)
			getMutations_single_multiple_nonORF_basepairs_NED(outputFolder, inputFile, desiredStart = args.Start , desiredEnd = args.End, numberOfMutations =args.Number, acceptedEDrestoringDiff = restoringThreshold, disruptiveEDLowestDiff = disruptingThreshold, verbose=args.verbose, smp=args.smp, quiet=args.quiet)
		else:
			# getMutations_single_multiple_nonORF_basepairs_Pstructure(args.Outputfolder, args.CTfile, disruptionFactor = disruptingThreshold, restorationFactor = restoringThreshold,  numberOfMutations =args.Number, verbose=args.verbose)
			getMutations_single_multiple_nonORF_basepairs_Pstructure(outputFolder, inputFile, disruptionFactor = disruptingThreshold, restorationFactor = restoringThreshold,  numberOfMutations =args.Number, verbose=args.verbose, quiet=args.quiet)
	else: #this is the ORF version
		orf = args.ORF.split("/")[-1]
		if(args.Prob == False):
			# getMutations_single_multiple_NED_altAAs(args.Outputfolder, args.CTfile, args.ORF, desiredStart = args.Start , desiredEnd = args.End, numberOfMutations =args.Number, acceptedEDrestoringDiff = restoringThreshold, disruptiveEDLowestDiff = disruptingThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose)
			# getMutations_single_multiple_NED_altAAs(outputFolder, inputFile, args.ORF, desiredStart = args.Start , desiredEnd = args.End, numberOfMutations =args.Number, acceptedEDrestoringDiff = restoringThreshold, disruptiveEDLowestDiff = disruptingThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose)
			getMutations_single_multiple_NED_altAAs(outputFolder, inputFile, orf, desiredStart = args.Start , desiredEnd = args.End, numberOfMutations =args.Number, acceptedEDrestoringDiff = restoringThreshold, disruptiveEDLowestDiff = disruptingThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose, smp=args.smp, quiet=args.quiet)
		else:
			# getMutations_single_multiple_Pstructure_altAAs(args.Outputfolder, args.CTfile, args.ORF, numberOfMutations =args.Number, disruptionFactor = disruptingThreshold, restorationFactor = restoringThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose)
			# getMutations_single_multiple_Pstructure_altAAs(outputFolder, inputFile, args.ORF, numberOfMutations =args.Number, disruptionFactor = disruptingThreshold, restorationFactor = restoringThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose)
			getMutations_single_multiple_Pstructure_altAAs(outputFolder, inputFile, orf, numberOfMutations =args.Number, disruptionFactor = disruptingThreshold, restorationFactor = restoringThreshold, blosumNum=args.blosum, checkAltAAs=args.alt, verbose=args.verbose, quiet=args.quiet)

"""
Notes:
mutate2test provides detailed reports on selected mutations. Ranked Mutations are outputted in the same directory as the inputput ct file. Each report file is in a space delimited .txt format.
* Using NED objective function, mutations are ranked in a descending order based on disrupting NED. Mutations with the same disruption NED, are ranked further in an ascending order using restoring NED. We also output another report file where we rank mutations in a descending order based on the total number of disrupted nucleotides in each disrupting mutation. 
* Similarly, using the probability objective function, we rank mutations based on disrupting probability in an ascending order and for groups with the same disrupting probability, we rank them further using restoring probability in a descending order. 

mutate2test outputs structure drawings for the top ranked mutation for each file report in the topHitDrawings folder created in the same directory as the input ct file.
For each mutation pair we draw the following structures:
* disrupting mutant predicted maximum expected accuracy structure (e.g. 1topDisruption_MEA.ps where the number represents the number of codons changed for disruption)
* disrupting mutant used to draw the model structure color-annotated using the base-pairing probabilites calculated using the mutant sequence (e.g. 1topDisruption_bpProbEffectOnModelStruct where the number represents the number of codons changed for disruption)
* restoring mutantpredicted maximum expected accuracy structure (e.g. 1topRestoration_MEA.ps where the number represents the number of codons changed for disruption)
* restoring mutant used to draw the model structure color-annotated using the base-pairing probabilites calculated using the mutant sequence (e.g. 1topRestoration_bpProbEffectOnModelStruct where the number represents the number of codons changed for disruption)

P.S. When using the NED function, there will be 2 subfolders in topHitDrawings: 
* totalNEDrank: this shows the top ranked mutation structure drawings based on the total disrupting NED
* nt-NEDrank: this shows the top ranked mutation structure drawings based on the number of disrupted nucleotides


mutate2test requires installing the following Python libraries: biopython, tqdm, and blosum (the 2.0.2 release). A package manager can be used to install them. For example, pip can be used using the following command: pip install biopython

mutate2test calls components of the RNAstructure package. Make sure to add the RNAstructure executables to your global path. The following steps show how to do so in a Mac terminal.
1- Use a command line-based text editor like nano to open up the bash profile. For example, use "nano ~/.bashrc"
2- Add this line to the file: export path/to/RNAstructure/exe/:${PATH}
3- Save changes and source the changed bashrc profile. You can use: source ~/.bashrc


The execution time and the number of mutations outputted by each metric might vary when using NED or Probability as the objective functions, therefore, it is important to adjust Accepting/Rejecting Mutation thresholds to best fit the user’s needs. Depending on the number of mutations outputted, the user may increase or decrease the disrupting and the restoring thresholds as needed.
* The NED disruption threshold ensures that the disrupting mutations are selected such that their NED has to be greater than or equal to the sum of the model NED and the disruption threshold.
* The NED restoration threshold means that the restoring mutations are only accepted if they have an NED lower than or equal to the sum of the model NED and the restoration threshold.
* The recommended start NED values are an NED disruption increase threshold of 0.15 and an NED restoration difference threshold of 0.05.

* The probability disruption threshold ensures that the disrupting mutations are only accepted if they have a probability lower than or equal to the model probability multiplied by the disruption factor.
* The probability restoration threshold determines that the restoring mutations are to be accepted only if they have a probability greater than or equal to the model probability divided by the restoration factor.
* I recommend starting with a probability factor threshold of 100 and a restoration factor threshold of 3.

There are two modes where mutatet2test can used to evaluate user-inputted mutations.
1- The exact base-paring mode (using -c, -C, --Calc) ensures that if a model includes a pair that cannot base pair because the mutated nucleotides do not allow pairing, the structure probability will be 0. This can be used to evaluate wildtype sequence and restoring mutants.
2- The flexible mode (using -f, -F, --flexible) allows flexibility in base pairing. This indicates that the probability will be calculated after any base pairs no longer pairing given the sequence are eliminated. This can be used to evaluate disrupting mutants where the sequence may have been mutated to eliminate base pairs in the model.

The output folder will be in the same directory as the input CT file.
Please make sure the ORF fasta file is in the same directory as the input CT file.
"""
def main():
	# start = timeit.default_timer()
	getCodonTable(os.path.dirname(__file__)+"/codonTable.txt")
	runWithFlags()
	# final_benchmarkRuns()
	# testRSV()
	# stop = timeit.default_timer()
	


if __name__ == "__main__":
    main() 

