#include "GeneticAlgorithm.h"

#define OUTPUT
#define ITERATIONS_PER_SAVE 1 // Save the state file after this many iterations.

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <cassert>
#include "../src/random.h"

using std::cout;
using std::cerr;
using std::string;
using std::vector;

bool all_digits(const string& s){
	bool ret = true;
	for(int i=0;i<s.size()-1;i++){
		if(!isdigit(s[i]))
			ret=false;
	}
	return s.size()==0?false:ret;
}

// Save the current state of the calculation, so it can be resumed later.
void saveStateFile(string filename, RNA_Container* r, int rnaCount) {
	if (!filename.empty()) toFile(filename, r, rnaCount);
}


int orega(std::string InputFile,
	std::string OutputFile,
	std::string savefile,    // file where current state should be saved so it can be resumed later.
	std::string restartFile, // previous savefile that should be opened to resume or extend calculation.
	int NumIterations,
	const char* const alphabet,
    int MutationStartNuc,//position at which to start nucleotide mutations
    int NumberOfNucs,
    double MutationRate,
    int RecombinationFrequency,
    double RecombinationRate,
    int NumberOfSequences,
	int objectiveFunction,
    int randomseed,
	int MutationSwitch,
	int limitG,
	int filterAUG,
	int filterCUG,
	int filteroligoA,
	double complexity_constant,
	double threshold
	)
{

	// if(argc!=12){
	// 	cerr << "Usage error. Please input: 1:[GeneticAlgorithm] 2:[sequence_file.seq] 3:[random_number_seed] 4:[RNA=1 or DNA=2] ";
	// 	cerr << "5:[Number_Of_Sequences] 6:[Nuc_#_for_mutation_start] 7:[Number_of_nucleotides] 8:[Number_of_iterations] ";
	// 	cerr << "9:[Recombinations_Frequency] 10:[Mutation_Rate] 11:[Recombination_Rate] 12:[Output_file_name]\n";
	// 	return 1;
	// }
// #ifdef OUTPUT
// 	cerr << "Reading user parameters...";
// #endif

	// MutationStartNuc=atoi(argv[5]);//read the start position for mutations
	// NumberOfNucs=atoi(argv[6]);//read the number of nucleotides to mutate
	// int RecombinationFrequency=atoi(argv[8]);
	// MutationRate=atof(argv[9]);//store the rate at which a nucleotide will be mutated
	// RecombinationRate=atof(argv[10]);//store the rate at which a nucleotide will be recombined
	// string savefile("");
	// bool resume = false;

	//Read if the sequence is an RNA or DNA
	// bool IsRNA;
	// if(atoi(argv[3])==1) IsRNA=true;
	// else if(atoi(argv[3])==2) IsRNA=false;
	// else{
	// 	cerr << "Usage error! Input #4 (RNA or DNA) is incorrect. Please input 1 for RNA or 2 for DNA.\n";
	// 	return 1;
	// }
	//END: RNA or DNA

	//seed the class to generate random number
	// string seedIN=argv[2];
	// for(int i=0;i<seedIN.size();++i){
	// 	if(!isdigit(seedIN[i])){
	// 		cerr << "Usage errror! Input #3 (random_number_seed) is incorrect. Please input numbers only.\n";
	// 		return 1;
	// 	}
	// }
	
	//Declare the instance of the class to generate random number
	randomnumber randomnum;
	//int seed=atoi(seedIN.c_str());
	randomnum.seed(randomseed);
	//srand(seed);
//	for(int i=0;i<100;i++) cout<<randomnum->roll_int(4,5)<<std::endl;
	//END: seed the class to genereate random numbers

// #ifdef OUTPUT
// 	cerr << "\t\tDONE\nCalculating wildtype probabilities...";    
// #endif

	//Initialize RNA to calculate the wild type probability RMSD
	RNA *rnaWT = new RNA(InputFile.c_str(), FILE_SEQ, alphabet);

	//Check for errors while reading the structure
	int error=rnaWT->GetErrorCode();
	if(error!=0){
		cerr << rnaWT->GetErrorMessage(error) << endl;
		delete rnaWT;
		return 1;
	}

	//Check if the number of nucleotides and starting position don't exceed the sequence length.
	if((NumberOfNucs+MutationStartNuc-1)>rnaWT->GetSequenceLength()){
		cerr << "\nThe start + length parameters specify a target segment that extends past the length of the given sequence.\n";
		return 1;
	}

	//Run the partition function on rnaWT
	//***mohammad	rnaWT->SetTemperature(t1);
	error=rnaWT->PartitionFunction();
	

	//Check for errors while reading the structure
	if(error!=0){
		cerr << rnaWT->GetErrorMessage(error) << endl;
		delete rnaWT;
		return 1;
	}    

//#ifdef OUTPUT
	cout << "Initializing mutant RNA sequences...";
//#endif

	//RNA_Container provides an array of RNA classes
	RNA_Container *rna;

	//Read the number of parallel sequences to mutate
	//int NumberOfSequences=atoi(argv[4]);

	//For each element in rna, the RNA_Container array, allocate the underlying RNA class.
	int NumberOfStructureElements=NumberOfSequences*2;

	

	if(restartFile.empty()){
		rna = fromSeqFile(InputFile.c_str(),NumberOfStructureElements,rnaWT);
	}
	else{//we are resuming a previous run
		rna = fromSaveFile(restartFile,rnaWT);
	}
	//Vector that stores the computed RMSDs
	//***mohammad	vector<vector<double> > CalculatedRMSDs(NumberOfStructureElements, vector<double>(3,0));
	vector<double> CalculatedFFTs(NumberOfStructureElements, 0);
	vector<double> mean_bp_prob(NumberOfStructureElements, 1);

#ifdef OUTPUT
	//cerr << "\tDONE\nCalculating RMSD at 39C...";
#endif


	//Compute the initial RMSD at 39C
	//Change the temperature to 39
	//rna[NumberOfSequences].Return_RNA()->SetTemperature(t2);
	//Run Partition Function
	rna[NumberOfSequences].Return_RNA()->PartitionFunction();
	//Compute the initial RMSD
	//double RMSD39=CalculateRMSD(rnaWT, rna[NumberOfSequences].Return_RNA());

	 double FFTinit=CalculateFFT(rna[NumberOfSequences].Return_RNA(), MutationStartNuc, NumberOfNucs, objectiveFunction,filterAUG, filterCUG, filteroligoA,complexity_constant,NumberOfNucs, threshold);
	 double meanBPinit= totalPairwiseProbability(rna[NumberOfSequences].Return_RNA(), MutationStartNuc, NumberOfNucs) / ((double)NumberOfNucs);
	//Store the initial RMSDs
	for(int i=0;i<NumberOfSequences;++i){
		//cerr << "\nRMSD39=" << RMSD39 << endl;
		//CalculatedRMSDs[i][1]=RMSD39;
		CalculatedFFTs[i]=FFTinit;
		mean_bp_prob[i] = meanBPinit;
	}

//#ifdef OUTPUT
	cout << "\t\tDONE\nRunning OREGA calculations:\n";
//#endif

	// THE FOR LOOP STARTS HERE
	for(int iterations=1;iterations<=NumIterations; ++iterations){

//#ifdef OUTPUT
		cout << "\tITERATION " << iterations <<std::endl;
//#endif
		//Run "Recombine" every predetermined number of steps
		if(iterations % RecombinationFrequency == 0){
			//Vector that stores the pair of sequences that will be recombined with each other
			vector<vector<int> > RecombinationPair (floor(NumberOfSequences/2),vector<int>(2,-1));
			//Assign a unique recombination pairs
			for(int n1=0;n1<RecombinationPair.size();++n1){
				for(int n2=0;n2<RecombinationPair[n1].size();++n2){
					//Keep track if the random number is a repeat
					bool Repeat=true;
					//Make sure that the sequence has not been marked for a recombination multiple times
					while(Repeat){
						double random=randomnum.roll();
						random=floor(random*NumberOfSequences);
						//cerr << "\n!!!random=" << random << endl;
						if(!CheckForRepeat(RecombinationPair,random)){
							Repeat=false;
							RecombinationPair[n1][n2]=random;
						}
					}

				}
			}
			//Recombine each pair of sequneces
			int storageindex=NumberOfSequences;//keeps track of indexing
			for(int pair=0;pair<RecombinationPair.size();++pair){

				Recombine(rna[RecombinationPair[pair][0]].Return_RNA(), rna[RecombinationPair[pair][1]].Return_RNA(), 
					rna[storageindex].Return_RNA(), rna[storageindex+1].Return_RNA(), NumberOfNucs, MutationStartNuc, RecombinationRate, &randomnum);
				storageindex=storageindex+2;      
			}
		}//End: run recombination

		//Else run the regular single nucleotide mutation
		else{
			//Mutate a nucleotide in every Sequence and calculate the RMSD.
			for(int i=NumberOfSequences;i<NumberOfStructureElements;++i){
			//Copy the permanent sequence to the temporary location
				CopySequence(rna[i].Return_RNA(),rna[i-NumberOfSequences].Return_RNA());
			//Mutate the nucleotide
			
			
				int MutateOut=MutateNuc(rna[i].Return_RNA(), MutationStartNuc, NumberOfNucs, MutationRate, MutationSwitch, limitG, &randomnum,filterAUG,filterCUG,filteroligoA);
				 //If the output of the 'MutateOut' function is 2="err", this means that an unknown nucleotide has been enrountered
				if(MutateOut==2){
					cerr << "\n!!!An unknown nucleotide has been encountered!!!\n";
					delete rna;
					delete rnaWT;
					return 1;        
				}
				}
			
		}
		//###################################
		// PARALLELL EXECUTION STARTS HERE 
		//##################################

//#pragma omp parallel for
		//Now calculate the RMSDs
		//Mutate a nucleotide in every Sequence and calculate the RMSD.
		for(int i=NumberOfSequences;i<NumberOfStructureElements;++i){
			//Change the temperature to 34
			//***mohammad	rna[i].Return_RNA()->SetTemperature(t1);
			//Run Partition Function
			rna[i].Return_RNA()->PartitionFunction();
			//Calculate and store the RMSD

			 CalculatedFFTs[i]=CalculateFFT(rna[i].Return_RNA(), MutationStartNuc, NumberOfNucs, objectiveFunction,filterAUG, filterCUG, filteroligoA,complexity_constant,NumberOfNucs,threshold);
			 mean_bp_prob[i] = totalPairwiseProbability(rna[i].Return_RNA(), MutationStartNuc, NumberOfNucs) / ((double)NumberOfNucs);
			//CalculatedRMSDs[i][0]=CalculateRMSD(rnaWT, rna[i].Return_RNA());
			//Change the temperature to 39
			//rna[i].Return_RNA()->SetTemperature(t2);
			//Run Partition Function
			//rna[i].Return_RNA()->PartitionFunction();
			//CalculatedRMSDs[i][1]=CalculateRMSD(rnaWT, rna[i].Return_RNA());
		}
		//##################################
		// PARALLELL EXECUTION ENDS HERE 
		//##################################

		

		//Initialize vector to store the order of the lowest RMSD
		vector<int> LowestRMSDorder (NumberOfSequences,0);
        vector<bool> KeepSequences (NumberOfStructureElements,false);

        for(int i=0;i<NumberOfStructureElements;++i){
            int numberBelow=0;
            int numberAbove=0;
            for(int j=0;j<NumberOfStructureElements;++j){
                if(j!=i){
                    if(CalculatedFFTs[i]>CalculatedFFTs[j]) numberBelow++;
                    else if(CalculatedFFTs[i]<CalculatedFFTs[j]) numberAbove++;
                    else if(CalculatedFFTs[i]==CalculatedFFTs[j]){
                        if(!KeepSequences[j]){
                            numberBelow++;
                            numberAbove++;
                        }
                    }
                }
            }
            //cerr << "\nNumberBelow=" << numberBelow << " NumberAbove " << numberAbove;
            if(numberBelow>=NumberOfSequences) KeepSequences[i]=true;
            //if(numberAbove>=NumberOfSequences) KeepSequences[i]=true;
        }

		//Perform the sequence swap
		int PositionToSwap1=0;
		int PositionToSwap2=NumberOfStructureElements-1;
		while(!AllSwapped(KeepSequences,NumberOfSequences)){

			while(KeepSequences[PositionToSwap1]) PositionToSwap1++;
			while(!KeepSequences[PositionToSwap2]) PositionToSwap2--;

			CopySequence(rna[PositionToSwap1].Return_RNA(),rna[PositionToSwap2].Return_RNA());
            CalculatedFFTs[PositionToSwap1]=CalculatedFFTs[PositionToSwap2];
			mean_bp_prob[PositionToSwap1]=mean_bp_prob[PositionToSwap2];
			KeepSequences[PositionToSwap1]=true;
			KeepSequences[PositionToSwap2]=false;
		}

		// //Output the structures
		// for(int i=0;i<NumberOfSequences;++i){
		// 	if(i==0){
		// 		//stringstream commentSS;
		// 		//commentSS << " ##->## RMSD@34C = " << CalculatedRMSDs[i][0] << " RMSD@39C = " << CalculatedRMSDs[i][1];
		// 		//rna[i].Return_RNA()->AddComment(commentSS.str().c_str());
		// 		rna[i].Return_RNA()->WriteCt(argv[11],false);
		// 	}
		// 	else rna[i].Return_RNA()->WriteCt(argv[11],true);
		// }
		//Output the structures
        for(int i=NumberOfSequences;i<NumberOfStructureElements;++i){
        	double complexity_value=0;
        	string small_seq=rna[i].Return_RNA()->GetSequence(MutationStartNuc, NumberOfNucs);
    		complexity_value=complexity_calc(small_seq,filterAUG, filterCUG, filteroligoA);
			double seg_bp_fraction = Basepair_fraction(rna[i].Return_RNA(),MutationStartNuc, NumberOfNucs,threshold) / ((double)NumberOfNucs);
    		

			//cerr << "ITERATION: " << iterations << " FFT" << i << "=" << CalculatedFFTs[i] << "	"<<"Complexity=	"<<complexity_value<<"	"<<small_seq<<endl;
            cout<<sfmt("Sequence %d:\t Score = %.4f \tMean Basepair Probability = %.4f\tComplexity = %.3f \t Basepair fraction: %.2f \t Optimized Segment : %s", i-NumberOfSequences+1, CalculatedFFTs[i], mean_bp_prob[i], complexity_value,seg_bp_fraction, small_seq.c_str())<<endl;
        }

		if (iterations % ITERATIONS_PER_SAVE == 0)
			saveStateFile(savefile, rna, NumberOfSequences);
    }

	// Save the final state file.
	saveStateFile(savefile, rna, NumberOfSequences);

	// Output the best sequence
	std::string seq = rna[0].Return_RNA()->GetSequence(MutationStartNuc, NumberOfNucs);
	double final_complexity_value=complexity_calc(seq,filterAUG, filterCUG, filteroligoA);
	double final_avg_bp = 0.0; //TODO: Calculate final average basepair probability for the optimized region
	double final_avg_fraction = 0.0;
	rna[0].Return_RNA()->PartitionFunction();
	final_avg_bp=totalPairwiseProbability(rna[0].Return_RNA(), MutationStartNuc, NumberOfNucs)/ ((double)NumberOfNucs);
    final_avg_fraction = Basepair_fraction(rna[0].Return_RNA(), MutationStartNuc, NumberOfNucs,threshold) / ((double)NumberOfNucs);
	const int FASTA_TYPE = 2;
	std::string label = rnaWT->GetCommentString() + sfmt("\tMean Basepair Probability: %.4f\t\tComplexity: %.3f\t\tBasepair fraction: %.2f", final_avg_bp, final_complexity_value, final_avg_fraction);
	//cout<<sfmt("\tMean Basepair Probabilty: %.4f\t\tComplexity: %.4f", final_avg_bp, final_complexity_value);
	rna[0].Return_RNA()->SetSequenceLabel(label);
	rna[0].Return_RNA()->GetStructure()->writeseq(OutputFile.c_str(), FASTA_TYPE, false);

		// for(int i=0;i<NumberOfSequences;++i){
		// 	if(i==0){
		// 		//stringstream commentSS;
		// 		//commentSS << " ##->## RMSD@34C = " << CalculatedRMSDs[i][0] << " RMSD@39C = " << CalculatedRMSDs[i][1];
		// 		//rna[i].Return_RNA()->AddComment(commentSS.str().c_str());
				
		// 	}
		// 	else rna[i].Return_RNA()->WriteCt(argv[11],true);
		// }	
    
    cout << "\nOutput:   "<<endl;

    //Output the structures
    for(int i=0;i< NumberOfSequences;++i){
    	rna[i].Return_RNA()->PartitionFunction();
		double complexity_value=0;
    	string small_seq=rna[i].Return_RNA()->GetSequence(MutationStartNuc, NumberOfNucs);
		complexity_value=complexity_calc(small_seq,filterAUG, filterCUG, filteroligoA);
		double mean_bp_prob=totalPairwiseProbability(rna[i].Return_RNA(), MutationStartNuc, NumberOfNucs)/((double) NumberOfNucs);
        double bp_fraction = Basepair_fraction(rna[i].Return_RNA(), MutationStartNuc, NumberOfNucs,threshold)/((double)NumberOfNucs);
		//cerr << "ITERATION: " << iterations << " FFT" << i << "=" << CalculatedFFTs[i] << "	"<<"Complexity=	"<<complexity_value<<"	"<<small_seq<<endl;
        cout<<sfmt("Final Sequence %d:\t Score = %.4f \tMean Basepair Probabilty = %.4f\tComplexity = %.3f \t Basepair fraction: %.2f \t Optimized Segment : %s", i+1, CalculatedFFTs[i], mean_bp_prob, complexity_value, bp_fraction, small_seq.c_str())<<endl;
        //cerr << "FINAL FFT" << i << "=" << CalculatedFFTs[i] << endl;
        //cout << "Final score " << i+1 << "=" << CalculatedFFTs[i] << endl;
    }

	delete[] rna;
	delete rnaWT;
	return 0;
}
//Calculate the linguistic complexity of a sequence. Reported value is between 0 and 1. 
double complexity_calc(const string& sequence, int filterAUG, int filterCUG, int filteroligoA){

	string input_segment=sequence;
	//Changing all nucleotides in the sequence to upper case letters
	std::transform(input_segment.begin(), input_segment.end(),input_segment.begin(), ::toupper);

	for (int i=0; i<input_segment.length(); i++){
		if (input_segment[i]=='T')
			input_segment[i]='U';
		if (input_segment[i]!='A' && input_segment[i]!='C' && input_segment[i]!='G' && input_segment[i]!='U')
			cerr<<endl<<"Only standard A,C,G,U/T nucleotides are allowed. Please check the input sequence."<<endl;
	}

	//Define the initial complexity parameter   C=U1*U2*U3*...*Uw
    double C=1;  

    int length=input_segment.length();
    //W parameter is the size of oligomers that are being evaluated. 
    //It depends on the length of the input segment 
    int W=0;
    if (length<18)
        W=3;
    else if (length<67)
        W=4;
    else if (length<260)
        W=5;
    else if (length<1029)
        W=6;
    else if (length<3910)
        W=7;
    else 
    	W=7;
    
    int size_olig[7];
    for (int i=1; i<=7; i++){
        int num1=pow(4, i);
        int num2=length-i+1;
        if(num1<num2)
            size_olig[i-1]=num1;
        else
            size_olig[i-1]=num2;
    }

    int olig_count[7];
    for(int i=0; i<7; i++)
        olig_count[i]=0;

    string **oligs;
    oligs = new string* [7];
    for (int i=0; i<7; i++)
        oligs[i]=new string[size_olig[i]];
    

    for(int i=1; i<=W; i++){
        for(int j=0; j<length-i+1; j++){
            string short_frag=input_segment.substr(j, i);
            bool uniq=true;
            int save_indx=0;
            int uniq_index=0;
            for(int k=0;k<size_olig[i-1]; k++){

                if(short_frag==oligs[i-1][k]){
                    uniq=false;
                    break;
                }
                if(oligs[i-1][k]!="")
                    uniq_index++;
            }
            if(uniq){
                olig_count[i-1]++;
                oligs[i-1][uniq_index]=short_frag;
            }                  
        }
    }

    for(int i=0; i<W; i++)
        C*=(double)olig_count[i]/size_olig[i];

    for (int i=0; i<7; i++)
        delete[] oligs[i];
    delete [] oligs;
	if(filteroligoA > 0){
		double polyA = 0;
		for(int i=1; i<=W; i++){
			for(int j=0; j<length-i+1; j++){
				string short_frag=input_segment.substr(j,i);
				if (short_frag =="A"){
					if (j != 0){
						if (input_segment.substr(j-1,i) == "A"){
							polyA++;
						}
						else if (polyA < filteroligoA)
						{
							polyA = 0;
						}
						else
						{
							polyA = polyA;
						}
						
					}
					else{
						polyA++;
					}
				}
			}
		}
		if (polyA >= filteroligoA){
			C = -100000;
		}
	}
    return C;
    
}

double totalPairwiseProbability(RNA* rna, int start, int length){
	double total = 0;
	for (int i=start;i<start+length;++i) {
		// All j before i
		for (int j=1; j<i; ++j)
			total+=rna->GetPairProbability(j, i);
			
		// All j after i
	    	for (int j=i+1; j<=rna->GetSequenceLength(); ++j)
    			total+=rna->GetPairProbability(i, j);
				
	}
	return total;
}

double Basepair_fraction(RNA* rna, int start, int length,double threshold){
	double total = 0;
	for (int i=start;i<start+length;++i) {
		double basepair_probability = 0;
		// All j before i
		for (int j=1; j<i; ++j)
			basepair_probability+=rna->GetPairProbability(j, i);
			
		// All j after i
	    	for (int j=i+1; j<=rna->GetSequenceLength(); ++j)
    			basepair_probability+=rna->GetPairProbability(i, j);
				
			if (basepair_probability >= threshold) {
				total += 1;
			}
	}
	return total;
}
double CalculateFFT(RNA* rna, int mutation_start, int number_of_nucs, int objectiveFunction, int filterAUG, int filterCUG, int filteroligoA,double complexity_constant,int total_length,double threshold){

    //double probs[2048]={0};

    //****for now.. come back later for this part;
    // int seq_len=512;
    // int number_of_5C_repeats=0;
    // int number_of_5A_repeats=0;
    // int number_of_5U_repeats=0;
    // int number_of_CA_repeats=0;

    // int number_of_3CA_repeats=0;
    // int number_of_CAA_repeats=0;

    double sum_53_utr_probs=0;
/*
    for (int i=mutation_start-10;i<=mutation_start+number_of_nucs;++i) {
    	for (int j=1221; j<(1221+106); j++){
    		sum_53_utr_probs+=rna->GetPairProbability(i, j);
    	}
    }
*/

	for (int i=1;i<=rna->GetSequenceLength();++i) {
		// All j before i
		for (int j=1; j<i; ++j)
			sum_53_utr_probs+=rna->GetPairProbability(j, i);
		// All j after i
    	for (int j=i+1; j<=rna->GetSequenceLength(); ++j)
    		sum_53_utr_probs+=rna->GetPairProbability(i, j);
	}

/*    
	for (int i=1;i<=106;++i) {
    	for (int j=mutation_start; j<mutation_start+number_of_nucs; j++){
    		sum_53_utr_probs+=rna->GetPairProbability(i, j);
    	}
	}
*/

    // for (int i=1;i<=rna->GetSequenceLength();++i) {
    //     char nuc=rna->GetNucleotide(i);

    //     PFPRECISION probability=0;

    //     //check all possible pairing partners
    //     for(int j=1;j<i;++j) {
    //     //  cout<<i<<"  "<<j<<" "<<ct->GetPairProbability(j,i)<<endl;   
    //         probability += rna->GetPairProbability(j,i);

    //     }
    //     for(int j=i+1;j<=rna->GetSequenceLength();++j) {
    //     //  cout<<i<<"  "<<j<<" "<<ct->GetPairProbability(i,j)<<endl;
    //         probability += rna->GetPairProbability(i,j);

    //     }
    //     probs[i-1]=probability;
    //     //write the probability to file:
    // }

    // double utr_prob_sum=0;

    // for (int i=mutation_start; i<(mutation_start+number_of_nucs); i++){
    //     char nuc=rna->GetNucleotide(i);
    //     if (nuc=='C'){
    //         int C_counter=1;
    //         int j=i;
    //         while(rna->GetNucleotide(j+1)=='C'){
    //             C_counter++;
    //             j++;
    //         }
    //         if (C_counter>=5)
    //             number_of_5C_repeats++;

    //         //if(rna->GetNucleotide(i+1)=='A')
    //         //	number_of_CA_repeats++;
    //     }
    //     if (nuc=='A'){
    //         int A_counter=1;
    //         int k=i;
    //         while(rna->GetNucleotide(k+1)=='A'){
    //             A_counter++;
    //             k++;
    //         }
    //         if (A_counter>=5)
    //             number_of_5A_repeats++;

    //         //if(rna->GetNucleotide(i+1)=='C')
    //         //	number_of_CA_repeats++;
    //     }
    //     if (nuc=='U' || nuc=='T'){
    //         int U_counter=1;
    //         int j=i;
    //         while(rna->GetNucleotide(j+1)=='U' || rna->GetNucleotide(j+1)=='T'){
    //             U_counter++;
    //             j++;
    //         }
    //         if (U_counter>=5)
    //             number_of_5U_repeats++;
    //     }
    //     utr_prob_sum+=probs[i];
    // }


    // for (int i=mutation_start; i<(mutation_start+number_of_nucs-6); i++){
    // 	string short_seq="";
    // 	for (int j=0; j<6; j++){
    // 		char nuc=rna->GetNucleotide(i+j);
    // 		short_seq+=nuc;
    // 	}

    // 	if(short_seq[0]=='C' && short_seq[1]=='A' && short_seq[2]=='C' && short_seq[3]=='A' && short_seq[4]=='C' && short_seq[5]=='A')
    // 		number_of_3CA_repeats++;
    //  }

    //  for (int i=mutation_start; i<(mutation_start+number_of_nucs-3); i++){
    // 	string short_seq="";
    // 	for (int j=0; j<3; j++){
    // 		char nuc=rna->GetNucleotide(i+j);
    // 		short_seq+=nuc;
    // 	}

    // 	if(short_seq[0]=='C' && short_seq[1]=='A' && short_seq[2]=='A')
    // 		number_of_CAA_repeats++;
    //  }
    //  double complexity_value;
    //  string small_seq="";
    //  for (int i=mutation_start; i<(mutation_start+number_of_nucs); i++){
    // 	char nuc=rna->GetNucleotide(i);
    // 	small_seq+=nuc;
    // }

    double segment_prob_sum = totalPairwiseProbability(rna, mutation_start, number_of_nucs);
	
	// objectiveFunction is one of the enum constants defined in OregaObjectiveFunction in GeneticAlgorithm.h
	switch(objectiveFunction) {
		case OREGA_SIMPLE:
			return 1.0 - segment_prob_sum/((double) number_of_nucs);
		case OREGA_COMPLEX: {
			const string segment = rna->GetSequence(mutation_start, number_of_nucs);
			const double complexity_value = complexity_calc(segment,filterAUG,filterCUG,filteroligoA);
			return complexity_constant*complexity_value + 1.0 - segment_prob_sum / ((double)number_of_nucs);
			
		}
		case OREGA_ALL: {
			const string segment = rna->GetSequence(mutation_start, number_of_nucs);
			const double complexity_value = complexity_calc(segment,filterAUG,filterCUG,filteroligoA);
			double segment_bp_fraction = total_length-Basepair_fraction(rna, mutation_start,total_length,threshold);
			return complexity_constant*complexity_value + 1.0 - segment_prob_sum / ((double)number_of_nucs)+segment_bp_fraction;
			
		}
		//cout<<"Complexity=	"<<complexity_value<<"	"<<small_seq<<endl;
		//return (-1*utr_prob_sum - number_of_5C_repeats - number_of_5A_repeats - number_of_5U_repeats- number_of_CA_repeats);
		//return (-1*utr_prob_sum - number_of_5C_repeats - number_of_5A_repeats - number_of_5U_repeats- number_of_3CA_repeats-number_of_CAA_repeats+5*complexity_value);
		//return (-1*utr_prob_sum);
		//return (-1*sum_53_utr_probs);
	}

	// We should have returned a value above. If we got here, there was an unknown value for objectiveFunction.
	throw "Unknown objective function.";
}


double CalculateRMSD(RNA* rnaWT, RNA* rna){
	double RMSD=0.0;//initialize the RMSD
	int NumberOfPossiblePairs=0;//Initialize the counter for number of possible pairs

	//Calculate the RMSD
	for(int i=1;i<=rnaWT->GetSequenceLength();++i){//For all 5'
		for(int j=i+4;j<=rnaWT->GetSequenceLength();++j){//for all 3'
			if(rnaWT->GetPairProbability(i,j)!=0&&rna->GetPairProbability(i,j)!=0){//if the pairing probability is not zero
				RMSD+=pow((rnaWT->GetPairProbability(i,j)-rna->GetPairProbability(i,j)),2);//square the of the probability difference
				NumberOfPossiblePairs++;//count the number of added members
			}
		}
	}
	//Divide the total sum of squared differences by the number of differences summed
	RMSD=RMSD/NumberOfPossiblePairs;
	
	//Take a square root of the RMSD
	return sqrt(RMSD);
}

//Function that mutates the nucleotides            
int MutateNuc(RNA* rna, int MutationStartNuc, int NumberOfNucs, double MutationRate, int MutationSwitch,int limitG, randomnumber *randomnum,int filterAUG,int filterCUG,int filteroligoA){

    //Declare the variable to store the errors
    int errorcode=0;
	int nucleotide;
	//string s1 ="ACCAC";
	//string s2 ="UCC";
	//list for filter
	string s4 = "AUG";
	string s5 = "CUG";
	string polyA;
    
		//creat filter based on filteroligoA number
		if(filteroligoA > 0){
		for (int i = 0; i < filteroligoA; ++i){
			polyA += 'A';
		}
		}

		//Mutated by mean base pair probability
		if (MutationSwitch == 1){
			double total_probab = totalPairwiseProbability(rna, MutationStartNuc, NumberOfNucs);
			double random = randomnum->roll();
			//initiate running probability check
			double running_sum = 0;
			for(int NucToMutate=MutationStartNuc;NucToMutate<=MutationStartNuc+NumberOfNucs;++NucToMutate){	
				//Calculate the position of the nucleotide to be mutated
				running_sum += totalPairwiseProbability(rna, NucToMutate,1)/total_probab;

				if(running_sum >= random){
				    int Design_done = 0;
					while (Design_done ==0){
					//check if there is any filter
					int filter_check = 0;
					//cerr << "\nrandom=" << random << " MutationStartNuc=" << MutationStartNuc << " NumberOfNucs=" << NumberOfNucs << endl;
						string choose_from_nucs="ACGU";
						int random_int;
						//choose mutation option
						if(limitG == 1){
							random_int=randomnum->roll_int(0,2);
							choose_from_nucs="ACU";
						}
						else{
							random_int=randomnum->roll_int(0,3);
							}
						
						//int dice=rand()%4;
						char random_nuc=choose_from_nucs[random_int];
						rna->GetStructure()->nucs[NucToMutate]=random_nuc;
						
						
						nucleotide=nuc2num(random_nuc);
						
						string input_segment = rna->GetSequence(MutationStartNuc, NumberOfNucs);
						//if (input_segment.find(s1,0) == input_segment.npos && input_segment.find(s2,0) == input_segment.npos && input_segment.find(s3,0) == input_segment.npos && input_segment.find(s4,0) == input_segment.npos && input_segment.find(s5,0) == input_segment.npos)Design_done=1;
					//} 
					if ( input_segment.find(polyA,0) != input_segment.npos && filteroligoA > 0)filter_check += 1;
					
					if ( input_segment.find(s4,0) != input_segment.npos && filterAUG > 0)filter_check += 1;
					
					if ( input_segment.find(s5,0) != input_segment.npos && filterCUG > 0)filter_check += 1;
					
					if (filter_check == 0) Design_done = 1;
						if(nucleotide==1||nucleotide==2||nucleotide==3||nucleotide==4) rna->GetStructure()->numseq[NucToMutate]=nucleotide;
						else{
							cerr << "\nNucleotide conversion error: unknown nucleotide is encountered\n";
							return 1;
						}
						break;
				}
				}
				
			}
			
			
			
			 
			
		}
		//Mutated by random
		else if (MutationSwitch == 0){
			for(int NucToMutate=1;NucToMutate<=NumberOfNucs;++NucToMutate){
				
				//Roll the random number generator
				double random = randomnum->roll();

				//If the random number is less then the set Mutation Rate, mutate this nucleotide
				if(random<MutationRate){
					int Design_done = 0;
					while (Design_done ==0){
						int filter_check = 0;
					//cerr << "\nrandom=" << random << " MutationStartNuc=" << MutationStartNuc << " NumberOfNucs=" << NumberOfNucs << endl;
						string choose_from_nucs="ACGU";
						int random_int;
						//choose mutation option
						if(limitG == 1){
							random_int=randomnum->roll_int(0,2);
							choose_from_nucs="ACU";
						}
						else{
							random_int=randomnum->roll_int(0,3);
							}
						
						//int dice=rand()%4;
						char random_nuc=choose_from_nucs[random_int];
						rna->GetStructure()->nucs[NucToMutate]=random_nuc;
						
						
						nucleotide=nuc2num(random_nuc);
						
						string input_segment = rna->GetSequence(MutationStartNuc, NumberOfNucs);
						//if (input_segment.find(s1,0) == input_segment.npos && input_segment.find(s2,0) == input_segment.npos && input_segment.find(s3,0) == input_segment.npos && input_segment.find(s4,0) == input_segment.npos && input_segment.find(s5,0) == input_segment.npos)Design_done=1;
					//}
					if ( input_segment.find(polyA,0) != input_segment.npos && filteroligoA > 0){
						filter_check += 1;
					}
					if ( input_segment.find(s4,0) != input_segment.npos && filterAUG > 0){
						filter_check += 1;
					}
					if ( input_segment.find(s5,0) != input_segment.npos && filterCUG > 0){
						filter_check += 1;
					}
					if (filter_check == 0) Design_done = 1;
					}
						if(nucleotide==1||nucleotide==2||nucleotide==3||nucleotide==4) rna->GetStructure()->numseq[NucToMutate]=nucleotide;
						else{
							cerr << "\nNucleotide conversion error: unknown nucleotide is encountered\n";
							return 1;
						}
						break;

				}
			}
		}
			
	
	
    return 0;
}
    



int nuc2num(char nuc){

	if (nuc=='A'||nuc=='a'){
		return 1;
	}
	else if (nuc=='C'||nuc=='c'){
		return 2;
	}
	else if (nuc=='G'||nuc=='g'){
		return 3;
	}
	else if (nuc=='T'||nuc=='t'||nuc=='U'||nuc=='u'){
		return 4;
	}
	//if the nucleotide is not ACGTU, return 5 as error
	return 5;
}


int Recombine(RNA* rna1, RNA* rna2, RNA* rnaStore1, RNA* rnaStore2, int NumberOfNucs, int MutationStartNuc, double RecombinationRate, randomnumber *randomnum){

	string sequence1;
	string sequence2;

	//vector<bool> recomID(rna1->GetSequenceLength(),0);
	vector<int> recomPosition;
	recomPosition.push_back(1);//Store the beginning of the sequence

	//Look for codons to recombine and store the positions of the start of the recombination
	for(int nuc=1;nuc<=NumberOfNucs; ++nuc){

		double random = randomnum->roll();
		if(random<RecombinationRate){
			recomPosition.push_back((nuc-1)*1+MutationStartNuc);
		}
	}

	recomPosition.push_back(rna1->GetSequenceLength()+1);//Store the end of the sequence


	bool Flip=true;//Keeps track which strand is being processed
	for(int recom=0;recom<recomPosition.size()-1;++recom){
		//cerr << "\nRunning from " << recomPosition[recom] << " to " << recomPosition[recom+1];
		for(int walk=recomPosition[recom]; walk < recomPosition[recom+1]; ++walk){
			//cerr << Flip << " ";
			if(Flip==true){
				sequence1+=rna1->GetNucleotide(walk);
				sequence2+=rna2->GetNucleotide(walk);
				//sequence1+='0';
				//sequence2+='1';
			}
			else{
				sequence1+=rna2->GetNucleotide(walk);
				sequence2+=rna1->GetNucleotide(walk);
				//sequence1+='1';
				//sequence2+='0';
			}
		}
		//cerr << "\n### FLIP=" << Flip;
		if(Flip)Flip=false;
		else Flip=true;
	}

	//cerr << "\nSequence1=" << sequence1 << endl << "Sequence2=" << sequence2 ;
	
	for(int nuc=1;nuc<=rna1->GetSequenceLength(); ++nuc){
		rnaStore1->GetStructure()->nucs[nuc]=sequence1[nuc-1];
		rnaStore1->GetStructure()->numseq[nuc]=nuc2num(sequence1[nuc-1]);
		rnaStore2->GetStructure()->nucs[nuc]=sequence2[nuc-1];
		rnaStore2->GetStructure()->numseq[nuc]=nuc2num(sequence2[nuc-1]);
	}
	return 0;
}

//Function that checks if the recombination position has been repeated
bool CheckForRepeat(vector<vector<int> > RecombinationPair, int position){

	bool Repeat=false;
	for(int n1=0;n1<RecombinationPair.size();++n1){
		for(int n2=0;n2<RecombinationPair[n1].size();++n2){
			if(RecombinationPair[n1][n2]==position){
				Repeat=true;
			}
		}
	}
	return Repeat;
}

bool AllSwapped(vector<bool> &KeepSequences,int NumberOfSequences){

	bool allswapped=true;

	for(int i=0;i<NumberOfSequences;++i){
		if(!KeepSequences[i]) allswapped=false;
	}
	return allswapped;
}


void CopySequence(RNA* rnaCopyTo, RNA* rnaCopyFrom){
	for(int nuc=1;nuc<=rnaCopyFrom->GetSequenceLength(); ++nuc){
		//cerr << rnaCopyTo->GetStructure()->nucs[nuc] << "-" << rnaCopyTo->GetStructure()->numseq[nuc] << endl;
		rnaCopyTo->GetStructure()->nucs[nuc]=rnaCopyFrom->GetStructure()->nucs[nuc];
		rnaCopyTo->GetStructure()->numseq[nuc]=nuc2num(rnaCopyFrom->GetStructure()->nucs[nuc]);
		//cerr << rnaCopyTo->GetStructure()->nucs[nuc] << "-" << rnaCopyTo->GetStructure()->numseq[nuc] << endl;
	}
}

double objective(double,double);

void ObjectiveFunction(vector<vector<double> > &CalculatedRMSDs, int i){
	CalculatedRMSDs[i][2]=objective(CalculatedRMSDs[i][0],CalculatedRMSDs[i][1]);
}

double objective(double RMSD_hi,double RMSD_low){
	return sqrt(pow(RMSD_hi,2)*15+pow(1-RMSD_low,2));
}

//new comment
