
#include "design.h"

//#define DEBUG
#undef DEBUG

//How close a slit on an exterior loop needs to be during decomposition
const double EXTERIORCLOSENESS = 0.10;

//How close a slit on an multibranch loop needs to be during decomposition
const double MULTICLOSENESS = 0.5;

//Minimum fragment length worth keeping
const int MINIMUMLENGTH = 6;


//Constructor:
design::design(const char filename[], const char* const alphabet):RNA(filename,FILE_CT,alphabet,true,true) {
	//For now, break any pseudoknots in the structure.
	BreakPseudoknot(false);
	MaxRedesign=30;
	MaxMutate=4;
	MaxLeafRedesign=3;
	defectweighted=true;
	numbering = 1;


}



//Design sequence:
//	pernucdefect:  IN: the required maximum per-nucleotide ensemble defect
//                 OUT: the per-nucleotide ensemble defect of the resulting (designed) structure. 
int design::design_sequence(double &pernucdefect, const bool random, const bool bias, const char prob_info[], const int maxdepth, bool heuristic, int MaxRedesignC, int MaxMutateC, int MaxLeafRedesignC, long randomSeed, bool allowisolated) {
	
	
	//This code handles a selection bias for nucleotides and pairs.  The default (bias=false) is a flat prob across the possibilities.
	if(bias){
		//move the following section to a separate function later, but for now:
		ifstream in;
		in.open(prob_info);
		string line;
		int stage=-1;
		int index;
		int count=1;
		char base, b0, b1;
		int b_ndx_0, b_ndx_1;
		double probability;

		//Now also hardwire the single-stranded nucs, give everything (but X) the same probability for now:
		//singles.resize(GetStructure()->GetThermodynamicDataTable()->alphabet.size());
		singlebias.resize(GetStructure()->GetThermodynamicDataTable()->alphabet.size());
		


		//hard wire some settings for now:

		//Decide how many pairs are allowed in the alphabet:
		//int totalpairs=0;

		//for (int i=0;i < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++i) {
		//	for (int j=0;j < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++j) {

		//		if (GetStructure()->GetThermodynamicDataTable()->pairing[i][j]) {
		//			++totalpairs;
		//		}

		//	}
		//}

		//pairs.resize(GetStructure()->GetThermodynamicDataTable()->alphabet.size() - GetStructure()->GetThermodynamicDataTable()->linker.size());
		pairbias.resize(GetStructure()->GetThermodynamicDataTable()->alphabet.size());

		for (int i = 0; i < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++i) {
			singlebias[i] = 0.0;//set x = 0

			pairbias[i].resize(GetStructure()->GetThermodynamicDataTable()->alphabet.size());
			for (int j = 0; j < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++j) {
				pairbias[i][j] = 0.0;
			}
			//totalpairs=0;//now keep track of how many pairs have been stored:
		}

		while(!in.eof()){
			in>>line;
			if(in.eof())
				break;

	        	//skip empty lines
			if (line=="")
				continue;

			//read line is "Singles" then stage=0,     
			if (line == "Singles")
				stage = 0;
			//read line is "Pairs" then stage=1,  
			else if (line == "Pairs")
				stage = 1;
			else if (line == "Bias-in-Leaf-Refinement") stage = 2;
			else{
				switch(stage){
					case 0:{
						in>>probability;
						base=line[0];
						index=GetStructure()->GetThermodynamicDataTable()->basetonum(base);
						//singles[count]=index;
						singlebias[index]=probability;
						count++;
					}
					break;

					case 1:{
						in>>b1>>probability;
						b0=line[0];
						pairbias[GetStructure()->GetThermodynamicDataTable()->basetonum(b0)][GetStructure()->GetThermodynamicDataTable()->basetonum(b1)] = probability/2;
						pairbias[GetStructure()->GetThermodynamicDataTable()->basetonum(b1)][GetStructure()->GetThermodynamicDataTable()->basetonum(b0)] = probability/2;
						/*b_ndx_0=GetStructure()->GetThermodynamicDataTable()->basetonum(b0);
						b_ndx_1=GetStructure()->GetThermodynamicDataTable()->basetonum(b1);
						switch(b_ndx_0){
							//for A-U basepair line in input file
							case 1:{
								pairs[0][0]=b_ndx_0; pairs[0][1]=b_ndx_1; pairbias[0]=probability;
								pairs[4][0]=b_ndx_1; pairs[4][1]=b_ndx_0; pairbias[4]=probability;

							}
							break;
							//for C-G basepair line in input file
							case 2:{
								pairs[1][0]=b_ndx_0; pairs[1][1]=b_ndx_1; pairbias[1]=probability;
								pairs[2][0]=b_ndx_1; pairs[2][1]=b_ndx_0; pairbias[2]=probability;
							}
							break;
							//for G-C or G-U basepair lines in input file
							case 3:{
								if(b_ndx_1==2){
								pairs[2][0]=b_ndx_0; pairs[2][1]=b_ndx_1; pairbias[2]=probability;

								pairs[1][0]=b_ndx_1; pairs[1][1]=b_ndx_0; pairbias[1]=probability;
								}
								else{
								pairs[3][0]=b_ndx_0; pairs[3][1]=b_ndx_1; pairbias[3]=probability;
								pairs[5][0]=b_ndx_1; pairs[5][1]=b_ndx_0; pairbias[5]=probability;
								}
							}
							break;
							//for U-A or U-G basepair lines in input file
							case 4:{
								if(b_ndx_1==1){
								pairs[4][0]=b_ndx_0; pairs[4][1]=b_ndx_1; pairbias[4]=probability;
								pairs[0][0]=b_ndx_1; pairs[0][1]=b_ndx_0; pairbias[0]=probability;
								}
								else{
								pairs[5][0]=b_ndx_0; pairs[5][1]=b_ndx_1; pairbias[5]=probability;
								pairs[3][0]=b_ndx_1; pairs[3][1]=b_ndx_0; pairbias[3]=probability;
								}
							}
							break;
						}*/
						//totalpairs+=2;				
					}
					break;

					case 3: {
						if (line == "Yes") bias_in_leaf_refinement = true;
						else bias_in_leaf_refinement = false;
					}
					break;


				}
			}
		}
		

		//print input probability values to screen
		cout << "\n\nBiases were specified from a bias file:" << endl;
		cout<<"Singles"<<endl;
		for(int i=0; i< GetStructure()->GetThermodynamicDataTable()->alphabet.size(); i++)
			if (singlebias[i]>0.0) cout<<GetStructure()->GetThermodynamicDataTable()->numtobase(i)<<"	"<<singlebias[i]<<endl;
		cout<<"Pairs"<<endl;
		//cout<<GetStructure()->GetThermodynamicDataTable()->numtobase(pairs[0][0])<<" "<<GetStructure()->GetThermodynamicDataTable()->numtobase(pairs[0][1])<<"	"<<pairbias[0]<<endl;
		for (int i = 0; i < GetStructure()->GetThermodynamicDataTable()->alphabet.size() ; ++i) {
			for (int j = 0; j < GetStructure()->GetThermodynamicDataTable()->alphabet.size() ; ++j) {
				if (pairbias[i][j]>0.0) cout << GetStructure()->GetThermodynamicDataTable()->numtobase(i) << "-" << GetStructure()->GetThermodynamicDataTable()->numtobase(j) << " =" << pairbias[i][j] << endl;
			}
		}
		if (bias_in_leaf_refinement) {
			cout << "Leaf refinement will use the specified bias." << endl;
		}
		else {
			cout << "Leaf refinement will not use the specified bias and will only use A, C, G, or U/T nucleotides and A-U/T or G-C pairs." << endl;
		}
		cout<<"\n\nResults"<<endl;
	}
	else bias_in_leaf_refinement = false; //no bias used, therefore none in leaf refinement

	int **tree;
	int i,j;
	MaxRedesign=MaxRedesignC;
	MaxMutate=MaxMutateC;
	MaxLeafRedesign=MaxLeafRedesignC;
	//Load the thermodynamic parameters:
	if (!VerifyThermodynamic()) return 5;

	//First decompose the structure
	//This decomposition breaks down the structure along the lines of multibranch loops
	//(This could be changed...)

	//Store the decomposition in **tree
	tree = new int *[maxdepth];
	for (i=0;i<maxdepth;i++) {
		tree[i] = new int [GetSequenceLength()+1];
	}

	//////////////////////
	/////Look for uninitialized values
	for (i=0;i<maxdepth;i++) {
		for (j=1;j<=GetSequenceLength();j++) {
			tree[i][j] = -99;
		}
	}
	decompose(1,GetSequenceLength(),0,maxdepth,tree);
#ifdef DEBUG
	//cerr << "DONE: structure decomposition\n";
	///////////////////////////
	//send debug info to disk:
	ofstream out;
	out.open("tree.out");

	for (i=1;i<=GetSequenceLength();i++) {
		out << i << "\t";
	}
	out<<"\n";

	for (i=1;i<=GetSequenceLength();i++) {
		if (GetPair(i)==0) out << "-";
		else if (GetPair(i)<i) out << "<";
		else out << ">";
		out << "\t";
	}
	out << "\n";

	for (i=0;i<maxdepth;i++) {
		for (j=1;j<=GetSequenceLength();j++) {
			
			out << tree[i][j]<<"\t";

		}
		out << "\n";
	}

	out.close();
#endif
	//cerr << "DONE: writing 'tree.out' file\n";

	//Now design elements  (and set pernucdefect to the resulting ensemble defect)
	if (!heuristic) 
		pernucdefect=SelectSequence(tree,random,bias,maxdepth,pernucdefect, randomSeed, allowisolated);
	else 
		pernucdefect=SelectSequenceHeuristic(tree,random,bias,maxdepth,pernucdefect, randomSeed,allowisolated);

	//cerr << "OUTSIDE" << endl;

	//Now cleanup:
	for (i=0;i<maxdepth;i++) {

		delete[] tree[i];

	}
	delete[] tree;

	//cerr << "DONE WITH OUTSIDE" << endl;
	return 0;

}



//Recursively perform a binary decomposition along multibranch loops
//Note that this only works for pseudoknot-free structures
void design::decompose(int nucstart,int nucend, int currentdepth, int maxdepth, int **tree, int missingstart, int missingend) {
	int i,j,branches,k;
	vector<int> stack;
	vector<int> shortstack;
	int beststart,bestend,currentlength;

	
	beststart = nucstart;
	bestend = nucend;

	
	int newmissingstart,newmissingend;

	//Find the optimal decomposition point for this level of the structure:

	


	currentlength = nucend - nucstart + 1;
	//check if this encomposses a missing fragment and correct:
	if (missingstart > nucstart) {

		currentlength-= (missingend - missingstart +1);

	}


	//Only decompose if sequence is long enough to have two fragments
	if (currentlength>2*MINIMUMLENGTH+1) {

		//Run between nucstart and nucend to find best division that splits sequence in two

		i=nucstart;

		//First probe the exterior loop
		while (i<=nucend) {
		
		
			//check if i is a good place to divide the sequence in two
			if (closeenoughtocut(i,nucend,nucstart,nucend,missingstart,missingend,EXTERIORCLOSENESS)) {
				//This is close enough to the middle

				//Mark the tree
				marktree(nucstart,i-1,nucstart,nucend,missingstart,missingend,currentdepth,tree);
				/*for (j=nucstart;j<i;j++) {
				  if (!(j>=missingstart&&j<=missingend)) {
				  tree[currentdepth][j] = numbering;
				  }

				  }
				  for (j=i;j<=nucend;j++) {
				  tree[currentdepth][j] = numbering + pow(two,currentdepth);
				  }*/

				if (currentdepth<maxdepth-1) {

					if (missingstart==0) {
						newmissingstart=0;
						newmissingend=0;
					}
					else if (missingstart<i-1) {
						newmissingstart=missingstart;
						newmissingend=missingend;

					}
					else {
						newmissingstart=0;
						newmissingend=0;

					}
					decompose(nucstart,i-1,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

					if (missingstart==0) {
						newmissingstart=0;
						newmissingend=0;
					}
					else if (missingstart>i-1) {
						newmissingstart=missingstart;
						newmissingend=missingend;

					}
					else {
						newmissingstart=0;
						newmissingend=0;
					}
					decompose(i,nucend,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

				}

				return;


			}//end if "close enough"
			else {
				//keep track of the best decomposition found
				bestdecomposition(nucstart,nucend,i,nucend,&beststart,&bestend,missingstart,missingend);

			}//end else

		
			if (GetPair(i)==0) {

				++i;
			

			}
			else {
				//add to the helix list
				stack.push_back(i);
				i = GetPair(i) + 1;

			}

			//check if i is entering part of the sequence that is missing
			if (i==missingstart) {
				i = missingend+1;
			}
		

		}//end while i<nucend


		//proceed while there are still structure fragments on the stack
		while (stack.size()>0) {
		
			//pull a fragment from the stack, i will be a nuc in a basepair
			i = stack.back();
			stack.erase (stack.begin() + stack.size() - 1);
		
			//Scan sequence structure fragment for multibranch loops
			while (GetPair(i+1)+1==GetPair(i)) {

				++i;
			}
		
			//interuption in pairing reached:

			//check if this is a multibranch loop:
			j = i+1;
			if (j==missingstart) j = missingend+1;
			branches = 0;
			while (j<GetPair(i)) {
				if (GetPair(j)>0) {
					//new branch
					++branches;

					//put this branch on the stack
					stack.push_back(j);
					j = GetPair(j);

				}
				

				++j;
				if (j==missingstart) {
					j = missingend+1;
					//count the missing fragment as a branch
					++branches;

				}
				

			}

			if (branches>1) {
				//This is a multibranch loop
				shortstack.clear();

				//Go through again looking for possible splits and adding fragments to stack

				j = i+1;
				if (j==missingstart) j = missingend+1;
				branches = 0;
				//check if i and its pair is a good place to split:
				//bestdecomposition(nucstart,nucend,i,GetPair(i),&beststart,&bestend,missingstart,missingend);



				while (j<GetPair(i)) {
					if (GetPair(j)>0) {
						//branch
					
						shortstack.push_back(j);

						//check all possible splits to see if they are good:
						for (k=0;k<shortstack.size();k++) {
							
							
							bestdecomposition(nucstart,nucend,shortstack[k],GetPair(j),&beststart,&bestend,missingstart,missingend);
							
						}

						j = GetPair(j);

					}//end if GetPair(j)>0
					else {

						//check all possible splits to see if they are good:
						for (k=0;k<shortstack.size();k++) {
							
							
							bestdecomposition(nucstart,nucend,shortstack[k],j,&beststart,&bestend,missingstart,missingend);
							
						}


					}
					

					++j;
					if (j==missingstart) {
						j = missingend+1;

						//put this edge on the shortstack:
						shortstack.push_back(j);
					}
					

				}//end while j<GetPair(i)


				//check to see if the decomposition is good enough
				if (closeenoughtocut(beststart,bestend,nucstart,nucend,missingstart,missingend,MULTICLOSENESS)) {
					//This is close enough to the middle

					

					//Mark the tree
					marktree(beststart,bestend,nucstart,nucend,missingstart,missingend,currentdepth,tree);

					if (currentdepth<maxdepth-1) {

						if (missingstart==0) {
							newmissingstart=0;
							newmissingend=0;
						}
						else {
							if (missingend==beststart-1) {

								newmissingend=0;
								newmissingstart=0;

							}
							else if (missingstart==bestend+1) {
								newmissingend=0;
								newmissingstart=0;
							}
							else {
								newmissingstart=missingstart;
								newmissingend=missingend;
							}

						}
						
						decompose(beststart,bestend,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

						if (missingend==beststart-1) {

							newmissingend=bestend;
							newmissingstart=missingstart;

						}
						else if (missingstart==bestend+1) {
							newmissingend=missingend;
							newmissingstart=beststart;
						}
						else {
							newmissingend = bestend;
							newmissingstart = beststart;
						}


						decompose(nucstart,nucend,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

					}

					return;


				}//end if 

			}//end if (branches >1)


		
		

		}//end while (stack.size()>0)

	}//end if (currentlength>2*MINIMUMLENGTH+1)

	//Decide if there is a decomposition that can be made
	if (beststart==nucstart&&bestend==nucend) {
		//There is no decomposition
		for (i=currentdepth; i<maxdepth; i++) {
			for (j=nucstart;j<=nucend;j++) {
				if (j==missingstart) j=missingend;
				else tree[i][j]=0;
			}
		}

	}
	else {
		//There is a decomposition possible


		//Mark the tree
		marktree(beststart,bestend,nucstart,nucend,missingstart,missingend,currentdepth,tree);

		if (currentdepth<maxdepth-1) {

			if (missingstart==0) {
				newmissingstart=0;
				newmissingend=0;
			}
			else {
				if (missingend==beststart-1) {

					newmissingend=0;
					newmissingstart=0;

				}
				else if (missingstart==bestend+1) {
					newmissingend=0;
					newmissingstart=0;
				}
				else {
					newmissingstart=missingstart;
					newmissingend=missingend;
				}

			}
						
			decompose(beststart,bestend,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

			if (missingend==beststart-1) {

				newmissingend=bestend;
				newmissingstart=missingstart;

			}
			else if (missingstart==bestend+1) {
				newmissingend=missingend;
				newmissingstart=beststart;
			}
			else {
				newmissingend = bestend;
				newmissingstart = beststart;
			}

			if (newmissingend==nucend) {

				decompose(nucstart,newmissingstart-1,currentdepth+1,maxdepth,tree,0,0);

			}

			else {

				decompose(nucstart,nucend,currentdepth+1,maxdepth,tree,newmissingstart,newmissingend);

			}
		}

	}//end else


}


bool design::closeenoughtocut(int i,int j,int nucstart,int nucend,int missingstart,int missingend,double CLOSENESS) {
	int targetlength,closeenoughlow,closeenoughhigh;
	int currentlength;

	//First make sure the cut is allowed.  Each Fragment can have at most one hole.  In other words, missingstart and missingend must be in the inside fragment.
	if (j!=nucend) {
		//if currentend==nucend, there is no missing fragment, so it is OK

		if (missingstart!=0) {

			if (missingstart<i||missingstart>j) {
				
				//This would not be a valid cut because it would make a fragment with two holes.
				return false;

			}

		}
	}

	//targetlength is the desired fragment length
	targetlength = nucend - nucstart+1;
	//correct for missing fragment in middle
	targetlength-= (missingend-missingstart+1);
	//now divide by two
	targetlength = targetlength/2;


	closeenoughlow = targetlength - (int)(((double)targetlength)*CLOSENESS);
	closeenoughhigh = targetlength + (int)(((double)targetlength)*CLOSENESS);

	currentlength = j-i;
	//correct for any missing pieces
	if (missingstart>i&&missingend<j) {
		currentlength-=(missingend-missingstart+1);

	}

	if (currentlength>closeenoughlow&&currentlength<closeenoughhigh) return true;
	else return false;


}

void design::marktree(int beststart,int bestend,int nucstart,int nucend,int missingstart,int missingend, int currentdepth,int **tree) {
	int j;

	for (j=nucstart;j<beststart;j++) {
		if (!(j>=missingstart&&j<=missingend)) {
			tree[currentdepth][j] = numbering;
		}
	}
	for (j=bestend+1;j<=nucend;j++) {
		if (!(j>=missingstart&&j<=missingend)) {
			tree[currentdepth][j] = numbering;
		}
	}
	numbering++;
	for (j=beststart;j<=bestend;j++) {
		if (!(j>=missingstart&&j<=missingend)) {
			//int twopoew = pow(two,currentdepth); 
			tree[currentdepth][j] = numbering;
			
		}
	}
	numbering++;

	

}

void design::bestdecomposition(int nucstart, int nucend, int currentstart, int currentend, int *beststart, int *bestend, int missingstart, int missingend) {
	
	int desiredlength;
	int bestsofarlength;
	int currentlength;

	

	//First make sure the cut is allowed.  Each Fragment can have at most one hole.  In other words, missingstart and missingend must be in the inside fragment.
	if (currentend!=nucend) {
		//if currentend==nucend, there is no missing fragment, so it is OK

		if (missingstart!=0) {

			if (missingend+1<currentstart||missingstart>currentend+1) {
				
				//This would not be a valid cut because it would make a fragment with two holes.
				return;

			}

		}
	}

	currentlength = currentend - currentstart + 1;
	//check if this encomposses a missing fragment and correct:
	if (missingstart > currentstart && missingstart < currentend) {

		currentlength -= (missingend - missingstart +1);

	}


	//Make sure the current length is longer than a minimum length
	if (currentlength<MINIMUMLENGTH||(nucend - nucstart + 1  - missingend + missingstart - currentlength)<MINIMUMLENGTH) return;


	desiredlength = (nucend - nucstart - missingend + missingstart)/2;

	bestsofarlength = *bestend - *beststart + 1;
	//check if this encomposses a missing fragment and correct:
	if (missingstart > *beststart && missingstart < *bestend) {

		bestsofarlength-= (missingend - missingstart +1);

	}

	
	
	if (abs(currentlength - desiredlength) < abs(bestsofarlength - desiredlength)) {

		//so this fragment is closer to the desired length
		*beststart = currentstart;
		*bestend = currentend;
	}
}

// Returns ensemble defect (per-nucleotide) if a sequence was designed.
double design::SelectSequence(int **tree, bool random, bool bias, int depth, const double bestpernucdefect, long seed, bool allowisolated) {
	//////////////////////HEADER/////////////////
	int i,j,k,l;
	RNA *fragment;
	//Thermodynamics *thermo;
	vector<int> stackstart, stackend, stackmissingstart, stackmissingend, stackfragmentdepth; //use a vector as a stack of fragments
	int start,end,missingstart,missingend,fragmentdepth;
	bool providesequence;
	//bool leafOpt;
	randomnumber dice;
	char *sequence;
	//char *sequenceBEST;
	char sequenceTEMP,sequenceTEMP2;
	double defect,defectOLD;
	double pernucdefect,temppernucdefect;
	int frag1,frag2,ident1;
	double roll;

#ifdef DEBUG
	ofstream out;
	out.open("stack-pulls.out");
#endif

	//initialize the random number generator
	dice.seed(seed);

	//initilaize the sequence to be large neough for any situation
	sequence = new char [GetSequenceLength()+9];

	//Read the thermodynamics from disk
	//thermo = new Thermodynamics(*this);
	//if (!thermo->VerifyThermodynamic()) ;
	///////////////////////END HEADER///////////////////////

	//start building the sequence from the bottom up, i.e. small fragments first
	//Do this by putting fragments on the stack, smallest fragments last so that they are handled first
	for (i=-1;i<depth;++i) {
		//For the current level, find each fragment, select a sequence if neccessary, and test that the correct structure forms
		//Identify fragments from left to right
		if (i>=0) FindFragments(tree,i,1,GetSequenceLength(),0,0,&stackstart, &stackend, &stackmissingstart, &stackmissingend, &stackfragmentdepth);
		else {
			//Place the whole sequence on the stack:
			stackstart.push_back(1);
			stackend.push_back(GetSequenceLength());
			stackmissingstart.push_back(0);
			stackmissingend.push_back(0);
			stackfragmentdepth.push_back(-1);
		}
	}//end for i
	//cerr << "depth=" << depth << " stackfragmentdepth.size()=" << stackfragmentdepth.size() << " stackfragmentdepth.back()=" << stackfragmentdepth.back()<< endl;
	double **EDtree;//Keep track of best ensemble defect for each tree decomposition level
	EDtree=new double *[depth+2];
	for(i=0;i<depth+1;i++){
		EDtree[i]=new double [GetSequenceLength()+1];
		for (j=0;j<=GetSequenceLength();++j) EDtree[i][j]=2.0; //initilaize to large value
	}
	EDtree=EDtree+1;

	char **SEQtree;//Keep track of the sequence with best ensemble defect for each tree decomposition level
	SEQtree=new char *[depth+2];
	for(i=0;i<depth+1;i++){
		SEQtree[i]=new char [GetSequenceLength()+1];
	}
	SEQtree=SEQtree+1;
	
	int **ReoptNum;//Keep track of the number of reoptimizations for each tree decomposition level
	ReoptNum=new int *[depth+2];
	for(i=0;i<depth+1;i++){
		ReoptNum[i]=new int [GetSequenceLength()+1];
		for (j=0;j<=GetSequenceLength();j++) ReoptNum[i][j] = 0;
	}
	ReoptNum=ReoptNum+1;


	vector<vector<string> > Helices (11,vector<string>());
	vector<vector<string> > Loops (11,vector<string>());
	//Read the helices and loops from the database if 'random==false'
	if(!random){
	
		//Read datapath.
		string hf;
		string lf;
		if(GetBackboneType()){//if GetBackbone returns true, read RNA sequence files
			hf="design.RNA.Helices.dat";//Holds the file name with helix sequences.
			lf="design.RNA.Loops.dat";//Holds the file name with loop sequences
		}
		else{//else, if GetBackbone returns false, read DNA files
			hf="design.DNA.Helices.dat";//Holds the file name with helix sequences.
			lf="design.DNA.Loops.dat";//Holds the file name with loop sequences
		}		
		string dp(getDataPath()); //Holds the path to data_tables.
		
		//Append the filename to the path 'dp' and store in 'fn'.
		hf=dp+"/"+hf;
		lf=dp+"/"+lf;
		string lineoftext;//String that temporarily holds strings of data from the read files.
		//Open file 
		ifstream readHfile(hf.c_str());
		ifstream readLfile(lf.c_str());
		//			cerr << "lf=" << lf.c_str() << endl;

		while(!readHfile.eof()){
			readHfile >> lineoftext;
			if(lineoftext.empty()) continue;
			Helices.at((lineoftext.size()-3)/2).push_back(lineoftext);
		}

		while(!readLfile.eof()){
			readLfile >> lineoftext;
			//				cerr << "lineoftext=" << lineoftext << endl;
			if(lineoftext.empty()) continue;
			Loops.at(lineoftext.size()).push_back(lineoftext);
		}
		readHfile.close();
		readLfile.close();
	}

	//cerr << "stackfragmentdepth.back()=" << stackfragmentdepth.back();
	//cerr << " stackfragmentdepth.back()=" << stackfragmentdepth.back();
	//Pull each fragment from the stack
	//pull a fragment from the stack
	while (stackstart.size()>0) {

	  //cerr << "STACK=" << stackstart.size() << endl;
		//		cerr << "TEST0" << endl;
		//cerr << "-----------------CHANGING STACK SIZE------------------\n";

		start = stackstart.back();
		stackstart.erase (stackstart.begin() + stackstart.size() - 1); 

		end = stackend.back();
		stackend.erase (stackend.begin() + stackend.size() - 1);

		missingstart = stackmissingstart.back();
		stackmissingstart.erase (stackmissingstart.begin() + stackmissingstart.size() - 1);

		missingend = stackmissingend.back();
		stackmissingend.erase (stackmissingend.begin() + stackmissingend.size() - 1);

		fragmentdepth = stackfragmentdepth.back();
		stackfragmentdepth.erase (stackfragmentdepth.begin() + stackfragmentdepth.size() - 1);
		//cerr << "fragmentdepth=" << fragmentdepth << " =" << stackfragmentdepth.back() << " =" << stackfragmentdepth.back()  << endl;
;		

#ifdef DEBUG
		out << start << "\t" << end << "\t"<<missingstart << "\t" << missingend << "\t" << fragmentdepth << "\n"<<flush;
#endif

		//cerr << " start=" << start << " end=" << end << " missingstart=" << missingstart << " missingend=" << missingend << endl;
		

		// cerr << " check1 fragmentdepth=" << fragmentdepth << " start=" << start << " sequenceLength=" << GetSequenceLength()<< endl;

		//Increment the counter for re-optimization
		++ReoptNum[fragmentdepth][start];
      
		//cerr << "STACK=" << stackstart.size() << " fragmentdepth=" << fragmentdepth << " start=" << start << " Re-Opt#=" << ReoptNum[fragmentdepth][start] << endl;


		//Provide a sequence if none has been provided for this fragment (because we are at the leaves of the tree)
		providesequence = false;
		if (fragmentdepth==depth-1) providesequence=true;
		if (fragmentdepth<depth-1) {
			if (tree[fragmentdepth+1][start]==0) providesequence=true;
		}
		
		//Calculate revised per nuc defect here:
		pernucdefect = bestpernucdefect;
		if (!providesequence) {
			//Find the subfragments, and determine their pernucdefect

			temppernucdefect=0;
			frag1=1;
			frag2=0;
			//temppernucdefect=EDtree[fragmentdepth+1][start];
			ident1=tree[fragmentdepth+1][start];
			for (i=start+1;i<=end;++i) {
				if (!(i>=missingstart&&i<=missingend)) {

					if (tree[fragmentdepth+1][i]==ident1) frag1++;
					else {
						frag2++;
						if (frag2==1) {
							temppernucdefect=EDtree[fragmentdepth+1][i];
						}

					}
					
					
				}
				
			}
			pernucdefect=(max(temppernucdefect,bestpernucdefect)*((double) frag2)+max(bestpernucdefect,EDtree[fragmentdepth+1][start])*((double) frag1))/(((double) frag1) + ((double) frag2));
			
				
#ifdef DEBUG
			if (pernucdefect>bestpernucdefect) {
				cerr <<"revised pernucdefect of " <<pernucdefect << "\n";
			}
#endif

			
		}

		//cerr << "TEST1" << endl;
		if (providesequence) {

			FillSequence(start,end,missingstart,missingend,random,bias,&dice,Helices,Loops);

		}//end if provide sequence
		//		cerr << "CompletedFillSequence" << endl;
		//Now check this fragment to see if it folds correctly
		//Populate sequence with the information
		k=0;

		for (j=start;j<=end;j++) {
			//iterate over each nucleotide
			if (j==missingstart) {
				j=missingend;
				for (l=0;l<6;l++) {
					sequence[k] = 'X';
					k++;
				}
			}
			else {
				sequence[k] = GetStructure()->nucs[j];
				k++;
			}//end else, j!=missingstart
		}//end for j
		sequence[k]='\0';//add the null terminator
	
		//cerr << "SIZEOF=" << sizeof(sequence) << " STRLEN=" << strlen(sequence) << " GetSequenceLength=" << GetSequenceLength()<< endl;

		//for(int f=0;f<=k;f++){
		//	cerr << "sequence[" << f << "]=" << sequence[f] << endl;
		//}

		//initialize the fragment. Copy thermodynamic parameters from *this* (a design:RNA instance that has already loaded the parameter tables).
		fragment = new RNA(sequence,SEQUENCE_STRING,this);

		//copy the thermodynamic parameters into fragment
		//fragment->CopyThermodynamic(thermo); -- RMW: no longer necessary. tables are copied in the constructor.

		//Perform the partition function calculation
		fragment->PartitionFunction("",-10.0,false,true,allowisolated);
		
		//Debug1(start,end,missingstart,missingend,sequence,fragment);

		//Now tabulate the ensemble defect
		defect = 0.0;
		
		vector<double> def (fragment->GetSequenceLength()+1,0);//stores ensemble defect for each nucleotide in a fragment
		GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'
		//If XXXXXX were used, subtract them from the total length of the fragment
		if(missingend!=0) defect= defect/((double)fragment->GetSequenceLength()-6);
		else defect= defect/((double)fragment->GetSequenceLength());

		//QUESTION: implement this function instead of the code in "testing new code" section
		//LeafOptimize(pernucdefect, &dice, &defect, &def, fragment, start, end, missingstart, missingend, sequence, thermo)

		////////////////////TESTING NEW CODE/////////////////////
		//bool redesign=true;
		double length;//stores the length of the fragment exclusing the 'xxxxxx'
		if(missingend!=0) length = (double)fragment->GetSequenceLength()-6;//set length if there is 'xxxxxx'
		else length = (double)fragment->GetSequenceLength();//set length if there is NO 'xxxxxx'
		int CountUnfavorable=0;//keeps track of unfavorable mutations
		vector <vector<bool> > Unfavorable (fragment->GetSequenceLength()+1,vector<bool> (GetStructure()->GetThermodynamicDataTable()->alphabet.size(),false));//vector keeps track of unfavorablel mutations
		//char sequenceTemp;//temporary stores the pre-mutated nucleotide in case of rejected mutation

		if(providesequence&&random){
			while(defect>pernucdefect&&CountUnfavorable<MaxMutate*length){



				//Find the most probable nucleotide to mutate 
				//choose the nucleotide
				i=0;//i is the position in the sequence string (0 indexed)
				double roll = dice.roll();
				double current = def[1]/(length*defect);//def, the per nuc defect array, is 1 indexed
			
				while (roll > current) {
					
					++i;
					//consider each nucleotide that is notthe XXXXXX linker
					if (sequence[i]!='X') {
						current+=def[i+1]/(length*defect);
					}
					else i+=5;
					if (i>= fragment->GetSequenceLength()) {
						cerr << "i ran out of bounds in leaf defect-weighted redesign\n";
					}
				}		
				//END: Find the most probable nucleotide to mutate

				char newnuc1, newnuc2;
			
				//now the i-th nucleotide should be mutated by 'change'
				if (bias_in_leaf_refinement) {
					if (GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)) > 0) {
						int newnuc1int, newnuc2int;
						randompair(&newnuc1int, &newnuc2int, &dice);
						newnuc1 = tonuc(newnuc1int);
						newnuc2 = tonuc(newnuc2int);

					}
					else {
						newnuc1 = tonuc(randomnuc(&dice));
					}

				}
				else {
					int change = ((int)(dice.roll() * 3.0)) + 1;
					newnuc1 = tonuc((toint(sequence[i]) + change - 1) % 4 + 1);
					if (GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)) > 0) {
						newnuc2 = tonuc(5 - toint(newnuc1));
					}
				}
				
				

				//If the mutation was already rejected before, return to the begining of 'while' loop and add 1 to CountUnfavorable
				if(Unfavorable[i][toint(newnuc1)]) CountUnfavorable++;
				else{//If the mutation was not rejected before
					sequenceTEMP=sequence[i];//save the i-th nuc so that it can be restored later
					defectOLD=defect;//save the ensemble defect of the fragment before the mutation
					vector<double> defTEMP (fragment->GetSequenceLength()+1,0);//save the ensemble defect per nucleotide before the fragment mutation

					sequence[i]= newnuc1;//mutate the i-th nucleotide
					
					//First check if the position is paired
					if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
						//if the position is paired mutate the paired nucleotide.
						sequenceTEMP2 = sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1];
						sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1] = newnuc2;
						
					}
					///////////////////////////
					//Evaluate the Ensemble Defect for the new fragment with mutated nucleotide
					///////////////////////////
					delete fragment;//delete the fragment with the old sequence
					//initialize the new fragment. Copy thermodynamic parameters from *this* (a design:RNA instance that has already loaded the parameter tables).
					fragment = new RNA(sequence,SEQUENCE_STRING,this);					
					fragment->PartitionFunction("",-10.0,false,true,allowisolated);//Perform the partition function calculation
					defTEMP = def;//temporarily store the 'def' in defTemp vector
					for (j=0;j<=fragment->GetSequenceLength();j++) def[j]=0;//restore 'def' vector to all 0s
					defect=0;//restore 'defect' to 0
					GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'
					defect= defect/length;
					//END: Evaluate Ensemble Defect for the new fragment with mutated nucleotide
					///////////////////////////
					
					//If the new defect is below the old defect, store the fragment
					if(defect<defectOLD){
						CountUnfavorable=0;
						defectOLD=defect;
						StoreMutation(start,end,missingstart,missingend,sequence);
						for(int p=0;p<=fragment->GetSequenceLength();p++){
							for(int q=0;q<=4;q++) Unfavorable[p][q]=false;
						}
					}
					//If the new defect is above the old defect
					else {
						sequence[i]=sequenceTEMP;//reverse the mutation
						defect = defectOLD;//restore the defect
						//also restore the pair if needed...
						if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
							sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=sequenceTEMP2;
						}												
						def=defTEMP;//restore the 'def' - vector that stores defects per nucleotide
						Unfavorable[i][toint(newnuc1)]=true;//store the rejected mutation
						CountUnfavorable++;//add 1 to CountUnfavorable to keep track of total rejected mutations
					}//END: if the defect is above the old defect
				}//END: if the mutation was not rejected before
			}//END: while(defect<pernucdefect && CountUnfavorable < MaxMutate*length)
		}
		////////////////END: TESTING NEW CODE/////////////////
		/*
		if (providesequence&&defect>pernucdefect&&random) {
			vector<double> defTEMP (fragment->GetSequenceLength()+1,0);//For temporary storage of 'def'
			//vector<double> defTOTAL (fragment->GetSequenceLength()+1,0);
			//vector<int> Mutated(fragment->GetSequenceLength()+1,0);//keeps track which nucleotide was mutated and how many times during the leaf optimization
			//cerr << "fragment_length=" << fragment->GetSequenceLength() <<" ";

			bool redesign=true;
			int changes = 0;
			defectOLD = defect;
			double length;

			if(missingend!=0) length = (double)fragment->GetSequenceLength()-6;
			else length = (double)fragment->GetSequenceLength();
#ifdef DEBUG
			int KeepMut=0;//Counter for the number of times the mutation was accepted, but needed more mutations
			int RejectMut=0;//Counter for the number of times the mutation was rejected
#endif
			vector <vector<int> > Unfavorable (length+1,vector<int> (5,0));//vector keeps track of unfavorablel mutations
			while (redesign) {


				//Find the most probable nucleotide to mutate 
				//choose the nucleotide
				i=0;//i is the position in the sequence string (0 indexed)
				double roll = dice.roll();
				double current = def[1]/(length*defect);//def, the per nuc defect array, is 1 indexed
			
				while (roll > current) {
					
					++i;
					//consider each nucleotide that is notthe XXXXXX linker
					if (sequence[i]!='X') {
						current+=def[i+1]/(length*defect);
					}
					else i+=5;
					if (i>= fragment->GetSequenceLength()) {
						cerr << "i ran out of bounds in leaf defect-weighted redesign\n";
					}
				}		
				//END: Find the most probable nucleotide to mutate


				//now the ith nucleotide should be mutated
				sequenceTEMP=sequence[i];//save the ith nuc so that it can be restored later
				int change = ((int)(dice.roll()*3.0))+1;
				//cerr << "BEFORE: sequence[" << i << "]=" << sequence[i];

				//cerr << "i=" << i << endl;
				//cerr << "start=" << start << " missingStart=" << missingstart << " missingEnd=" << missingend << " end=" << end << endl;
				//cerr << "MapFragmenttoNuc=" << MapFragmenttoNuc(i+1,start,missingstart,missingend) << endl;

				sequence[i]= tonuc((toint(sequence[i])+change-1)%4+1);

				//cerr << " AFTER: sequence[" << i << "]=" << sequence[i]  << endl;
				

				//cerr << "CHECK1: sizeof=" << sizeof(sequence) << " strlen=" << strlen(sequence) << " GetSequenceLength=" << GetSequenceLength()<< endl;
				//now check if the position is paired
				if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {

					//	cerr << "MapFragmenttoNuc=" << MapFragmenttoNuc(i+1,start,missingstart,missingend) << " GetPair(MapFrToNuc)=" << GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend));

					//cerr << " MapNucToFragm(GetPair(MapFrToNuc))=" << MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend) << endl;

// cerr << "BEFORE: sequence[" << MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1 << "]=" << sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1] << endl;		
	
//cerr << "sequence[i]=" << sequence[i] << " toint(sequence[i])=" << toint(sequence[i]) << " tonuc(5-toint(sequence[i]))=" << tonuc(5-toint(sequence[i])) << endl;


					sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));

					//cerr << "AFTER: sequence[" << MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1 << "]=" << sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1] << endl;					



				}
				//cerr << "CHECK2: sizeof=" << sizeof(sequence) << " strlen=" << strlen(sequence) << " GetSequenceLength=" << GetSequenceLength()<< endl;
				delete fragment;
				
				//cerr << "Just deleted the fragment" << endl;
				//cerr << "sequence[" << i << "]=" << sequence[i] << " GetBackboneType=" << GetBackboneType()<< endl;
				

				//initialize the fragment
				fragment = new RNA(sequence,GetBackboneType());
				//cerr << "CHECK3: sizeof=" << sizeof(sequence) << " strlen=" << strlen(sequence) << " GetSequenceLength=" << GetSequenceLength()<< endl;
				//cerr << "  FRAGMENT2 size = " << fragment->GetSequenceLength() << endl;

				//for(int p=0;p<=fragment->GetSequenceLength();p++) cerr << "fragment[" << p << "]=" << fragment->GetNucleotide(p) << endl;

				//copy the thermodynamic parameters into fragment
				fragment->CopyThermodynamic(thermo);
				//Perform the partition function calculation
				fragment->PartitionFunction();
				//	cerr << "fragment=" << fragment->GetSequenceLength() << " sequence=" << strlen(sequence) << endl;

				defTEMP = def;


				//reset the values in defect
				for (j=0;j<=fragment->GetSequenceLength();j++) {

					def[j]=0;
				}
				defect=0;
				GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'
				if(missingend!=0) defect= defect/((double)fragment->GetSequenceLength()-6);
				else defect= defect/((double)fragment->GetSequenceLength());


				if (defect<pernucdefect) {
					StoreMutation(start,end,missingstart,missingend,sequence);
					redesign = false;
#ifdef DEBUG
					cerr << " ACCEPTING THE FRAGMENT. Accepted: " << KeepMut << " | Rejected: " << RejectMut << " \n\n";
#endif
				}

				else {
					if (changes>=MaxMutate*fragment->GetSequenceLength()) {
						//give up...
						redesign = false;
#ifdef DEBUG
						cerr << " REJECTING THE FRAGMENT. Accepted: " << KeepMut << " | Rejected: " << RejectMut << "\n\n";
#endif
					}
					else {

						if (defect<defectOLD) {
							//defect improved, but is still not quite as good as desired.
							
							//cerr << "defect=" << defect << endl;
							defectOLD=defect;
#ifdef DEBUG
							
							//cerr << "mutate " << i << " Keeping the mutation, but need to mutate more\n";
							KeepMut++;
#endif
						}
						else {

							//defect got larger with change, so restore the original.
#ifdef DEBUG
							//cerr << "mutate " << i << " Restoring the mutation\n";
							RejectMut++;
#endif
							sequence[i]=sequenceTEMP;

							defect = defectOLD;

							//also restore the pair if needed...
							//now check if the position is paired
							if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
								sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));

							}

							def=defTEMP;

						}

					}
					changes++;

				}


			}

		}
		*/
		//Store best sequence found so far, and best defect found so far:
		if(defect<EDtree[fragmentdepth][start]){
			
			//store this defect
			EDtree[fragmentdepth][start]=defect;
#ifdef DEBUG
			cerr << " defect=" << defect << endl;
#endif
			//Now store the sequence:
			for (i=start;i<=end;++i) {
				if (!(i>=missingstart&&i<=missingend)) {

					SEQtree[fragmentdepth][i] = GetStructure()->nucs[i];
					
				}
			}

		}
		
		//cerr << "test1" << endl;
		if (defect>pernucdefect) {
			int Maximum;
			if (providesequence) Maximum = MaxLeafRedesign;
			else  Maximum = MaxRedesign;
			if (ReoptNum[fragmentdepth][start]<Maximum) {

				//Redesign this fragment by adding all the components back on the stack
				for (j=fragmentdepth;j<depth;j++) {


					
					
					if (j>=0) {
						
						FindFragments(tree,j,start,end,missingstart,missingend,&stackstart, &stackend, &stackmissingstart, &stackmissingend, &stackfragmentdepth);


						////Only redesign half the sequence.
						if (defectweighted&&j>fragmentdepth) {

							if (stackfragmentdepth.size()>1) {
								if (stackfragmentdepth.back()==stackfragmentdepth.at(stackfragmentdepth.size()-2)&&stackfragmentdepth.back()>fragmentdepth) {
									//two fragments at the same level were placed on the stack, remove one
									roll=dice.roll();

									if (roll<(EDtree[j][stackstart.at(stackstart.size()-2)]/(EDtree[j][stackstart.at(stackstart.size()-2)]+EDtree[j][stackstart.at(stackstart.size()-1)]))) {
										//remove the first fragment

								

										stackstart.erase (stackstart.begin() + stackstart.size() - 2); 
										stackend.erase (stackend.begin() + stackend.size() - 2);
										stackmissingstart.erase (stackmissingstart.begin() + stackmissingstart.size() - 2);
										stackmissingend.erase (stackmissingend.begin() + stackmissingend.size() - 2);
										stackfragmentdepth.erase (stackfragmentdepth.begin() + stackfragmentdepth.size() - 2);

									}
									else {
										//remove the second fragment
										stackstart.erase (stackstart.begin() + stackstart.size() - 1); 
										stackend.erase (stackend.begin() + stackend.size() - 1);
										stackmissingstart.erase (stackmissingstart.begin() + stackmissingstart.size() - 1);
										stackmissingend.erase (stackmissingend.begin() + stackmissingend.size() - 1);
										stackfragmentdepth.erase (stackfragmentdepth.begin() + stackfragmentdepth.size() - 1);

									}

									//Now narrow the area to be redesign as we ascend the tree towards leaves

									start = stackstart.back();
									end = stackend.back();
									missingstart = stackmissingstart.back();
									missingend = stackmissingend.back();
							


								}

							}//end if (stackfragmentdepth.size()>=2)
						}//end if (defectweighted&&j>fragmentdepth)
					}//end j<=0

					else {
						
						//Place the whole sequence on the stack:
						stackstart.push_back(1);
						stackend.push_back(GetSequenceLength());
						stackmissingstart.push_back(0);
						stackmissingend.push_back(0);
						stackfragmentdepth.push_back(-1);
					}
										


				}
			
			}//end if ReoptNum[fragmentdepth][start]<MaxRedesign
			else {
				//Will be going up the tree, so reset re-optimization counter
				ReoptNum[fragmentdepth][start]=0;

#ifdef DEBUG
				cerr << "failure to redesign with multiple tries\n";
#endif

				
				//Failure at getting a re-design < pernucdefect, so choose the best encountered sequence:
				for (i=start;i<=end;++i) {
					if (!(i>=missingstart&&i<=missingend)) {

						GetStructure()->nucs[i] = SEQtree[fragmentdepth][i];
						if (SEQtree[fragmentdepth][i]=='A') {
							GetStructure()->numseq[i]=1;

						}
						else if (SEQtree[fragmentdepth][i]=='C') {
							GetStructure()->numseq[i]=2;

						}
						else if (SEQtree[fragmentdepth][i]=='G') {
							GetStructure()->numseq[i]=3;

						}
						else if (SEQtree[fragmentdepth][i]=='T'||SEQtree[fragmentdepth][i]=='U') {
							GetStructure()->numseq[i]=4;

						}
						else {

							GetStructure()->numseq[i]=0;
						}
					
					}
				}


			}
		}//end if defect>pernucdefect
		else { ///defect<=pernucdefect
			//success at this level, reset re-optimization counter

			ReoptNum[fragmentdepth][start]=0;

		}//end else  (defect<=pernucdefect)
		//		cerr << "completed" << endl;
		//clean up memory use
		delete fragment;
		//cerr << "CHECK2" << endl;
		//cerr << "TEST5" << endl;
	}//end while stackstart.size()>0
	
	//#ifdef DEBUG
	//cout << "NED= " << EDtree[-1][1] << endl << flush; 
	//#endif
	//cerr << "DONE" << endl;

	double ned = EDtree[-1][1];

	//clean up memory use
	delete[] sequence;
	//delete[] sequenceTEMP;
	//delete thermo;
	EDtree=EDtree-1;
	SEQtree=SEQtree-1;
	ReoptNum=ReoptNum-1;
	for(i=0;i<depth+1;i++){
		delete[] EDtree[i];
		delete[] SEQtree[i];
		delete[] ReoptNum[i];
	}
	delete[] EDtree;
	delete[] SEQtree;
	delete[] ReoptNum;




#ifdef DEBUG
	out.close();
	
#endif
	//cerr << "DONE2" << endl;
	return ned;
}


double design::SelectSequenceHeuristic(int **tree, bool random, bool bias, int depth, const double bestpernucdefect, long seed, bool allowisolated) {
	//////////////////////HEADER/////////////////
	int i,j,k,l,m;//,j,k,l;
	RNA *fragment;
	Thermodynamics *thermo;
	vector<int> stackstart, stackend, stackmissingstart, stackmissingend, stackfragmentdepth; //use a vector as a stack of fragments
	//int start,end,missingstart,missingend,fragmentdepth;
	bool providesequence;
	//bool leafOpt;
	randomnumber dice;
	
	//char *sequenceBEST;
	//char sequenceTEMP;
	double defect;//,defectOLD;
	//double pernucdefect;//,temppernucdefect;
	//int frag1,frag2,ident1;
	//double roll;

	char *sequence;

#ifdef DEBUG
	ofstream out;
	out.open("stack-pulls.out");
#endif

	//initialize the random number generator
	dice.seed(seed);

	
	//sequenceTEMP = new char [GetSequenceLength()+9];
	//Read the thermodynamics from disk
	thermo = new Thermodynamics(GetBackboneType());
	thermo->ReadThermodynamic();
	///////////////////////END HEADER///////////////////////

	//allocate sequence large enough to handle anything
	sequence = new char [GetSequenceLength()+13];


	//Find all fragments.
	for (i=-1;i<depth;++i) {
		//For the current level, find each fragment, select a sequence if neccessary, and test that the correct structure forms
		//Identify fragments from left to right
		if (i>=0) FindFragments(tree,i,1,GetSequenceLength(),0,0,&stackstart, &stackend, &stackmissingstart, &stackmissingend, &stackfragmentdepth);
		else {
			//Place the whole sequence on the stack:
			stackstart.push_back(1);
			stackend.push_back(GetSequenceLength());
			stackmissingstart.push_back(0);
			stackmissingend.push_back(0);
			stackfragmentdepth.push_back(-1);
		}
	}//end for i


	//Throw out everything but the leaves
	for (i=0;i<stackstart.size();i++) {

		//Provide a sequence if none has been provided for this fragment (because we are at the leaves of the tree)
		providesequence = false;
		if (stackfragmentdepth[i]==depth-1) providesequence=true;
		if (stackfragmentdepth[i]<depth-1) {
			if (tree[stackfragmentdepth[i]+1][stackstart[i]]==0) providesequence=true;
		}

		if (!providesequence) {


			stackstart.erase(stackstart.begin()+i);
			stackend.erase(stackend.begin()+i);
			stackmissingstart.erase(stackmissingstart.begin()+i);
			stackmissingend.erase(stackmissingend.begin()+i);
			stackfragmentdepth.erase(stackfragmentdepth.begin()+i);
			--i;

		}

	}

#ifdef DEBUG
	//output the contents of the stacks for debugging
	for (i=0;i<stackstart.size();i++) {

		out << stackstart[i]<<"\t"<<stackend[i]<<"\t"<<stackmissingstart[i]<<"\t"<<stackmissingend[i]<<"\t"<<stackfragmentdepth[i]<<"\n"<<flush;

	}
#endif

	

	double *EDtree;//Keep track of best ensemble defect for each tree decomposition level
	EDtree=new double [GetSequenceLength()+1];
	for (i=0;i<=GetSequenceLength();++i) EDtree[i]=2.0; //initilaize to large value


	char *SEQtree;//Keep track of the sequence with best ensemble defect for each tree decomposition level
	SEQtree=new char [GetSequenceLength()+1];
	
	int *ReoptNum;//Keep track of the number of reoptimizations for each tree decomposition level
	ReoptNum=new int [GetSequenceLength()+1];
	for (i=0;i<=GetSequenceLength();++i) ReoptNum[i]=0;	


	vector<vector<string> > Helices (11,vector<string>());
	vector<vector<string> > Loops (11,vector<string>());
	//Read the helices and loops from the database if 'random==false'
	if(!random){
			
		//Read datapath.
		string hf;
		string lf;
		if(GetBackboneType()){//if GetBackbone returns true, read RNA sequence files
			hf="design.Helices.dat";//Holds the file name with helix sequences.
			lf="design.Loops.dat";//Holds the file name with loop sequences
		}
		else{//else, if GetBackbone returns false, read DNA files
			hf="design.DNA.Helices.dat";//Holds the file name with helix sequences.
			lf="design.DNA.Loops.dat";//Holds the file name with loop sequences
		}		
		string dp(getDataPath());//Holds the path to data_tables.
			
		//Append the filename to the path 'dp' and store in 'fn'.
		hf=dp+"/"+hf;
		lf=dp+"/"+lf;

		string lineoftext;//String that temporarily holds strings of data from the read files.
		//Open file 
		ifstream readHfile(hf.c_str());
		ifstream readLfile(lf.c_str());
		//			cerr << "lf=" << lf.c_str() << endl;
		while(!readHfile.eof()){
			readHfile >> lineoftext;
			if(lineoftext.empty()) continue;
			Helices.at((lineoftext.size()-3)/2).push_back(lineoftext);
		}
		while(!readLfile.eof()){
			readLfile >> lineoftext;
			//				cerr << "lineoftext=" << lineoftext << endl;
			if(lineoftext.empty()) continue;
			Loops.at(lineoftext.size()).push_back(lineoftext);
		}
		readHfile.close();
		readLfile.close();
	}


	//pernucdefect = bestpernucdefect;
	//design a sequence for each leaf
	for (int pos=0;pos<stackstart.size();++pos) {
		
		EDtree[stackstart[pos]] = leafdesign(stackstart[pos],stackend[pos],stackmissingstart[pos],stackmissingend[pos],random,bias,&dice,&Helices,&Loops,bestpernucdefect,thermo,allowisolated);

		

	}//end for (pos=0;pos<stackstart.size();++pos) 

	

	//Now check fragments against each other
	for (i=0;i<stackstart.size()-1;++i) {
		for (j=i+1;j<stackstart.size();++j) {


			bool repeat = true;
			int repeats = 0;
			double bestdefect=2.0;

			while (repeat) {
				//Now check these fragments to see if they fold correctly
				//Populate sequence with the information
				k=0;
				int length1;

				for (l=stackstart[i];l<=stackend[i];l++) {
					//iterate over each nucleotide
					if (l==stackmissingstart[i]) {
						l=stackmissingend[i];
						for (m=0;m<6;m++) {
							sequence[k] = 'X';
							k++;
						}
					}
					else {
						sequence[k] = GetStructure()->nucs[l];
						k++;
					}//end else, l!=missingstart
				}//end for l

				length1=k;
				for (m=0;m<6;m++) {
					sequence[k] = 'X';
					k++;
				}

				for (l=stackstart[j];l<=stackend[j];l++) {
					//iterate over each nucleotide
					if (l==stackmissingstart[j]) {
						l=stackmissingend[j];
						for (m=0;m<6;m++) {
							sequence[k] = 'X';
							k++;
						}
					}
					else {
						sequence[k] = GetStructure()->nucs[l];
						k++;
					}//end else, l!=missingstart
				}//end for l
				sequence[k]='\0';//add the null terminator
		
				//initialize the fragment
				fragment = new RNA(sequence,SEQUENCE_STRING,this);
				//copy the thermodynamic parameters into fragment
				//fragment->CopyThermodynamic(thermo);
				//Perform the partition function calculation
				fragment->PartitionFunction("",-10.0,false,true,allowisolated);

				
				//Now tabulate the ensemble defect
				defect = 0.0;


				//check the first fragment
				for (l=1;l<=length1;++l) {

					if (fragment->GetNucleotide(l)!='X') {
						if (GetPair(MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]))==0) {
							//Nucleotide l should be unpaired

							for (m=1;m<=fragment->GetSequenceLength();m++) {
					
								if (m>l) {
							
									defect+=fragment->GetPairProbability(l,m);
								}
					
								else if (l>m) {
						
									defect+=fragment->GetPairProbability(m,l);
								}
							}//end for m


						}
						else {
							//Nucleotide l should be paired
							int p5,p3;
							if (GetPair(MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]))>MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i])) {
								p5 = MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]);
								p3 = GetPair(MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]));

							}
							else {
								p3 = MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]);
								p5 = GetPair(MapFragmenttoNuc(l,stackstart[i],stackmissingstart[i],stackmissingend[i]));

							}
							defect+=GetPairProbability(p5,p3);

						}

					}


				}//end for l
				//check the second fragment:
				for (l=length1+6;l<=fragment->GetSequenceLength();++l) {

					if (fragment->GetNucleotide(l)!='X') {
						if (GetPair(MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]))==0) {
							//Nucleotide l should be unpaired

							for (m=1;m<=fragment->GetSequenceLength();m++) {
					
								if (m>l) {
							
									defect+=fragment->GetPairProbability(l,m);
								}
					
								else if (l>m) {
						
									defect+=fragment->GetPairProbability(m,l);
								}
							}//end for m


						}
						else {
							//Nucleotide l should be paired
							int p5,p3;
							if (GetPair(MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]))>MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j])) {
								p5 = MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]);
								p3 = GetPair(MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]));

							}
							else {
								p3 = MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]);
								p5 = GetPair(MapFragmenttoNuc(l-6-length1,stackstart[j],stackmissingstart[j],stackmissingend[j]));

							}
							defect+=GetPairProbability(p5,p3);

						}

					}


				}//end for l, check the second fragment


				//normalize defect by length
				int totallength = fragment->GetSequenceLength()-6;
				if (stackmissingstart[i]>0) totallength-=6;
				if (stackmissingstart[j]>0) totallength-=6;

				defect = defect / ((double) totallength);

			
				//calculate the target defect:
				int l1,l2;
				l1 = stackend[i]-stackstart[i]+1;
				if (stackmissingstart[i]>0) l1-=(stackmissingend[i]-stackmissingstart[i]+1);
				l2 = stackend[j]-stackstart[j]+1;
				if (stackmissingstart[j]>0) l2-=(stackmissingend[j]-stackmissingstart[j]+1);
				double targetdefect = (max(EDtree[stackstart[i]],bestpernucdefect)*((double)l1)+max(EDtree[stackstart[j]],bestpernucdefect)*((double)l2)) / ((double)(l1+l2));

				if (defect>targetdefect) {
					//This isn't yet as good as we want


					++repeats;
					if (defect<bestdefect) {
						//this is the best seen yet, store this sequence

						bestdefect = defect;
						for (l=stackstart[j];l<=stackend[j];++l) {
							if (l==stackmissingstart[j]) l=stackmissingend[j];
							else {

								SEQtree[l]=GetStructure()->nucs[l];

							}

						}


					}
					else {
						//This is not the best, restore the best seen yet
						for (l=stackstart[j];l<=stackend[j];++l) {
							if (l==stackmissingstart[j]) l=stackmissingend[j];
							else {

								GetStructure()->nucs[l]=SEQtree[l];

							}

						}

					}

					if (repeats>MaxRedesign) {
						//Too many redesigns, stop now

						repeat = false;
					}
					else {

						EDtree[stackstart[j]] = leafdesign(stackstart[j],stackend[j],stackmissingstart[j],stackmissingend[j],random,bias,&dice,&Helices,&Loops,bestpernucdefect,thermo,allowisolated);

					}


					


				}

				else {
					//This defect is better than we expected, stop redesigning
					repeat = false;

				}

				//clean up
				delete fragment;

			}//end while repeat
		}//end for j
	}//end for i


	//TODO: Check this -- We want the ensemble deffect of the resulting sequence.
	double ned = EDtree[0];
	
	delete thermo;

	delete[] EDtree;
	delete[] SEQtree;
	delete[] ReoptNum;

	delete[] sequence;




#ifdef DEBUG
	out.close();
	
#endif
	//cerr << "DONE2" << endl;
	return ned;

}


//used by SelectSequenceHeuristic to design a sequence for a leaf
//This is largely code duplicated from SelectSequence and this should later be unified in one function
//Returns the defect for the leaf

double design::leafdesign(int start,int end,int missingstart,int missingend, bool random, bool bias, randomnumber *dice, vector<vector<string> > *Helices,vector<vector<string> > *Loops, double pernucdefect, Thermodynamics *thermo, bool allowisolated) {
	char *sequence,*bestseq;
	RNA *fragment;
	double defect,defectOLD,bestdefect;
	int reopt=0;
	char sequenceTEMP, sequenceTEMP2;
	bool godesign;
	int i,j,k,l;

	//initilaize the sequence to be large neough for any situation
	sequence = new char [GetSequenceLength()+9];
	bestseq = new char [GetSequenceLength()+9]; 
	
	
	//++ReoptNum[start];
	godesign=true;
	bestdefect=2;
	while(godesign) {	

		++reopt;
	
		FillSequence(start,end,missingstart,missingend,random,bias,dice,*Helices,*Loops);
	
		
		//Now check this fragment to see if it folds correctly
		//Populate sequence with the information
		k=0;

		for (j=start;j<=end;j++) {
			//iterate over each nucleotide
			if (j==missingstart) {
				j=missingend;
				for (l=0;l<6;l++) {
					sequence[k] = 'X';
					k++;
				}
			}
			else {
				sequence[k] = GetStructure()->nucs[j];
				k++;
			}//end else, j!=missingstart
		}//end for j
		sequence[k]='\0';//add the null terminator
		
		//initialize the fragment
		fragment = new RNA(sequence,SEQUENCE_STRING,this);
		//copy the thermodynamic parameters into fragment
		//fragment->CopyThermodynamic(thermo);
		//Perform the partition function calculation
		fragment->PartitionFunction("",-10.0,false,true,allowisolated);

				
		//Now tabulate the ensemble defect
		defect = 0.0;
		
		vector<double> def (fragment->GetSequenceLength()+1,0);//stores ensemble defect for each nucleotide in a fragment
		
		
		GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'

		//cerr << "check1.1" << endl;
		//If XXXXXX were used, subtract them from the total length of the fragment
		if(missingend!=0) defect= defect/((double)fragment->GetSequenceLength()-6);
		else defect= defect/((double)fragment->GetSequenceLength());

		if (defect>pernucdefect&&random) {
			vector<double> defTEMP (fragment->GetSequenceLength()+1,0);//For temporary storage of 'def'
			//vector<double> defTOTAL (fragment->GetSequenceLength()+1,0);
			//vector<int> Mutated(fragment->GetSequenceLength()+1,0);//keeps track which nucleotide was mutated and how many times during the leaf optimization
			//cerr << "fragment_length=" << fragment->GetSequenceLength() <<" ";

			bool redesign=true;
			int changes = 0;
			defectOLD = defect;
			double length;

			if(missingend!=0) length = (double)fragment->GetSequenceLength()-6;
			else length = (double)fragment->GetSequenceLength();

			while (redesign) {

				//choose the nucleotide
				i=0;//i is the position in the sequence string (0 indexed)
				double roll = dice->roll();
				double current = def[1]/(length*defect);//def, the per nuc defect array, is 1 indexed

				while (roll > current) {
					
					++i;
					//if (MapFragmenttoNuc(i+1,start,missingstart,missingend)<missingstart||MapFragmenttoNuc(i+1,start,missingstart,missingend)>missingend||missingend==0) {
					//consider each nucleotide that is notthe XXXXXX linker
					if (sequence[i]!='X') {
						current+=def[i+1]/(length*defect);
					}
#ifdef DEBUG
					if (i>= fragment->GetSequenceLength()) {
						cerr << "i ran out of bounds in leaf defect-weighted redesign\n";
					}
#endif


				}		

				//char newnuc1, newnuc2;
				//now the ith nucleotide should be mutated
				sequenceTEMP=sequence[i];//save the ith nuc so that it can be restored later
				//int change = ((int)(dice->roll()*3.0))+1;
				//sequence[i]= tonuc((toint(sequence[i])+change-1)%4+1);

				//now the i-th nucleotide should be mutated by 'change'
				if (bias_in_leaf_refinement) {
					if (GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)) > 0) {
						int newnuc1int, newnuc2int;
						randompair(&newnuc1int, &newnuc2int, dice);
						sequence[i] = tonuc(newnuc1int);
						sequenceTEMP2 = sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1];
						sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1] = tonuc(newnuc2int);

					}
					else {
						sequence[i] = tonuc(randomnuc(dice));
					}

				}
				else {
					int change = ((int)(dice->roll() * 3.0)) + 1;
					sequence[i] = tonuc((toint(sequence[i]) + change - 1) % 4 + 1);
					if (GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)) > 0) {
						sequenceTEMP2 = sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1];
						sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1] = tonuc(5 - toint(sequence[i]));
					}
				}

				//now check if the position is paired
				//if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
				//	sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));

				//}

				delete fragment;
					
				//initialize the fragment
				fragment = new RNA(sequence,SEQUENCE_STRING,this);
				//copy the thermodynamic parameters into fragment
				//fragment->CopyThermodynamic(thermo);
				//Perform the partition function calculation
				fragment->PartitionFunction("",-10.0,false,true,allowisolated);

				defTEMP = def;


				//reset the values in defect
				for (j=0;j<=fragment->GetSequenceLength();j++) {

					def[j]=0;
				}
				defect=0;
				GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'
				if(missingend!=0) defect= defect/((double)fragment->GetSequenceLength()-6);
				else defect= defect/((double)fragment->GetSequenceLength());


				if (defect<pernucdefect) {
					StoreMutation(start,end,missingstart,missingend,sequence);
					redesign = false;

				}
				else {
					
					if (changes>=MaxMutate*fragment->GetSequenceLength()) {
						//give up...
						redesign = false;

					}
					else {

						if (defect<defectOLD) {
							//defect improved, but is still not quite as good as desired.

							defectOLD=defect;
							

						}
						else {

							//defect got larger with change, so restore the original.

							sequence[i]=sequenceTEMP;

							defect = defectOLD;

							//also restore the pair if needed...
							//now check if the position is paired
							//if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
							//	sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));

							//}

							if (GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)) > 0) {
								sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i + 1, start, missingstart, missingend)), start, missingstart, missingend) - 1] = sequenceTEMP2;
							}

							def=defTEMP;

						}

					}
					changes++;

				}


			}

		}//end leaf redesign -- defect>pernucdefect&&random
		
		
		if (defect<=pernucdefect) {
			godesign=false;

		}
		else if (reopt>MaxLeafRedesign) {

			godesign=false;

		}
		

		

		

		//Store best sequence found so far, and best defect found so far:
		if(defect<bestdefect){
			
			//store this defect
			bestdefect=defect;
			j=0;
			//Now store the sequence:
			for (i=start;i<=end;++i) {
				if (!(i>=missingstart&&i<=missingend)) {

					bestseq[j] = GetStructure()->nucs[i];
					++j;
					
				}
			}

		}//end defect<EDtree[start]
		else {
			//restore the best sequence so far

			j=0;
			for (i=start;i<=end;++i) {
				if (!(i>=missingstart&&i<=missingend)) {

					GetStructure()->nucs[i] = bestseq[j];
					GetStructure()->numseq[i] = toint(GetStructure()->nucs[i]);
					++j;
					
				}
			}

			defect = bestdefect;


		}
		
		//clean up memory use
		delete fragment;

	}//end while (godesign)
	//cerr << "test1" << endl;

	//		cerr << "completed" << endl;
	

	//clean up memory use
	delete[] sequence;
	delete[] bestseq;

	return defect;

}


void design::FindFragments(int **tree,int level, int start, int stop, int missingstart, int missingstop, vector<int> *stackstart, vector<int> *stackend, vector<int> *stackmissingstart, vector<int> *stackmissingend, vector<int> *stackfragmentdepth) {
	int i,j;
	int a,b,c;
	double defect1,defect2;

	i = start;
	while (i<=stop) {
		//fast forward through zeros (these are fragments not split at this level
		if (tree[level][i]==0) i++;

		else {
			//define the fragment starting at i
			a=i;
			while (tree[level][a]==tree[level][i]) {

				++a;
				if (a>stop) {
					//stop running forward if there is no room left
					break;
				}

			}
			--a;

			//now a ends the first (maybe only fragment)

			//look for any other fragment with same label
			b=a+1;
			if (b<=stop) {
				while (tree[level][b]!=tree[level][i]) {
					++b;
					if (b>stop) {
						//stop running forward if there is no room left
						break;
					}
				}
			}
			//b now marks another fragment or is past stop
			if (b<=stop) {
				c=b+1;
				if (c<=stop) {
					while (c<=stop&&tree[level][c]==tree[level][i]) {
						++c;
						if (c>stop) {
							//stop running forward if there is no room left
							break;
						}
					}
					
				}
				--c;

				stackstart->push_back(i);
				stackend->push_back(c);
				stackmissingstart->push_back(a+1);
				stackmissingend->push_back(b-1);
				stackfragmentdepth->push_back(level);

				//put the middle fragment on the stack unless it coincides with the missing region
				
				if (a+1==missingstart&&b-1!=missingstop) {
					//one end matches a missing region, so shorten
					FindFragments(tree,level,missingstop+1,b-1,0,0,stackstart,stackend,stackmissingstart,stackmissingend,stackfragmentdepth);


				}
				else if (b-1==missingstop&&a+1!=missingstop) {
					//other end matches a missing region, so shorten
					FindFragments(tree,level,a+1,missingstart-1,0,0,stackstart,stackend,stackmissingstart,stackmissingend,stackfragmentdepth);

				}
				else if (!(a+1==missingstart&&b-1==missingstop)) {
					//Both ends do not match missing region
					FindFragments(tree,level,a+1,b-1,missingstart,missingstop,stackstart,stackend,stackmissingstart,stackmissingend,stackfragmentdepth);

				}
				i=c+1;
				 

			}
			else {
				stackstart->push_back(i);
				stackend->push_back(a);
				stackmissingstart->push_back(0);
				stackmissingend->push_back(0);
				stackfragmentdepth->push_back(level);

				i = a +1;

			}
			
			
		}

	}

	return;
}

void design::FillSequence(int start, int end, int missingstart, int missingend, bool random, bool bias,  randomnumber *dice, vector<vector<string> > &Helices, vector<vector<string> > &Loops) {
	double roll;

	//Fill in the sequence
	for (int j=start;j<=end;j++) {

		//iterate over each nucleotide
		if (j==missingstart) j=missingend;
		else {
			if(bias){
				if (GetPair(j)==0) {
					GetStructure()->numseq[j]=randomnuc(dice);

					GetStructure()->nucs[j]=data->numtobase(GetStructure()->numseq[j]);
				}
				else {

					int p5,p3;
					randompair(&p5,&p3,dice);
					GetStructure()->numseq[j]=p5;
					GetStructure()->nucs[j]=data->numtobase(p5);
					GetStructure()->numseq[GetPair(j)]=p3;
					GetStructure()->nucs[GetPair(j)]=data->numtobase(p3);
				}

			}

			else if (random && !bias) {
				//USE a random selection for the sequence

				roll = dice->roll();
				if (GetPair(j)==0) {
					//j is unpaired, pull a nuc from a flat distribution
					if (roll<0.25) {
						//A
						GetStructure()->numseq[j]=1;
						GetStructure()->nucs[j]='A';
					}
					else if (roll<0.5) {
						//C
						GetStructure()->numseq[j]=2;
						GetStructure()->nucs[j]='C';
					}
					else if (roll<0.75) {
						//G
						GetStructure()->numseq[j]=3;
						GetStructure()->nucs[j]='G';
					}
					else {
						//U/T
						GetStructure()->numseq[j]=4;
						if (GetBackboneType()) GetStructure()->nucs[j]='U';
						else GetStructure()->nucs[j]='T';
					}
								
				}//end if GetPair(j)==0
				else if (GetPair(j)>j) {
					//j is paired to a 3' partner, choose the sequence
					if (roll<0.25) {
						//A-U
						GetStructure()->numseq[j]=1;
						GetStructure()->nucs[j]='A';
						GetStructure()->numseq[GetPair(j)]=4;
						if (GetBackboneType()) GetStructure()->nucs[GetPair(j)]='U';
						else GetStructure()->nucs[GetPair(j)]='T';

					}
					else if (roll<0.5) {
						//C-G
						GetStructure()->numseq[j]=2;
						GetStructure()->nucs[j]='C';
						GetStructure()->numseq[GetPair(j)]=3;
						GetStructure()->nucs[GetPair(j)]='G';
					}
					else if (roll<0.75) {
						//G-C
						GetStructure()->numseq[j]=3;
						GetStructure()->nucs[j]='G';
						GetStructure()->numseq[GetPair(j)]=2;
						GetStructure()->nucs[GetPair(j)]='C';
					}
					else {
						//U-A
						GetStructure()->numseq[j]=4;
						if (GetBackboneType()) GetStructure()->nucs[j]='U';
						else GetStructure()->nucs[j]='T';
						GetStructure()->numseq[GetPair(j)]=1;
						GetStructure()->nucs[GetPair(j)]='A';
					}

				}//end if GetPair(j)>j

			}//end if random
			//IF NOT random:
			else if(!random && !bias) {
				//cerr << "test1" << endl;
				//cerr << "InLoop" << " GetPair(" << j << ")=" << GetPair(j)<< endl;
				if(GetPair(j)!=0&&j<GetPair(j)){//If it is a 5' part of the helix
					int Hcount=1;//set counter of helix length to 1
					vector<int> HelixSplit;//initiate a vector 'HelixSplit' to hold the start positions of the helix
					HelixSplit.push_back(j);//store the start position of the helix
					while(GetPair(j)-1==(GetPair(j+1))&&j!=missingstart){//while there is a stack without any interruptions
						Hcount++;//add 1 to the helix length
						j++;//go to the next nucleotide
					}
					if(Hcount>10){//if the helix is longer then 10 base pairs
						while(Hcount>20){//while the helix is longer then 20 base pairs
							Hcount-=10;//trim the helix by 10 base pairs
							HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+10);//store the location where the starting location of the next helix will be
						}

						int split=((int)floor((double)Hcount/2));//'split' splits the helix that is at most 20 base pairs long into two strands of largest size
						//while(split>10||Hcount-split>10)split=((int)floor((double)Hcount/2));//if any of the strands are longer then 10 base pairs, re-split the helix
						HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+split);//store the start location of the second split helix

						Hcount-=split;
						HelixSplit.push_back(HelixSplit[HelixSplit.size()-1]+Hcount);//store the location of the end of the second split helix +1

					}//END if(Hcount>10)
					else HelixSplit.push_back(j+1);//if the helix length is less then or equal to 10, store the location of the end of the helix +1
		
					vector<string> seqlist;
					for(int y=0;y<HelixSplit.size()-1;y++){
						roll = dice->roll();
						seqlist.push_back(Helices[HelixSplit[y+1]-HelixSplit[y]][(int)floor((double)Helices[HelixSplit[y+1]-HelixSplit[y]].size()*roll)]);//should I cast int instead of floor?
					}

					//const char* seq=Helices[floor((double)Helices[Hcount].size()*roll)].c_str();//should I cast int instead of floor?
						
					int k=0;
					for(int y=0;y<HelixSplit.size()-1;y++){
						//cerr << "Helices:" << HelixSplit[y] << "-" << HelixSplit[y+1] << endl;
						Hcount=HelixSplit[y+1]-HelixSplit[y];
						int start5=HelixSplit[y];//store the biginning of the helix
						int end5=HelixSplit[y+1];//store the end of the helix. end is not included in the helix
						const char* seq=seqlist[y].c_str();
						for(int i=start5;i<end5;++i){
							GetStructure()->nucs[i]=seq[k];
							if(seq[k]=='A') GetStructure()->numseq[i]=1;
							else if(seq[k]=='C') GetStructure()->numseq[i]=2;
							else if(seq[k]=='G') GetStructure()->numseq[i]=3;
							else if(seq[k]=='U') GetStructure()->numseq[i]=4;
							else if(seq[k]=='T') GetStructure()->numseq[i]=4;
							GetStructure()->nucs[GetPair(i)]=seq[Hcount*2+2-k];//TAKING INTO ACCOUNT THE 3 SEPERATORS BETWEEN HELIX STRANDS
							if(seq[Hcount*2+2-k]=='A') GetStructure()->numseq[GetPair(i)]=1;
							else if(seq[Hcount*2+2-k]=='C') GetStructure()->numseq[GetPair(i)]=2;
							else if(seq[Hcount*2+2-k]=='G') GetStructure()->numseq[GetPair(i)]=3;
							else if(seq[Hcount*2+2-k]=='U') GetStructure()->numseq[GetPair(i)]=4;
							else if(seq[Hcount*2+2-k]=='T') GetStructure()->numseq[GetPair(i)]=4;
							k++;
						}//END for(int i=start5;i<end5;++i){
						k=0;
					}// END for(int y=0;y<HelixSplit.size()-1;y++){

				}//END if(GetPair(j)!=0&&j<GetPair(j)){

				else if(GetPair(j)==0){
					int Lcount=1;//set counter of loop length to 1
					vector<int> LoopSplit;//initiate a vector 'LoopSplit' to hold the start and end positions of the loop
					LoopSplit.push_back(j);//store the start position of the loop
					while(GetPair(j+1)==0&&j!=missingstart&&j!=GetSequenceLength()){//If the next nucleotide is unpaired, and if j is still a part of the fragment(not equal to missingstart), and if j is not at the end of the sequence
						Lcount++;
						j++;
					}
					if(j==GetSequenceLength())j++;//If j is at the end of the sequence, count it once
					if(Lcount>10){//if the loop is longer then 10 base pairs
						while(Lcount>20){//while the loop is longer then 20 base pairs
							Lcount-=10;//trim the loop by 10 base pairs
							LoopSplit.push_back(LoopSplit[LoopSplit.size()-1]+10);//store the location where the starting location of the next loop will be
						}
						int split=((int)floor((double)Lcount/2));//'split' splits the loop that is at most 20 base pairs long into two strands of largest size
						//while(split>10||Lcount-split>10)split=((int)floor((double)Lcount/2));//if any of the strands are longer then 10 base pairs, re-split the loop
						LoopSplit.push_back(LoopSplit[LoopSplit.size()-1]+split);//store the start location of the second split loop
						Lcount-=split;
						LoopSplit.push_back(LoopSplit[LoopSplit.size()-1]+Lcount);//store the location of the end of the second split loop +1
					}//END if(Lcount>10)
					else LoopSplit.push_back(j+1);//if the loop length is less then or equal to 10, store the location of the end of the loop +1
					//cerr << "test2" << endl;
					vector<string> seqlist;
					for(int y=0;y<LoopSplit.size()-1;y++){
						roll = dice->roll();
						//cerr << "roll=" << roll << " y=" << y << " LoopSplit.size=" << LoopSplit.size() /*<< "(double)Loops[LoopSplit[y+1]-LoopSplit[y]].size()=" << (double)Loops[LoopSplit[y+1]-LoopSplit[y]].size()*/ << " LoopSplit[y+1]=" << LoopSplit[y+1] << " LoopSplit[y]=" << LoopSplit[y] << " Loops[LoopSplit[y+1]-LoopSplit[y]].size()=" << Loops[LoopSplit[y+1] - LoopSplit[y]].size() << " (double)( Loops[LoopSplit[y+1]-LoopSplit[y]].size())/roll=" << (double)( Loops[ LoopSplit[y+1] - LoopSplit[y] ].size() )*roll<< endl;
 seqlist.push_back(Loops[ LoopSplit[y+1] - LoopSplit[y] ][(int) floor((double)( Loops[ LoopSplit[y+1] - LoopSplit[y] ].size() )*roll) ]);//should I cast int instead of floor?
					}
					//const char* seq=Loops[floor((double)Loops[Lcount].size()*roll)].c_str();//should I cast int instead of floor?
					int k=0;
					for(int y=0;y<LoopSplit.size()-1;y++){
						//cerr << "Loops:" << LoopSplit[y] << "-" << LoopSplit[y+1] << endl;
						int start5=LoopSplit[y];//store the biginning of the loop
						int end5=LoopSplit[y+1];//store the end of the loop. end is not included in the loop
						const char* seq=seqlist[y].c_str();
						for(int i=start5;i<end5;++i){
							GetStructure()->nucs[i]=seq[k];
							if(seq[k]=='A') GetStructure()->numseq[i]=1;
							else if(seq[k]=='C') GetStructure()->numseq[i]=2;
							else if(seq[k]=='G') GetStructure()->numseq[i]=3;
							else if(seq[k]=='U') GetStructure()->numseq[i]=4;
							else if(seq[k]=='T') GetStructure()->numseq[i]=4;
							k++;
						}
						k=0;
					}//END for(int y=0;y<LoopSplit.size()-1;y++){
					//cerr << "test3" << endl;
				}//END else if(GetPair(j)==0){
			}//END else, !random
		}//END else, (j!=missingstart)
	}//end for j
}

int design::MapNuctoFragment(int j,int start,int missingstart,int missingend) {
	int i,k;


	//correct for the start location
	i = j;
	//cerr << "IN MapNucToFragment: i=" << i << " start=" << start << " missingstart=" << missingstart << " missingend=" << missingend ;
	i-=start;
	++i;
	//cerr << "\n MapNuctoNuc:  j=" << j << " start=" << start ;

	if (j>missingend&&missingend!=0) {
		//k=missingend-missingstart-1;
		i-=(missingend-missingstart+1);
		//i=i-k;
		i+=6;//add back the six X's

	}
	//cerr  << " FINAL: i=" << i << '\n';
	return i;
	//cerr << " missingstart=" << missingstart << " missingend=" << missingend << " i=" << i << '\n';
}

int design::MapFragmenttoNuc(int j,int start,int missingstart,int missingend) {
	int i,k;

	
	//correct for the start location
	i = j;
	//cerr << "IN MapFragmentToNuc: i=" << i << " start=" << start << " missingstart=" << missingstart << " missingend=" << missingend ;
	//cerr << "i=" << i << " j=" << j;
	i+=start;
	//cerr << " i+=start=" << i;
	--i;
	//cerr << "i=" << i << endl;
	//	cerr << "\n MapFregmenttoNuc:  j=" << j << " start=" << start ;

	if (i>=missingstart&&missingend!=0) {
		//k=missingend-missingstart-1;
		i+=(missingend-missingstart+1);
		//cerr << " InMissingEnd: i=" << i << endl;
		//i=i+k;
		i-=6;//subtract the six X's
		
	}
	//cerr  << " FINAL: i=" << i << '\n';
	return i;

}
/*
void design::Mutation(int maxDefPos, int start, int missingstart, int missingend, char* sequence, vector<int>& Mutated) {

	if(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend))==0){
		if(sequence[maxDefPos-1]=='A') sequence[maxDefPos-1]='C';
		else if(sequence[maxDefPos-1]=='C') sequence[maxDefPos-1]='G';
		else if(sequence[maxDefPos-1]=='G'){ 
			if(GetBackboneType()) sequence[maxDefPos-1]='U';
			else sequence[maxDefPos-1]='T';
		}
		else if(sequence[maxDefPos-1]=='T'||sequence[maxDefPos-1]=='U') sequence[maxDefPos-1]='A';
		Mutated.at(maxDefPos)++;
	}
	else {
		if(sequence[maxDefPos-1]=='A'){
			sequence[maxDefPos-1]='C';
			sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='G';
		}
		else if(sequence[maxDefPos-1]=='C'){
			sequence[maxDefPos-1]='G';
			sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='C';
		}
		else if(sequence[maxDefPos-1]=='G'){ 
			if(GetBackboneType()){
				sequence[maxDefPos-1]='U';
				sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='A';
			}
			else{
				sequence[maxDefPos-1]='T';
				sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='A';
			}
		}
		else if(sequence[maxDefPos-1]=='T'){
			sequence[maxDefPos-1]='A';
			sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='T';
		}
		else if(sequence[maxDefPos-1]=='U'){
			sequence[maxDefPos-1]='A';
			sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(maxDefPos,start,missingstart,missingend)),start,missingstart,missingend)-1]='U';
		}
		Mutated.at(maxDefPos)++;
	}
}
*/

//def is an array that stores the per nucleotide ensemble defect as a function of position.  This needs to be zeroed before passing to design.
//defect will store the total defect and this needs to be zeroed before calling this function.
void design::GetDefect(int start, int end, int missingstart, int missingend, vector<double>& def, double& defect, RNA* fragment) {
	int j,k,l;
	//cerr << "InGD" << endl;
	for (j=start;j<=end;j++) {
		
		//iterate over each nucleotide
		if (j==missingstart) j=missingend;
		else {
			if (GetPair(j)==0) {
				
				//j is unpaired
				k = MapNuctoFragment(j,start,missingstart,missingend);
				
				//Get the probability that k is paired in fragment, that is the defect
				for (l=1;l<=fragment->GetSequenceLength();l++) {
					
					if (l>k) {
						def.at(k)+=fragment->GetPairProbability(k,l);
						defect+=fragment->GetPairProbability(k,l);
					}
					
					else if (k>l) {
						def.at(k)+=fragment->GetPairProbability(l,k);
						defect+=fragment->GetPairProbability(l,k);
					}
				}//end for l
			}//end if GetPair(j)==0
			
			else if (GetPair(j)>j) {
				//cerr << "j=" << j << " GetPair(j)=" << GetPair(j) << " start=" << start << " end=" << end << " missingstart=" << missingstart << " missingend=" << missingend << endl;
				//j is paired to a 3' partner
				defect+=2*(1-fragment->GetPairProbability(MapNuctoFragment(j,start,missingstart,missingend),MapNuctoFragment(GetPair(j),start,missingstart,missingend)));

				def.at(MapNuctoFragment(j,start,missingstart,missingend))=1-fragment->GetPairProbability(MapNuctoFragment(j,start,missingstart,missingend),MapNuctoFragment(GetPair(j),start,missingstart,missingend));
				//cerr << "before_def.at 1=" << MapNuctoFragment(GetPair(j),start,missingstart,missingend) << endl;
				def.at(MapNuctoFragment(GetPair(j),start,missingstart,missingend))=1-fragment->GetPairProbability(MapNuctoFragment(j,start,missingstart,missingend),MapNuctoFragment(GetPair(j),start,missingstart,missingend));
				//cerr << "after_def.at" << endl;
			}//end if GetPair(j)>j
		}//!end else, (j!=missingstart)
	}//end for j
}


void design::Debug1(int start, int end, int missingstart, int missingend, char* sequence, RNA* fragment) {
	cerr << "start=" << start << " missingstart=" << missingstart << " missingend=" << missingend << " end=" << end << endl;
	for(int q=0;q<fragment->GetSequenceLength();q++){
		cerr << sequence[q] << '\t';}
	cerr << endl;
	for(int q=start;q<=end;q++){
		if(q==missingstart){
			q=missingend;
			cerr << "X\tX\tX\tX\tX\tX\t";
		}
		else cerr << GetStructure()->nucs[q] << '\t';
	}
	cerr << endl;
	for(int q=0;q<fragment->GetSequenceLength();q++){
		cerr << q+1 << '\t';}
	cerr << endl;
	//for(int q=0;q<fragment->GetSequenceLength();q++){
	//	cerr << fragment->GetNucleotide(q+1)<< '\t';}
	//cerr << "\n";
	for(int q=1;q<=fragment->GetSequenceLength();q++){
		if(q<missingstart||q>missingstart+5||missingend==0){
			if(GetPair(MapFragmenttoNuc(q,start,missingstart,missingend))==0) cerr << "-\t";
			else cerr << "P\t";
		}
		else cerr << '\t';
	}
	cerr << "\n";
	for(int q=1;q<=fragment->GetSequenceLength();q++){
		if(q<missingstart||q>missingstart+5||missingend==0){
			cerr << MapFragmenttoNuc(q,start,missingstart,missingend) << "\t";
		}
		else cerr << '\t';
	}				
	cerr << endl;
}

void design::SpecifyRedesignLimits(int leaf, int parent, int mutate) {

	MaxRedesign=parent;
	MaxMutate=mutate;
	MaxLeafRedesign=leaf;


}

void design::SpecifyWeightedDefect(bool DefectWeighted) {

	defectweighted =  DefectWeighted;
}

void design::PlaceSeqOnStack(vector<int> *stackstart, vector<int> *stackend, vector<int> *stackmissingstart, vector<int> *stackmissingend, vector<int> *stackfragmentdepth) {
	//Place the whole sequence on the stack:
	stackstart->push_back(1);
	stackend->push_back(GetSequenceLength());
	stackmissingstart->push_back(0);
	stackmissingend->push_back(0);
	stackfragmentdepth->push_back(-1);
}


void design::StoreMutation(int start, int end, int missingstart, int missingend, char* sequence) {
	int k=0;
	for (int j=start;j<=end;j++) {
		//iterate over each nucleotide
		if (j==missingstart){
			j=missingend;
			k+=6;
		}
		else {
			if (GetPair(j)==0) {
				GetStructure()->nucs[j]=sequence[k];
				//cerr << " sequence[" << k << "-1]=" << sequence[k];
				if(sequence[k]=='A') GetStructure()->numseq[j]=1;
				else if(sequence[k]=='C') GetStructure()->numseq[j]=2;
				else if(sequence[k]=='G') GetStructure()->numseq[j]=3;
				else if(sequence[k]=='U'||sequence[k]=='T') GetStructure()->numseq[j]=4;
				k++;
			}
			else if (GetPair(j)!=0) {
				//cerr << " sequence[" << k << "-1]=" << sequence[k];
				GetStructure()->nucs[j]=sequence[k];
				//GetStructure()->nucs[GetPair(j)]=sequence[GetPair(j)-1];
				if(sequence[k]=='A'){
					GetStructure()->numseq[j]=1;
					//GetStructure()->numseq[GetPair(j)]=4;
				}
				else if(sequence[k]=='C'){
					GetStructure()->numseq[j]=2;
					//GetStructure()->numseq[GetPair(j)]=3;
				}
				else if(sequence[k]=='G'){
					GetStructure()->numseq[j]=3;
					//GetStructure()->numseq[GetPair(j)]=2;
				}
				else if(sequence[k]=='U'||sequence[k]=='T'){
					GetStructure()->numseq[j]=4;
					//GetStructure()->numseq[GetPair(j)]=1;
				}
				k++;
			}
		}
	}
}

void design::StoreBestSequence(int start, int end, int missingstart, int missingend, char** sequence, int fragmentdepth) {
	int k=0;
	for (int j=start;j<=end;j++) {
		//iterate over each nucleotide
		if (j==missingstart){
			j=missingend;
			for(int l=k;l<k+6;++l) sequence[fragmentdepth][l]='X';
			k+=6;
		}
		else {
			sequence[fragmentdepth][k]=GetStructure()->nucs[j];
			k++;
		}
	}
}

char design::tonuc(int i){
	/*if (i==1) return 'A';
	else if (i==2) return 'C';
	else if (i==3) return 'G';
	else if (i==4) {

		if (GetBackboneType()) return 'U';
		else return 'T';

	}
	else return 'X';*/

	//update to extended alphabet code
	return GetStructure()->GetThermodynamicDataTable()->numtobase(i);

}

int design::toint(char i){
	/*if (i=='A') return 1;
	else if (i=='C') return 2;
	else if (i=='G') return 3;
	else if (i=='T'||i=='U') return 4;
	else return 0;*/


	//update to extended alphabet code
	return GetStructure()->GetThermodynamicDataTable()->basetonum(i);

}


//This whole function is unused at this time:
void design::LeafOptimize(const double pernucdefect, randomnumber &dice, double& defect, vector<double>& def, RNA* fragment, int start, int end, int missingstart, int missingend, char* sequence, Thermodynamics *thermo, bool allowisolated){
	
	bool redesign=true;
	double length;//stores the length of the fragment exclusing the 'xxxxxx'
	if(missingend!=0) length = (double)fragment->GetSequenceLength()-6;//set length if there is 'xxxxxx'
	else length = (double)fragment->GetSequenceLength();//set length if there is NO 'xxxxxx'
	int CountUnfavorable=0;//keeps track of unfavorable mutations
	vector <vector<bool> > Unfavorable (fragment->GetSequenceLength()+1,vector<bool> (GetStructure()->GetThermodynamicDataTable()->alphabet.size(),false));//vector keeps track of unfavorablel mutations
	char sequenceTEMP;//temporary stores the pre-mutated nucleotide in case of rejected mutation
	double defectOLD;//temporary stores the pre-mutated defect in case of rejected mutation
	while(defect>pernucdefect&&CountUnfavorable<MaxMutate*length){
		
		//Find the most probable nucleotide to mutate 
		//choose the nucleotide
		int i=0;//i is the position in the sequence string (0 indexed)
		double roll = dice.roll();
		double current = def[1]/(length*defect);//def, the per nuc defect array, is 1 indexed
		
		while (roll > current) {
			
			++i;
			//consider each nucleotide that is notthe XXXXXX linker
			if (sequence[i]!='X') {
				current+=def[i+1]/(length*defect);
			}
			else i+=5;
			if (i>= fragment->GetSequenceLength()) {
				cerr << "i ran out of bounds in leaf defect-weighted redesign\n";
			}
		}		
		//END: Find the most probable nucleotide to mutate
		
		//now the i-th nucleotide should be mutated by 'change'
		int change = ((int)(dice.roll()*3.0))+1;
		
		//If the mutation was already rejected before, return to the begining of 'while' loop and add 1 to CountUnfavorable
		if(Unfavorable[i][(toint(sequence[i])+change-1)%4+1]) CountUnfavorable++;
		else{//If the mutation was not rejected before
			sequenceTEMP=sequence[i];//save the i-th nuc so that it can be restored later
			defectOLD=defect;//save the ensemble defect of the fragment before the mutation
			vector<double> defTEMP (fragment->GetSequenceLength()+1,0);//save the ensemble defect per nucleotide before the fragment mutation
			
			sequence[i]= tonuc((toint(sequence[i])+change-1)%4+1);//mutate the i-th nucleotide
			
			//First check if the position is paired
			if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
				//if the position is paired mutate the paired nucleotide.
				sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));
			}
			///////////////////////////
			//Evaluate the Ensemble Defect for the new fragment with mutated nucleotide
			///////////////////////////
			delete fragment;//delete the fragment with the old sequence
			fragment = new RNA(sequence,SEQUENCE_STRING,this);//initialize the fragment
			//fragment->CopyThermodynamic(thermo);//copy the thermodynamic parameters into fragment
			fragment->PartitionFunction("",-10.0,false,true,allowisolated);//Perform the partition function calculation
			defTEMP = def;//temporarily store the 'def' in defTemp vector
			for (int j=0;j<=fragment->GetSequenceLength();j++) def[j]=0;//restore 'def' vector to all 0s
			defect=0;//restore 'defect' to 0
			GetDefect(start,end,missingstart,missingend,def,defect,fragment);//Calculate defect of the fragment and store defect contribution of each nucleotide in 'def'
			defect= defect/length;
			//END: Evaluate Ensemble Defect for the new fragment with mutated nucleotide
			///////////////////////////
			
			//If the new defect is below the old defect, store the fragment
			if(defect<defectOLD){
				CountUnfavorable=0;
				defectOLD=defect;
				StoreMutation(start,end,missingstart,missingend,sequence);
				for(int p=0;p<=fragment->GetSequenceLength();p++){
					for(int q=0;q<=4;q++) Unfavorable[p][q]=false;
				}
			}
			//If the new defect is above the old defect
			else {
				sequence[i]=sequenceTEMP;//reverse the mutation
				defect = defectOLD;//restore the defect
				//also restore the pair if needed...
				if (GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend))>0) {
					sequence[MapNuctoFragment(GetPair(MapFragmenttoNuc(i+1,start,missingstart,missingend)),start,missingstart,missingend)-1]=tonuc(5-toint(sequence[i]));
				}												
				def=defTEMP;//restore the 'def' - vector that stores defects per nucleotide
				Unfavorable[i][(toint(sequence[i])+change-1)%4+1]=true;//store the rejected mutation
				CountUnfavorable++;//add 1 to CountUnfavorable to keep track of total rejected mutations
			}//END: if the defect is above the old defect
		}//END: if the mutation was not rejected before
	}//END: while(defect<pernucdefect && CountUnfavorable < MaxMutate*length)
}


//return the integer for a random nucleotide choice
//This is used when a bias is being used.
int design::randomnuc(randomnumber *dice){
	double roll = dice->roll();
	double accumulate = 0;



	for (int count=0;count < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++count) {
		accumulate+=singlebias[count];
		if (roll<accumulate) return count;
	}

	//If the cumulative percentage does not equal exactly 1, then return last option...
	for (int count = 0; count < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++count) {
		
		if (singlebias[count]>0.0) return count;
	}


}


//set integer sequence representations for a random pair choice
//This is used when a bias is being used.
void design::randompair(int *i, int *j, randomnumber *dice) {
	double roll = dice->roll();

	double accumulate = 0;

	for (int a=0;a < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++a) {
		for (int b=0; b < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++b) {
			accumulate+=pairbias[a][b];
			if (roll < accumulate) {
				*i = a;
				*j = b;
				return;

			}

		}
	}

	//If the cumulative percentage does not equal exactly 1, then return the first option...
	for (int a = 0; a < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++a) {
		for (int b = 0; b < GetStructure()->GetThermodynamicDataTable()->alphabet.size(); ++b) {
			if (pairbias[a][b] > 0.0) {
				*i = a;
				*j = b;
				return;

			}
		}
	}
	
	return; 


}

