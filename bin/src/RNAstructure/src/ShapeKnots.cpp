// Read description in ShapeKnots.h

#include "ShapeKnots.h"
#include "outputconstraints.h" // HERE

//Flags for output to the screen
//#undef OUTPUT_TO_SCREEN
#define OUTPUT_TO_SCREEN

// Flags for debug features
#undef DEBUG_MODE
//#define DEBUG_MODE

void addtoAggregate(structure * folded, structure * psa, int r){
	//Function that takes a folded structure (after call to dynamic) and adds it to the final list of folded structures
	//folded is the newly folded structure
	//psa is a pointer to the final structure and is equal to pseudoStructAggregate.
	//NOTE!!! only structure r is added to the list, no matter how many tracebacks there are

	psa->AddStructure();	 
	psa->SetCtLabel(folded->GetCtLabel(1),psa->GetNumberofStructures());
	psa->SetEnergy(psa->GetNumberofStructures(),folded->GetEnergy(r));
	for (int i=1;i<=folded->GetSequenceLength();i++){ 
		if (folded->GetPair(i,r)>i) psa->SetPair(i,folded->GetPair(i,r),psa->GetNumberofStructures());
	}
}

void pseudoknotFold(pkHelix &pknot, RNA * st, RNA * psa, int energyPrune, datatable *data, int maxtracebacks, int percent, int window, int &numstructures, double P1, double P2, double Ss, double Si, string DMSFile, string SHAPEFile, double DSs, string DSHAPEFile, string doubleOffsetFile)
//Function that removes a helix (forces all nucleotides in that helix to be single stranded) from an RNA, folds the RNA, and then 
//	adds the helix back IF the addition of the helix to the folded RNA lowers the overall energy of the RNA below the energy of the
//	pseudoknot-free minimum free energy structure. 
//pknot is a pointer to the helix to be removed from the RNA (forced to be single stranded).
//st is a pointer to the RNA class object that will be folded.
//psa is a pointer to the RNA class object that will hold the pseudoknot free minimum free energy structure and 
//	all folded RNA structures that contain a pseudoknot.  This is equivalent to pseudoStructAggregate.
//energyPrune is the total energy of the minimum free energy structure
//data is a pointer to a datatable object (defined by DHM, needed for the call to Dynamic)
//maxtracebacks defines number of sub-optimal structures to be included in the dynamic function
//percent defines how close in energy sub-optimal structures can be
//numstructures indicates how many RNAs are in psa

{
    
	pknot.removeHelixFromStructure(st->GetStructure());	   //force all nts in helix to be single stranded


	while (st->GetStructureNumber()>0) st->GetStructure()->RemoveLastStructure();	//Remove all previous structures from st
	dynamic(st->GetStructure(), data, maxtracebacks, percent, window);	  //fold the rest of the structure

    // remove single-stranded constraints imposed by removeHelixFromStructure
    pknot.reverseRemoveHelix(st->GetStructure());

	//determine if the energy of structure is less than that of the pseudoknot-free MFE	
	for(int r=1; r<=st->GetStructureNumber(); r++){
		pknot.addHelixToStructure(st->GetStructure(),r);	//add the helix back to the folded structure
        
        //Store a copy of base pairings
        vector<int> rnacopy2= RNAcopy(st, r);
        FillMismatch(st,r);//fill single and tandem mismatches
        RemoveIsolatedPairs(st,r);//remove isloated pairs
        


		//If there is a pseudoknot:
        if(st->ContainsPseudoknot(r)){
			
            RestoreStructure(st,rnacopy2,r);//restore the broken pairs back
			//	if((st->energy[r]+pknot->getEnergy())<(energyPrune)){ 
			addtoAggregate(st->GetStructure(), psa->GetStructure(), r);			//add the entire structure to the final list of pseudoknotted structures
			//		}		
		
			
		}
	}
}

//Print the file with helices
void printhelixListtoFile(vector<pkHelix> pkhelixList){
	//Function that prints the list of 40 helices to a text file. Useful for debugging how the list is built
	ofstream checkfile;
	checkfile.open("HelixEnergyFile.txt");
	checkfile<<"helix energy  helix_length	  helix nucleotides"<<endl;

	for(int i=0;i<pkhelixList.size();i++){
		checkfile<<setw(5)<<setiosflags(ios::right)<<pkhelixList[i].getEnergy()<<setw(15)<<pkhelixList[i].getSize()/2<<"		   ";
		for (int j=0;j<pkhelixList[i].getSize();j++){
			checkfile<<pkhelixList[i].getelement(j)<<" ";
		}
		checkfile<<endl;
	}
}


void printPseudoknotList(RNA* pseudoStructAggregate){
	//Function that checks through each of the final folded structures for the presence of pseuodoknots.  Checks by iterating through
	//	all helices looking for the condition that satisfies i<i'<j<j', where i is paired to j, and i' is paired to j'.
	//This algorithm is not exhaustive in that a helix might have interleaved basepairs with more than one helix.  Therefore,
	//	manual inspection of any structure which contains interleaved basepairs will have to be examined.
	//pseudoStructAggregate is the final list of structures, each structure is individually checked.
	ofstream cck;
	cck.open("PseudoknotOutput.txt");
	bool pseudoknot = false;
	for(int l=1; l<=pseudoStructAggregate->GetStructureNumber(); l++){
		pseudoknot = false;
		for( int r=1;r<=pseudoStructAggregate->GetSequenceLength();r++){
			if(pseudoStructAggregate->GetPair(r,l)!=0 && pseudoknot == false){
				for(int x=1;x<pseudoStructAggregate->GetSequenceLength();x++){
					if(pseudoStructAggregate->GetPair(x,l)!=0 && pseudoknot == false){
						if( ( (r< x) &&(x < pseudoStructAggregate->GetPair(r,l)) &&( pseudoStructAggregate->GetPair(r,l) < pseudoStructAggregate->GetPair(x,l)) )||
						    ( (x< r) &&(r < pseudoStructAggregate->GetPair(x,l))&& (pseudoStructAggregate->GetPair(x,l)<pseudoStructAggregate->GetPair(r,l)) )	){
							pseudoknot = true;
							//cout<<"In structure "<<l<<", pseudoknot between helix starting with basepair "<<r<<"-"<<pseudoStructAggregate->basepr[l][r]<<endl<<
							//	" and helix starting with basepair "<<x<<"-"<<pseudoStructAggregate->basepr[l][x]<<endl;
							//cck<<"In structure "<<l<<", pseudoknot between helix starting with basepair "<<r<<"-"<<pseudoStructAggregate->basepr[l][r]<<endl<<
							//	" and helix starting with basepair "<<x<<"-"<<pseudoStructAggregate->basepr[l][x]<<endl;
						}			
					}
				}
			}
		}
	}
}


void convertToHundredHelices(vector<pkHelix> &pkHelixList, int finallistSize){
	//Function that "trims" the list of possible pseudoknots to the 100 (or 'finallistSize') most energetically stable (lowest energy).
	// Initial list can be a few to several thousand possibilities depending on the size of the RNA.
	//pkHelixList is a pointer to the list of helices, a list of pkHelix objects.
	//pkHelisListSize is the current size of the list of helices
		
	//if number of possible pseudoknots is less than this, return
	//QUESTION: Should I changed the '<' to '<='
	if(pkHelixList.size()<finallistSize){
#ifdef DEBUG_MODE
		printhelixListtoFile(pkHelixList);
#endif
		return;
	}

	//start at position finallistSize (default = 100)
	for (int i=finallistSize; i<pkHelixList.size(); i++){
		//examine the energy of the first finallistSize (defalut = 100) helices compared to each helix after finallistSize
		for (int j=0;j<finallistSize; j++){			
			//if the energy of the ith helix is less than the energy of the jth helix, swap 
			if ((pkHelixList[i].getEnergy())<(pkHelixList[j].getEnergy())){
				pkHelix temppkhelix (pkHelixList[i]);
				pkHelixList[i].setEqual(pkHelixList[j]);
				pkHelixList[j].setEqual(temppkhelix);
			}
		}
	}
	pkHelixList.erase (pkHelixList.begin()+finallistSize,pkHelixList.end());
	//pkHelixList.size()=finallistSize;
#ifdef DEBUG_MODE
	printhelixListtoFile(pkHelixList);
#endif
}

void findhelix(short **arrayPointer, vector<pkHelix>  &pkHelixList, int x, int y, datatable *data, structure *ct){
	//Function to find helices in the 2-d array generated fill.	 
	//arrayPointer is the pointer to the 2-d array.
	//pkHelixList is a pointer to the list of helices, a list of pkHelix objects.
	//pkHelisListSize is the current size of the list of helices
	//x and y are coordinates in the 2-d array
	//data is a pointer to a data object defined by DHM, passed to the pkHelix constructor and ultimately used to calculate the energy of the new helix
	//ct is a pointer to a structure object defined by DHM, passed to the pkHelix constructor and ultimately used to calculate the energy of the new helix
	//This function is a little awkward.  Takes a cell (x,y) identified to hold a base pair in the defineHelixList function and looks to see if there are more
	//	basepairs on the diagonal with the same energy, indicating a helix.	 Only helices greater than 2nts and less than 10 nts are considered to be possible pseudoknots.
	short helixLength = 1;
	
	//make sure that the next cell on the diagonal exists first
	//by subtracting the Y value from the X value...if at the end, the result will be less than 2
	//return without changing anything.
	//check to see if the next element on the diagonal is equal to calling cell
	while(	(((x-y)>2)||((x-y)<-2)) && (arrayPointer[x][y] == arrayPointer[x-1][y+1] ) ){
		arrayPointer[x][y] = 28001;	  //reset the current cell so it won't be found again 
		x--;						  //decrement x to be tested next while loop	
		y++;						  //increment y to be tested next while loop	
		helixLength++;
	}
	//make sure that the calling helix is  greater than two nts and less than 10 nts
	if (helixLength>2){
		arrayPointer[x][y] = 28001; //set the last element of the helix to 28001
		//create an array of shorts to hold numbers involved in helix
		//use standard array notation 0 to n-1
		short *temphelix;
		//temphelix = new short [2*helixLength-1];//DELETED
		temphelix = new short [2*helixLength];//ADDED
		//fill the temporary array of shorts backward since variables hold last elements of helix
		for (int i =2*helixLength-1; i>0;i--){//ADDED
			temphelix[i]=x;
			temphelix[i-1]=y;
			x++;
			y--;
			i--; //filling two elements at a time, so decrement i again
		}

		//create a temporary pkHelix object to hold new helix
		//pkHelix temppkhelix(data, ct, (2*helixLength), temphelix);
		//assign temporary pkHelix object to next empty position in list
		//QUESTION: Memory Leak Here:		
		//pkHelixList.push_back(temppkhelix);
		pkHelixList.push_back(pkHelix(data, ct, (2*helixLength), temphelix));
		delete[] temphelix;	
	}
}

void defineHelixList(short **arrayPointer, vector<pkHelix> &pkHelixList, int & maxsize, int lowvalue, datatable *data, structure *ct){
	//Function that iterates through the 2-d array generated by the fill function, and looks for structures with total energy within 25% of the pseudoknot-free minimum free energy structure.
	//arrayPointer is the pointer to the 2-d array.
	//pkHelixList is a pointer to the list of helices, a list of pkHelix objects.
	//pkHelisListSize is the current size of the list of helices
	//data is a pointer to a data object defined by DHM, passed to findhelix function and then the pkHelix constructor and ultimately used to calculate the energy of the new helix
	//ct is a pointer to a structure object defined by DHM, passed to findhelix function and then the pkHelix constructor and ultimately used to calculate the energy of the new helix
	//max size is the maximum size of the list of helices -	 resized as needed
	//lowvalue is the energy of the pseudoknot-free minimum free energy structure.	This was calculated in fill and stored when the 2-d array was built
//<<<<<<< ShapeKnots.cpp
//	int upperEnergyValue;
//	if(ct->GetSequenceLength()>100)
//=======
	double upperEnergyValue;
	if(ct->GetSequenceLength()>100)
//>>>>>>> 1.11
		upperEnergyValue =( lowvalue-(lowvalue*0.25));		//defines upper end of energy range when looking through helices in other structures
	else
		upperEnergyValue =( lowvalue-(lowvalue*0.999));
	//set upper energy range to be withing 25% of pseudoknot-free MFE structure
	//now iterate through the array looking for matrix positions that might be good starting points for the pseudoknot
	for (int i=1;i<=ct->GetSequenceLength(); i++){
		for(int j=i;j<=ct->GetSequenceLength(); j++){
			//checks the following in order:
			//are the bases less than 600 nt apart?
			//is the energy greater than pseudoknot-free MFE, this helps prevent false positives, 
			//position greate and make sure the position is not right at the diagonal 
			if((j-i<=600)&&(arrayPointer[j][i]<=upperEnergyValue)&& (arrayPointer[j][i]>lowvalue) && ( ((i-j)>2) ||( (i-j)<-2) ) ){
				//once a position is found, pass the findhelix function the 2-d array, the 3-d array to be populated and the positions of the cell
				//REMEMBER using j to iterate along the X-axis and i to iterate along Y-axis
				
				findhelix(arrayPointer, pkHelixList, j, i, data, ct);
			}
		}
	}		
}

void convertStructureToHelixList(structure *st, datatable *data, vector<pkHelix> &structureHelices){
	//Function that converts the pseudoknot-free MFE (OR ANY STRUCTURE) to a list of helices.  This function uses a pointer to a list of pkHelix objects 
	//	to hold the list.  This list will use this to compare MFE helices to proposed pseudoknotted structure to look for false positives
	//	Only looks at lowest energy structure of each traceback, will only report first structure anyway
	//st is a pointer to a structure object, holds all data associated with an RNA structure 
	//data is a pointer to datatable object 
	//structureHelices is a pointer to the list of helices in the pk-free MFE
	//structureHelicesSize is the number of helices in the list
   
	//temporary object to hold new list
	int thlength = 0;
	short *temphelix;

	temphelix = new short [st->GetSequenceLength()];
	//look through the array that holds the list of basepairs in a structure, look for helices.
	for (int i=1; i<=st->GetSequenceLength(); i++){
		//look for a base pair and only pick it up once by only taking it if is greater than position
		if( st->GetPair(i)>i ){
			(temphelix[thlength++]=i); 
			(temphelix[thlength++]=st->GetPair(i));
						
			//statement that checks that helix is closed properly
			if(	 st->GetPair(i+1) != st->GetPair(i)-1 ){
				structureHelices.push_back(pkHelix(data,st,thlength,temphelix));
				thlength=0;
			}//end of end of helix - inside loop finding basepair >i
		}//end of find basepair >i
	}//end of for loop
	delete[] temphelix;
}//end of function

void comparepkHelixListToMFE(vector<pkHelix> &pkHelixList, structure *ct, vector<pkHelix> &MFEhelixList){
	//Function that removes any helices from the list of possible pseudoknots that will significantly repair the pseudknot-free MFE.
	//this is simply an empirically derived method to reduce the number of false positives in the prediction
	//pkHelixList is a pointer to the list of helices, a list of pkHelix objects. tHIS IS THE LIST OF POSSIBLE PSEUDOKNOTS
	//pkHelisListSize is the current size of the list of possible pseudoknots
	//ct is a pointer to a structure object, holds all data associated with an RNA structure (programmed by DHM)
	//data is a pointer to datatable object (programmed by DHM)
	//MFEhelixlist is a pointer to a list of helices, a list of pkHelix objects.  THIS IS THE LIST OF HELICES IN THE PK-FREE MFE.
	//MFEhelixListSize is the number of helices in the pk-free MFE
	bool same =false;
	//first iterate through the list of helices in the pk-free MFE and remove all helices that are less than 3 basepairs.
	//But keep all helices that have 1nt bulges
	//Again, just an empirical rule that allows repairing of short helices
	for(int y=0; y<MFEhelixList.size(); y++){
		if((MFEhelixList[y].getSize()<=6)&&((ct->GetPair((MFEhelixList[y].get5primeLast())+2))+1!=(ct->GetPair((MFEhelixList[y].get5primeLast()))))
		   &&((ct->GetPair((MFEhelixList[y].get3primeLast())+2))+1!=(ct->GetPair((MFEhelixList[y].get3primeLast()))))  
		   &&((ct->GetPair((MFEhelixList[y].get5primeFirst())-2))-1!=(ct->GetPair((MFEhelixList[y].get5primeFirst()))))
		   &&((ct->GetPair((MFEhelixList[y].get3primeFirst())-2))-1!=(ct->GetPair((MFEhelixList[y].get5primeFirst())))) ){
			for(int z=0; z<MFEhelixList[y].getSize();z+=2){
				ct->RemovePair(MFEhelixList[y].getelement(z));
				//ct->RemovePair(MFEhelixList[y].getelement(z+1));//RemovePair function takes care of this automatically.
			}
		}
	}
	int removefirst= 0;
	int removesecond=0;//indicates the number of times nucleotides in each half of possible pseudoknotted helix overlaps MFE helices
	//If more than half of nucleotides in both halves of possible pseudoknotted helix are involved in helices in pk-free MFE, 
	//the helix significantly repairs MFE, thereforeremove helix from list
	
	//iterate through the list of possible pseudoknots
	for (int i=0; i<pkHelixList.size(); i++){
		//compare against all helices in pseudoknot free MFE structure
		for(int j=0; j<pkHelixList[i].getSize(); j+=2){
			if(ct->GetPair((pkHelixList[i].getelement(j)))!=0){
				removefirst++;
				
			}
			if(ct->GetPair(pkHelixList[i].getelement(j+1))!=0){
				removesecond++;
				
			}
		}//end of j
		//check if more than half of nucleotides in both 5' and 3' halves of possible pseudoknot are basepaired in pk-free MFE
		//note that checking the possible pseudoknotted helix, not the helix in the pk-free MFE
		if((removefirst>pkHelixList[i].getSize()/4)&&(removesecond>pkHelixList[i].getSize()/4)){
			same=false;
			
			for(int f=0; f<MFEhelixList.size();f++){
				if(pkHelixList[i].isEqual(MFEhelixList[f]))
					same=true;
			}
			if(same==false){
				pkHelixList.erase(pkHelixList.begin()+i);
				i--;
			}
		}
		removefirst=removesecond=0;
		
	}//end of i

}


//Check for duplicate structures in a structure class and remove them.
void checkForDuplicates(structure* &pseudoStructAggregate){
	//function that checks the final folded list of structures to ensure that no structure is duplicated
	//pseudostructAggregate is a pointer to the final list of folded structures.
	bool same;
	for(int l=1; l<pseudoStructAggregate->GetNumberofStructures(); l++){

		same=true;
		for(int r=l+1;r<pseudoStructAggregate->GetNumberofStructures();r++){
			//only check structures which have the same energy
			//if (pseudoStructAggregate->energy[l]==pseudoStructAggregate->energy[r]){
			same=true;
			for (int h=1; h<=pseudoStructAggregate->GetSequenceLength();h++){		
				//if the basepair is not equal at any point, then structures are not equivalent
				//can stop comparing these two structures and break back to r
				if(pseudoStructAggregate->GetPair(h,l)!=pseudoStructAggregate->GetPair(h,r)){
					same = false;						
					break;
				}//end of if
			
			}//end of j
			//will either come out of j with same equal to true or false
			//if same equals true want to exit r and go back to i same is equal to false if the break occurs
			//if same is still true after going all the way through any of them, want to skip this i position and go to the next
			if(same==true){
			
				//Use INIFINTE_ENERGY as a flag for duplicates
				pseudoStructAggregate->SetEnergy(l,INFINITE_ENERGY);
			
				break;	//break out of r
			}
			//	}//end of if energy is the same
			
		}	
		
	}//end of R	


	//Use the flag of high folding free energy to put all duplicates at the end
	//of the list of structures by sorting and then remove them.
	pseudoStructAggregate->sort();
	for (int i=pseudoStructAggregate->GetNumberofStructures();i>0;--i) {
		if (pseudoStructAggregate->GetEnergy(i)==INFINITE_ENERGY) {
			pseudoStructAggregate->RemoveLastStructure();
		}
		else {
			break;
		}
	}
}		
void copyStructure(structure * ct, structure* &pseudoStructAggregate, string DMSFile, string SHAPEFile){
	//Function that initializes pseudoStructAggregate to have some of the same attributes as the original structure ct
	//ct is allocated in main, and reads in the initial data
	//pseudoStructAggregate is a pointer to a structure that holds the final list of folded structures

	if (ct->DistanceLimited()) pseudoStructAggregate->SetPairingDistance(ct->GetPairingDistanceLimit());		//only allow base pairs closer than 600
	
    // copy constraints; NOTE that NMR and microArray constraints are not copied!!!
    // TONY 
    for (int i=0; i<ct->GetNumberofDoubles(); i++)
        pseudoStructAggregate->AddDouble( ct->GetDouble(i) );
    for (int i=0; i<ct->GetNumberofSingles(); i++)
        pseudoStructAggregate->AddSingle( ct->GetSingle(i) );
    for (int i=0;i<ct->GetNumberofModified(); i++)
        pseudoStructAggregate->AddModified( ct->GetModified(i) );
    for (int i=0;i<ct->GetNumberofPairs(); i++)
        pseudoStructAggregate->AddPair( ct->GetPair5(i), ct->GetPair3(i) );
    for (int i=0;i<ct->GetNumberofGU(); i++)
        pseudoStructAggregate->AddGUPair( ct->GetGUpair(i) );
    for (int i=0;i<ct->GetNumberofForbiddenPairs();i++)
        pseudoStructAggregate->AddForbiddenPair( ct->GetForbiddenPair5(i), ct->GetForbiddenPair3(i) );



	pseudoStructAggregate->allocate(ct->GetSequenceLength()); //indicates the size of the new structure
	if(!SHAPEFile.empty()||!DMSFile.empty()){
		pseudoStructAggregate->SHAPE = new double [2*ct->GetSequenceLength()+1];
		pseudoStructAggregate->SHAPEss = new double [2*ct->GetSequenceLength()+1];
		
		pseudoStructAggregate->SHAPEss_region = new short int *[pseudoStructAggregate->GetSequenceLength() + 1];
		for (int i = 1; i <= pseudoStructAggregate->GetSequenceLength(); i++) {
			pseudoStructAggregate->SHAPEss_region[i] = new short int [i];
			for (int j=0;j<i;++j) pseudoStructAggregate->SHAPEss_region[i][j]=0;
		}

	}
	
	if(!SHAPEFile.empty()||!DMSFile.empty()) pseudoStructAggregate->shaped=true;
	else pseudoStructAggregate->shaped=false;

    if (ct->EX!=NULL) {
        pseudoStructAggregate->experimentalPairBonusExists = true;
        pseudoStructAggregate->EX = new double *[ 2*ct->GetSequenceLength()+1 ];
        for (int i = 0; i < 2*ct->GetSequenceLength()+1; i++) {
            pseudoStructAggregate->EX[i] = new double [ 2*ct->GetSequenceLength()+1 ];
            for (int j=0; j < 2*ct->GetSequenceLength()+1; j++){
                pseudoStructAggregate->EX[i][j] = ct->EX[i][j];
            }
        }
    }

	//set the identifier for the structure 
	for (int q=1;q<=ct->GetSequenceLength();q++){	
		pseudoStructAggregate->numseq[q] = ct->numseq[q];	//set to equivalent position in original structure
		pseudoStructAggregate->hnumber[q]= ct->hnumber[q];
		pseudoStructAggregate->nucs[q] = ct->nucs[q];
		if(!SHAPEFile.empty()||!DMSFile.empty()){
			pseudoStructAggregate->SHAPE[q]=ct->SHAPE[q];		
			pseudoStructAggregate->SHAPEss[q]=ct->SHAPEss[q];
		}
	}
	if(!SHAPEFile.empty()||!DMSFile.empty()){
		pseudoStructAggregate->SHAPEslope=ct->SHAPEslope;
		pseudoStructAggregate->SHAPEintercept=ct->SHAPEintercept;
	}
}


//	Function gets the names of data files to open
void getdat(char *fspec, char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA, bool isEnthalpy)

{

	//copy the path to each name
  strcpy(fspec,datapath);
  strcpy (loop,datapath);
  strcpy (stackf,datapath);
  strcpy (tstackh,datapath);
  strcpy (tstacki,datapath);
  strcpy (tloop,datapath);
  strcpy (miscloop,datapath);
  strcpy (danglef,datapath);
  strcpy (int22,datapath);
  strcpy (int21,datapath);
  strcpy (triloop,datapath);
  strcpy (coax,datapath);
  strcpy (tstackcoax,datapath);
  strcpy (coaxstack,datapath);
  strcpy (tstack,datapath);
  strcpy (tstackm,datapath);
  strcpy (int11,datapath);
  strcpy (hexaloop,datapath);
  strcpy (tstacki23,datapath);
  strcpy (tstacki1n,datapath);
  


  if( !isRNA) {
	  //these are dna parameters and so they need to start with "dna"
	  strcat (fspec,"dna.");
	strcat (loop,"dna.");
	strcat (stackf,"dna.");
	  strcat (tstackh,"dna.");
	  strcat (tstacki,"dna.");
	  strcat (tloop,"dna.");
	  strcat (miscloop,"dna.");
	  strcat (danglef,"dna.");
	  strcat (int22,"dna.");
	  strcat (int21,"dna.");
	  strcat (triloop,"dna.");
	  strcat (coax,"dna.");
	  strcat (tstackcoax,"dna.");
	  strcat (coaxstack,"dna.");
	  strcat (tstack,"dna.");
	  strcat (tstackm,"dna.");
	  strcat (int11,"dna.");
	  strcat (hexaloop,"dna.");
	  strcat (tstacki23,"dna.");
	  strcat (tstacki1n,"dna.");  
	 
  }
  else { //Right now, this means that the files are for RNA
	  //these are dna parameters and so they need to start with "dna"
	  strcat (fspec,"rna.");
	strcat (loop,"rna.");
	strcat (stackf,"rna.");
	  strcat (tstackh,"rna.");
	  strcat (tstacki,"rna.");
	  strcat (tloop,"rna.");
	  strcat (miscloop,"rna.");
	  strcat (danglef,"rna.");
	  strcat (int22,"rna.");
	  strcat (int21,"rna.");
	  strcat (triloop,"rna.");
	  strcat (coax,"rna.");
	  strcat (tstackcoax,"rna.");
	  strcat (coaxstack,"rna.");
	  strcat (tstack,"rna.");
	  strcat (tstackm,"rna.");
	  strcat (int11,"rna.");
	  strcat (hexaloop,"rna.");
	  strcat (tstacki23,"rna.");
	  strcat (tstacki1n,"rna.");  
	 
  }
	 
	//append the actual file name
  strcat(fspec, "specification.dat");//include full extension here, it is the same for DG and DH
  strcat (loop,"loop.");
  strcat (stackf,"stack.");
  strcat (tstackh,"tstackh.");
  strcat (tstacki,"tstacki.");
  strcat (tloop,"tloop.");
  strcat (miscloop,"miscloop.");
  strcat (danglef,"dangle.");
  strcat (int22,"int22.");
  strcat (int21,"int21.");
  strcat (triloop,"triloop.");
  strcat (coax,"coaxial.");
  strcat (tstackcoax,"tstackcoax.");
  strcat (coaxstack,"coaxstack.");
  strcat (tstack,"tstack.");
  strcat (tstackm,"tstackm.");
  strcat (int11,"int11.");
  strcat (hexaloop,"hexaloop.");
  strcat (tstacki23,"tstacki23.");
  strcat (tstacki1n,"tstacki1n.");
  
  if (isEnthalpy) {
	  //thse are enthalpy partameters so they need to end in .dh
	strcat (loop,"dh");
	strcat (stackf,"dh");
	  strcat (tstackh,"dh");
	  strcat (tstacki,"dh");
	  strcat (tloop,"dh");
	  strcat (miscloop,"dh");
	  strcat (danglef,"dh");
	  strcat (int22,"dh");
	  strcat (int21,"dh");
	  strcat (triloop,"dh");
	  strcat (coax,"dh");
	  strcat (tstackcoax,"dh");
	  strcat (coaxstack,"dh");
	  strcat (tstack,"dh");
	  strcat (tstackm,"dh");
	  strcat (int11,"dh");
	  strcat (hexaloop,"dh");
	  strcat (tstacki23,"dh");
	  strcat (tstacki1n,"dh");  
  }
  else {
	  //these are free energy parameters and the files end in .dat
	strcat (loop,"dg");
	strcat (stackf,"dg");
	  strcat (tstackh,"dg");
	  strcat (tstacki,"dg");
	  strcat (tloop,"dg");
	  strcat (miscloop,"dg");
	  strcat (danglef,"dg");
	  strcat (int22,"dg");
	  strcat (int21,"dg");
	  strcat (triloop,"dg");
	  strcat (coax,"dg");
	  strcat (tstackcoax,"dg");
	  strcat (coaxstack,"dg");
	  strcat (tstack,"dg");
	  strcat (tstackm,"dg");
	  strcat (int11,"dg");
	  strcat (hexaloop,"dg");
	  strcat (tstacki23,"dg");
	  strcat (tstacki1n,"dg");  


  }
}

void pseudoknot(RNA* rnaCT, datatable *data, int maxStructures, int percent, int windowSize, string  ctFile, double P1, double P2, double Ss, double Si, string DMSFile, string SHAPEFile, double DSs, string DSHAPEFile, string doubleOffsetFile, int OutPercent, int OutWindowSize, int OutMaxStructures, bool ifWindowOptions, int finallistSize)
//Function that drives the rest of the code.  This function is called in main after the sequence file, shape file and shape parameters are
//	are read.
//rnaCT is a pointer to a RNA class object, holds all data associated with an RNA structure (programmed by DHM)
//data is a pointer to datatable object (programmed by DHM)
//maxtracebacks defines number of sub-optimal structures to be included in the dynamic function
//percent defines how close in energy sub-optimal structures can be
//ctFile is the name of the file where results are printed
//SHAPEslope and SHAPEintercept are the SHAPE parameters
{
	integersize *w5,*w3;
	bool *lfce,*mod;
	//int crit = 0;//DELETED because unised
	int i = 0;
	int j = 0;
	int vmin;
	DynProgArray<integersize> *w2,*wmb2;
	DynProgArray<integersize> *w, *v, *wmb;
	forceclass *fce;
	w = new DynProgArray<integersize>(rnaCT->GetSequenceLength());
	v = new DynProgArray<integersize>(rnaCT->GetSequenceLength());
	wmb = new DynProgArray<integersize>(rnaCT->GetSequenceLength());
	fce = new forceclass(rnaCT->GetSequenceLength());

	lfce = new bool [2*rnaCT->GetSequenceLength()+1];
	mod = new bool [2*rnaCT->GetSequenceLength()+1];

	for (i=0;i<=2*rnaCT->GetSequenceLength();i++) {
		lfce[i] = false;
		mod[i] = false;
	}

	w5 = new integersize [rnaCT->GetSequenceLength()+1];
	w3 = new integersize [rnaCT->GetSequenceLength()+2];

	for (i=0;i<=rnaCT->GetSequenceLength();i++) {
		w5[i] = 0;
		w3[i] = 0;

	}
	w3[rnaCT->GetSequenceLength()+1] = 0;

	if (rnaCT->GetStructure()->intermolecular) {
		w2 = new DynProgArray<integersize>(rnaCT->GetSequenceLength());
		wmb2 = new DynProgArray<integersize>(rnaCT->GetSequenceLength());
	}
	else {
		w2 = NULL;
		wmb2 = NULL;
	}
	

	force(rnaCT->GetStructure(),fce,lfce);
	vmin=DYNALIGN_INFINITY;
#ifdef OUTPUT_TO_SCREEN
	cout << "\t\t\t\tDONE\nGenerating energy dotplot..." << flush;
#endif
	//perform the fill steps:(i.e. fill arrays v and w.)
	fill(rnaCT->GetStructure(), *v, *w, *wmb, *fce, vmin,lfce, mod,w5, w3, false, data, w2, wmb2);
#ifdef OUTPUT_TO_SCREEN
	cout << "\t\t\tDONE\n" << flush;
#endif
	//create the 2d array using the first dimension as an array of pointers to short
	//the second dimension is then an array of shorts which grows with each allocation +1 to create an unequal or triangular shape
	
	short **arrayPointer;	 //allocate the pointer to the 2d array for the dot plot
	arrayPointer = new short *[rnaCT->GetSequenceLength()+1];
	for (int t =1; t<=rnaCT->GetSequenceLength();t++)
		arrayPointer[t]=new short [t+1];

	/*This algorithm calculates the lowest dG structure containing the i-j pair such that i is less than j
	  is based on calculating the energy of the helix internal to the pair (i,j) (v.f(i,j)
	  and the energy for the part of the helix external to the pair (i,j) as v.f(j,i+ct->GetSequenceLength())
	  this function indexes the first argument as the lower
	  therefore iterate through the 2d array and fill it with this calculation:
	*/
	int lowvalue = 0;	//holds energy of the pseudoknot free MFE
	//i MUST START AT ONE 
	for(i=1; i<=(rnaCT->GetSequenceLength()); i++){
		
		for(j=i;j<=rnaCT->GetSequenceLength();j++){
			arrayPointer[j][i] = (v->f(i,j)+ (v->f(j,i+rnaCT->GetSequenceLength())) );
			if (arrayPointer[j][i]<lowvalue)
				lowvalue=arrayPointer[j][i];
		}
			
	}

	//set up the list that will hold all helices that may be involved in a pseudoknot

	int maxsize=10;
	vector<pkHelix> pkhelixList;

	//this function creates a 2-d array that is a mirror of that created in fill function, and iterates through it
	//looking for starting points for helices
	defineHelixList(arrayPointer, pkhelixList, maxsize, lowvalue, data, rnaCT->GetStructure());
  
	for(j=1;j<=rnaCT->GetSequenceLength();j++) delete[] arrayPointer[j];
	delete []arrayPointer;
	//pkhelisList is now full WITH ALL POSSIBLE HELICES GREATER THAN 2NTS AND LESS THAN 10 NTS 
	//this list will be trimmed down to a manageble size once the pseudoknot-free MFE is folded

	//structure that will hold the final list of pseudoknotted helices
	//all pseudoknotted structures will be added to this structure
	//structure* pseudoStructAggregate= new structure;
    RNA* pseudoStructAggregateCT=new RNA;
    structure* pseudoStructAggregate=pseudoStructAggregateCT->GetStructure();
	copyStructure(rnaCT->GetStructure(), pseudoStructAggregate, DMSFile, SHAPEFile);
    

#ifdef OUTPUT_TO_SCREEN
	cout << "Folding the MFE structure..." << flush;
#endif
	//fold the pseudoknot free MFE structure
	dynamic(rnaCT->GetStructure(),data,maxStructures,percent,windowSize);
#ifdef OUTPUT_TO_SCREEN
	cout << "\t\t\tDONE\n" << flush;
#endif
    

	//will add all structures indicated by the user with maxStructures

	for (int h=1; h<=rnaCT->GetStructureNumber();h++)
		addtoAggregate(rnaCT->GetStructure(), pseudoStructAggregate,h);

	//Now want to "trim" the list of all possible helices to a manageble size
	//First will remove any helices that significantly repair the pseudoknot free MFE
	//only compare list of possible pseudoknots to pseudoknot-free MFE if RNA is long enough to generate more than 10 possibilities 
	if(pkhelixList.size()>10){
		//2-d list to store the pseudoknot free structure as a list of helices -  will be used later to test for psuedoknot false positives
		//pkHelix *MFEhelixList= new pkHelix[ct->GetSequenceLength()];
		vector<pkHelix> MFEhelixList(0);
		//		int MFEhelixListSize = 0;
		convertStructureToHelixList(rnaCT->GetStructure(), data, MFEhelixList);
		comparepkHelixListToMFE(pkhelixList, rnaCT->GetStructure(), MFEhelixList);
	}

	//Convert the list of helices to 100, this is final list of helices to be tested 
	//These helices will be tested for the ability to form pseudoknots by forcing all nucleotides in each 
	//to be single stranded.

	convertToHundredHelices(pkhelixList, finallistSize);

	//now test each helix for its ability to form a pseudoknot
	int numstructures=0;

	
	//iterate through the list of possible pseudoknots
	for (i=0;i<pkhelixList.size();i++){
#ifdef OUTPUT_TO_SCREEN
		cout << '\r' << "Folding modified structure "<<(i+1)<<" of "<<pkhelixList.size()<< flush;
#endif
		//fold the structure while forcing each nt in the current helix to be single stranded, and then add the helix back 
		pseudoknotFold(pkhelixList[i], rnaCT, pseudoStructAggregateCT, lowvalue, data, maxStructures, percent, windowSize, numstructures, P1, P2, Ss, Si, DMSFile, SHAPEFile, DSs, DSHAPEFile, doubleOffsetFile);

		

	}
#ifdef OUTPUT_TO_SCREEN
	cout << "\t\tDONE\nChecking for duplicate structures..." << flush;
#endif
	



	//check list of folded structures for duplicates
	checkForDuplicates(pseudoStructAggregate);




#ifdef OUTPUT_TO_SCREEN
	cout << "\t\tDONE\n" << flush;
#endif



    double pseud_param[16];//Array that holds pseudoknot penalty calculation constants.
    //Read pseudoknot parameters from the datapath
    ReadPseudoParam(pseud_param);

#ifdef DEBUG_MODE
 	int pseudnum=0;
#endif

	//SuboptimaLpseudoStructAggregate->ReadOffset(NULL,DSO);//read the offset data

	int numofstr = pseudoStructAggregateCT->GetStructureNumber();

	//For every structure in pseudoStructAggregate execute the following:
	for(int strnum=1;strnum<=pseudoStructAggregateCT->GetStructureNumber();++strnum){
#ifdef OUTPUT_TO_SCREEN
		if(strnum%100==0||strnum==numofstr) cout << "\rCalculating energy of structure " << strnum << " of " << numofstr << flush;
#endif

        vector<int> rnacopy2=RNAcopy(pseudoStructAggregateCT,strnum);
        FillMismatch(pseudoStructAggregateCT,strnum);
        RemoveIsolatedPairs(pseudoStructAggregateCT,strnum);
		//If there is a pseudoknot in rna:
		//if(rna->ContainsPseudoknot(1)){ 
        if(pseudoStructAggregateCT->ContainsPseudoknot(strnum)){
//<<<<<<< ShapeKnots.cpp
			//Copy the suboptimal structure number 'strnum' from 'pseudoStructAggregate' to SuboptimaLpseudoStructAggregate'
//			for(int i=0;i<pseudoStructAggregate->GetSequenceLength();++i){
//				SuboptimaLpseudoStructAggregate->SetPair(i,pseudoStructAggregate->GetPair(i,strnum)); 
//			}
//=======
//>>>>>>> 1.11

#ifdef DEBUG_MODE
			++pseudnum; 
#endif
            //Calculate pseudoknot energy penalty
			//double energyPenalty = pkEnergyCalc(rna, pseud_param, P1, P2);
			//std::stringstream ctfilename;
			//ctfilename << "pseudoStructAggregateCT_" << strnum << ".ct";
			//pseudoStructAggregateCT->WriteCt(ctfilename.str().c_str());
            double energyPenalty = pkEnergyCalc(pseudoStructAggregateCT, pseud_param, P1, P2, strnum);
            RestoreStructure(pseudoStructAggregateCT,rnacopy2,strnum);
			//Since energies stored in 'structure' class are in tenth precision, energyPenalty needs to be converted to tenth precision.
            
			pseudoStructAggregate->SetEnergy(strnum,pseudoStructAggregate->GetEnergy(strnum)+(energyPenalty*conversionfactor));
            
		}
	}
	

#ifdef DEBUG_MODE
	cerr << '\n'<< pseudnum << " structures contain pseudoknot \n";
	cerr << pseudoStructAggregateCT->GetStructureNumber()-pseudnum << " structures don't contain pseudoknot \n" ;
#endif

	



	//Sort the generated structures based on their energy
	pseudoStructAggregate->sort();
#ifdef OUTPUT_TO_SCREEN
	cout << "\tDONE\nSorting structures by their free energy..." << flush;
#endif

	//If the OutWindowSize option was not specified by the user on the command line
	if(!ifWindowOptions){
		//Scale OutWindowSize based on the length of the sequence
		if(pseudoStructAggregateCT->GetSequenceLength()>1200) OutWindowSize=20;
		else if(pseudoStructAggregateCT->GetSequenceLength()>800) OutWindowSize=15;
		else if(pseudoStructAggregateCT->GetSequenceLength()>500) OutWindowSize=11;
		else if(pseudoStructAggregateCT->GetSequenceLength()>300) OutWindowSize=7;
		else if(pseudoStructAggregateCT->GetSequenceLength()>120) OutWindowSize=5;
		else if(pseudoStructAggregateCT->GetSequenceLength()>50) OutWindowSize=3;
		else OutWindowSize=2;
	}

	//Filter the output to a more reasonable number of structures
	filter(pseudoStructAggregate, OutPercent, OutMaxStructures, OutWindowSize);

	//Store the path to output CT file
	pseudoStructAggregate->ctout(ctFile.c_str());

#ifdef OUTPUT_TO_SCREEN
	cout << "\tDONE\n" << flush;
	cout << "############################################\n";
	cout << "                    DONE\n";
	cout << "############################################\n" << flush;
#endif

#ifdef DEBUG_MODE
	//Output the pseudoknot list
	printPseudoknotList(pseudoStructAggregateCT);
#endif

	delete[] lfce;
	delete[] mod;
 
	delete w;
	delete v;
	delete wmb;
	delete fce;
 

	delete[] w5;
	delete[] w3;

	if (rnaCT->GetStructure()->intermolecular) {
		delete w2;
		delete wmb2;
	}
	delete pseudoStructAggregateCT;
}
