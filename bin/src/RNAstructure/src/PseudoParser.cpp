/*
 * Pseudoparser functions that are used by the pkEnergyCalc cunction to
 * calculate the pseudoenergy of the pseudoknot. Used in ShapeKnots.
 *
 * (c) 2013 Stanislav Bellaousov
 * Mathews Lab, University of Rochester Medical Center
 */

#include "PseudoParser.h"
#include <iostream>
#include <iomanip>
#include <sstream>

void FillMismatch(RNA *rna, int strnum)
{
	/////////////////////////////////////////////////////
	//FILLING IN PAIRS IN SINGLE OR TANDEM MISMATCHES
	/////////////////////////////////////////////////////

	int error;//Stores errorcode.
	int j;

	//Loop over over all nucleotides in the strucutre and stop eight nucleotides to the end
	//of the structure. Basepairs closer then eight nucleotides to the end of the structure
	//cannot be a part of single or tandem mismatch.
	for(int i=1;i<=rna->GetSequenceLength()-8;++i) {
		//Only move 5' to 3' and ignore unpaired nucleotides. 
		if(rna->GetPair(i,strnum)>i){
			//Setting 'j' equal to the pairing partner of 'i'.
			j=rna->GetPair(i,strnum);

			//Start single mismatch filling
			//If there is a single mismatch...
			if(rna->GetPair(i+1,strnum)==0&&rna->GetPair(rna->GetPair(i,strnum)-1,strnum)==0&&rna->GetPair(i+2,strnum)==j-2&&isCanonical(rna->GetNucleotide(i+1),rna->GetNucleotide(rna->GetPair(i,strnum)-1))){
                //...fill in the mismatch.
				error=rna->SpecifyPair(i+1,j-1,strnum);
				//If there was an error output it to the screen and exit.
				if (error!=0) {
					cerr << rna->GetErrorMessage(error);
					delete rna;
				}
			}//Finish single mismatch filling.

			//Start double mismatch filling.
			//If there is a tandem mismatch...
			if(rna->GetPair(rna->GetPair(i,strnum)-1,strnum)==0&&rna->GetPair(rna->GetPair(i,strnum)-2,strnum)==0&&rna->GetPair(i+1,strnum)==0&&rna->GetPair(i+2,strnum)==0&&rna->GetPair(i+3,strnum)==j-3&&isCanonical(rna->GetNucleotide(i+1),rna->GetNucleotide(rna->GetPair(i,strnum)-1))&&isCanonical(rna->GetNucleotide(i+2),rna->GetNucleotide(rna->GetPair(i,strnum)-2))){ 
				//...fill in the first mismatched pair.
				error = rna->SpecifyPair(i+1,j-1,strnum);
				//If there was an error output it to the screen and exit.
				if (error!=0) {
					cerr << rna->GetErrorMessage(error);
					delete rna;
				}
				//And fill in the second mismatched pair.
				error = rna->SpecifyPair(i+2,j-2,strnum);
				//If there was an error output it to the screen and exit.
				if (error!=0) {
					cerr << rna->GetErrorMessage(error);
					delete rna;
				}
			}//Finish double mismatch filling.
		}
	}//Close "FILLING IN PAIRS IN SINGLE OR TANDEM MISMATCHES".


}


void RemoveIsolatedPairs(RNA *rna,int strnum)
{
	/////////////////////////////////////////////////////
	//FINDING AND REMOVING ISOLATED PAIRS
	/////////////////////////////////////////////////////  

	int pairs;//Helix length counter for removing isolated pairs. 
	int j;
	int error;//Stores errorcode.

	//For every nucleotide in the structure.
	for(int i=1;i<=rna->GetSequenceLength();++i) {
		//Execute the following 5' to 3' ignoring unpaired nucleotides. 
		if(rna->GetPair(i,strnum)>i){
			//Start with the helix counter 'pairs' equal to 1.	
			pairs=1;	  
			//Set 'j' equal to the pairing partner of 'i'.
			j=rna->GetPair(i,strnum);
	  
			//FINDING isolated pairs.
			//If there is a stacked helix (including slippage) enter the 'while' loop. 
			while (rna->GetPair(i+1,strnum)==j-1||rna->GetPair(i+2,strnum)==j-1||rna->GetPair(i+1,strnum)==j-2){
				//In the case of non-slipped stacked helix...
				if (rna->GetPair(i+1,strnum)==j-1){
					++i;//...move 'i' counter one nucleotide up,
					--j;//...move 'j' counter one nucleotide down,
					++pairs;//... and add ONE to the helix counter. 
				}
				//In the case of stacked helix slipped on 5' end...
				else if (rna->GetPair(i+2,strnum)==j-1){
					//...first check for an isolated pair between the slipped nucleotides.
					if(rna->GetPair(i+1,strnum)!=0){
						//If there is one, remove it.
						rna->RemoveBasePair(i+1,strnum);
					}
					i=i+2;//Then move 'i' counter two nucleotides up,...
					--j;//...move 'j' counter one nucleotide down,...
					++pairs;//... and add ONE to the helix counter.
				}
				//In the case of stacked helix slipped on 3' end... 
				else {
					//...first check for an isolated pair between the slipped nucleotides.
					if(rna->GetPair(j-1,strnum)!=0){
						//If there is one, remove it.
						rna->RemoveBasePair(j-1,strnum);
					}
					++i;//Then move 'i' counter one nucleotides up,...
					j=j-2;//...move 'j' counter two nucleotide down,...
					++pairs;//... and add ONE to the helix counter.
				}
			}//Finish loop "FINDING isolated pairs".
	  
			//REMOVING isopated pairs.
			//If helix length counter 'pairs' is equal to 1, it means that there was an isolated
			//pair is encountered. 

			if (pairs==1){
				//Remove this isolated pair
				error = rna->RemoveBasePair(i,strnum);
				//If encountered error, output error to the screen.
				if (error!=0) {
					cerr << rna->GetErrorMessage(error);
					delete rna;
				}
			}//Finish loop "REMOVING isolated pairs".
		}
	}//Close "FINDING AND REMOVING ISOLATED PAIRS".
}


vector<int> RNAcopy(RNA *rna, int strnum)
{
	/////////////////////////////////////////////////////
	//CREATING rnacopy
	/////////////////////////////////////////////////////  
	vector <int> rnacopy;//Vector that holds a copy of the inputted structure.

	//Store copy of the structure that is in rna class after the preprocessing (filled single and tandem mismatches
	//and removed isolated pairs) in a vector 'rnacopy'.
	//Initialize 'rnacopy':
	rnacopy = vector<int>(rna->GetSequenceLength()+1,0);

	//Copy all pairs from the RNA class 'rna' to 'rnacopy.
	for(int i=1;i<=rna->GetSequenceLength();++i) {
		rnacopy[i]=rna->GetPair(i,strnum);
	}
	return rnacopy;
}


void BreakPseudoknot(RNA *rna, int strnum)
{
	/////////////////////////////////////////////////////
	//BREAKING PSEUDOKNOTS .
	/////////////////////////////////////////////////////  
    int tempenergy= rna->GetStructure()->GetEnergy(strnum);


	//Then using 'BreakPseudoknot' function break all pseudoknots that were in the inputed structure which is stored in the RNA class 'rna'.
	int error=rna->BreakPseudoknot(false, strnum);

	//Check for errors.
	if (error!=0){
		//There was an error.
		cerr << rna->GetErrorMessage(error)<< '\n';
		delete rna;
	}
    rna->GetStructure()->SetEnergy(strnum,tempenergy);
}


vector<int> PairsBroken(RNA *rna, vector<int> &rnacopy, int strnum)
{
	/////////////////////////////////////////////////////
	//RETURNING VECTOR WITH BROKEN PAIRS 'pairsBroken'.
	/////////////////////////////////////////////////////  

	vector<int> pairsBroken=vector<int>(rna->GetSequenceLength()+1,0);
	
	
	//First assign values in 'pairsBroken' to be same as pairs in the rna structure of interest.
	for (int i=1;i<=rna->GetSequenceLength();++i) pairsBroken[i]=rna->GetPair(i,strnum);
	
    //Break the pseudoknot
    BreakPseudoknot(rna, strnum);
    
	//Then change the 'pairsBroken' values to '0' for every case when the value in 'pairsBroken' 
	//matches the value in rna class (which now only holds nonpseudoknotted pairs). This will leave 'pairsBroken' 
	//with non zero values only for nucleotides that were broken after BreakPseudoknot was applied to rna class.
	for(int i=1; i<=rna->GetSequenceLength(); ++i){
		if(pairsBroken[i]==rna->GetPair(i,strnum)) pairsBroken[i]=0;
	}//Now pairsBroken only holds values for nucleotides that were broken when the pseudoknot was broken.

    return pairsBroken;
}

vector<Pseudoknot> PseudoParser(RNA *rna, vector<int> &pairsBroken, int strnum)
{
	/////////////////////////////////////////////////////
	//CREATING LIST OF PSEUDOKNOTS
	/////////////////////////////////////////////////////  

  	int i, k, thepk, itnum;
	unsigned int pknum;
	bool newpk=true;
	//Declare vector of pseudoknots 'pk' that holds a class 'pseudoknot' of vectors of classes.
	vector <Pseudoknot> pk;
	
	//Start loop over all nicleotides that were broken (that are in pairsBroken).
	for(i=1;i<=rna->GetSequenceLength();++i){
		//Move 5' to 3' direction skipping unpaired nucleotides.
		if(pairsBroken[i]!=0&&i<pairsBroken[i]){
			//Start loop over all pseudoknots 'pk'.
			for(pknum=0;pknum<pk.size();++pknum){
				//Start loop over all intact pairs in the 'pknum' instance of pseudoknot 'pk'.
				for(itnum=0;itnum<pk[pknum].GetNumberIntact();++itnum){
					//Start loop over all pairs within the boundaried formed by the helix i-GetPair(i).
					for (k=i+1;k<pairsBroken[i];++k){
						//Look for pairs that form a pseudoknot with pair i-GetPair(i).
						if(rna->GetPair(k,strnum)!=0&&(rna->GetPair(k,strnum)<i||rna->GetPair(k,strnum)>pairsBroken[i])){
							//If at least one of the pairs that form a pseudoknot with pair i has already been stored
							//then the new instance of a pseudoknot doesn't have to be created, and 'newpk' is set to 
							//'false', and the instance of a pseudoknot where this pair has been stored is stored in 'thepk'.
							if(k==pk[pknum].GetIntact(itnum)){  
								newpk=false;
								thepk=pknum;
							}
						}
					}//Close loop over all pairs that form pseudoknot with pair i,GetPair(i).
				}//Close loop over all intact pairs in pseudoknot 'pk'.
			}//Close loop over all instances of pseudoknot 'pk'.
	  
			//If no new pseudoknots need to be created and 'newpk' is false, the broken pair 'i' is then stored
			//in the pseudoknot 'thepk'.
			if(!newpk){
				pk[thepk].SetBroken(i,pairsBroken[i]);
				//Start Looking over all nucleotides within the pseudoknot 'thepk'
				for(k=i+1;k<pairsBroken[i];++k){
					//If the pair involving nucleotide 'k' is a pseudoknotted pair
					if(rna->GetPair(k,strnum)!=0&&(rna->GetPair(k,strnum)<i||rna->GetPair(k,strnum)>pairsBroken[i])){
						bool newIntact=true;//keep track of new intact pairs
						//Start loop over all intact pairs in the 'pknum' instance of pseudoknot 'pk'.
						for(itnum=0;itnum<pk[thepk].GetNumberIntact();++itnum){
							//If the pair is already stored, set 'newIntact' to false
							if(k==pk[thepk].GetIntact(itnum)) newIntact=false;
						}
						//If the 'newIntact' was never set to false, it means that the pair have not been stored yet, then store it
						if(newIntact==true){
							pk[thepk].SetIntact(k,rna->GetPair(k,strnum));
							pk[thepk].SetIntact(rna->GetPair(k,strnum),k);
						}
					}
				}
				newpk=true;
			}
			//If 'newpk' is true, this means that the new instance of a pseudoknot has to be created.
			else {
				//New instance of the pseudoknot 'pk' is created.
				pk.push_back(Pseudoknot());
				//Broken pair 'i' is stored in a newly created pseudoknot 'pk'...
				pk[pk.size()-1].SetBroken(i,pairsBroken[i]);
				//...and all the pairs that form a pseudoknot with the pair 'i' are also stored in this newly created
				//pseudoknot 'pk' instance.
				//Start search over all pairs that form a pseudoknot with pair i,GetPair(i).
				for(k=i+1;k<pairsBroken[i];++k){
					//Find pairs that form a pseudoknot with 'i'.
					if(rna->GetPair(k,strnum)!=0&&(rna->GetPair(k,strnum)<i||rna->GetPair(k,strnum)>pairsBroken[i])){
						//Assigning pairs to the new pseudoknot 'pk' instance that form a pseudoknot with 'i' in such
						//a way that the 5' nucleotide is always smaller then 3' as required by 'basepair'.
						if(k<rna->GetPair(k,strnum)) pk[pk.size()-1].SetIntact(k,rna->GetPair(k,strnum));
						else pk[pk.size()-1].SetIntact(rna->GetPair(k,strnum),k);
					}
				}//Close loop over all pairs that form pseudoknot with pair i,GetPair(i).
			}//Close the condition when the new pseudoknot 'pk' instance needs to be created.
		}//Close loops over all nucleotides that were broken.
	}//Close loops over all nucleotides that were broken.
	//Sort all pairs in all pseudoknots 'pk' from 5' to 3'.
	for(pknum=0;pknum<pk.size();++pknum) pk[pknum].Sort();

	return pk;
}


void DistanceCounter(RNA *rna, vector<Pseudoknot> &pk, vector<int> &rnacopy, vector<int> &UNcount, vector<int> &IBcount, vector<vector<int> > &SHLcount, int strnum)
{
	/////////////////////////////////////////////////////////////////////////////////////////////
	//FINDING # OF UNPAIRED NUCLEOTIDES, # OF INTERNAL BRANCHES, AND SPANNING HELIX LENGTHS
	/////////////////////////////////////////////////////////////////////////////////////////////
	//If there are pseudoknots in this structure start the search .

	int tempSHL;//Holds temporary SHLcount values.
	int itnum, b, e;
	
	//Start loop over all pseudoknots 'pk'.
	for(unsigned int pknum=0;pknum<pk.size();++pknum){

		//Creating vector 'pairsInPseudoknot'
		//Vector 'pairsInPseudoknot' temporarily stores pairs that comprise a given pseudoknot.
		//This allowes to quickly see if a pair in question is a part of a given psuedoknot.
		//Initialize all values in 'pairInPseudoknot' vector to 0.
		vector<int> pairsInPseudoknot = vector<int> (rna->GetSequenceLength()+1,0);
		//Copy all broken pseudoknotted pairs in a given pseudoknot to vector 'pairsInPseudoknot'.

		for(itnum=0;itnum<pk[pknum].GetNumberBroken();++itnum){
			//First copy the 3' side,...
			pairsInPseudoknot[pk[pknum].GetBroken(itnum).Get5pNucleotide()]=pk[pknum].GetBroken(itnum).Get3pNucleotide(); 
			//...then copy the 5' side.
			pairsInPseudoknot[pk[pknum].GetBroken(itnum).Get3pNucleotide()]=pk[pknum].GetBroken(itnum).Get5pNucleotide();
		}

		//Copy all intact pseudoknotted pairs in a given pseudoknot to vector 'pairsInPseudoknot'.
		for(itnum=0;itnum<pk[pknum].GetNumberIntact();++itnum){
			//First copy the 3' side,...
			pairsInPseudoknot[pk[pknum].GetIntact(itnum).Get5pNucleotide()]=pk[pknum].GetIntact(itnum).Get3pNucleotide(); 
			//...then copy the 5' side
			pairsInPseudoknot[pk[pknum].GetIntact(itnum).Get3pNucleotide()]=pk[pknum].GetIntact(itnum).Get5pNucleotide();
		}//Finished "Creating vector 'pairsInPseudoknot'".

		//Find the beginning and the end of the search.
		//Set 'b' be the LAST nucleotide of the FIRST 5'-most stacked helix of the given pseudoknot, ignoring 
		//single bulges. 
		//Preset 'b' value to 1 so that the search starts at the beginning of the sequence.
		b=1;
		//Move 'b' to the first nucleotide of the pseudoknot.
		while(pairsInPseudoknot[b]==0&&b<=rna->GetSequenceLength()) ++b;		 
	  
		//Set 'e' be the FIRST nucleotide of the LAST 3'-most stacked helix of the given pseudoknot, ignoring 
		//single bulges. 
		//Preset 'e' value to the last nucleotide of the sequence.
		e=rna->GetSequenceLength();
		//Move 'e' to the last nucleotide of the pseudoknot.
		while(pairsInPseudoknot[e]==0&&e>0)--e;
	  
		//Move the 'b' value to the end of the first 5'-most stacked helix of the given pseudoknot, ignoring single bulges.
		while(rnacopy[b+1]==rnacopy[b]-1||(rnacopy[b+1]==rnacopy[b]-2&&rnacopy[rnacopy[b]-1]==0)||(rnacopy[b+2]==rnacopy[b]-1&&rnacopy[b+1]==0)){
			//If there is no bulge, advance 'b' by 1.
			if(rnacopy[b+1]==rnacopy[b]-1||rnacopy[b+1]==rnacopy[b]-2) ++b;
			//If there is a bulge, advance 'b' by 2.
			else b=b+2;
		}
		//Move starting nucleotide 'b' by one nucleotide to leave the first stacked helix of the pseudoknot.
		b=b+1;

		//Move the 'e' value to the beginning of the last 3'-most stacked helix of the given pseudoknot, ignoring single bulges.
		while(rnacopy[e-1]==rnacopy[e]+1||(rnacopy[e-1]==rnacopy[e]+2&&rnacopy[rnacopy[e]+1]==0)||(rnacopy[e-2]==rnacopy[e]+1&&rnacopy[e-1]==0)){
			//If there is no bulge, advance 'e' by 1.
			if(rnacopy[e-1]==rnacopy[e]+1||rnacopy[e-1]==rnacopy[e]+2) --e;
			//If there is a bulge, advance 'e' by 2.
			else e=e-2;
		}//Finished "Find the beginning and the end of the search."

		//Create vector 'SHLrepeat' to hold nucleotides that are paired to those that were counted by SHLcount.  
		//This will allow to identify kissing loop pseudoknots and NOT count them twice.
		vector<int> SHLrepeat;
		bool SHLrepeatBOOL=false;//Bool that keeps track of repeats. "true" when there is a repeat.

		//Start counting UN, IB, and SHL for every pseudoknot 'pk' in the range from 'b' to 'e'.
		while(b<e){
		
			//UNcount: If nucleotide 'b' is not paired, +1 is added to UNcount.
			if(rnacopy[b]==0) {
				++UNcount[pknum];
				++b;
			}
			//If nicleotide 'b' is a part of the removed internal pseudoknot, +1 is added to UNcount.
			else if(rnacopy[b]!=0&&rnacopy[b]!=pairsInPseudoknot[b]&&rna->GetPair(b,strnum)==0){
				++b;
				++UNcount[pknum];
			}

			//IBcount: If the nucleotide 'b' is paired and not a part of any pseudoknotted helix (not part of 
			//'pairsInPseudoknot'), which means that it is a part of internal branch or a part of internal pseudoknot, 
			//skip the nucleotide counter 'b' to the pair of nucleotide 'b', and add 1 to the IBconter.
			else if(rnacopy[b]!=0&&rnacopy[b]!=pairsInPseudoknot[b]&&rna->GetPair(b,strnum)!=0){
				b=rnacopy[b]+1;
				++IBcount[pknum];
			}

			//SHLcount: If the nucleotide 'b' is paired and is also a part of pseudoknotted helix (part of pairsInPseudoknot),
			//start determining what size of spanning helix it is.
			else if(rnacopy[b]!=0&&rnacopy[b]==pairsInPseudoknot[b]){
							
				//Reset the temporary SHL counter to 1.
				tempSHL=0;
				//Start walking down till the end of the stacked pseudoknotted helix, ignoring single bulges.
				while((rnacopy[b+1]==rnacopy[b]-1&&rnacopy[b]-1!=0)||(rnacopy[b+1]==rnacopy[b]-2&&rnacopy[rnacopy[b]-1]==0&&rnacopy[b]-1!=0&&rnacopy[b]-2!=0)||(rnacopy[b+2]==rnacopy[b]-1&&rnacopy[b+1]==0&&rnacopy[b]-1!=0)){

					//Check if the nucleotide 'b' was already counted by SHLcount
					for(int i=0;i<SHLrepeat.size();++i){
						if(b==SHLrepeat[i]) SHLrepeatBOOL=true;
					}
					//Add the nucleotide that pairs to 'b' to SHLrepeat. Nucleotide that pairs to 'b' is a candidate for SHLcount. 					
					SHLrepeat.push_back(rnacopy[b]);

					//If 'b' was counted, advance 'b' depending on the bulge presence/location.
					if(SHLrepeatBOOL){
						//If there is no bulge in the stack:
						if(rnacopy[b+1]==rnacopy[b]-1) ++b;
						//If there is a bulge adjacent to rnacopy[b]:
						else if((rnacopy[b+1]==rnacopy[b]-2)&&(rnacopy[b]-2!=0)) ++b;
						//If there is a bulge adjacent to b:
						else if(rnacopy[b+2]==rnacopy[b]-1) b=b+2;
					}
					
					//If 'b' was not counted before, count:
					else {
						//Count how many stacked pairs there are in the given helix.
						++tempSHL;
						
						//If there is no bulge in the stack:
						if(rnacopy[b+1]==rnacopy[b]-1){
							//If there is no bulge, advance 'b' only by 1.
							++b;
						}
						
						//If there is a bulge adjacent to rnacopy[b].
						else if((rnacopy[b+1]==rnacopy[b]-2)&&(rnacopy[b]-2!=0)){
							//If there is a bulge on the 3' end, advance 'b' only by 1.
							++b;
						}
						
						//If there is a bulge adjacent b.
						else if(rnacopy[b+2]==rnacopy[b]-1){

							//If there is a bulge on 5' end, advance 'b' by 2.
							b=b+2;
						}
						
					}
					SHLrepeatBOOL=false;
				}	

				//HERE IS THE POINT FOR THE LAST NUCLEOTIDE IN THE SPANNING HELIX
				//Add the LAST nucleotide in the stacked pseudoknotted helix that pairs to 'b' to SHLrepeat.  					
				SHLrepeat.push_back(rnacopy[b]);
				
				//Check if the LAST nucleotide 'b' was already counted by SHLcount
				for(int i=0;i<SHLrepeat.size();++i){
					if(b==SHLrepeat[i]) SHLrepeatBOOL=true;
				}
				//If 'b' was NOT counted, add 1 to tempSHL
				if(!SHLrepeatBOOL) ++tempSHL;
				SHLrepeatBOOL=false;
				//Advance 'b' by ONE at the end of the count.
				++b;
				//Store the size of the helix in the proper location in 'SHLcount' vector.
				//If the size of the helix is smaller then 15 nucleotides, store every size instance separately.
				if(tempSHL<15) ++SHLcount[pknum][tempSHL];
				//If the size of the helix larger or equal to 15 ucleotides, store all together.
				else ++SHLcount[pknum][15];
			}//Finished counting SHLcount.
		}//Finished counting UN, IB, and SHL.
	}
}


int ReadPseudoParam(double *pseudo_params)
{
	// fn holds the file name with pseudoknot penalty calculation constants. The file is usually data_tables/pseudconst.dat
	string fn(getDataPath()); fn+="/pseudconst.dat";
	string lineoftext;//String that temporarily holds strings of data from 'pseudconst.dat'.

	//Open file with address 'fn' that holds pseudoknot penalty calculation constants.
	ifstream pseudoConst(fn.c_str());
	//If no file, output to user.
	if(!pseudoConst.good()){
		cerr << "File \"pseudconst.dat\" with pseudoknot penalty constants was not found in the DATAPATH directory." << endl;
		cerr << "Please check the DATAPATH environment variable and verify the existence of this file." << endl;
		return 33; // Error opening pseudoknot penalty constants file
	}
	//Copy data from file to string 'lineoftext'.
	pseudoConst >> lineoftext;
	//There are 16 constants that are stored in the file 'pseudconst.dat'. For all 16 constants iterate over i.
	for(int i=0;i<=15;++i){
		//Keep reading file till '-->' is reached. 
		while(strcmp(lineoftext.c_str(),"-->")!=0) pseudoConst>>lineoftext;
		//When '-->' is reached, store it in the constant.
		pseudoConst>>pseudo_params[i];
		//Read next entry from file.
		pseudoConst >> lineoftext;
	}
	//Finish "Read datapath".
    return 0;
}

//Function that restores removed pseudoknotted helix back to the structure
void RestorePseudoknot(RNA *rna, vector<int> &pairsBroken, int strnum)
{
    //Restore the broken pairs from 'rna
    for(int i=1;i<=rna->GetSequenceLength();++i){
        if(pairsBroken[i]!=0) rna->SpecifyPair(i,pairsBroken[i],strnum);
    }

    //Finished
    return;
}

//Function that restores any changes made to the structure
void RestoreStructure(RNA *rna, vector<int> &pairsChanged, int strnum)
{

    //Restore the changed pairs from 'rna
    for(int i=1;i<=rna->GetSequenceLength();++i){

        if(pairsChanged[i]!=rna->GetPair(i,strnum)){
            if(pairsChanged[i]!=0) rna->SpecifyPair(i,pairsChanged[i],strnum);
            else rna->RemoveBasePair(i,strnum);
        }
    }

    //Finished
    return;
}

//Function that returns the pseudoknot statistics to the screen
void RetToScreen (vector<Pseudoknot> &pk, vector<int> &UNcount, vector<int> &IBcount, vector<vector<int> > &SHLcount)
{
	/////////////////////////////////////////////////////////////////////////////////////////////
	//RETURN INFORMATION TO THE SCREEN
	/////////////////////////////////////////////////////////////////////////////////////////////

	int i, itnum;		  
			
    cerr << "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";

	//If there are no pseudoknots in this structure let user know and exit.
	if(pk.size()==0){
		cout << "THERE ARE NO PSEUDOKNOTS IN THIS STRUCTURE\n";
	}

	//If there are pseudoknots in this structure, proceed.
	else{
		/*
		//If there is only one pseudoknot report that.
		if(pk.size()==1){
		cout << ":::::::::::::::::::::::::::::::::::::::::::::::" << endl;
		cout << "	THERE IS ONLY ONE PSEUDOKNOT IN THIS STRUCTURE" << endl;
		cout << ":::::::::::::::::::::::::::::::::::::::::::::::" << endl << endl;
		}

		//If there are more then one pseudoknots, report the number of pseudoknots.
		if(pk.size()>1){
		cout << ":::::::::::::::::::::::::::::::::::::::::::::::" << endl;
		cout << "	THERE ARE " << pk.size() << " PSEUDOKNOTS IN THIS STRUCTURE" << endl;
		cout << ":::::::::::::::::::::::::::::::::::::::::::::::" << endl << endl;
		}
		*/
		//Start looping over all pseudoknots.
		for(unsigned int pknum=0;pknum<pk.size();++pknum){

			//Report UN, IB, SHL for every pseudonot.
			cout << "		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			//			cout << "					 FOR PSEUDOKNOT #" << pknum+1 << ": " << endl  << "		  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
			cout << "Unpaired Nucleotide Count is:	   " << UNcount[pknum] << " nucleotide(s)" << endl;
			cout << "Internal Branch Count is:		   " << IBcount[pknum] << " branch(s)" << endl;
			cout << "Spanning Helix Lengths are:	   " << endl;
			for(i=1;i<=15;++i){
				if(SHLcount[pknum][i]!=0&&i<15){
					cout << "								   " << SHLcount[pknum][i] << " count(s) of "<< i << "nt spanning helixe(s)" << endl;
				}
				else if(SHLcount[pknum][i]!=0&&i==15){
					cout << "								   " << SHLcount[pknum][i] << " count(s) of " << "spanning helixe(s) with length longer then 15nt" << endl;
				}
			}
			cout << endl << endl;

			//Report all pairs that make pseudoknots for every pseudoknot
			cout << "Pseudoknot # " << pknum+1 << " consists of helixes with the following pairs:\n"<< '\n';
			//Loop over all broken pseudoknotted pairs in every pseudoknot 'pk'
			for(itnum=0;itnum<pk[pknum].GetNumberBroken();++itnum){
				//Report the 5' pair
				cout << "		   "<< pk[pknum].GetBroken(itnum).Get5pNucleotide(); 
				cout << " is paired to ";
				//Report the 3' pair		
				cout << pk[pknum].GetBroken(itnum).Get3pNucleotide();
				cout << '\n';
			}
			cout << '\n' << "				   AND" << '\n' << '\n';
			//Loop over all intact pseudoknotted pairs in every pseudoknot 'pk'
			for(itnum=0;itnum<pk[pknum].GetNumberIntact();++itnum){
				//Report 5' pair
				cout << "		   " << pk[pknum].GetIntact(itnum).Get5pNucleotide(); 
				cout << " is paired to ";
				//Report 3' pair		
				cout << pk[pknum].GetIntact(itnum).Get3pNucleotide();
				cout << '\n';
			}
			cout << '\n' << '\n';
		}
	}
}


//Function that calculates the energy penalty for a pseudoknot.
//Input: RNA class - instance of RNA class that holds the structure
//'a' - holds the pseudoknot experimental distances
//P1 and P2 - pseudoknot constants
//strnum - structure number.
double pkEnergyCalc(RNA *rna, double *a, double P1, double P2, int strnum)
{

	double energyPenalty=0;//Stores energy penalty for predicting pseudoknots in the structure. 

	//Check if 'rna' still contains a pseudoknot after removing isolated pairs. If it does, execute:
	if(rna->ContainsPseudoknot(strnum)){
        
        //Create a vector 'rnacopy' that stores a copy of all pairs from rna structure 'rna'
		vector<int> rnacopy = RNAcopy(rna, strnum);
        
        //PairsBroken outputs vector that stores all broken pairs
        //IMPORTANT: it also breaks the pseudoknot in 'rna'
		vector<int> pairsBroken = PairsBroken(rna, rnacopy, strnum);
        
		//Declare vector of pseudoknots 'pk' that holds a class 'pseudoknot' of vectors of classes.
		//Parse through the structure to separete all pseudoknots and store pairing values for all
		//pseudoknotted pairs in the structure 'rna'
		vector <Pseudoknot> pk=PseudoParser(rna, pairsBroken, strnum);
        //Initialise UNcount, IBcount, and SHLcount for use in DistanceCounter.
		//Vector UNcount keeps track of Unpaired Nucleotides for each pseudoknot. Unpaired Nucleotides 
		//are not calculated when they are a part of a single nucleotide bulge of a pseudoknotted helix. 
		//They are also not calculated when they is inside of an internal branch.
		//Initialize vector 'UNcount' and set its length equal to the number or pseudoknots in the 'rna' structure 
		//so there is an UNcount entry for every pseudoknot. And set all values of 'UNcount' to 0 by default.
		vector<int> UNcount = vector<int>(pk.size(),0);
		//Vector IBcount keeps track of number of Internal Branches for each pseudoknot. Internal Branch
		//is a helix that is inside a specific pseudoknot, but is not a part of that pseudoknot. 
		//It can be either a pseudoknotted or non-pseudoknotted helix. 
		//Initialize vector 'IBcount' and set its length equal to the number or pseudoknots in the 'rna' structure
		//so there is a IBcount entry for every pseudoknot. And set all velues of 'IBcount' to 0 by default.
		vector<int> IBcount = vector<int> (pk.size(),0);
		//Vector SHLcount is a 2D vector that keeps track of the Spanning Helix Length for each pseudoknot. 
		//Spanning Helix Length is length of pseudoknotted helixes that are a part of the given pseudoknot 
		//excluding the first and last stacked helixes. Internal pseudoknots are ignored. 
		//Initialize 2D vector 'SHLcount' and set its length equal to the number or pseudoknots in this structure 
		//in one dimention (rows), and to 16 in the other dimention(columnts), so there is 'SHLcount' entry 
		//for every pseudoknot and for every spanning helix length up to 14 nucleotides. Any helix with length
		//longer then 14 nucleotides stored in position 15.
		vector<vector<int> > SHLcount;
		for(unsigned int l=0;l<pk.size();++l) SHLcount.push_back(vector<int>(16,0));
		//Finish "Initialise UNcount, IBcount, and SHLcount."

		//Calculate the number of unpaired nucleotides, number of internal branches, and length and type 
		//of spanning helixes for each pseudoknot from vector 'pk'.
		DistanceCounter(rna, pk, rnacopy, UNcount, IBcount, SHLcount, strnum);

		//RetToScreen (pk, UNcount, IBcount, SHLcount);
        
        //Restore the broken helix
        RestorePseudoknot(rna, pairsBroken, strnum);

		//Calculate Energy Penalty for the structure in 'rna' class.
		
		double SHLcountTOTAL=0;//Holds the sum of products of spanning helix lengths and its corresponding constants.
		//Calculate 'energyPenalty' for all pseudoknots in the structure.
		for(unsigned int k=0;k<pk.size();++k){
			//Calculate the sum of prodictes of spanning helix lengths and its corresponding constants.
			for(int j=2;j<=15;++j){
				SHLcountTOTAL = SHLcountTOTAL + a[j]*a[j]*(double)SHLcount[k][j];
			}
            //If only single base-pair spanning helices are present, set SHLcountTOTAL equal to 1 so that the "P2" term 
            //of the energy penalty equation equals to 0
            if(SHLcountTOTAL==0)SHLcountTOTAL=1;

			//Calculate the energy penalty.
            if(UNcount[k]==0&&IBcount[k]==0){
                energyPenalty = energyPenalty + (double)(((P1+P2*log(SHLcountTOTAL))));
            }

            else{
                energyPenalty = energyPenalty + (double)(((P1*log(a[0]*a[0]*((double)UNcount[k])+a[1]*a[1]*((double)IBcount[k]))+P2*log(SHLcountTOTAL))));
            }
        }
	}
    return energyPenalty;
}

//Function that calculates the energy of the removed pseudoknotted helix.
//It requires the following input:
//RNA class, structure number, data tables
double CalculatePKhelixEnergy(RNA* rna, int strnum, datatable* data) { return CalculatePKhelixEnergy(rna->GetStructure(), strnum, data);}
double CalculatePKhelixEnergy(structure* ct, int strnum, datatable* data)
{
  double energy=0;//store energy in 'energy'
    for(int i=1;i<=ct->GetSequenceLength()-7;++i){//for every nucleotide up to 7 nucleotides till the end of the sequence
        if(i<ct->GetPair(i,strnum)&&ct->GetPair(i,strnum)!=0){//if it is a 5' to 3' pair
            //If there is a stack
            if(ct->GetPair(i+1,strnum)==ct->GetPair(i,strnum)-1){
                //get the energy of the stack
                energy+=erg1(i, ct->GetPair(i,strnum), i+1, ct->GetPair(i,strnum)-1, ct, data);

                //calculate AU penalty if it is opening or closing pair
                //Check if it is the first nucleotide
                if(i==1) energy+=penalty(i,ct->GetPair(i,strnum),ct,data);
                //check if it is the last nucleotide
                else if(ct->GetPair(i,strnum)==ct->GetSequenceLength()) energy+=penalty(i,ct->GetPair(i,strnum),ct,data);
                else if((ct->GetPair(i-1,strnum)==0&&ct->GetPair(ct->GetPair(i,strnum)+1,strnum)==0)||
                   (ct->GetPair(i+1,strnum)==0&&ct->GetPair(ct->GetPair(i,strnum)-1,strnum)==0)){
                    energy+=penalty(i,ct->GetPair(i,strnum),ct,data);
                }

            }
            //If there is a bulge or internal loop
            /*  
            else{
                //store current i
                int Ci=i;
                //store current j
                int Cj=ct->GetPair(i,strnum);
                //initialize the i prime and j prime
                int Pi, Pj;
                
                //find the i prime. 
                while(ct->GetPair(i+1,strnum)==0) i++;
                //If the i+1 is still 5 prime
                if(i+1 < ct->GetPair(i+1,strnum)){
                    //Store i prime
                    Pi=i+1;
                    //Store j prime
                    Pj=ct->GetPair(i+1,strnum);
                    if((Pi-Ci)<3&&(Cj-Pj)<3){
                        energy+=erg2(Ci, Cj, Pi, Pj, ct, data, 0, 0);
                    }
                }
            }
            */
        }
    }
    return (double)energy/(double)conversionfactor;
}

double CalculateFreeEnergywithPseudoknot(RNA* rna, int strnum, double P1, double P2, double* pseudo_params, datatable* data, bool UseSimpleMBLoopRules, bool fillMismatch, bool remIsolatedPairs, ofstream* writeDetails)
{
    
    //double pseudo_params[16];//Array that holds pseudoknot penalty calculation constants.
    //Read pseudoknot parameters from the datapath
    //ReadPseudoParam(pseudo_params);

    //Stores the total energy
    double TotalEnergy;
    vector<int> rnacopy2;

    if(fillMismatch||remIsolatedPairs){
    //Copy the rna pairings into a vector so we can restore them later
    rnacopy2= RNAcopy(rna, strnum);
    }

    //If the 'fillMismatch' option was specified, proceed
    if(fillMismatch){
        FillMismatch(rna,strnum);

    }
    //If the 'remIsolatedPairs' option was specified, proceed
    if(remIsolatedPairs){
        RemoveIsolatedPairs(rna,strnum);
    }

    //Check if there is still a pseudoknot in the 'rep' structure
    if(rna->ContainsPseudoknot(strnum)){
        //Calculate the pseudoknot penalty
        double PseudPenalty=pkEnergyCalc(rna, pseudo_params, P1, P2, strnum);
        //Calculate the energy of the non pseudoknotted structure and the broken helices

        //Make a vector that would hold all pairs from rnaCT2
        vector<int> rnacopy = RNAcopy(rna, strnum);
		// RMW TODO: use more efficient pseudoknot code, e.g.:
		//  vector<int> normal(rnacopy.size()), knots(rnacopy.size());
		//  rna->GetStructure()->FindPseudoknots(strnum, &knots, &normal);

        //Break the pseudoknot
        BreakPseudoknot(rna, strnum);
        
        //Set the energy to 0, since it will be recalculated
        rna->GetStructure()->SetEnergy(strnum,0);

        //Calculate the energy of the structure with the broken pseudoknot (Energy 1)
        double Energy1=rna->CalculateFreeEnergy(strnum,UseSimpleMBLoopRules);
        //double Energy1=rna->CalculateFreeEnergy(strnum);

        //Only specify the broken pseudoknotted helix in structure
        for(int i=1;i<= rna->GetSequenceLength();++i){
            if(rnacopy[i]== rna->GetPair(i,strnum)&rnacopy[i]<i){
                rna->RemoveBasePair(i,strnum);
            }
            else if(rnacopy[i]<i){
                rna->SpecifyPair(i,rnacopy[i],strnum);
            }
        }

        double Energy2=0;//Holds the energy of the broken pseudoknot
        //Do not support very complicated pseudoknots for now
        if( rna->ContainsPseudoknot(strnum)){
            Energy2=999;
        }

        //If the pseudoknot is supported, 
        else{
            //Calculate the energy of the structure with only the broken helix (Energy 2)
            //Energy2=((rna[rep].Return_RNA()->CalculateFreeEnergy(1)-4.09)*conversionfactor;
            Energy2=CalculatePKhelixEnergy(rna,strnum,data);
        }

        TotalEnergy = Energy1+Energy2+PseudPenalty;

		if (writeDetails!=NULL) {
			// Output should be compatible with efn2's thermodynamic detail output.
			vector<int> knots(rnacopy.size()), normal(rnacopy.size());
			for(int i=1;i<= rna->GetSequenceLength();++i)
				((rnacopy[i]==rna->GetPair(i,strnum))?knots[i]:normal[i])=rnacopy[i];
			ostream &out = *writeDetails;
			out << "Pseudoknot Energy Details for Structure " << strnum << ":" << endl;
			out << "\tKnotted Pairs:";
			int pairCount=0;
			for(int n=1;n<=rna->GetSequenceLength();n++)
				if (knots[n]!=0 && n<knots[n]) {
					if ((pairCount++%8)==0) out << endl << "\t    ";
					else out << " ";
					out << n << ":" << knots[n];
				}
			out << endl;
			out << "\tEnergy of Broken Structure: " << fixed << setprecision(1) << Energy1 << endl;
			out << "\tEnergy of Pseudo-Helices:   " << fixed << setprecision(1) << Energy2 << endl;
			out << "\tPseudoknot Energy Penalty:  " << fixed << setprecision(1) << PseudPenalty << endl;
			out << "\tTotal Free Energy:          " << fixed << setprecision(1) << TotalEnergy << endl;
		}

        RestorePseudoknot(rna, rnacopy, strnum);

    }
    
    //if there is no pseudoknotted pairs in the structure calculate free energy
    else{


        TotalEnergy = rna->CalculateFreeEnergy(strnum,UseSimpleMBLoopRules);
        
    }

    if(fillMismatch||remIsolatedPairs){
        RestoreStructure(rna, rnacopy2, strnum);
    }
	rna->GetStructure()->SetEnergy(strnum,(int)(TotalEnergy*conversionfactor));
    return TotalEnergy;
}
