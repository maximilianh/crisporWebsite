/*=======================================================================
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
intermolecular_test.cpp fold the whole sequence, save partition function in a file 
and reuse it.
stochastic sampling method were also embedded.

They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Aug. 6, 2006
															John(Zhi Lu)



Changes made by DHM: 7/12/09:
Removed folding when break local structure is used.

=======================================================================*/
//#define debugmode true
#undef debugmode
#define TESTMODE false   //read the save file instead of calculating if in TEST mode
#include "intermolecular.h"
#include "common_utils.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

// Verifies that the ct->datatable matches the datatable passed to dynamic.
#define SAFE_DYNAMIC(CT,DATA,...) do { \
			if (CT->GetThermodynamicDataTable()!=DATA) { \
				cerr << endl << "RMW: CT/datatable Mismatch!! " \
				<< "CT:" << CT->GetThermodynamicDataTable() << " DATA:" << DATA \
				<< " line " << __LINE__ << endl; \
				CT->SetThermodynamicDataTable(DATA); } \
			dynamic(CT,DATA,__VA_ARGS__); \
		break; } while (true)
// #define SAFE_DYNAMIC(CT,DATA,...) dynamic(CT,DATA,__VA_ARGS__)

//define the parameters for alltrace()
#define percent 100	//100% of the suboptimal structure windows
#define delta 6	//only the suboptimal structures having energy less than 2.8kcal/mol higher than optimal one
                    //not store structures in ct file for alltrace() to save memory
// #define MaxStructure 10000000 (10M) is defined in alltrace.h

/*
 olig is the backend function for oligo walk

 table -- A table (2D array) that is filled with the thermodynamic data 
 	for each oligo. The first dimension is an array of rows, 
		each representing an oligo bound to a specific site on the target.
 	The second dimension is an array with the following indices:
 	[0] - overall DG
 	[1] - duplex DG
 	[2] - free energy of breaking target structure
 	[3] - intramolecular oligo free energy
 	[4] - intermolecular oligo free energy
 	[5] - Tm of duplex (x10)
 
 mode -- How the target structure should be handled.
 	1 - break local target structure to bind oligo 
 	2 - refold target RNA after oligo binding
 	3 - no target structure considered
 
 suboptimal - How suboptimal structures should be handled.
 	0 - only consider optimal structure
 	1 - like choice 3,using suboptimal structures,but the whole set from alltarce() function prediction
 	2 - using partition funcion considering every possible structure  of target
 		suboptimal 2 can only used with mold 2
 	3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
 	4 - using stochasit sampling method to sample 1000 structures
 
 prefilter  -- Whether criteria should be used to prefill functional siRNA
 	0 - no prefilter
 	1 - use prefilter
 
 foldsize -- Only fold a fragment with size=foldsize+binding_length, 
 	which is centered on the siRNA binding region.
 	when foldsize > 1, only mode=2 and suboptimal=2 are valid options.
 
 distance  -- limit the maximum distance between nucleotides that can pair
 
 shapefile -- specify a SHAPE datafile (set to "" to not use SHAPE data)
 
 TESTS -- run tests, set to -1 for no tests.

 writesav - write sav files to save time in test mode
*/
	void olig(const bool isdna, const int mode, structure *ct, const int length, const double conc, const int suboptimal,
		const int start, const int stop, int foldsize, const int distance, 
		int **table, int **numofsubstructures, const char *const shapefile, int *TEST, const bool writesav,
		datatable& data, datatable& ddata, thermo* helixstack, rddata *hybriddata, 
		siPREFILTER *prefilter, ProgressHandler *update) {
	
	int i,j,k;
	int ip,jp;
	int foldstart, foldstop;
	int dh,ds,dgeff;
	int *energyarray, *ctenergy;//ctenergy is the energy array to store structure energies from alltrace()
	int numofstructures;
	long double k1b,k1u;
	long double fn,sn,sum;
	char savefile[250],pfsfile[250],pos[250];
	int energy=INFINITE_ENERGY;//to store the free energy of lowest free energy structure
	int *temp,**temp2;//temp will store basepairing info in mode == 1
	

	PFPRECISION Q,Qc,Qolig,pftemp=310.15;//store partion function without and with constrains
	PFPRECISION rescaleinrefill;//rescaling facor for Qc relative to Q
	pfdatatable *pfdata,*dpfdata;  //store the data tables for partionfunction
	OligoPclass *interoligo,*intraoligo;
	OligoPclass *target,*targetcopy,*targettemp;
	structure *oligo,*oligo1;//to store the structural info for the oligo
						   //oligo is intermolecular structure; oligo1 is intramolecular
	structure *fracct;// to store folded fraction of target centered at siRNA binding site
	

	//Infrastructure for 
	//define shape information for the whole sequnence (ct)
	if (!is_blank(shapefile)) {
		ct->SHAPEslope = 32;
		ct->SHAPEintercept = -10;
		ct->ReadSHAPE(shapefile) ;
		if (suboptimal==2) pfshape(ct,pftemp);
	}


	
	//----------------------------------------------------------------------------------------------
	//allocate space in oligo for the sequence
	//difine for inter oligo sequence, oligo structure can be used both by intra and inter molecular
	//oligo:  intermolecular structure   aagucXXXggcaa
	//oligo1: intramolecular structure   aaguc
	oligo = new structure();
	oligo1 = new structure();

	oligo->allocate(2*length+3);
	oligo1->allocate(length);

	ct->SetThermodynamicDataTable(&data); // always RNA
	if (isdna) {
		oligo->SetThermodynamicDataTable(&ddata);
		oligo1->SetThermodynamicDataTable(&ddata);
	} else {
		oligo->SetThermodynamicDataTable(&data);
		oligo1->SetThermodynamicDataTable(&data);
	}

	// Show complement sequence.
	// cout << "SEQ Complement:\t"; for(int i=1;i<=ct->GetSequenceLength();i++) cout << complement(i, ct); cout << endl;

	string label;
	label = "Oligo_inter";

	oligo->SetSequenceLabel(label);

	label = "Oligo_intra";

	oligo1->SetSequenceLabel(label);

	
	for (j=1;j<=3;j++) {		oligo->inter[j-1] = length + j;	}

	foldsize= foldsize/2*2; // a trick to make foldsize even
	//----------------------------------------------------------------------------------------------
	//define the base pair distance constraint
	if (distance > 0) {
		ct->SetPairingDistance(distance);
		
	}
	//define the size of fracct, which is the region to be folded
	if (foldsize >0) {
		fracct = new structure();
		fracct->allocate(foldsize+length);
		fracct->SetThermodynamicDataTable(&data);
		label="fraction_of_target";
		fracct->SetSequenceLabel(label);
		
		//define the base pair distance constraint
		if (distance > 0) {
			fracct->SetPairingDistance(distance);
		}
	}	
	
	//----------------------------------------------------------------------------------------------
	//some variables declared for each mode 
	if (mode ==1) {//Break local structure mode
		if (suboptimal==0)	{//no suboptimals considered
			temp = new int [length];//allocate an array for storing base pairing information
			//if (foldsize == 0) {	
			ct->RemoveConstraints();
				//SAFE_DYNAMIC(ct,&data,1,10,0,0);
			efn2(&data,ct,1);
		   	energy = ct->GetEnergy(1);
			//}
		}
		else if (suboptimal==3) {//standard mfold heuristic for suboptimals
			//if (foldsize == 0) {	
			ct->RemoveConstraints();
				//SAFE_DYNAMIC(ct,&data,1000,10,0);
			efn2(&data,ct);
			energyarray = new int [ct->GetNumberofStructures()+1];
			for (k=1;k<=ct->GetNumberofStructures();k++) energyarray[k] = ct->GetEnergy(k);
			
			temp2 = new int *[ct->GetNumberofStructures()+1];
			for (i=1;i<=ct->GetNumberofStructures();i++)   		temp2[i] = new int [ct->GetSequenceLength()];
			
			numofstructures=ct->GetNumberofStructures();
			//}
      	}
		
	}
	else if (mode==2) {//refold for each oligo
   		temp = new int [ct->GetSequenceLength()+1];
		for (j=1;j<=ct->GetSequenceLength();++j) {
      		temp[j] = ct->GetPair(j);
		}
		if (suboptimal ==0 ) {//no suboptimals
			//only interested in lowest free energy structure, reduce the number of structures to 1
			for (int last = ct->GetNumberofStructures();last>1;--last) ct->RemoveLastStructure();

			if (foldsize == 0) {
				ct->RemoveConstraints();
				//ct->nnopair=0;
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s0.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&energy);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//write the save file information so that the fold need not to be done again
					
					// RMW_QUESTION: Should the following call to `dynamic` use RNA or DNA params when isdna is true?
					SAFE_DYNAMIC(ct,&data,1000,10,0,0,true);

					energy = ct->GetEnergy(1);
					if (writesav) {
						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&energy);
						sav.close();
					}
				}
			}
		}
		else if	(suboptimal==1)	 {	//use AllSub to generate suboptimals
			ctenergy=new int [MaxStructure+1];
			if (foldsize ==0) {	
				ct->RemoveConstraints();
				//ct->nnopair=0;
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s1.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//calculate the whole set of suboptimal structures
					alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ctenergy[j];}
					if (writesav) {
 						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if	(suboptimal==3)	 {//use the mfold heuristic to generate suboptimals	
			if (foldsize ==0) {	
				//ct->nnopair=0;
				ct->RemoveConstraints();
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s3.sav");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				if(TESTMODE && sav ) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					SAFE_DYNAMIC(ct,&data,1000,10,0);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->GetEnergy(j);}
					if (writesav) {
 						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if (suboptimal==2) {//use partition funtion to generate ensemble free energies
			if (isdna) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,dpfdata);
				interoligo = new OligoPclass(oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,pfdata);
				interoligo = new OligoPclass(oligo,pfdata);
			}
			//calculate partion function for the whole target without constrain
			//ct->nnopair=0;
			ct->RemoveConstraints();
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(savefile, ct->GetSequenceLabel().c_str());
				if (savefile[strlen(savefile)-1]=='\n') savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				//strcpy(savefile, ct->ctlabel[1]);
				//savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s2.sav");
				strcpy(pfsfile, ct->GetSequenceLabel().c_str());
				if (pfsfile[strlen(savefile)-1]=='\n') pfsfile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream sav(savefile,std::ios::binary);
				std::ifstream pfs(savefile,std::ios::binary);
				if (TESTMODE && sav && pfs  ){
					pfs.close();
					//read information from file if it was found
					read(&sav,&Q);
					for (i=0;i<=ct->GetSequenceLength();++i) {
						read(&sav,&(target->copyw5[i]));
						for (j=0;j<=ct->GetSequenceLength();++j) {
							read(&sav,&(target->copyv->f(i,j)));
							read(&sav,&(target->copyw->f(i,j)));
							read(&sav,&(target->copywca[i][j]));
							read(&sav,&(target->copywmb->f(i,j)));
							read(&sav,&(target->copywl->f(i,j)));
							read(&sav,&(target->copywmbl->f(i,j)));
							read(&sav,&(target->copywcoax->f(i,j)));
						}	
					}
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					pfs.close();
					//calculate the partition function if no file found
					if (writesav) {
						target->partition4refill(&Q,pfsfile);
					}
					else target->partition4refill(&Q);

				
					//write the save file information so that the partition function can be re-folded,
					if (writesav) {
						std::ofstream sav(savefile,std::ios::binary);
						write(&sav,&Q);
						for (i=0;i<=ct->GetSequenceLength();++i) {
							write(&sav,&(target->copyw5[i]));
							for (j=0;j<=ct->GetSequenceLength();++j) {
								write(&sav,&(target->copyv->f(i,j)));
								write(&sav,&(target->copyw->f(i,j)));
								write(&sav,&(target->copywca[i][j]));
								write(&sav,&(target->copywmb->f(i,j)));
								write(&sav,&(target->copywl->f(i,j)));
								write(&sav,&(target->copywmbl->f(i,j)));
								write(&sav,&(target->copywcoax->f(i,j)));
							}		
						}
						sav.close();
					}
						
				}
			}
			else if(prefilter->useit) {//using prefilter and fold different region each time
				target = new OligoPclass(fracct,pfdata);
			}
		}
		else if (suboptimal==4) {//use stochastic to generate a sample
			if (isdna) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,dpfdata);
				interoligo = new OligoPclass(oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(oligo1,pfdata);
				interoligo = new OligoPclass(oligo,pfdata);
			}
			//calculate partion function for the whole target without constraints
			ct->RemoveConstraints();
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(pfsfile, ct->GetSequenceLabel().c_str());
				if (pfsfile[strlen(savefile)-1]=='\n') pfsfile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				std::ifstream pfs(pfsfile,std::ios::binary);
				if(TESTMODE && pfs) {//read information from file if it was found
					pfs.close();				
				}
				else{
					//close the readonly file sav
					pfs.close();
					//calculate the partition function if no file found
					target->partition4refill(&Q,pfsfile); 
					
				}
				stochastic(ct,pfsfile,1000,5);
				efn2(&data,ct);
				//ofstream sout("RL_sample_oligo->out");
				//for (j=1;j<=1000;j++)	{ sout << ct->energy[j]<<"\n";	}
				//sout.close();
				numofstructures=ct->GetNumberofStructures();
				energyarray = new int [numofstructures+1];
				for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->GetEnergy(j);}
			}
			else if(prefilter->useit) {//using prefilter and fold different region each time
				//allocate target only, if prefilter==0, allocate both target and targetcopy later
				target = new OligoPclass(fracct,pfdata);

			}
		}
		
	}
	//having single-strand constrained for the calculation of constrained energy
	//ct->nnopair=length;
	//ct->checknopair();	

	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//begin to scan the target from start to stop
	for (i=start; i<=stop; ++i) {

		
		//communicate progress
		if (update!=NULL) {
   			update->update(i-start, stop-start);
		}

		//Remove previous constraints that were imposed:
		ct->RemoveConstraints();

		//define the oligo sequence
		//oligo1->intermolecular = false;
		//oligo1->GetSequenceLength() = length;
		//set sequence
		 for (j=1;j<=length;j++) {
			oligo1->numseq[j] = complement(i+length-j,ct);
			oligo->numseq[j] = complement(i+length-j,ct);
		}
   		//-------------------------------------------------------------------------------------------
		//prefiltering the functional siRNA
		if (prefilter->useit) {		
			prefilter->count(oligo1,i,TEST[i]);
			if (prefilter->score[i] < FILTER_PASS) {
				continue;
			}
		}
	
		//-------------------------------------------------------------------------------------------
		//-------------------------------------------------------------------------------------------
		//calculate the free energy of breaking target structure

		//Option 1: break local structure only
		if (mode==1) {
			if (foldsize ==0) {
				//-------------------------------------------------------------------------------------------
				//mode 1 + consider only the first structure:
				if (suboptimal ==0){
					//store basepairing info
		     		for (j=0;j<length;j++) {
						temp[j] = ct->GetPair(i+j);
		        		if (ct->GetPair(i+j) > 0) {
							ct->RemovePair(i+j);
		        			
			       		}
		     		}
	         	
					efn2(&data,ct);
		    		table[i][2] = ct->GetEnergy(1) - energy;
	         	
					//restore basepairing
		    		for (j=0;j<length;j++) {
		      			if (temp[j]>0) {
							ct->SetPair(i+j,temp[j]);
		  
		        		}
		    		}
				}	
				//-------------------------------------------------------------------------------------------
				//mode 1 + consider all subotimal structures:
				else if (suboptimal==3) { 
					sn = 0;
					sum = 0;
				  //store the basepairing info and force the region of hybridization
				 //to have no structure
				  for (k=1;k<=ct->GetNumberofStructures();++k) {
		          		for (j=0;j<length;j++) {//store basepairing info
		     				temp2[k][j] = ct->GetPair(i+j,k);
		           			if (ct->GetPair(i+j,k) != 0) {
								ct->RemovePair(i+j,k);
		         				
		          			}
		      			}

					}
					efn2(&data,ct);
				  for (k=1;k<=ct->GetNumberofStructures();k++) {

		          		fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
					   sn = sn + fn;
		   				sum = sum+fn*(((long double)(ct->GetEnergy(k)- energyarray[k])));
				   }
					//restore the basepairing
					for (k=1;k<=ct->GetNumberofStructures();k++) {
			       		for (j=0;j<length;j++) {
							if (temp2[k][j]>0) {
								ct->SetPair(i+j,temp2[k][j],k);
           						
           					}
        				}
					}
	         		
					table[i][2] = (int)(sum/sn);
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 

				}
			}
			if (foldsize > 0 ) {
				//define the region for refolding 
				if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->GetSequenceLength() ) {
					for (j=1;j<=foldsize+length;j++) {
      					fracct->numseq[j] = ct->numseq[i-foldsize/2+j-1];
					}
				}
				else if ( (i-foldsize/2)<=1 ) 
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[j];
				}
				else if( (i-1+length+foldsize/2)>=ct->GetSequenceLength() )
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[(ct->GetSequenceLength())+j-foldsize-length];
				}
				fracct->RemoveConstraints();
				//fracct->nnopair=0;
				//calculate the strucutre without constrains
				if (suboptimal==0)	{
					//Remove structures, after the first because only need 1 structure
					//for (int structurenum=fracct->GetNumberofStructures();structurenum>1;--structurenum) {
					//	fracct->RemoveLastStructure();

					//}
					SAFE_DYNAMIC(fracct,&data,1,10,0,0);
					efn2(&data,fracct);
		   			energy = fracct->GetEnergy(1);
				}
				else if (suboptimal==3) {
					SAFE_DYNAMIC(fracct,&data,1000,10,0);
					efn2(&data,fracct);
					energyarray = new int [fracct->GetNumberofStructures()+1];
					for (k=1;k<=fracct->GetNumberofStructures();k++) energyarray[k] = fracct->GetEnergy(k);
					temp2 = new int *[fracct->GetNumberofStructures()+1];
					for (k=1;k<=fracct->GetNumberofStructures();k++)   		temp2[k] = new int [fracct->GetSequenceLength()];
					numofstructures=fracct->GetNumberofStructures();
				
      			}
				//set the constrained position
				if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
					foldstart=i;
					foldstop=i+length-1;
				} 
				else if(i+length-1+foldsize/2>= (ct->GetSequenceLength()) ) {
					foldstart=foldsize+length-(ct->GetSequenceLength()) +i;
					foldstop=foldsize+length+length-1-(ct->GetSequenceLength()) +i;
				}
				else {//folded region begin to move, constrained position will be always in the middle of this region
					foldstart=foldsize/2+1;
					foldstop=foldsize/2+length;
				}
				
				//recalc the energy again
				if (suboptimal==0) {
					//store basepairing info
		     		for (j=foldstart;j<=foldstop;j++) {
						temp[j-foldstart] = fracct->GetPair(j);
		        		if (fracct->GetPair(j) > 0) {
							fracct->RemovePair(j);
			       		}
		     		}
	         	
					efn2(&data,fracct);
		    		table[i][2] = fracct->GetEnergy(1) - energy;
	         	
					//restore basepairing
		    		for (j=foldstart;j<=foldstop;j++) {
		      			if (temp[j-foldstart]>0) {
							fracct->SetPair(j,temp[j-foldstart]);
		        		}
		    		}
				}
				else if (suboptimal==3) {
					sn = 0;
					sum = 0;
					//store the basepairing info and force the region of hybridization
					//to have no structure
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
		          		for (j=foldstart;j<=foldstop;j++) {//store basepairing info
		     				temp2[k][j] = fracct->GetPair(j,k);
		           			if (fracct->GetPair(j,k) != 0) {
								fracct->RemovePair(j,k);
		         				
		          			}
		      			}

					}
					efn2(&data,fracct);
					for (k=1;k<=fracct->GetNumberofStructures();k++) {

		          		fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						sn = sn + fn;
		   				sum = sum+fn*(((long double)(fracct->GetEnergy(k)- energyarray[k])));
					}
					//restore the basepairing
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
			       		for (j=foldstart;j<=foldstop;j++) {
							if (temp2[k][j]>0) {
								fracct->SetPair(j,temp2[k][j],k);
           						
           					}
        				}
					}
				
					table[i][2] = (int)(sum/sn);
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 

					delete[] energyarray;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {	delete[] temp2[k];	}
					delete[] temp2;

				}
			}
		}//end mode 1, breaking local structure


		//-------------------------------------------------------------------------------------------
		//mode 2: refold  rna
		else if (mode == 2) {
			//not scan: refolding the whole sequence
			if (foldsize ==0) {
				//reset the nopair constrains
				ct->RemoveConstraints();
				for (j=0;j<length;j++)		{\
					ct->AddSingle(i+j);
				}		
				
				//mode 2: consider only the first suboptimal structure
				if (suboptimal ==0) {
					SAFE_DYNAMIC(ct,&data,1000,10,0,0,true);
				
					table[i][2] = ct->GetEnergy(1) - energy;
				}
				//mode 2: consider all suboptimal structures with ensemble energy
				else if (suboptimal==1) {
					//sum of unconstrained energy from energyarray[]
					Q = (PFPRECISION) 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Q = Q + fn;
         			}
									
					//calculate the whole set of suboptimal structures with constrain
					alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
					//sum of constrained energy now
					Qc = (PFPRECISION) 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ctenergy[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Qc = Qc + fn;
         			}
					//final free energy difference					
					table[i][2] = (int)( conversionfactor*RT_37C*log( (long double)Q/(long double)Qc ) );

					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 
				}
				//mode 2: consider heuristic suboptimal structures with average energy
				else if (suboptimal==3) {
					sum = 0;
					sn = 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)energyarray[k]);
					}
					table[i][2] = (int)(sum/sn);
								
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(savefile,pos);
					strcat(savefile,"_s3_0_constrain.sav");
					//fold the constrained structure
					SAFE_DYNAMIC(ct,&data,1000,10,0);
					if (writesav) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = ct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=ct->GetNumberofStructures();j++)		{
							localint = ct->GetEnergy(j);
							write(&sav, & (localint));
						}
						sav.close();
					}
					sum = 0;
					sn = 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ct->GetEnergy(k)-ct->GetEnergy(1))/(RT_37C*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)ct->GetEnergy(k));
						//cout<< ct->energy[k]<<"\n";
					}
					table[i][2] = (int)(sum/sn) - table[i][2] ;
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 
		   		}
				//mode 2: partionfunction
				else if (suboptimal ==2) {
					//save the energies in a save file
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(pfsfile,pos);
					strcat(pfsfile,"_s2_0_constrain.pfs");
					if (writesav) {
						target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
					}
					else target->refill(ct,&Qc,i, i+length-1,rescaleinrefill);
					
					//table[i][2] = (int)(conversionfactor*RT_37C* 
					//	(PF_LOG_PFPRECISION(PROD(POWER(rescaleinrefill,ct->GetSequenceLength()),DIV(Q, Qc)) ))); //LOG_LINEAR_TAG

					table[i][2] = (int)(conversionfactor*RT_37C* 
						( (long double)(ct->GetSequenceLength())*PF_LOG_PFPRECISION(rescaleinrefill) 
						+ PF_LOG_PFPRECISION(DIV(Q,Qc) ) ));
					
				}
				//mode 2: stochastic
				else if (suboptimal ==4) {
					sum = 0;
					for (k=1;k<=numofstructures;k++) {
	          			sum = sum + energyarray[k];
					}
					table[i][2] = (int)(sum/numofstructures);
				
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d",i);
					strcat(pfsfile,pos);
					strcat(pfsfile,"_s2_0_constrain.pfs");
					std::ifstream pfs(pfsfile,std::ios::binary);
					if(TESTMODE && pfs) {
						pfs.close();
					}
					else{
						pfs.close();
						target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
						
					}
					stochastic(ct,pfsfile,1000,5);
					efn2(&data,ct);
					sum = 0;
					for (k=1;k<=ct->GetNumberofStructures();k++) {
	          			sum = sum + ct->GetEnergy(k);
					}
					table[i][2] = (int)(sum/ct->GetNumberofStructures()) - table[i][2] ;
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=ct->GetNumberofStructures(); 			
				}
				
			}
			//----------------------------------------------------------
			//----------------------------------------------------------
			//scan: only folding the region close to the binding site
			else if (foldsize >0) {
				
				//define the region for refolding 
				if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->GetSequenceLength() ) {
					for (j=1;j<=foldsize+length;j++) {
      					fracct->numseq[j] = ct->numseq[i-foldsize/2+j-1];
					}
				}
				else if ( (i-foldsize/2)<=1 ) 
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[j];
				}
				else if( (i-1+length+foldsize/2)>=ct->GetSequenceLength() )
				for (j=1;j<=foldsize+length;j++) {
      				fracct->numseq[j] = ct->numseq[(ct->GetSequenceLength())+j-foldsize-length];
				}
				fracct->RemoveConstraints();

				//----------------------------------------------------------
				//fold the scanned region without any constraint:
				//refold for different suboptimal options
				if (suboptimal ==0) {//no suboptimals
					SAFE_DYNAMIC(fracct,&data,1,10,0,0,true);
					//fracct->GetNumberofStructures() = 1;//only interested in lowest free energy structure
   					energy = fracct->GetEnergy(1);
					
				}
				else if(suboptimal==1){//use allsub
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s1_%d.sav",i,foldsize+length);
					strcat(savefile,pos);
					std::ifstream sav(savefile,std::ios::binary);
					if (TESTMODE && sav) {
						//read information from file if it was found
						read(&sav,&numofstructures);
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
						sav.close();
					}
					else {
						sav.close();
						alltrace(false,ctenergy,fracct,&data, percent, delta, NULL,NULL);
						numofstructures=fracct->GetNumberofStructures();
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)		energyarray[j]=ctenergy[j];
					}

					
   				}
				else if(suboptimal==3){//use mfold heuristic
					SAFE_DYNAMIC(fracct,&data,1000,10,0);
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s3_%d.sav",i,foldsize+length);
					strcat(savefile,pos);
					if (writesav) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = fracct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=fracct->GetNumberofStructures();j++) {
							int localint = fracct->GetEnergy(j);
							write(&sav, & (localint));

						}
						sav.close();
					}
					numofstructures=fracct->GetNumberofStructures();
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct->GetEnergy(j);
   				}
				else if(suboptimal==2) {//use partition function
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
			
					//not using prefilter, so arrays can be reused when region move to the right
					//fold the first region, the folded region begin to move to right in the middle of target
					if (prefilter->useit == 0) {
						if (i==start) {
							target=new OligoPclass(fracct,pfdata);
							targetcopy=new OligoPclass(fracct,pfdata);
							target->partition(true,&Q,NULL);
							//char *report1="report1.out";				
							//target->partition(true,&Q,NULL,report1);
						}
						//when folded region begin moving,reuse some arrays overlapped expect for those on the edges
						else if((i-foldsize/2)>1 && (i+length-1+foldsize/2)<= (ct->GetSequenceLength()) ){
				   		//char *report2="report2.out";		
						target->scanfill(fracct,&Q,0);
						}
						//copy the arrays to be reused for next scan region
						//folded region is not moving at two ends, so copy the array outside for next folding without constrain
						if ( (i-foldsize/2)<1 || (i+length-1+foldsize/2)>=(ct->GetSequenceLength()) ) {
							scancopyend(target,targetcopy);
						}
						//folded region begin moving next, copy the overlapped region outside with different index
						else	scancopy(target,targetcopy);
				
					}
					//using prefilter, not reuse arrays,refold the new region every time
					else {
						
						//refold the new region without using any information from previous folding 
						target->reset4oligo(fracct);
						if (writesav) {
							target->partition(true,&Q,NULL,pfsfile);
							
						}
						else target->partition(true,&Q,NULL);
					}
				}
				else if(suboptimal==4) {//use stochastic sampling
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					//not using prefilter, so arrays can be reused when region move to the right
					//fold the first region, the folded region begin to move to right in the middle of target
					if (prefilter->useit == 0) {
						std::cout <<"Must use prefilter!\n";
						return;
					}
					//using prefilter, not reuse arrays,refold the new region every time
					else {
						std::ifstream pfs(pfsfile,std::ios::binary);
						if(TESTMODE && pfs) {
							pfs.close();
						}
						else{
							pfs.close();
							//refold the new region without using any information from previous folding 
							target->reset4oligo(fracct);
							target->partition(true,&Q,NULL,pfsfile);
							
						}
						stochastic(fracct,pfsfile,1000,5);
						efn2(&data,fracct);
					/*	ofstream sout("RL_sample_oligo_scan.out");
						for (j=1;j<=1000;j++)	{ sout << fracct->energy[j]<<"\n";	}
						sout.close();
					*/
						numofstructures=fracct->GetNumberofStructures();
						energyarray = new int [numofstructures+1];
						for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct->GetEnergy(j);
					}
				}

				
				
				
				//set the single-stranded constrain 
				//fracct->nnopair=length;
				//fracct->checknopair();
				fracct->RemoveConstraints();
				if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
					for (j=0;j<length;j++) {
						fracct->AddSingle(i+j); 
					}
					foldstart=i;
					foldstop=i+length-1;
				} 
				else if(i+length-1+foldsize/2>= (ct->GetSequenceLength()) ) {
					for (j=0;j<length;j++) {
						fracct->AddSingle(foldsize+length-(ct->GetSequenceLength())+i+j);
						//fracct->nopair[j+1] = foldsize+length-(ct->GetSequenceLength())+i+j;
					}
					foldstart=foldsize+length-(ct->GetSequenceLength()) +i;
					foldstop=foldsize+length+length-1-(ct->GetSequenceLength()) +i;
				}
				else {//folded region begin to move, constrained position will be always in the middle of this region
					for (j=0;j<length;j++) {
						fracct->AddSingle(foldsize/2+1+j);
						//fracct->nopair[j+1] = foldsize/2+1+j;
					}
					foldstart=foldsize/2+1;
					foldstop=foldsize/2+length;
				}

				//---------------------------------------------------------------------
				//refold with constraint:
				//mode 2 + refold + consider only the first suboptimal structure
				if (suboptimal==0) {
					SAFE_DYNAMIC(fracct,&data,1000,10,0,0,true);
					table[i][2] = fracct->GetEnergy(1) - energy;			
				}
				//---------------------------------------------------------------------------------------------
				//mode 2 + refold + consider all suboptimal structures ensemble energy
				else if (suboptimal==1) {
					
					//sum of unconstrained energy from energyarray[]
					Q = (PFPRECISION) 0;
					for (k=1;k<=numofstructures;k++) {
	          			fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Q = Q + fn;
         			}
					//calculate the whole set of suboptimal structures with constrain
					alltrace(false,ctenergy,fracct,&data, percent, delta, NULL,NULL);
					//sum of constrained energy now
					Qc = (PFPRECISION) 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          			fn = exp(-(((long double)ctenergy[k]-energyarray[1])/(RT_37C*conversionfactor)));
						Qc = Qc + fn;
         			}
					//final free energy difference					
					table[i][2] = (int)( conversionfactor*RT_37C*log( (long double)Q/(long double)Qc ) );

					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 
									      			
					delete[] energyarray;
				}
				//---------------------------------------------------------------------------------------------
				//mode 2 + refold + consider heuristic suboptimal structures average free energy
				else if (suboptimal==3) {
      				sum = 0;
					sn = 0;
					for (k=1;k<=numofstructures;k++) {
	          				fn = exp(-(((long double)energyarray[k]-energyarray[1])/(RT_37C*conversionfactor)));
							sn = sn + fn;
         					sum = sum + fn*((long double)(energyarray[k]));
					}
					table[i][2] = (int)(sum/sn);
					delete[] energyarray;

					//SAFE_DYNAMIC(fracct,&data,1000,10,0,0,true);
      				SAFE_DYNAMIC(fracct,&data,1000,10,0);
					//save the energies in a save file
					strcpy(savefile, ct->GetCtLabel(1).c_str());
					savefile[strlen(savefile)-1]='\0'; 
					sprintf(pos,"_%d_s3_%d_constrain.sav",i,foldsize+length);
					strcat(savefile,pos);
					if (writesav) {
						std::ofstream sav(savefile,std::ios::binary);
						int localint = fracct->GetNumberofStructures();
						write(&sav,&(localint));
						for (j=1;j<=fracct->GetNumberofStructures();j++)		{
							int localint=fracct->GetEnergy(j);
							write(&sav, & (localint));
						}
						sav.close();
					}
					
   					sum = 0;
					sn = 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          				fn = exp(-(((long double)fracct->GetEnergy(k)-fracct->GetEnergy(1))/(RT_37C*conversionfactor)));
							sn = sn + fn;
         					sum = sum + fn*((long double)fracct->GetEnergy(k));
							//cout<< fracct->energy[k]<<"\n";
					}
					table[i][2] = (int)(sum/sn) - table[i][2];
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 
				}
				//---------------------------------------------------------------------------------------------
				//mode 2 + refold + partionfunction
				else if (suboptimal==2) {
					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					//reuse the arrays filled without constrained 
					if (writesav) {
						target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);

					}
					else target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill);
					//exchange arrays to be used for next folding site without constrain
					if (prefilter->useit==0) {//not using targetcopy when prefilter is used
					targettemp=target;
					target=targetcopy;
					targetcopy=targettemp;
					}
					table[i][2] = (int)(conversionfactor*RT_37C*
					( (long double)(fracct->GetSequenceLength())*PF_LOG_PFPRECISION(rescaleinrefill) 
						+ PF_LOG_PFPRECISION(DIV(Q,Qc) ) ));
				}
				//---------------------------------------------------------------------------------------------
				//mode 2 + refold + stochastic sample
				else if (suboptimal==4) {			
					
   					sum = 0;
					for (k=1;k<=numofstructures;k++) {
	          				sum = sum + energyarray[k];
					}
					table[i][2] = (int)(sum/numofstructures);
					delete[] energyarray;

					strcpy(pfsfile, ct->GetCtLabel(1).c_str());
					pfsfile[strlen(pfsfile)-1]='\0'; 
					sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
					strcat(pfsfile,pos);
					std::ifstream pfs(pfsfile,std::ios::binary);
					if(TESTMODE && pfs) {
						pfs.close();
					}
					else{
						pfs.close();
						//reuse the arrays filled without constrained 
						target->scanconstrain(fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);
						
					}
					stochastic(fracct,pfsfile,1000,5);
					efn2(&data,fracct);
					sum = 0;
					for (k=1;k<=fracct->GetNumberofStructures();k++) {
	          			sum = sum + fracct->GetEnergy(k);
							
					}
					table[i][2] = (int)(sum/fracct->GetNumberofStructures()) - table[i][2];
					numofsubstructures[i][0]=numofstructures; //report in report() function
					numofsubstructures[i][1]=fracct->GetNumberofStructures(); 

					//exchange arrays to be used for next folding site without constrain
					if (prefilter->useit==0) {//not using targetcopy when prefilter is used
					/*targettemp=target;
					target=targetcopy;
					targetcopy=targettemp;
					*/
						std::cout << "must use prefilter for stochasitic folding\n";
						return;
					}
				}
			
			}
		}//end mode 2

		//-----------------------------------------------------------------------------------------------------
		//mode 3: no local structure considered

		else {
			table[i][2] = 0;
		}

		
	
		//-----------------------------------------------------------------------------------------------------
		//-----------------------------------------------------------------------------------------------------
		//calculate the stability of the hybrid duplex
		
		//oligo is DNA:
		if (isdna) {
			
			table[i][1] = hybriddata->init;//initiation
			for (j=0;j<(length-1);j++) {
			
				table[i][1] += 
					hybriddata->stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			}
		}

		//oligo is RNA:
		else {
			table[i][1] = data.init;//initiation
      		for (j=0;j<(length-1);j++) {
          		table[i][1] += //LOG_LINEAR_TAG
					data.stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			 }
			//consider AU end effects for RNA/RNA duplexes
			if (ct->numseq[i]==1||ct->numseq[i]==4) table[i][1]+=data.auend;  //LOG_LINEAR_TAG
			if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) table[i][1]+=data.auend;  //LOG_LINEAR_TAG

		}

		//calculate the Tm of the duplex
		ds = helixstack->dsi;
		dh = helixstack->dhi;
		for (j=0;j<(length-1);j++) {
	  		dh +=helixstack->dh[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
			ds +=helixstack->ds[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
   		}
		if (ct->numseq[i]==1||ct->numseq[i]==4) {
			dh = dh + helixstack->dha;
			ds = ds + helixstack->dsa;
		}
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
			dh = dh + helixstack->dha;
			ds = ds + helixstack->dsa;
		}
	      
		table[i][5] = (int) ((conversionfactor*((double)dh*1000)/
								((double)ds+conversionfactor*Rgas*log(conc)))-273.15*conversionfactor);  //LOG_LINEAR_TAG



		//-----------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------
		//calculate the free energy of intramolecular folding of oligo

	    
		if (suboptimal==2) { //-2 is unavailabe, so not using partition function as interoligo cannot use it
			//reuse some arrays when scanning along the sequence for intramolecule
			//only calculate the partition function for the first one
			//cannot reuse array when prefilter is used
			if (prefilter==0) {
				if (i==start) {	
					intraoligo->reset4oligo(oligo1);
					intraoligo->partition(true,&Qolig);
				
				}
				//reuse some arrays overlapped, only change the index of left binding site
				//be careful that the oligo is complementary to the target, copy direction is reversed
				else {
				
					for (ip=length-1;ip>2;ip--) {
						for (jp=length-1;jp>=ip;jp--) {
								
							intraoligo->wca[ip][jp]=intraoligo->wca[ip-1][jp-1];
							intraoligo->w->f(ip,jp)=intraoligo->w->f(ip-1,jp-1);
							intraoligo->v->f(ip,jp)=intraoligo->v->f(ip-1,jp-1);
							intraoligo->wmb->f(ip,jp)=intraoligo->wmb->f(ip-1,jp-1);
							intraoligo->wl->f(ip,jp)=intraoligo->wl->f(ip-1,jp-1);
							intraoligo->wmbl->f(ip,jp)=intraoligo->wmbl->f(ip-1,jp-1);
							intraoligo->wcoax->f(ip,jp)=intraoligo->wcoax->f(ip-1,jp-1);
						}
					}
					//oligo is complementary to target, set reverse=1 for scanfill()
					intraoligo->scanfill(oligo1,&Qolig,1);

				}
			}
			else {
				intraoligo->reset4oligo(oligo1);
				intraoligo->partition(true,&Qolig);
			}
			oligo1->SetEnergy(1,(int)( conversionfactor*RT_37C*
				((long double)length*PF_LOG_PFPRECISION( (intraoligo->data)->scaling) - PF_LOG_PFPRECISION(Qolig)) ));//LOG_LINEAR_TAG

		}
		else {//use the lowest free energy for suboptimal 0 and 1 since oligo is small
    		if (isdna)
    			SAFE_DYNAMIC(oligo1,&ddata,1000,10,0,0,true);
			else
				SAFE_DYNAMIC(oligo1,&data,1000,10,0,0,true);
		}
	    	
		table[i][3] = oligo1->GetEnergy(1); //LOG_LINEAR_TAG
		//if the intra structure is unfavorable , not consider it in the total energy
		if (table[i][3]>0) {
			table[i][3]=0;
		}

		//-----------------------------------------------------------------------------------------
		//calculate the free energy of intermolecular folding of oligo
		oligo->intermolecular = true;
		//oligo->GetSequenceLength() = 2*length + 3;
		//if (i==start) oligo->allocate(2*length + 3);
		for (j=1;j<=length;j++) {
			oligo->numseq[j+length+3] = oligo->numseq[j];
		}
		for (j=1;j<=3;j++) {
			oligo->numseq[length+j] = 5;
		}
	    
		if (suboptimal==2) {	//-2 is unavailable, as partition function is not available for inter oligo	
				
			interoligo->reset4oligo(oligo);
			interoligo->partition(true,&Qolig);
					
			oligo->SetEnergy(1,(int)( conversionfactor*RT_37C*
				( (long double)(oligo->GetSequenceLength())*PF_LOG_PFPRECISION((interoligo->data)->scaling)  - PF_LOG_PFPRECISION(Qolig) ) )) ;
		}
		else {
			if (isdna)
				SAFE_DYNAMIC(oligo,&ddata,1000,10,0,0,true);
			else
				SAFE_DYNAMIC(oligo,&data,1000,10,0,0,true);
		}
			      	
		table[i][4] = oligo->GetEnergy(1); //LOG_LINEAR_TAG
		//if the intra structure is unfavorable , do not consider it in the total energy
		if (table[i][4]>0) {
			table[i][4]=0;
		}
		


		//-----------------------------------------------------------------------------------------
		//-----------------------------------------------------------------------------------------
		//calculate the overall free energy for the binding
		k1b = exp(((-1)*(long double)(table[i][4]))/((1.9872)*(3.1)));
		k1u = exp(((-1)*(long double)(table[i][3]))/((1.9872)*(3.1)));

		//I don't know what's going on here.:)  ---John , Nov.9,2005
		//if (i==867) {
		//	table[0][0]=56;
		//}

		if ((k1u/(conc*k1b))<100) {

      		dgeff = (int)(-((1.9872)*(3.1)*log(((4.0*k1b*(long double)(conc))/(-1-k1u+sqrt(pow((long double)1.0+k1u,(long double) 2.0)
							+8.0*k1b*(long double)(conc))))-1.0)));
		}
		else {
			dgeff = table[i][3];

		}
		//this line converts the free energy to the convention explained in
		//the written version of the algorithm:
		table[i][2] = -table[i][2]; //LOG_LINEAR_TAG

		if ((table[i][2]<0)&&(dgeff<0)) {
	   		table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log( (exp(-(long double)(table[i][2])/ //LOG_LINEAR_TAG
						  (RT_37C*conversionfactor))+1.0)*  (exp(-(long double)(dgeff)/(RT_37C*conversionfactor))+1.0) )+0.5);
		}
		else if (table[i][2]<0) {
			table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log( (exp(-(long double)(table[i][2])/
						  (RT_37C*conversionfactor))+1.0))+0.5);
		}
		else {
			table[i][0] = table[i][1] + (int)((RT_37C*conversionfactor)*log(  (exp(-(long double)(dgeff)/
						  (RT_37C*conversionfactor))+1.0) )+0.5);
		}

		//table[i][0] = table[i][1] + table[i][2] - dgeff;

	}//end of main loop over all oligonucleotides
	//The scan is finished now
	//-------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------  



	//clean up the ct file
	ct->SetEnergy(1,energy) ;

	//clean up memory use
	if (mode==1) {
		if (suboptimal==3 && foldsize==0) {
			ct->RemoveConstraints();//nnopair=0;
      		for (i=1;i<=ct->GetNumberofStructures();i++) delete[] temp2[i];
			delete[] temp2;
			delete[] energyarray;
		}
		else if (suboptimal==0) delete[] temp;
		
	}
    else if (mode ==2) {
   		ct->RemoveConstraints();//nnopair=0;
   		for (j=1;j<=ct->GetSequenceLength();j++) {
    		ct->SetPair(j,temp[j],1);//basepr[1][j] = temp[j];
		}
		delete[] temp;
		if (suboptimal==1)	delete[] ctenergy;
		if ( (suboptimal==1||suboptimal==3||suboptimal==4) && foldsize==0 )	delete[] energyarray;
		else if (suboptimal==2 || suboptimal ==4) {
			if(isdna)	delete dpfdata;
			delete pfdata;
			delete target;
			if (foldsize>0 && prefilter->useit == 0)	delete targetcopy;
			if (suboptimal==2) {
				delete intraoligo ;
				delete interoligo ;
			}
		}
	}

	if (foldsize > 0) delete fracct;

	delete oligo;
	delete oligo1;
//cout <<"\n"<< numofstructures<<"\n";	
}


//=======================================================================
int readrd (rddata* data, const string &dnarnafile) {
	int count,i,k,j,l;
	std::ifstream dr(dnarnafile.c_str());
	char lineoftext[100];

	//make sure the file exists
	if (!dr.good()) return 0;

	/* Read info from stackdr */
	//add to the stack table the case where X (represented as 0) is looked up:
	for (count=1;count<=2;count++) dr >> lineoftext;//get past text in file
	dr >> lineoftext;
	data->init =(int)floor(conversionfactor*(atof(lineoftext)));

	for (count=1;count<=42;count++) dr >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (i!=0) for (count=1;count<=60;count++) dr >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->stack[i][j][k][l]=0;
					}
					else {
						dr >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->stack[i][j][k][l] =(int)floor(conversionfactor*(atof(lineoftext))+.5);
						}
						else data->stack[i][j][k][l] = INFINITE_ENERGY;
					}
					//cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->stack[i][j][k][l]<<"\n";
				}
				//cin >> m;
			}
		}
	}
	return 1;
}



//=======================================================================
//return the numerical equivalent of the complementary base to nucleotide i
int complement(const int i, structure *ct) {
	if (ct->numseq[i] == 0) return 0;
	return 5 - ct->numseq[i];
}

char numtobase(const int basenum, structure *ct, const bool isDNA) {
	char base = ct->GetThermodynamicDataTable()->numtobase(basenum);
	if (base=='U' && isDNA) return 'T';
	return base;
}

//=======================================================================
//This is a older post-filter of siRNA 
//use siRNA selection criteria to filter the output
//mask will contain true for those that meet the criteria
void filterbysirna ( structure *ct, int **table, int length, datatable *data, 
				     bool *mask, double asuf, double tofe, double fnnfe) {
	
	int iasuf,itofe,ifnnfe;
	int i,j,*k;
	k = new int [length];

	iasuf = (int) (conversionfactor*asuf);
	itofe = (int) (conversionfactor*tofe);
	ifnnfe = (int) (conversionfactor*fnnfe);


	for (i=1;i<=(ct->GetSequenceLength()-length+1); i++) {
		//filter all oligos
		mask[i]=true;
		if (iasuf>table[i][3]) {
			mask[i]=false;
		}
		if (itofe>table[i][2]) {
			mask[i]=false;
		}

		//filter out those that have a repeat of A, G, or U of more than 4. 
		if (length>4) {
			for (j=length-1;j>=0;j--) {
   				k[j] = complement(i+j,ct);
      		
			}
			for (j=0;j<length-3;j++) {
				if (k[j]==1) {
					if (k[j+1]==1&&k[j+2]==1&&k[j+3]==1)	mask[i]=false;
				}
				else if (k[j]==3) {
					if (k[j+1]==3&&k[j+2]==3&&k[j+3]==3)	mask[i]=false;
				}
				else if (k[j]==4) {
					if (k[j+1]==4&&k[j+2]==4&&k[j+3]==4)	mask[i]=false;
				}
			}
		}
		table[i][6]= data->stack[complement(i+length-1,ct)][ct->numseq[i+length-1]]
								[complement(i+length-2,ct)][ct->numseq[i+length-2]];
		//account for change in AU end if necessary:
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
			if (ct->numseq[i+length-1]==2||ct->numseq[i+length-1]==3) {
				table[i][6]+=data->auend;
			}
		}
		else {
			if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
				table[i][6]-=data->auend;
			}
		}
		if (ifnnfe>table[i][6])
			mask[i]=false;

	}

	delete[] k;
}

/*
 This function writes a tab delimited (or HTML) file with the oligowalk data
 Parameters:
 out			-- The output stream (e.g. to an output file).
 ct			-- Contains the sequence
 table		-- The table filled by OligoWalk
 numofsubstructure	-- 
 length		-- The length of the oligonucleotides
 isdna		-- Whether the oligos are dna (true = DNA, false = RNA)
 conc		-- The oligonucleotide concentration
 suboptimal	-- Whether suboptimal structures were used
 start		-- The nucleotide postition of the start of the walk on the target sequence
 stop		-- The nucleotide postition of the end of the walk on the target sequence
 prefilter	-- John Lu's siRNA prefilter
 foldsize	-- The size of the folding region that is centered around the oligo
 mask		-- A bool array in which mask[i] is true if the oligo at table[i] has passed 
 			   the filter criteria (i.e. asuf, tofe, fnnfe).
 asuf		-- Minimum Antisense strand unimolecular folding free energy
 tofe		-- Minimum Target opening free energy
 fnnfe		-- Minimum First nearest neighbor free energy
 isHTML		-- Whether this is HTML or tab delimited text (false = tab delimited, true = HTML)
 writeHeader	-- Whether header information should be written in the report.
 writeBody	-- Whether the report body should be written (as opposed to just the header).
*/
void report(ostream& out, 
	const bool isdna, const int mode, structure *ct, const int length, const double conc, const int suboptimal,
	const int start, const int stop, const int foldsize, const int distance,
	int **table, int **numofsubstructures, const char * const shapefile,
	const siPREFILTER *prefilter,
	const bool *mask, const double asuf, const double tofe, const double fnnfe,
	const bool isHTML, const bool writeHeader, const bool writeBody) {
	
	int i,j,k;
	const double scale = conversionfactor;

	FixWindowsExponentFormat();

	//define some constants based on whether or not html output is desired.
	const char *H1, *H1c, *H3, *H3c, *TR, *TRc, *TD, *TDc, *HR, *BR, *TABLE, *TABLEc, *INDENT; // XXc is a closing tag for XX;
	if (isHTML) {
		H1="<h1>";H1c="</h1>"; H3="<h3>";H3c="</h3>";
		TR="<tr>";TRc="</tr>"; TD="<td>";TDc="</td>";
		TABLE="<table>";TABLEc="</table>";
		HR="<hr>";BR="<br>";
		INDENT="&nbsp;&nbsp;&nbsp;&nbsp;";
	} else {
		H1=H1c=H3=H3c="";
		TR="";TRc="";TD="";TDc="\t";
		TABLE=TABLEc="";
		HR="--------------------------------------------------------------------------------";
		BR="";
		INDENT="    ";
	}

	if (writeHeader) {
		out << H1 << "RNAstructure OligoWalk calculation:" << H1c << endl;
		out << HR << endl;
		out << "Target sequence: " << ct->GetSequenceLabel() << BR << endl;
		out << "Total size of the target: " << ct->GetSequenceLength() << " nucleotides" << BR << endl;
		out << "Scanned position on target: " << start << " to " << stop << " (nt)" << BR << endl;
		out << BR << endl;
		out << "Oligonucleotides:" << BR << endl;
		out << INDENT << "Type:  " << (isdna?"DNA":"RNA") << BR << endl;
		out << INDENT << "Length: " << length << BR << endl;
		out << INDENT << "Concentration: " << conc << " M" << BR << endl;
		out << BR << endl;
		out << "Method options:" << BR << endl;
		out << INDENT << "Mode: ";
		switch(mode) {
			case 1: out << "Break local structure of the target to bind oligo."; break;
			case 2: out << "Refold target structure after oligo binding."; break;
			case 3: out << "Not considering the structure of target."; break;
		}
		out << BR << endl;
		if (mode ==2 || mode==1) {	
			out << INDENT << "Folding region size: ";
			if (foldsize==0) 
				out << "(Global Target Region)";
			else
				out << (foldsize+length) << " nt";
		}
		out << BR << endl;
		out << INDENT << "Suboptimal: ";
		switch(suboptimal) {
			case 0: out << "None -- Only one (optimal) structure was considered."; break;
			case 1: out << "All Suboptimal structures within an energy difference window were considered."; break;
			case 2: out << "All possible structures were considered using the Partition Function."; break; 
			case 3: out << "Heuristic suboptimal structures were considered for both oligo-free and olig-bound target." 
						<< BR << endl << INDENT << INDENT 
						<< "There are at most 1000 suboptimal structures (free energy within 10% of optimal structure)."; break; 
			case 4: out << "Considered 1000 stochastically sampled structuresfor both oligo-free and olig-bound target."; break; 
		}
		out << BR << endl;

		if (distance >0)			out << "The base pairs with distance larger than " << distance << " nt are not allowed." << BR << endl;
		if (prefilter->useit != 0)  out << "Prefilter was used." << BR << endl;
		if (!is_blank(shapefile))  	out << "The shape information was used: " << shapefile << "." << BR << endl;

		if (mask!=NULL) {
			out << "Oligonucleotides filtered by siRNA selection criteria: " << BR << endl;
			out << INDENT << "Antisense strand unimolecular folding free energy >= "<< asuf << BR << endl;
			out << INDENT << "Target opening free energy >= "<< tofe << BR << endl;
			out << INDENT << "First nearest neighbor free energy >= "<< fnnfe << BR << endl;
			out << INDENT << "Sequences with more than three G's, U's, or A's in a row removed." << BR << endl;
			out << "Antisense strand shown 5' to 3'." << BR << endl;
		}
		out << HR << endl;	
	}

	if (writeBody) {
		//-------------Table Header-------------
		out << BR << endl;
	 	out << H3 << "Energy table:" << H3c << BR << endl;
		out << TABLE << endl;
		out << TR; //start a new row
		out << TD << "Pos." << TDc;
		out << TD << "Oligo(5\'->3\')" << TDc;
		out << TD << "Overall (kcal/mol)" << TDc;
		out << TD << "Duplex (kcal/mol)" << TDc;
		out << TD << "Tm-Dup (degC)" << TDc;
		out << TD << "Break-Target (kcal/mol)" << TDc;
		out << TD << "Intra-oligo (kcal/mol)" << TDc;
		out << TD << "Inter-oligo (kcal/mol)" << TDc;
		if (prefilter->useit) {
			out << TD << "End Diff (kcal/mol)" << TDc;
			out << TD << "Prefilter Score" << TDc;
		}
		if (suboptimal==1||suboptimal==3||suboptimal==4) {
			out << TD << "Structures #" << TDc;
			out << TD << "Constrained structures #" << TDc;
		}
		if (mask!=NULL) out << TD << "First NN Opening (kcal/mol)" << TDc;
		out << TRc << endl; //end of row

		//-------------Table Rows-------------
		for (i=start;i<=stop; i++) {
			// Skip if it doesn't pass the prefilter.
			if (prefilter->useit && prefilter->score[i] < FILTER_PASS)
				continue;
			// Skip if it doesn't pass the selection criteria.
			if (mask!=NULL && !mask[i]) 
				continue;

			out << TR; //start a new row
			// row index
			out << TD << i << TDc;
			// oligo sequence
			out << TD;
			for (j=length-1;j>=0;j--) out << numtobase(complement(i+j,ct),ct, isdna);
			out << TDc;

			// calculated energies
			int *row=table[i];
			out << TD << row[0]/scale << TDc; // overall
			out << TD << row[1]/scale << TDc; // duplex
			out << TD << row[5]/scale << TDc; // Tm
			out << TD << row[2]/scale << TDc; // breaking target structure
			out << TD << row[3]/scale << TDc; // intramolecular oligo free energy
			out << TD << row[4]/scale << TDc; // inter-molecular oligo free energy
			
			// prefilter
			if (prefilter->useit) {
				out << TD <<  prefilter->enddiff[i] << TDc;
				out << TD <<  prefilter->score[i] << TDc;
			}
			// suboptimal
			if (suboptimal==1||suboptimal==3||suboptimal==4)	 {
				out << TD << numofsubstructures[i][0]  << TDc;
				out << TD << numofsubstructures[i][1]  << TDc;
			}
			// First NN Opening
			if (mask!=NULL) {
				out << TD << row[6]/scale <<  TDc;
			}
			out << TRc << endl;
		}
		out << TABLEc << endl;
	}
}


//=================================================================================================
//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	//only copy the overlapped region which is in the middle of sequence
	for (i=2;i<=(copyregion->number)-2;i++) {
		for (j=i;j<=(copyregion->number)-2;j++) {
			//w5 is not copied as it will be recalc when i==1						
			copyregion->wca[i][j]=region->wca[i+1][j+1];
			copyregion->w->f(i,j)=region->w->f(i+1,j+1);
			copyregion->v->f(i,j)=region->v->f(i+1,j+1);
			copyregion->wmb->f(i,j)=region->wmb->f(i+1,j+1);
			copyregion->wl->f(i,j)=region->wl->f(i+1,j+1);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i+1,j+1);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i+1,j+1);
			
			

		}
	}
}
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	for (i=1;i<=copyregion->number;i++) {
		for (j=i;j<=copyregion->number;j++) {
						
			if(i==1)	copyregion->w5[j]=region->w5[j];
			copyregion->wca[i][j]=region->wca[i][j];
			copyregion->w->f(i,j)=region->w->f(i,j);
			copyregion->v->f(i,j)=region->v->f(i,j);
			copyregion->wmb->f(i,j)=region->wmb->f(i,j);
			copyregion->wl->f(i,j)=region->wl->f(i,j);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i,j);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i,j);
		
		}
	}
}

