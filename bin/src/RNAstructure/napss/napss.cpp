/*
 * NAPSS, a program that predicts RNA secondary structures with pseudoknots with the aid of NMR constraints data.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by James Hart, David Mathews, Jon Chen, Stanislav Bellaousov
 *
 * Modified from the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 * Modified to parse command line interface by Jessica Reuter and Stanislav Bellaousov 2012
 *
 * Modified to read Triplet Constraints by Stanislav Bellaousov in August 2012
 *
 * Updated to new structure class interface by Richard M. Watson in 2017
 *     Also split napss.cpp into a text interface (NAPSS_Interface.cpp) and core (src/napss.cpp)
 *                          
 */

#include "napss.h"
#include "../src/algorithm.h"

// Flags for text output
#undef VERBOSE_MODE
//#define VERBOSE_MODE

// Flags that output dotplots
//#define DOTPLOTS_OUT
#undef DOTPLOTS_OUT

// Flags for output associated with Triplet Constraints
#undef TRIPLET_VERBOSE_MODE
//#define TRIPLET_VERBOSE_MODE

// Flags for debug features
#undef DEBUG_MODE
//#define DEBUG_MODE
#undef DEBUG2_MODE
//#define DEBUG2_MODE
#undef ALREADY_USED_CHECK
//#define ALREADY_USED_CHECK

//#define OUTPUT_MATCHES
#undef OUTPUT_MATCHES


//Define the ProbKnot mode, where the pairs are populated using ProbKnot method instead of dynamic
//#define PROBKNOT_MODE
#undef PROBKNOT_MODE

//Define to enable removing matches with broken pairs after "RemoveIsolatedPairs"
#define RemoveBrokenMatch 
//#undef RemoveBrokenMatch 

//Define to enables helix extension before the folding
#define ExtendBefore
//#undef ExtendBefore

//Defined to enable helix extension after the folding
#define ExtensionAfter
//#undef ExtensionAfter

//Define to allow matched pairs during 'dynamic'
#define AllowMatchedPairs
//#undef AllowMatchedPairs

//Define to enable debug execution of the code, which doesn't sort structures or doesn't remove
//duplicate structures. It also outputs the matches that produce lowest energy structures
//#define DEBUG_MATCHES
#undef DEBUG_MATCHES

int napss(RNA *rnaCT, const string inNMRconstraints, 
		/* pass in a NULL or uninitialized structure* it will be set to the address of the output structure on success. The caller must then delete the returned structure. */
		structure* &outputStructureResults, 
		const int windowsize, int cutoff, const int percent, 
		const int maxtracebacks, const int warningLimit,
		double pseud_param[16], const double P1, const double P2, const bool pseudoknotFree,
		const string inSHAPEfile, const double slope, const double intercept,
		const string constraintFile,
		const bool useSimpleMBRules) {
	
	structure *ct   = rnaCT->GetStructure();
	datatable *data = rnaCT->GetDatatable();
	const char* const sequence = rnaCT->GetSequence();

	int i,j,k,iter,iter2;
	string a;
	int errorcode = 0;
	bool ifwarningMessage=false;//keeps track if the warning message was printed. false - warning message was not printed yet, so should be printed when warning limit is reached; true - warning message has already been printed, so no need to print it again.

#if defined(VERBOSE_MODE)
	// Display basic info
	cout << "\n\n";
	for (i = 1; i <= ct->GetSequenceLength(); i++) cout << ct->nucs[i];
	cout << "\n";
	for (i = 1; i <= ct->GetSequenceLength(); i++) cout << ct->numseq[i];
	cout << "\nLength of sequence: " << ct->GetSequenceLength() << "\n";
#endif

	// Create a basepair type lookup matrix (AU = 5, GC = 6, GU = 7)
	short bpLookup[5][5] = {{0,0,0,0,0},
	                        {0,0,0,0,5},
	                        {0,0,0,6,0},
	                        {0,0,6,0,7},
	                        {0,5,0,7,0}};

	// Create 2D matrices to store converted dotplot, DeltaG, and vmb/vext values
	short** convertedDotplot;
	short** dgArray;
	short** mbDotplot;
	convertedDotplot = new short* [(ct->GetSequenceLength())+1];
	dgArray = new short* [(ct->GetSequenceLength())+1];
	mbDotplot = new short* [(ct->GetSequenceLength())+1];	
	for(i = 0; i <= ct->GetSequenceLength(); i++) {
		convertedDotplot[i] = new short[(ct->GetSequenceLength())+1];
		dgArray[i] = new short[(ct->GetSequenceLength())+1];
		mbDotplot[i] = new short[(ct->GetSequenceLength())+1];
	}
	
	// Clear the values stored in the matrices
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		for (j = 0; j <= ct->GetSequenceLength(); j++) {
			convertedDotplot[i][j] = 0;
			dgArray[i][j] = 0;
			mbDotplot[i][j] = 0;
		}
	}

	//allocate space for v, vmb, and vext arrays:
	DynProgArray<integersize> v(ct->GetSequenceLength());
	DynProgArray<integersize> vmb(ct->GetSequenceLength());
	DynProgArray<integersize> vext(ct->GetSequenceLength());
	double pairs;
	double mbpairs;
	double extpairs;
	int dgMin = 0;

	// Load NMR constraints
	short maxConLength = 0, numOfCon = 0, totalConLength = 0;
	
	ifstream inNMRConFile(inNMRconstraints.c_str());
	if (!inNMRConFile.good()) return NAPSS_BAD_CONSTRAINT_FILE;
	string s;

	// Loop through each line of the file to get count of constraints and maximum length of a constraint
	while (getline(inNMRConFile,s)) {
		if(!s.empty()&&s.length() < 2) return NAPSS_ERR_CONSTRAINT_LEN;
		if(!s.empty()){//make sure not to read empty lines
			numOfCon++;
			//Counting the constraints omitting the constraints in parentheses
			const char *TempRead = s.c_str();//Read the string into character array
			short TempConLength=0;//Re-Set the constraint length counter
			while(*TempRead != '\0'){//Check if it is the end of array. If not, execute while loop
				if(*TempRead++ != '(')++TempConLength;//if not the parenthesis, count the constraints
				else{//if find parenthesis
					while(*TempRead++ != ')');//skip till the closing parenthesis
				}
			}
#if defined(TRIPLET_VERBOSE_MODE)
			cerr << "-------------\nTempConLength = " << TempConLength << endl;	//debug	
#endif
			totalConLength+=TempConLength;//add to the total constraint length counter
			if(TempConLength > maxConLength) maxConLength = TempConLength;//find the longest constraint
			
		}//END:Making sure not to read empty lines
	}//END:while loop that reads file line by line
	inNMRConFile.close();
	inNMRConFile.clear();
	
#if defined(VERBOSE_MODE)
	cout << numOfCon << " NMR constraints\n";
	cout << "Total length of NMR constraints: " << totalConLength << "\n";
	cout << "Length of largest NMR constraint: " << maxConLength << "\n\n" << flush;
#endif


	// Allocate 2D array for storing constraints (first column contains length of each constraint)
	// Also allocate two more copies of this array which will be used later for storing temporary
	// dotplot coordinates
	char**** tripletArray;//Initialize a 4D array to hold the triplet constraints. 1D = non-triplet constraint number;
	                      //2D = constraint position; 3D = # of triplet constraint (could be 1 or 2);
	                      //4D = the triplet constraint, ex: +RGY
	short** conArray;
	short** xCoords;
	short** yCoords;
	tripletArray = new char*** [numOfCon+1];
	conArray = new short* [numOfCon+1];
	xCoords = new short* [numOfCon+1];
	yCoords = new short* [numOfCon+1];
	for(i = 0; i <= numOfCon; i++){
		tripletArray[i] = new char**[(maxConLength+1)];
		for(int j=0;j<=maxConLength;j++){
			tripletArray[i][j] = new char*[2];
			for(int k=0;k<2;k++){
				tripletArray[i][j][k] = new char[4];
				for(int m=0;m<4;m++) tripletArray[i][j][k][m] = 0;
			}
		}
		*(conArray + i) = new short[(maxConLength+1)];
		*(xCoords + i) = new short[(maxConLength+1)];
		*(yCoords + i) = new short[(maxConLength+1)];
	}

	// Re-open NMR constraints file and fill the 2D arrays (only fill Coords arrays with 0's)
	i = 0;
	inNMRConFile.open(inNMRconstraints.c_str());
    
    while (getline(inNMRConFile,s)) {

		if(!s.empty()){//make sure that the line is not empty

			while(s.find('(') != string::npos){//check if there are any open parenthesis in the string

                if(s[s.find('(')+1]!=')'){//check if the constrains are not empty

                    for(int tr=0;tr<4;tr++){//for the four positions of the triplet constraint
                        tripletArray[i][s.find('(')][0][tr]=s.at(s.find('(')+tr+1);//store the triplet constraint in the 'tripletArray' array
                    }
                
                    if(s.at(s.find('(')+5)=='/'){//check if there is a second constraint separated by '/'
                        for(int tr=0;tr<4;tr++){//for the four positions of the triplet constraint
                            tripletArray[i][s.find('(')][1][tr]=s.at(s.find('(')+tr+6);//store the triplet constraint in the 'tripletArray' array
                        }
                        s.erase(s.find('('),s.find(')')-s.find('(')+1);
                    }//END: if there is a second constraint
                    else(s.erase(s.find('('),s.find(')')-s.find('(')+1));//erase the triplet constraint from the string
                }//END: while loop that is looking for triplet constraints

                else{//if the constraints are empty
                    s.erase(s.find('('),s.find(')')-s.find('(')+1);
                }
            }
		
			conArray[i][0] = s.size();
			xCoords[i][0] = 0;
			yCoords[i][0] = 0;
			
			for (j = 1; j <= s.size(); j++) {
				conArray[i][j] = atoi(s.substr(j-1,1).c_str());
				xCoords[i][j] = 0;
				yCoords[i][j] = 0;
			}
			i++;		
		}
	}
    inNMRConFile.close();
	inNMRConFile.clear();

#if defined(VERBOSE_MODE)
	// Display the constraints array
	cout << "NMR Constraints:\n";
	for (int im = 0; im < numOfCon; im++) {
		for (int jm = 1; jm <= conArray[im][0]; jm++) {
			cout << conArray[im][jm];
			for (int km = 0; km <2; km++){
				for(int m=0;m<4;m++) cout << tripletArray[im][jm][km][m];
			}
		}
		cout << "\n";
	}
#endif

	// Create 1D bool array for storing the nucleotides that are being considered for a match
	bool* alreadyUsed = new bool[ct->GetSequenceLength()+1];
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		alreadyUsed[i] = false;
	}

	// Create storage container for matches
	vector<conMatch> matchVector;

	// Create 1D bool array for storing information about which constraints are symmetric
	bool* isSymmetric = new bool[numOfCon];
	for (i = 0; i < numOfCon; i++) {
		isSymmetric[i] = false;
	}


	// Test to see which constraints are symmetric
	bool tempEquality;
	for (i = 0; i < numOfCon; i++) {
		j = 0;
		tempEquality = true;
		while (j < conArray[i][0] && tempEquality) {
			if (conArray[i][conArray[i][0]-j] != conArray[i][j+1]) {
				tempEquality = false;
				continue;
			}
			j++;
			if (j == conArray[i][0]) {
				isSymmetric[i] = true;

#if defined(VERBOSE_MODE)
				cout << "Symmetric constraint: " << i << "\n";
#endif

			}
		}
	}
	short currConNum = 0;
	short currConPos = 1;

	cout << "\t\t\t\tDONE\nRunning dynamic programming algorithm..." << flush;
	
	// Fill arrays 'v', 'vmb', and 'vext' with thermodynamic data
	dynamic(ct,data,maxtracebacks,cutoff,windowsize,&v,&vmb,&vext);
    
	cout << "\tDONE\nGenerating dotplot..." << flush;

    // Calculate the dot plot information that tracks MB loops and exterior loops
	for (i=1;i<=ct->GetSequenceLength();i++) {
		for (j=(i+1);j<=ct->GetSequenceLength();j++) {
			if ((v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))<0) {
				pairs= ((double)(v.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0);
				mbpairs=min( ((double)(vmb.f(i,j)+v.f(j,i+ct->GetSequenceLength()))/10.0), ((double)(v.f(i,j)+vmb.f(j,i+ct->GetSequenceLength()))/10.0));
				extpairs= ((double)(v.f(i,j)+vext.f(j,i+ct->GetSequenceLength()))/10.0);
				//cerr << i << "\t" << j << "\t" << pairs<<"\t"<< mbpairs<< "\t"<<extpairs<<"\n"; 
				//             if (instream >> x         >> y            >> z         >> vmb           >> vext) {

				convertedDotplot[i][j] = bpLookup[ct->numseq[i]][ct->numseq[j]];
				dgArray[i][j] = short(pairs*10);
				mbDotplot[i][j] = min(short(mbpairs*10),short(extpairs*10));
				// Copy these values to the other half of the matrices
				convertedDotplot[j][i] = bpLookup[ct->numseq[i]][ct->numseq[j]];
				dgArray[j][i] = short(pairs*10);
				mbDotplot[j][i] = min(short(mbpairs*10),short(extpairs*10));
				if(pairs*10 < dgMin) dgMin = short(pairs*10); // Finds and stores minimum free energy value from dotplot
				//cerr << i << "\t" << j << "\t" << dgArray[i][j] << "\t" << mbDotplot[i][j] << "\t" << convertedDotplot[i][j] << endl;
			}
		}
	}

#if defined (DOTPLOTS_OUT)
    cerr <<endl << "mbDotplot:\n";;
    for(int pi=1; pi<=ct->GetSequenceLength();pi++){
        cerr << pi << ":\t";
        for(int pj=1;pj<=ct->GetSequenceLength();pj++){
            cerr << mbDotplot[pi][pj] << " ";
        }
        cerr << '\n';
    }

    cerr <<endl << "dgArrayDotplot:\n";;
    for(int pi=1; pi<=ct->GetSequenceLength();pi++){
        cerr << pi << ":\t";
        for(int pj=1;pj<=ct->GetSequenceLength();pj++){
            cerr << dgArray[pi][pj] << " ";
        }
        cerr << '\n';
    }
#endif
	// Store cutoff DG value
	int dgCutoff = dgMin*(100-percent)/100;
    //cerr << "dgMin " << dgMin << "\ndgCutoff " << dgCutoff << "\nDelta " << dgCutoff - dgMin << "\n";
	// Loop back through convertedDotplot and remove those values that have a DG value greater than dgCutoff
	for (i = 0; i <= ct->GetSequenceLength(); i++) {
		for (j = 0; j <= ct->GetSequenceLength(); j++) {
			if (dgArray[i][j] > dgCutoff) convertedDotplot[i][j] = 0;
		}
	}

	cout << "\t\t\t\tDONE\nLooking for Matches: " << flush;
	
#if defined(DOTPLOTS_OUT)
	// Display the converted dotplot
	cout << "\nConverted Dotplot:\n";
	for(i = 1; i <= ct->GetSequenceLength(); i++) {
        cout << i << ":" << tobase(ct->numseq[i]) <<  '\t';
		for(j = 1; j <= ct->GetSequenceLength(); j++){
            if(convertedDotplot[i][j]==0) cout << "-";
            else cout << convertedDotplot[i][j];
        }
		cout << "\n";
	}
	cout << "\n";
#endif
    //cerr << "\nifwm=" << ifwarningMessage << endl;
    //if(*ifwarningMessage==true) return;
	// Ready to begin recursive constraint matching
	firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,&matchVector,
	              &currConNum,&currConPos,convertedDotplot,mbDotplot,dgCutoff,ct->GetSequenceLength(),ct,tripletArray, warningLimit, &ifwarningMessage);
    //if(ifwarningMessage==1){
        //cerr << "!!!!!!!!!RETURNED!!!!!!!!1\n";
    //    return;
    //}

    // cerr << endl<< endl;

    

	// When finished with constraint matching
	cout << "\rLooking for matches: " << matchVector.size() << " matches found\tDONE\n" << flush;

#if defined(VERBOSE_MODE)
    cerr  << "number of constraints: " << totalConLength << endl;
    for(int foo=0;foo<matchVector.size();foo++){
        cerr << "--- MATCH " << foo << " ---"<< endl;
        for(int faa=0;faa<matchVector[foo].size()-1;faa=faa+2){
            cerr << matchVector[foo][faa] << "-" << matchVector[foo][faa+1]<< endl;
        }
    }

	// Display the matches
	//conMatch tempMatchVector;
	//for (i = 0; i < matchVector.size(); i++) {
	//	tempMatchVector = matchVector[i];
	//	for (j = 0; j < totalConLength; j++) {
    //        //cout << tempMatchVector.xCoords[j*2] << "," << tempMatchVector.xCoords[j*2+1] << " ";
	//	}
	//	cout << "\n";
	//}
#endif

	// Determine if refolding should occur
	if (matchVector.size() == 0) {
		cout << "No possible constraint matches - NAPSS will terminate.\n";
        // Insert steps to clean up allocated memory
        matchVector.clear();
        delete[] alreadyUsed;
        delete[] isSymmetric;
	
        for (i=0;i<=ct->GetSequenceLength();i++){
            delete[]convertedDotplot[i];
            delete[]dgArray[i];
            delete[]mbDotplot[i];
        }
        delete[]convertedDotplot;
        delete[]dgArray;
        delete[]mbDotplot;

        for(i=0;i<=numOfCon;i++){
            for(j=0;j<=maxConLength;j++){
                for(int k=0;k<2;k++){
                    delete[]tripletArray[i][j][k];
                }
            }
        }
        for(i=0;i<=numOfCon;i++){
            for(j=0;j<=maxConLength;j++){
                delete[]tripletArray[i][j];
            }
        }
        for(i=0;i<=numOfCon;i++){
            delete[]tripletArray[i];
        }
        delete[]tripletArray;

        for (i=0;i<=numOfCon;i++){
            delete[]conArray[i];
            delete[]xCoords[i];
            delete[]yCoords[i];
        }
        delete[]conArray;
        delete[]xCoords;
        delete[]yCoords;
        return NAPSS_ERR_NO_MATCHES; //RMW: TODO: Convert to error number/message
	}

    //Create another vector that holds matches
    vector<conMatch> matchVectorExtended;
    //Depending on the #define flag 'matchVectorExtended' will either mimic the 'matchVector' 
    //or will hold matches with extensions made by HelicalExtension
#ifdef ExtendBefore

	//Initialize RNA class with RNA sequence.
	RNA *rnaMATCH=new RNA(sequence,SEQUENCE_STRING,rnaCT);
	//Read the SHAPE data from the disk
	if(!inSHAPEfile.empty()) rnaMATCH->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);
	//Read the constraint data from the disk
	if(!constraintFile.empty()) rnaMATCH->ReadConstraints(constraintFile.c_str());

	if (rnaMATCH->GetErrorCode() != 0) return NAPSS_BAD_RNA_COPY;

    //structure class will point to the RNA class
    structure *ctMATCH = rnaMATCH->GetStructure();

	// Also create a 1D array for temporarily storing possible extensions and a control switch for that loop
	int* helixExtend2 = new int[ct->GetSequenceLength()+1];
    //cerr  << "number of constraints: " << totalConLength << endl;

    //for every match extend the matched pairs using helicalexetnsion function
    for(int matchnum=0;matchnum<matchVector.size();matchnum++){
        ctMATCH->AddStructure(); // when created with a sequence, no structures are created by default.
        //Store matches in the rnaMATCH structure
        for(int matchpos=0;matchpos<matchVector[matchnum].size()-1;matchpos=matchpos+2){
			rnaMATCH->SpecifyPair(matchVector[matchnum][matchpos],matchVector[matchnum][matchpos+1]);
            //cerr << matchVector[matchnum][matchpos] << "-" << matchVector[matchnum][matchpos+1]<< endl;
        }
        //rnaMATCH->WriteCt("before.ct");
        //assign structure class member 'ctMATCH' to point to rnaMATCH
        //ctMATCH=rnaMATCH->GetStructure();
        //Extend the matched helices using helical extension function
        HelicalExtension(ctMATCH, convertedDotplot, helixExtend2, dgArray);
        //rnaMATCH->WriteCt("after1.ct");

        //tempMatch vector will temporarily store the pairs after they were extended
        conMatch tempMatch;

        //go through the structure rnaMATCH and populate tempMatch vector with all the pairs
        for(int index=1;index<rnaMATCH->GetSequenceLength();index++){
            if(ctMATCH->GetPair(index)>index){
                tempMatch.push_back(index);
                tempMatch.push_back(rnaMATCH->GetPair(index));
                //cerr << "pair " << index << "-" << rnaMATCH->GetPair(index) << endl;
            }
        }
        //store the tempMatch vector in the 'matchVectorExtended' for every match set
        matchVectorExtended.push_back(tempMatch);
        //rnaMATCH->WriteCt("after2.ct");
        rnaMATCH->RemovePairs();
    }
    
    delete rnaMATCH;
    delete[] helixExtend2;

#else
    //if the extension was not defined, set 'matchVectorExtended' equal to the original 'matchVector'
    matchVectorExtended=matchVector;
#endif
    
#ifdef OUTPUT_MATCHES
    cout << endl;
    for(int pi=0;pi<matchVectorExtended.size();++pi){
        cout << "MATCH #" << pi << " ";
        for(int pp=0;pp<matchVectorExtended[pi].size();pp=pp+2){
            cerr <<matchVectorExtended[pi][pp] << "-" << matchVectorExtended[pi][pp+1] << " ";
        }
        cout << endl;
    }
#endif


	// Initialize ct's that will be used for storing all refolded structures and for calculating pseudoknot energies
	structure *ct2;
	//structure *ctbroken;

	//Initialize RNA class with RNA sequence.
	// There should be no errors in the next few lines that create a copy of rnaCT, 
	// since presumably the caller has already performed these steps on the original 
	// RNA and would not have called napss if there was an error.

	RNA *rnaCT2=new RNA(sequence,SEQUENCE_STRING,rnaCT);
	//RNA *rnaCTBROKEN=new RNA(inseq.c_str(),FILE_SEQ);

	//Read the SHAPE data from the disk
	if(!inSHAPEfile.empty()) rnaCT2->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);
	//if(!inSHAPEfile.empty()) rnaCTBROKEN->ReadSHAPE(inSHAPEfile.c_str(),slope,intercept);
	//Read the constraint data from the disk
	if(!constraintFile.empty()) rnaCT2->ReadConstraints(constraintFile.c_str());
	//if(!constraintFile.empty()) rnaCTBROKEN->ReadConstraints(constraintFile.c_str());

	if (rnaCT2->GetErrorCode() != 0) return NAPSS_BAD_RNA_COPY;

	//Set pointer 'ct2' to point to 'structure' in 'rnaCT2'
	ct2=rnaCT2->GetStructure();
	//Set pointer 'ctbroken' to point to 'structure' in 'rnaCTBROKEN'
	//ctbroken=rnaCTBROKEN->GetStructure();

	int start = 0;
	int count = 0;
	// Also create a 1D array for temporarily storing possible extensions and a control switch for that loop
	int* helixExtend = new int[ct->GetSequenceLength()+1];
	bool extensionAdded = true;

	// Create 2D boolean array (triangular matrix with sides of length ct->GetSequenceLength()+1) to store dotplot data
	ct->allocatetem();

    //cerr << "Num_Of_Matches " << matchVectorExtended.size() << '\n';

	// Loop over all matches
	for (iter = 0; iter < matchVectorExtended.size(); iter++) {

		cout << "\rFolding structure " << iter+1 << " of " << matchVectorExtended.size() << flush;
		//if (iter>10) continue;

		if(!pseudoknotFree){// If pseudoknot-free mode is NOT specified (default behavior)

			// Initialize all positions to false:
			for (i = 0; i < ct->GetSequenceLength(); i++) {
				for (j = i+1; j <= ct->GetSequenceLength(); j++) {
					ct->tem[j][i] = false;
				}
			}
			
#if defined(VERBOSE_MODE)
			cout << "Non-conflicting match " << iter+1 << ":\n";
#endif

			// Copy original dotplot into tem** - change all non-zero convertedDotplot values to true
			for (i = 0; i < ct->GetSequenceLength(); i++) {
				for (j = i+1; j <= ct->GetSequenceLength(); j++) {
					if (convertedDotplot[j][i] != 0){
                        ct->tem[j][i] = true;
                        //cerr << '\n' << "ct->tem i=" << i << " j=" << j << '\n';
                    }
				}
			}

			// Create new array to temporarily hold matched basepairs (initialize all to 0)
			int* tempbasepr = new int[ct->GetSequenceLength()+1];
			for (i=0; i<=ct->GetSequenceLength(); i++) tempbasepr[i]=0;
		
			// Trim the dotplot according to the base pairs of the current match
			int trimming_i, trimming_j;

			// Loop over all base pairs in the match combination
			for (iter2 = 0; iter2 < matchVectorExtended[iter].size(); iter2=iter2+2) {
				trimming_j = matchVectorExtended[iter][(iter2)];
				trimming_i = matchVectorExtended[iter][(iter2)+1];
                
				// Store coordinates in temporary basepair array
				tempbasepr[trimming_i] = trimming_j;
				tempbasepr[trimming_j] = trimming_i;
                
				// Change all tem values to false where (i or j) = (trimming_i or trimming_j)
				for (j = 0; j <= trimming_i-1; j++) {
					ct->tem[trimming_i][j] = false;
                }
				for (i = trimming_i + 1; i <= ct->GetSequenceLength(); i++) {
					ct->tem[i][trimming_i] = false;
                }
				for (j = 0; j <= trimming_j-1; j++) {
					ct->tem[trimming_j][j] = false;
				}
				for (i = trimming_j + 1; i <=ct->GetSequenceLength(); i++) {
					ct->tem[i][trimming_j] = false;
                } 
#ifdef AllowMatchedPairs
                if(trimming_i<trimming_j) ct->tem[trimming_j][trimming_i] = true;
                else ct->tem[trimming_i][trimming_j] = true;
#endif
            }

			// Search for regions of dotplot that can be excluded
            //STAS: removed the pseudodptrim because probably don't need it. Probably was needed when the previous energy model (DP) was used.
			//pseudodptrim(ct, tempbasepr, &count);
			delete[] tempbasepr;


            //STAS: calculate PROBKNOT

#ifdef PROBKNOT_MODE
            int iterations=10;
            int MinHelixLength=2;
            rnaCT->PartitionFunction();
            rnaCT->ProbKnot(iterations,MinHelixLength);
            
#else

            /////END CALCULATE PROBKNOT

			// Do the structure prediction considering only the allowed basepairs from the trimmed dotplot
			//RMW-DBG: cout << "Structures before: " << ct->GetNumberofStructures() << endl;
			for(int sn=ct->GetNumberofStructures(); sn>0; sn--)
				ct->RemoveStructure(sn);
			//RMW-DBG: cout << "Structures before: " << ct->GetNumberofStructures() << endl;
			dynamic(ct,data,maxtracebacks,cutoff,windowsize);
			//RMW-DBG: cout << "Structures after dynamic: " << ct->GetNumberofStructures() << endl;

			// cerr << "Verifying After Dynamic" << endl;
			// verifyCTAll(ct);			
#endif

            /*OUTPUT FILE FOR DEBUGGING
              ostringstream ss;
              ss << iter;
              string s = ss.str();
              string name="file_";
              name.append(s);
              ctout(ct,name.c_str());
            */
#if defined(VERBOSE_MODE)
			cout << "Refolding yields " << ct->GetNumberofStructures() << " structures.\n\n";
#endif

            //rnaCT->WriteCt("after3.ct");

			// Reinsert the basepairs that match the NMR constraints
			for (i = 1; i <= ct->GetNumberofStructures(); i++) {
                //cerr << "\n#########ct->basepr[10]=" << ct->GetPair(10,i) << " ###########\n";
				for (iter2 = 0; iter2 < matchVectorExtended[iter].size(); iter2+=2) {
					trimming_j = matchVectorExtended[iter][iter2];
					trimming_i = matchVectorExtended[iter][iter2+1];
					ct->RemovePair(trimming_i,i);
					ct->RemovePair(trimming_j,i);
					ct->SetPair(trimming_i,trimming_j,i);
                    //cerr << "\ntrimmingj=" << trimming_j << " trimming_i=" << trimming_i << endl;
				}
			}
			// cerr << "Verifying After trimming" << endl;
			// verifyCTAll(ct);			
            //cerr << "\n!!!!!!!!! num of structures = " << ct->GetNumberofStructures() << " !!!!!!!!!!!!!!\n";
			// Check for helical extension possibilities
            //for(int foo=1;foo<=ct->GetNumberofStructures();++foo){
            //if(ct->GetPair(10,foo)==23) cerr << "\nBEFORE!!!!!!!!! 10=23 at" << foo << " !!!!!!!!!!!\n";
            //}
            //rnaCT->WriteCt("after4.ct");

            //STAS: add the "HelicalExtension" back
#ifdef ExtensionAfter
			HelicalExtension(ct, convertedDotplot, helixExtend, dgArray);
			// cerr << "Verifying After HelicalExtension" << endl;
			// verifyCTAll(ct);				
#endif
            //rnaCT->WriteCt("after5.ct");

            //for(int foo=1;foo<=ct->GetNumberofStructures();++foo){
            //if(ct->GetPair(10,foo)==23) cerr << "\nAFTER!!!!!!!!! 10=23 at" << foo << " !!!!!!!!!!!\n";
            //}
		}
        
		else{// If pseudoknot-free mode is specified
#if defined(VERBOSE_MODE)//output forsed paires
            cerr << "\nForced Pairs in Pseudoknot-Free mode:\n";
#endif
            // Loop over all base pairs in the match combination
            for (iter2 = 0; iter2 < matchVectorExtended[iter].size(); iter2=iter2+2) {
				if(matchVectorExtended[iter][(iter2)]<matchVectorExtended[iter][(iter2)+1]){//if the 2*iter2 position is a 5' position
                    ct->AddPair(matchVector[iter][(2*iter2)],matchVector[iter][(2*iter2)+1]);
					//ct->pair[(iter2)/2+1][0] = matchVectorExtended[iter][(iter2)];//force paired 5'
                    //ct->pair[(iter2)/2+1][1] = matchVectorExtended[iter][(iter2)+1];//force paired 3'
                    //cerr << endl << matchVectorExtended[iter][(iter2)] << " - " << matchVectorExtended[iter][(iter2)+1] << " pair[" << (iter2)/2+1 << "][0]=" << ct->pair[(iter2)/2+1][0] << endl;
                }
                else{//else if the 2*iter2 position is a 3' position
					ct->AddPair(matchVector[iter][(2*iter2)+1],matchVector[iter][(2*iter2)]);
					//ct->pair[(iter2)/2+1][1] = matchVectorExtended[iter][(iter2)];//force paired 3'
                    //ct->pair[(iter2)/2+1][0] = matchVectorExtended[iter][(iter2)+1];//force paired 5'
                    //cerr << endl << matchVectorExtended[iter][(iter2)+1] << " - " << matchVectorExtended[iter][(iter2)] << " pair[" << (iter2)/2+1 << "][0]=" << ct->pair[(iter2)/2+1][0] <<  endl;
				}
#if defined(VERBOSE_MODE)//output forsed paires
				cerr << "pair\t" << matchVectorExtended[iter][(2*iter2)] << "-" << matchVectorExtended[iter][(2*iter2)+1] << endl;
#endif
                
			}
#if defined(VERBOSE_MODE)//output forsed paires
			cerr << "------------------------------------------\n";
#endif
			//cout << "Structures before: " << ct->GetNumberofStructures() << endl;
			dynamic(ct,data,maxtracebacks,cutoff,windowsize);
			//cout << "Structures after dynamic: " << ct->GetNumberofStructures() << endl;
			// cout << "Verifying After Dynamic" << endl;
			// verifyCTAll(ct);
		}
        
		vector<int> removed, repeated, kept;
        //rnaCT->WriteCt("after10.ct");
        //ctout(ct2,"file_0");
		// Copy current ct to new ct file that is used for appending the results from all matches
		//cerr << "998 Structures: " << ct->GetNumberofStructures() << endl;
		for (i = 1; i <= ct->GetNumberofStructures(); i++) {
			// cout << "Verifying " << i << endl;
			// // verifyBP(rnaCT, i);
			//cerr << "Structure: " << i << endl;
            //IMPORTANT: We should uncomment the fillmismatch and removeisolatedpairs
            //Fill in single or tandem mismatches in rna structure 'rna'
            FillMismatch(rnaCT, i);

			// cout << "Verifying FillMismatch" << i << endl;
			// // verifyBP(rnaCT, i);            

            //Remove isolated pairs in rna structure 'rna'
            RemoveIsolatedPairs(rnaCT, i);

			// cout << "Verifying RemoveIsolatedPairs" << i << endl;
			// // verifyBP(rnaCT, i);

            bool removedMatch=false;

#ifdef RemoveBrokenMatch //if "RemoveBrokenMatch" is defined, run the following code

            //For every match in the match list
            for(int iter2=0;iter2<matchVectorExtended[iter].size();iter2=iter2+2){
                
                //check if the match is present in the structure
                if(rnaCT->GetPair(matchVectorExtended[iter][iter2], i)!=matchVectorExtended[iter][iter2+1]){
                    //if it is not present, set 'removedMatch' to true
					removed.push_back(i);
					removedMatch=true;
					break;
                    //cerr << "\nMVpair " << matchVectorExtended[iter][iter2] << "-" << matchVectorExtended[iter][iter2+1] << "  GetPair(" << matchVectorExtended[iter][iter2] << ")=" << rnaCT->GetPair(matchVectorExtended[iter][iter2]) << endl;
                }
            }
#endif
            // Check if the structure is a repeat
            //int copiedStr=0;//int that keeps track of the number of copied structures
            bool repeat=false;
            //cerr << "\n###############################\n" << "i=" << i << endl;
            if(ct2->GetNumberofStructures()!=0){
                for(int l=1;l<=ct2->GetNumberofStructures();++l){
                    //cerr << "\n###############################\n" << "l=" << l << endl;
                    int rep=1;
                    while(rep<=ct2->GetSequenceLength() && ct->GetPair(rep,i)==ct2->GetPair(rep,l)){
                        //cerr << "\nrep=" << rep << " " << ct->GetPair(rep,i) << " = " << ct2->GetPair(rep,l) << endl;
                        rep++;
                    }
                    //cerr << "\nct=" << i << " ct2=" << l << " rep=" << rep << " = " << ct2->GetSequenceLength()+1 ;
                    if(rep==ct2->GetSequenceLength()+1){
                        //cerr << " THERE WAS A REPEAT\n";
                        repeated.push_back(i * 100000 + l);
						repeat=true;
                        break;
                    }
                }
            }
            // If the structure is not a repeat, store it in ct2
#ifndef DEBUG_MATCHES
            if(!repeat&&!removedMatch){
#endif 				
				kept.push_back(i);
				//copiedStr++;
                ct2->AddStructure();
                //ctbroken->AddStructure();
				int strnum = ct2->GetNumberofStructures();
                for (j = 1; j <= ct->GetSequenceLength(); j++) {
                    ct2->SetPair(j,ct->GetPair(j,i), strnum);
                }
                //cerr << " Label=" << ct->ctlabel[i] << endl;
				ct2->SetCtLabel(ct->GetCtLabel(i),strnum);
				// cout << "copied to ct2" << endl;
				// verifyCT(ct2, strnum);
#ifndef DEBUG_MATCHES
            }
#endif
            //cerr << "\n#########i=" << i << "  ct2->basepr[10]=" << ct->GetPair(10,i) << " " << ct->GetPair(23,i) << " ct2->basepr[10]=" << ct2->GetPair(10,start+i) << " " << ct2->GetPair(23,start+i) << " ###########\n";
		}
		//start += ct->GetNumberofStructures();
		//RMW-DBG: showVec("Removed",removed);
		//RMW-DBG: showVec("Repeated",repeated);
		//RMW-DBG: showVec("Kept",kept);
	}
	
    //ctout(ct2,"file_0");
	cout << "\t\t\tDONE\n" << flush;

#if defined(VERBOSE_MODE)
	cout << "Refolding yields " << ct2->GetNumberofStructures() << " total structures.\n\n" << flush;
#endif
	//QUESTIONS QUESTIONS:	
	delete[] helixExtend;

    

    // cerr << "Num_Of_Structures " << ct2->GetNumberofStructures() << '\n';

	// Re-calculate free energy with modified efn2
	for (i=1;i<=rnaCT2->GetStructureNumber();i++){
        cout << "\rRecalculating energies " << i << " of " << ct2->GetNumberofStructures() << flush;
		CalculateFreeEnergywithPseudoknot(rnaCT2, i, P1, P2, pseud_param, data, useSimpleMBRules, false, false);
	}
	cout << "\t\tDONE\nSorting structures..." << flush;
    
#ifdef DEBUG_MATCHES
    int lowestenergy=0;
    for(int ii=1;ii<=rnaCT2->GetStructureNumber();++ii){
        if((int)(rnaCT2->GetFreeEnergy(ii)*10)<lowestenergy){
            lowestenergy=rnaCT2->GetFreeEnergy(ii)*10;
        }
    }
    
    for(int ii=1;ii<=rnaCT2->GetStructureNumber();++ii){
        if((int)(rnaCT2->GetFreeEnergy(ii)*10)==lowestenergy){
            cerr << "\nmatch_number=" << ii << " energy=" << rnaCT2->GetFreeEnergy(ii) << " match=";
            for(int pp=0;pp<matchVectorExtended[ii-1].size();pp=pp+2){
                cerr << matchVectorExtended[ii-1][pp] << "-" << matchVectorExtended[ii-1][pp+1] << " ";
            }
            cerr << endl;
        }
    }
#endif

	// Re-sort the structures according to efn2mod-calculated energy
#ifndef DEBUG_MATCHES
    ct2->sort();
#endif
	cout << "\t\t\t\tDONE\n" << flush;

	// Calculate value of cutoff for output structures
	if (cutoff !=0) filterByEnergy(ct2, (100-cutoff)*ct2->GetEnergy(1)/100); // remove structures with energy above cutoff.

	// Re-output ct file
	cout << ct2->GetNumberofStructures() << " structures are within the specified cutoff percentage.\n";
	outputStructureResults = ct2;

	// Insert steps to clean up allocated memory
	matchVector.clear();
    matchVectorExtended.clear();
	delete[] alreadyUsed;
	delete[] isSymmetric;
	
	for (i=0;i<=ct->GetSequenceLength();i++){
		delete[]convertedDotplot[i];
		delete[]dgArray[i];
		delete[]mbDotplot[i];
	}
	delete[]convertedDotplot;
	delete[]dgArray;
	delete[]mbDotplot;

	for(i=0;i<=numOfCon;i++){
		for(j=0;j<=maxConLength;j++){
			for(int k=0;k<2;k++){
				delete[]tripletArray[i][j][k];
			}
		}
	}
 	for(i=0;i<=numOfCon;i++){
		for(j=0;j<=maxConLength;j++){
			delete[]tripletArray[i][j];
		}
	}

	for(i=0;i<=numOfCon;i++){
		delete[]tripletArray[i];
	}
	delete[]tripletArray;

	for (i=0;i<=numOfCon;i++){
		delete[]conArray[i];
		delete[]xCoords[i];
		delete[]yCoords[i];
	}

	delete[]conArray;
	delete[]xCoords;
	delete[]yCoords;
	//delete rnaCT2;
	//delete rnaCTBROKEN;

	return 0; // success
}

// Remove structures with energy above a specified cutoff.
int filterByEnergy(structure *ct, double energyCutoff) {
	for(int i = ct->GetNumberofStructures(); i > 0; i--)
		if (ct->GetEnergy(i) > energyCutoff)
			ct->RemoveStructure(i);
	return ct->GetNumberofStructures();
}


void errmsg2(int err,int erri) {

	if (err==30) {
		cout << "End Reached at traceback #"<<erri<<"\n";
		exit(1);
	}
	if (err==100) {
		cout << "error # "<<erri;
		exit(1);
	}
	switch (err) {
	case 1:
		cout << "Could not allocate enough memory";
		break;
	case 2:
		cout << "Too many possible base pairs";
		break;
	case 3:
		cout << "Too many helixes in multibranch loop";
	case 4:
		cout << "Too many structures in CT file";
	default:
		cout << "Unknown error";
	}
	cin >> err;
	exit(1);
	return;
}


void firstDotMatch(short** conArray, bool* alreadyUsed, bool* isSymmetric, short totalConLength,
                   short numOfCon, short** xCoords, short** yCoords, vector<conMatch>* matchVector,
                   short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
                   int dgCutoff, short sequenceLength,structure *ct,char**** tripletArray, int warningLimit, bool* ifwarningMessage) {

	short i, j;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);

#if defined(DEBUG_MODE)
	cout << "currConNum="<< currConNum << " currConPos=" << currConPos << "\n";
	cin >> i;
#endif
    //    cerr << "1: iwm=" << *ifwarningMessage << endl;
	if(matchVector->size()>warningLimit && *ifwarningMessage==false){
		cout << "\r                                                   ";
		cout << "\n!! There are too many matches. This will make the program run very long time.\n"<< flush;
		cout << "!! To speed up the calculation try lowering the -d and -m parameters (see help),\n"<< flush;
		cout << "!! or try removing shortest constraints from the constraint file.\n\n" << flush;
		*ifwarningMessage=true;
        cerr << "\nNum_Of_Matches " << matchVector->size();
        cerr << "\nTOO MANY MATCHES1 FOUND\n";
        cerr << "\nTo proceed enter 'P'; to exit enter 'E'\n";
        char usermessage;
        cin >> usermessage;
       
        if(usermessage=='E'){
            cerr << "You have entered 'E', the program is going to exit now\n";
            exit(1);
        }
	}
	//else if(matchVector->size()%1000 == 0){
	//	cout << "\rLooking for matches: " << matchVector->size() << " matches found" << flush;
	//}
    //cerr << "\nNum_Of_Matches " << matchVector->size();
    //cerr << "\nCHECK1\nifm=" << *ifwarningMessage << endl;
    //if(*ifwarningMessage==1){
        //cerr << "!!!!!!!!!RETURNED!!!!!!!!1\n";
    //    exit(1);
    //}
    //cerr << "\nCHECK1.1\n";
	// Loop through x,y of dotplot to find potential match for first position of the current constraint
	// Skip all values that have already been used by matches for previous constraints
	for (j = 1; j < sequenceLength+1; j++) {
		if (alreadyUsed[j] == true) continue;

		for (i = 1; i < sequenceLength+1; i++) {
			if (alreadyUsed[i] == true) continue;

            //Stas QUESTION: do we want the symmetric constraint even when we have the triplet constr.
            //cerr << "\nCHECK1.2\n";

			// If constraint is symmetric, only consider matches from upper-right half of dotplot
            //STAS: I commented this code because of triplet constraint symmetry.
			if (isSymmetric[currConNum]==false || (isSymmetric[currConNum]==true && i<j)) {
               

            // Check the first value only here - if this value matches,
            // then call the recursive method, passing it the proper parameters
            if (conArray[currConNum][1] == convertedDotplot[i][j]) {
				//cerr << "\nfirstDotMatch: " << i << "-" << j << endl;
                xCoords[currConNum][1] = i;
                yCoords[currConNum][1] = j;
                alreadyUsed[i] = true;
                alreadyUsed[j] = true;
                *(pCurrConPos) = 2;
                
#if defined(DEBUG2_MODE) 
                cerr << "\n####################\n#RECURSIVE_MATCH_IS_CALLED\n########################\nfirst match is found at " << i << "-" << j << endl;;
#endif
                recursiveMatch(conArray, alreadyUsed, totalConLength, numOfCon, xCoords, yCoords, 
                               matchVector, pCurrConNum, pCurrConPos, convertedDotplot, mbDotplot, dgCutoff, sequenceLength, 
                               isSymmetric,ct,tripletArray,warningLimit,ifwarningMessage);
                currConNum = *(pCurrConNum);
                currConPos = *(pCurrConPos);
            }
            } 
		}
	}
    //cerr << "\nCHECK2\n";
	// Once we've looped through all coordinates for this base pair, step back values for constraint 
	// number and position, flip switches for last base pair, and return
	currConNum--;
	if (currConNum > -1) {
		currConPos = conArray[currConNum][0];
		alreadyUsed[xCoords[currConNum][currConPos]] = false;
		alreadyUsed[yCoords[currConNum][currConPos]] = false;
	}
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;
}

void recursiveMatch(short** conArray, bool* alreadyUsed, short totalConLength, short numOfCon,
                    short** xCoords, short** yCoords, vector<conMatch>* matchVector, 
                    short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
                    int dgCutoff, short sequenceLength, bool* isSymmetric,structure *ct,char**** tripletArray, int warningLimit, bool* ifwarningMessage) {

	short i,j,k,xIter,yIter,start,stop;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);
	//conMatch* tempMatch;
	string a;
	bool searchX;  // Switch to indicate in which direction we're searching
	bool searchExtCoax; // Switch to indicate whether we're currently searching for a coaxial stack in 
	// the external direction
	
 subroutine: // Had to use this so that part of the code could return here from a nested "for" loop

	// This function will loop until we reach the end of the current constraint, 
	// or until no more matches are found using the current nucleotides for the
	// first base pair of this constraint
	while (currConPos > 1) {
		
#if defined DEBUG2_MODE 
		cout << currConNum << " " << currConPos << " " << xCoords[currConNum][currConPos-1] << 
			" " << yCoords[currConNum][currConPos-1] << " " << xCoords[currConNum][currConPos] <<
			" " << yCoords[currConNum][currConPos] << " First\n";
        //	cin >> i;

        cerr << "\ncurrConPos=" << currConPos << " conArray[" << currConPos << "][0]+1=" << (conArray[currConNum][0]+1) << endl;
#endif
		// Check if we've matched all the base pairs for the current constraint
		if (currConPos == conArray[currConNum][0]+1) {
#if defined DEBUG2_MODE
            cerr << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            cerr << "currConNum=" << currConNum << " == numOfCon-1=" << (numOfCon-1) << endl;
#endif
 			// Check if this is the final constraint
			if (currConNum == numOfCon-1) {
				// If so, we've completed a match - compile 1D short array of dotplot 
				// coordinates and append this to matchVector

				//				conMatch tempMatch;
				/*
                  tempMatch = new conMatch;
                  tempMatch->coords = new short[totalConLength*2];
				*/
				//				tempMatch.coords=vector<short>(totalConLength*2);
				conMatch tempMatch(totalConLength*2);

				i = 0; // current base pair to be written
				j = 0; // takes place of currConNum
				k = 1; // takes place of currConPos
				
				while (i < totalConLength) {
					if (k > conArray[j][0]) {j++; k=1; continue;} // reached end of current constraint
					tempMatch[i*2] = xCoords[j][k];
					tempMatch[i*2+1] = yCoords[j][k];
					i++;
					k++;
				}

#if defined(OUTPUT_MATCHES)
				cout << "\nMatch found!: ";
				for (i = 0; i < totalConLength*2; i=i+2) {
					cout << tempMatch[i] << "-" << tempMatch[i+1] << " ";
				}
				cout << "\n";
#endif

#if defined(ALREADY_USED_CHECK)
				// Verify that the number of "true" alreadyUsed values = twice the number of basepairs in all constraints
				j = 0;
				for (i = 1; i<=sequenceLength; i++){
					if (alreadyUsed[i]) j++;
				}
				if (j!=totalConLength*2){cerr<<"ERROR: Discrepancy in alreadyUsed array.\n";}
#endif

				matchVector->push_back(tempMatch);

#if defined(VERBOSE_MODE)
				if (matchVector->size() % 50000 == 0) {
					// Display a current count of matches
					cout<<matchVector->size()<<" ";
					for (i=0; i<numOfCon; i++){
						cout <<xCoords[i][1]<<","<<yCoords[i][1]<<" ";
					}
					cout << "\n";
				}
#endif
                //return;
				// Now step back two base pairs, clear final position in the Coords arrays,
				// flip the alreadyUsed values of last 2 positions, and continue searching
				currConPos=currConPos-2;
				alreadyUsed[xCoords[currConNum][currConPos+1]] = false;
				alreadyUsed[yCoords[currConNum][currConPos+1]] = false;
				xCoords[currConNum][currConPos+1] = 0;
				yCoords[currConNum][currConPos+1] = 0;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If this is not the last constraint, we need to move on and match the first position 
			// of the next constraint
			currConNum++;
			currConPos = 1;
			*(pCurrConNum) = currConNum;
			*(pCurrConPos) = currConPos;
			firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,matchVector,
			              pCurrConNum,pCurrConPos,convertedDotplot,mbDotplot,dgCutoff,sequenceLength,ct,tripletArray, warningLimit, ifwarningMessage);
            //if(*ifwarningMessage==1){
                //cerr << "!!!!!!!!!RETURNED!!!!!!!!1\n";
            //    return;
            //}
			// Once this returns, need to continue the search for matches to this constraint
			currConNum = *(pCurrConNum);
			currConPos = *(pCurrConPos);
			continue;
		}
		
		// Check if we're still in-bounds on the next position across and down.
		// If we're not, step back, flip switches, and continue
		if (xCoords[currConNum][currConPos-1]-1 < 1 || yCoords[currConNum][currConPos-1]+1 > sequenceLength) {
#if defined DEBUG2_MODE
            cerr << "\nINSIDE FIRST IF STATEMENT\n";
#endif
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;

		}
#if defined DEBUG2_MODE
        cerr << "\nCHECK7\n";
        cerr << "\n\n#################################################\n\n";
        cerr << "\ncurrConPos=" << currConPos << " currConNum=" << currConNum << " conArray[" << currConNum << "][" << currConPos << "]=" << conArray[currConNum][currConPos] << " xCoords[" << currConNum << "][" << (currConPos-1) << "]=" << xCoords[currConNum][currConPos-1] << " yCoords[" << currConNum << "][" << (currConPos-1) << "]=" << yCoords[currConNum][currConPos-1] << " xCoords[" << currConNum << "][" << (currConPos) << "]=" << xCoords[currConNum][currConPos] << " yCoords[" << currConNum << "][" << (currConPos) << "]=" << yCoords[currConNum][currConPos]  << endl;
#endif
		// If this is the first or last step, or if the previous basepair is unstable as a closing pair for 
		// external and mb loops, only look for non-bulged/non-stacked match on current step
        //		if (/*currConPos==2 || currConPos==(conArray[currConNum][0]) || */
        //  mbDotplot[xCoords[currConNum][currConPos-1]][yCoords[currConNum][currConPos-1]] >= dgCutoff /*||
		//	                            xCoords[currConNum][currConPos-2]-xCoords[currConNum][currConPos-1]>1 ||
		//	                            yCoords[currConNum][currConPos-1]-yCoords[currConNum][currConPos-2]>1*/) {
#if defined DEBUG2_MODE
        //cerr << "\nCHECK8\n";
        //if(xCoords[currConNum][currConPos-1]==43)
        //cerr << "currConPos==2==" << currConPos << " currConPos==" << currConPos << "==conArray[" << currConNum << "][0]==" << conArray[currConNum][0] << " dgCutoff==" << dgCutoff << "<=" << mbDotplot[xCoords[currConNum][currConPos-1]][yCoords[currConNum][currConPos-1]] << endl;
 
        //cerr << "xCoords[" << currConNum << "][" << currConPos << "]=" << xCoords[currConNum][currConPos] << " yCoords=" << yCoords[currConNum][currConPos] << endl;
#endif
        // If this is the first or last step,
        //if ((currConPos==2 || currConPos==(conArray[currConNum][0])){
        //&&(xCoords[currConNum][currConPos] != 0 || yCoords[currConNum][currConPos] !=0)) {

        //if(yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
            //xCoords[currConNum][currConPos] <= xCoords[currConNum][currConPos-1]-1){

        //STAS QUESTION: this used to be "!=". I changed it to ">" because the program was not switching the search from X direction to Y direction
        //if (xCoords[currConNum][currConPos] != 0 || yCoords[currConNum][currConPos] !=0) {
	/*		
        if (xCoords[currConNum][currConPos] > 0 || yCoords[currConNum][currConPos] > 0) {
#if defined DEBUG2_MODE
            cerr << "\nTHE PAIR HAS BEEN CHECKED - DELETING POSITION\n";
#endif
            xCoords[currConNum][currConPos] = 0;
            yCoords[currConNum][currConPos] = 0;
            currConPos--;
            alreadyUsed[xCoords[currConNum][currConPos]] = false;
            alreadyUsed[yCoords[currConNum][currConPos]] = false;
            continue;
            //}
            // This would indicate that this basepair has already been checked
            // Reset coords, step back, flip switches, and continue
            //else{
            //   xCoords[currConNum][currConPos] = 0;
            //   yCoords[currConNum][currConPos] = 0;
            //   currConPos--;
            //   alreadyUsed[xCoords[currConNum][currConPos]] = false;
            //   alreadyUsed[yCoords[currConNum][currConPos]] = false;
                //continue;
            //}
        
        }
*/		
#if defined DEBUG2_MODE
        cerr << "\ncurrConPos=" << currConPos << " currConNum=" << currConNum << " conArray[" << currConNum << "][" << currConPos << "]=" << conArray[currConNum][currConPos] << " xCoords[" << currConNum << "][" << (currConPos-1) << "]=" << xCoords[currConNum][currConPos-1] << " yCoords[" << currConNum << "][" << (currConPos-1) << "]=" << yCoords[currConNum][currConPos-1] << " xCoords[" << currConNum << "][" << (currConPos) << "]=" << xCoords[currConNum][currConPos] << " yCoords[" << currConNum << "][" << (currConPos) << "]=" << yCoords[currConNum][currConPos]  << endl;
#endif
		/*
        if (conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1][yCoords[currConNum][currConPos-1]+1] 
            && alreadyUsed[xCoords[currConNum][currConPos-1]-1] == false
            && alreadyUsed[yCoords[currConNum][currConPos-1]+1] == false) {
            //cerr << "\nCHECK10\n";
            // Match found at lastX-1,lastY+1
            xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
            yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
            cerr << "Location:A Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
                 << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
                 << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
                 << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
#if defined DEBUG2_MODE
			cerr << "\ncurrConPos=" << currConPos << " currConNum=" << currConNum << " conArray[" << currConNum << "][" << currConPos << "]=" << conArray[currConNum][currConPos] << " xCoords[" << currConNum << "][" << (currConPos-1) << "]=" << xCoords[currConNum][currConPos-1] << " yCoords[" << currConNum << "][" << (currConPos-1) << "]=" << yCoords[currConNum][currConPos-1] << " xCoords[" << currConNum << "][" << (currConPos) << "]=" << xCoords[currConNum][currConPos] << " yCoords[" << currConNum << "][" << (currConPos) << "]=" << yCoords[currConNum][currConPos]  << endl;
#endif
            //Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
            //It also checks if there are any bulges or loops between the last 3 matched positions.
            //If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
            //constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
            //Otherwise it by default returns 'true'. 
            if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
#if defined DEBUG2_MODE
                cerr << "\nTRIPLET CONSTRAINT MATCHED\n";
#endif
                alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
                alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
                currConPos++;
                continue;
            }
            //If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
            else{
                // If not a match, step back, flip switches, and continue
                // Match found at lastX-1,lastY+1
                xCoords[currConNum][currConPos] = 0;
                yCoords[currConNum][currConPos] = 0;
                currConPos--;
                alreadyUsed[xCoords[currConNum][currConPos]] = false;
                alreadyUsed[yCoords[currConNum][currConPos]] = false;
                continue;
            }
        }
        
        // If not a match, step back, flip switches, and continue
		//	currConPos--;
		//	alreadyUsed[xCoords[currConNum][currConPos]] = false;
		//	alreadyUsed[yCoords[currConNum][currConPos]] = false;
		//	continue;
		//}
        */

        /*

/////////////////
//STAS: fixing the problem of non-matching constraint in a short hairpin structure with a bulge
/////////////////
//Look for bulge on non-first/last nucleotide
else if((conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-2][yCoords[currConNum][currConPos-1]+1]
|| conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1][yCoords[currConNum][currConPos-1]+2]) 
&& alreadyUsed[xCoords[currConNum][currConPos-1]-2] == false
&& alreadyUsed[yCoords[currConNum][currConPos-1]+1] == false){
                
cerr << "\nCHECK1\n";
//Looking for bulge in x-direction
if(conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-2][yCoords[currConNum][currConPos-1]+1]){
// Match found at lastX-2,lastY+1
xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-2;
yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
                    
#if defined(TRIPLET_VERBOSE_MODE)
cerr << "Location:AA Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
<< " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
<< " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
<< "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
//It also checks if there are any bulges or loops between the last 3 matched positions.
//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
//Otherwise it by default returns 'true'. 
if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
//cerr << "\nCHECK11\n";
alreadyUsed[xCoords[currConNum][currConPos-1]-2] = true;
alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
currConPos++;
continue;
}
//If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
else{
//cerr << "\nCHECK12\n";
// If not a match, step back, flip switches, and continue
// Match found at lastX-2,lastY+1
xCoords[currConNum][currConPos] = 0;
yCoords[currConNum][currConPos] = 0;
currConPos--;
alreadyUsed[xCoords[currConNum][currConPos]] = false;
alreadyUsed[yCoords[currConNum][currConPos]] = false;
continue;
}
}
                
//Looking for bulge in y-direction
else if(conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1][yCoords[currConNum][currConPos-1]+2]){
// Match found at lastX-1,lastY+1
xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+2;
                    
#if defined(TRIPLET_VERBOSE_MODE)
cerr << "Location:A Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
<< " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
<< " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
<< "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
//It also checks if there are any bulges or loops between the last 3 matched positions.
//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
//Otherwise it by default returns 'true'. 
if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
//cerr << "\nCHECK11\n";
alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
alreadyUsed[yCoords[currConNum][currConPos-1]+2] = true;
currConPos++;
continue;
}
//If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
else{
//cerr << "\nCHECK12\n";
// If not a match, step back, flip switches, and continue
// Match found at lastX-1,lastY+1
xCoords[currConNum][currConPos] = 0;
yCoords[currConNum][currConPos] = 0;
currConPos--;
alreadyUsed[xCoords[currConNum][currConPos]] = false;
alreadyUsed[yCoords[currConNum][currConPos]] = false;
continue;
}
}
}




/////////////////
//END STAS: fixing the problem of non-matching constraint in a short hairpin structure with a bulge
/////////////////
            
*/


        /*
        else{   
            // If not a match, step back, flip switches, and continue
            currConPos--;
            alreadyUsed[xCoords[currConNum][currConPos]] = false;
            alreadyUsed[yCoords[currConNum][currConPos]] = false;
            continue;
        }
        */
#ifdef DEBUG2_MODE
        cerr << "\nCHECK14\n";
#endif
        /*        
/////////////////
//STAS: fixing the problem of non-matching constraint in a short hairpin structure with a bulge
/////////////////
//Look for bulge on non-first/last nucleotide
else if((conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-2][yCoords[currConNum][currConPos-1]+1]
|| conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1][yCoords[currConNum][currConPos-1]+2]) 
&& alreadyUsed[xCoords[currConNum][currConPos-1]-2] == false
&& alreadyUsed[yCoords[currConNum][currConPos-1]+1] == false){
            
cerr << "\nCHECK1\n";
//Looking for bulge in x-direction
if(conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-2][yCoords[currConNum][currConPos-1]+1]){
// Match found at lastX-2,lastY+1
xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-2;
yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
                
#if defined(TRIPLET_VERBOSE_MODE)
cerr << "Location:AA Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
<< " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
<< " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
<< "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
//It also checks if there are any bulges or loops between the last 3 matched positions.
//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
//Otherwise it by default returns 'true'. 
if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
//cerr << "\nCHECK11\n";
alreadyUsed[xCoords[currConNum][currConPos-1]-2] = true;
alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
currConPos++;
continue;
}
//If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
else{
//cerr << "\nCHECK12\n";
// If not a match, step back, flip switches, and continue
// Match found at lastX-2,lastY+1
xCoords[currConNum][currConPos] = 0;
yCoords[currConNum][currConPos] = 0;
currConPos--;
alreadyUsed[xCoords[currConNum][currConPos]] = false;
alreadyUsed[yCoords[currConNum][currConPos]] = false;
continue;
}
}
            
//Looking for bulge in y-direction
else if(conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1][yCoords[currConNum][currConPos-1]+2]){
// Match found at lastX-1,lastY+1
xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+2;
                
#if defined(TRIPLET_VERBOSE_MODE)
cerr << "Location:A Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
<< " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
<< " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
<< "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
//It also checks if there are any bulges or loops between the last 3 matched positions.
//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
//Otherwise it by default returns 'true'. 
if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
//cerr << "\nCHECK11\n";
alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
alreadyUsed[yCoords[currConNum][currConPos-1]+2] = true;
currConPos++;
continue;
}
//If the triplet constraint hasn't matched, return the current position xCoords and yCoords to previous state	
else{
//cerr << "\nCHECK12\n";
// If not a match, step back, flip switches, and continue
// Match found at lastX-1,lastY+1
xCoords[currConNum][currConPos] = 0;
yCoords[currConNum][currConPos] = 0;
currConPos--;
alreadyUsed[xCoords[currConNum][currConPos]] = false;
alreadyUsed[yCoords[currConNum][currConPos]] = false;
continue;
}
}
}
        
        
        
        
/////////////////
//END STAS: fixing the problem of non-matching constraint in a short hairpin structure with a bulge
/////////////////
*/
    
#if defined DEBUG2_MODE
        cerr << "\nPASSED First Part\n";
#endif 
        // Next step is in-bounds - find out where the search left off
		// Next step is in-bounds - find out where the search left off
	
		if (xCoords[currConNum][currConPos] == 0 && yCoords[currConNum][currConPos] == 0){
			// This is the first time searching forward from this position
			// Set coords values to -1,-1 to act as placeholder
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			searchX = true;
			searchExtCoax = false;
			start = xCoords[currConNum][currConPos-1] - 1;
		}

		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
		         xCoords[currConNum][currConPos] <= xCoords[currConNum][currConPos-1]-1) {
			// We left off searching in the X direction

			if (xCoords[currConNum][currConPos-1] -1 < 1){ 
				// Last found a match at left edge of dot plot - start searching in Y direction
				xCoords[currConNum][currConPos] = -1;
				yCoords[currConNum][currConPos] = -1;
				continue;
			}

			// Otherwise, continue searching in X direction
			searchX = true;
			searchExtCoax = false;
			start = xCoords[currConNum][currConPos] - 1;
		}

		else if (yCoords[currConNum][currConPos] == -1) {
			// Search in X direction was completed

			if (yCoords[currConNum][currConPos-1] + 2 > sequenceLength) {
				// Cannot search in Y direction because last match is near bottom edge of dot plot
				// Time to search for external coaxial stacks
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
				xCoords[currConNum][currConPos] = sequenceLength+1;
				continue;
			}
			
			// Otherwise, start at +2 because +1 was already covered in searchX
			searchX = false;
			searchExtCoax = false;
			start = yCoords[currConNum][currConPos-1] + 2;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
		         yCoords[currConNum][currConPos] >= yCoords[currConNum][currConPos-1]+1) {
			// We left off searching in the Y direction

			if (yCoords[currConNum][currConPos] == sequenceLength){ 
				// Last found a match at bottom edge of dot plot - start searching for external coaxial stacks
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
				xCoords[currConNum][currConPos] = sequenceLength+1;
				continue;
			}

			// Otherwise, continue searching in Y direction
			searchX = false;
			searchExtCoax = false;
			start = yCoords[currConNum][currConPos] + 1;
		}

		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
		         xCoords[currConNum][currConPos] > yCoords[currConNum][currConPos]) {
			// We're currently searching for a external coaxial stack in the X direction
			searchX = true;
			searchExtCoax = true;
			start = xCoords[currConNum][currConPos] - 1;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
		         yCoords[currConNum][currConPos] != -1 && 
		         yCoords[currConNum][currConPos] < xCoords[currConNum][currConPos]) {
			// We're currently searching for a external coaxial stack in the Y direction
			searchX = false;
			searchExtCoax = true;
			start = yCoords[currConNum][currConPos] + 1;
		}


		else {
			cerr<<"ERROR: Unclassified looping condition in constraint matching.\n"
			    <<xCoords[currConNum][currConPos]<<" "<<yCoords[currConNum][currConPos]<<"\n";			
			cin>>i;break;
		}

		// Enter this loop if we need to search in the X direction
		if (searchX && !searchExtCoax) {

			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = -1;
				yCoords[currConNum][currConPos] = -1;
				continue;
			}
			
			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,lastY + 1}, y = lastY+ 1
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (yCoords[currConNum][currConPos-1]);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,1}, y = lastY + 1
			else stop = 1;

#if defined(DEBUG_MODE)
			cout << "start="<< start << " stop=" << stop << " Second\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter] ||
				    (xCoords[currConNum][currConPos-1]-(start-xIter) > 2 &&
				     xCoords[currConNum][currConPos-1]-(start-xIter) < 8)) continue;

#if defined(DEBUG_MODE)
				cout << "start-xIter=" << start-xIter << " Third";
#endif

				if (conArray[currConNum][currConPos]==convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1]
                    ){
                    //&& 
                    //(start-xIter == xCoords[currConNum][currConPos-1] - 1 ||
                    // mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff)) {
                    // Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = start-xIter;
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:B Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all X values at this Y coordinate
			// Time to continue the search in the Y direction
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			continue;
		}

		// Enter this loop if we need to search in the y direction 
		else if (!searchX && !searchExtCoax) {

			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used
				
				// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
				if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
					xCoords[currConNum][currConPos] = sequenceLength + 1;
					continue;
				}
				
				// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,lastX - 1}
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (xCoords[currConNum][currConPos-1]-1);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,sequenceLength}
			else stop = sequenceLength;

#if defined(DEBUG_MODE)
			cout << start << " " << stop << " Fourth\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter] ||
				    ((start+yIter)-yCoords[currConNum][currConPos-1] > 2 &&
				     (start+yIter)-yCoords[currConNum][currConPos-1] < 8)) continue;

#if defined(DEBUG_MODE)
				cout << "start+yIter=" << start+yIter << " Fifth";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter]
                    ){
                    //&&
                    //mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {
					
					// Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
					yCoords[currConNum][currConPos] = start+yIter;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:C Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}

			}

			// If we reach this point, we've iterated through all x and y values for this base pair

			// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
			if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
				xCoords[currConNum][currConPos] = sequenceLength + 1;
				continue;
			}

			// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}

		// Enter this loop if we're searching for an external coaxial stack in the X direction
		else if (searchX && searchExtCoax) {

			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
				yCoords[currConNum][currConPos] = 0;
				continue;
			}
			
			stop = (yCoords[currConNum][currConPos-1]+1);
			
#if defined(DEBUG_MODE)
			cout << "start=" << start << " stop=" << stop << " Sixth\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter]) continue;

#if defined(DEBUG_MODE)
				cout << "start-xIter=" << start-xIter << " Seventh";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1]
                    ){
                    //&& 
                    //mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff) {
                    // Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = start-xIter;
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:D Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all x values at this y coordinate
			// Time to continue the external coaxial stack search in the y direction
			xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
			yCoords[currConNum][currConPos] = 0;
			continue;
		}

		// Enter this loop if we need to search for external coaxial stacks in the y direction 
		else if (!searchX && searchExtCoax) {

			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used, time to reset, step back, flip, and continue 
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}
			
			stop = (xCoords[currConNum][currConPos-1]-1);

#if defined(DEBUG_MODE)
			cout << "start=" << start << " stop=" << stop << " Eighths\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter]) continue;

#if defined(DEBUG_MODE)
				cout << "start+yIter="<< start+yIter << " Nineth";
#endif

				if (conArray[currConNum][currConPos]==
				    convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter]
                    ){
                    //&&
                    //mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {
                    
                    //Match found: store coords, flip switches, step forward, and continue
					xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
					yCoords[currConNum][currConPos] = start+yIter;

#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "Location:E Constraint_Number:"<< currConNum << " Constraints_From/To:" << conArray[currConNum][currConPos-1] << "->" << conArray[currConNum][currConPos] 
					     << " Constraint_Position_From/To:" << currConPos-1 << "->" << currConPos 
					     << " Seqence_Position_From/To: xCoords[" << xCoords[currConNum][currConPos-1]<< "->" << xCoords[currConNum][currConPos] 
					     << "] Ycoord[" << yCoords[currConNum][currConPos-1]<< "->" << yCoords[currConNum][currConPos] << "]\n";
#endif
					//Run the 'TripletMatch' function. This function test if there is a triplet constraint for the current location on the sequence.
					//It also checks if there are any bulges or loops between the last 3 matched positions.
					//If there is a triplet constraint and there are no loops or bulges separating the constraint, it checks if the triplet
					//constraint matches the current matched helix, and returns 'true' for a MATCH, and 'false' for NON-MATCH.
					//Otherwise it by default returns 'true'. 
					
					if(TripletMatch(ct,tripletArray,xCoords,yCoords,currConNum,currConPos)){
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
					}
				}
			}

			// If we reach this point, we've iterated through all x and y values for this base pair
			// Time to reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}
	}

	// The current position has gone back to the first base pair of a constraint
	// Flip switches and return
	alreadyUsed[xCoords[currConNum][currConPos]] = false;
	alreadyUsed[yCoords[currConNum][currConPos]] = false;
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;

}

// Search dotplot for complicated pseudoknot folds to exclude
void pseudodptrim(structure *ct, int *tempbp, int * count) {
	int i, iter, iter2, iter3, init=1, a=0, b=0, c=0, d=0, e=0, f=0, j=0;
	bool rightSingle, findCF, leftSingle, findAD;
	while(init <= ct->GetSequenceLength()) { 
		if(tempbp[init] > init) {
			for(i = 1; i < init; i++) {
				if(tempbp[i] > init && tempbp[i] < tempbp[init]) {
					// Define stems of pseudoknotted region
					j = tempbp[init];
					e = tempbp[i];
					b = tempbp[j];

					rightSingle = false;
					findCF = false;
					f = e+1; // Initial guess
					while( findCF == false && rightSingle == false) {
						if( tempbp[f] != 0 && tempbp[f] < e) {
							c = tempbp[f];
							findCF = true;
						}
						else {
							f++;
							if( f >= j) {
								rightSingle = true;
								f = j;
								c = b;
							}
						}
					}

					leftSingle = false;
					findAD = false;
					a = b - 1; // Initial guess
					while( findAD == false && leftSingle == false) {
						if( tempbp[a] != 0 && tempbp[a] >  c) {
							d = tempbp[a];			
							findAD = true;
						}
						else {
							a--;
							if( a <= i) {
								leftSingle = true;
								a = i;
								d = e;
							}
						}
					}
					
#if defined(VERBOSE_MODE)
					cout << "Pseudoknot detected: " << 
						i << "-" << e << "," << a << "-" << d << " ; " <<
						c << "-" << f << "," << b << "-" << j << "\n";
#endif


					// Remove all possible base pairs that would involved pairing between the 2 gap regions
					int count = 0;
					for (iter = a+1; iter < b; iter++) {
						for (iter2 = e+1; iter2 < f; iter2++) {
							if (ct->tem[iter2][iter] == true) count++;
							ct->tem[iter2][iter] = false;
						}
					}

#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Excluded pairings between gaps: " << count << "\n";
#endif
				
					// Scan gaps for contained secondary structures - forbid this region from pairing outside the gap
					// a-b gap:
					count = 0;
					iter = a+1;
					while (iter < b ) {
						if (tempbp[iter] > iter && tempbp[iter] > a && tempbp[iter] < b) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= a; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = b; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in a-b gap: " << count << "\n";
#endif
					
					// c-d gap:
					count = 0;
					iter = c+1;
					while (iter < d ) {
						if (tempbp[iter] > iter && tempbp[iter] > c && tempbp[iter] < d) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= c; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = d; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in c-d gap: " << count << "\n";
#endif
				
					// e-f gap:
					count = 0;
					iter = e+1;
					while (iter < f ) {
						if (tempbp[iter] > iter && tempbp[iter] > e && tempbp[iter] < f) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= e; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = f; iter3 <= ct->GetSequenceLength(); iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
					
#if defined(VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in e-f gap: " << count << "\n";
#endif
					
					// Continue loop at 3' end of pseudoknot
					init = j;
					continue;
				}
			}				
		}
		// If no pseudoknot detected, step up one position and recheck
		init++;
		continue;
	}
	return;
}

int pairout(structure *ct, const char* pairsoutfile) {
	int iter, iter2, a, b, c, d, bpCount;
	ofstream outFile(pairsoutfile);
	if (!outFile.good())
		return 34;

	for (iter = 1; iter <= ct->GetNumberofStructures(); iter++) {
		//if (cutoff !=0 && ct->GetEnergy(iter) > cutoff) break;
		outFile << "#  Structure #" << iter << "; ENERGY = " << double(ct->GetEnergy(iter))/10 << "; " << 
			ct->GetCtLabel(iter) << endl;
		outFile << "1" << endl;
		for (iter2 = 1; iter2 <= ct->GetSequenceLength(); iter2++) outFile << ct->nucs[iter2];
		outFile << endl << endl;
		outFile << "Positions paired." << endl;
		iter2 = 1;
		while (iter2 <= ct->GetSequenceLength()) {
			if (ct->GetPair(iter2,iter) < iter2) {
				iter2++;
				continue;
			}
			else {
				a = iter2;
				b = ct->GetPair(iter2,iter);
				bpCount = 1;
				while (ct->GetPair(iter2+bpCount,iter) == (b-bpCount)) bpCount++;
				c = iter2 + bpCount - 1;
				d = ct->GetPair(iter2+(bpCount-1),iter);
				outFile << a << "-" << c << "; " << d << "-" << b << endl;
				iter2 = c + 1;
				continue;
			}
		}
		outFile << endl << "---------------------------------------" << endl << endl;
	}
	if (!outFile.good()) return 35;
	outFile.close();
	return 0;
}

bool TripletMatch(structure *ct, char**** tripletArray, short** xCoords, short** yCoords, int currConNum, int currConPos) {
	
	bool TripConCheck=true;
		
	//Check if the nucleotides in the constraint are not separated by bulges are loops 
	if(currConPos>2&&xCoords[currConNum][currConPos-2]-xCoords[currConNum][currConPos]==2&&yCoords[currConNum][currConPos]-yCoords[currConNum][currConPos-2]==2){
		//if there is a triplet constraint in the PREVIOUS position, this means that there is an active constraint
		if(tripletArray[currConNum][currConPos-1][0][0]!=0){
			//if the constraint is on XCoord side
			if(NucCompare(tripletArray[currConNum][currConPos-1][0][2],ct->numseq[xCoords[currConNum][currConPos-1]],0)){
#if defined(TRIPLET_VERBOSE_MODE)
				cerr << "    Triplet Constraint found on xCoord:\n    ";
#endif
				//if the nucleotides matchs the constraint
				if((NucCompare(tripletArray[currConNum][currConPos-1][0][1],ct->numseq[xCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][0][0])) &&
				   (NucCompare(tripletArray[currConNum][currConPos-1][0][3],ct->numseq[xCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][0][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n       Match at (x+1)[" << xCoords[currConNum][currConPos] << "]=" << ct->numseq[xCoords[currConNum][currConPos]] << " and (x-1)[" << xCoords[currConNum][currConPos-2] << "]=" << ct->numseq[xCoords[currConNum][currConPos-2]] << endl;
					//QUESTION: temporary line
					cerr << " foo " << tripletArray[currConNum][currConPos-1][0][0] << " " << xCoords[currConNum][currConPos] << "-" << xCoords[currConNum][currConPos-2] << endl;
#endif
					//if there is a second constraint for that match
					if(tripletArray[currConNum][currConPos-1][1][0]!=0){
                        //Check if both constraints are not "+" constraints
                        if((tripletArray[currConNum][currConPos-1][0][0]=='+'&&tripletArray[currConNum][currConPos-1][1][0]=='+')||(tripletArray[currConNum][currConPos-1][0][0]=='+'&&tripletArray[currConNum][currConPos-1][1][0]=='-')||(tripletArray[currConNum][currConPos-1][0][0]=='-'&&tripletArray[currConNum][currConPos-1][1][0]=='+')){
                            cout << "\nERORR in constraints file:\nIf triplet constraints come as a pair, both MUST be negative constraints.\nPlease check the constraints file.";
                            exit(1);
                        }

						//if the second constraint doesnt' match
						if(!(NucCompare(tripletArray[currConNum][currConPos-1][1][1],ct->numseq[xCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][1][0])) &&
						   !(NucCompare(tripletArray[currConNum][currConPos-1][1][3],ct->numseq[xCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][1][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
							cerr << "\n             !!! Second Constraint Did Not Match !!!";
#endif
							TripConCheck=false;//set TripConCheck to FALSE if the second constraint didn't match
						}//END: if the second constraint doesn't match
					}//END: if there is a second constraint
				}
				else{
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n                     !!! First Constraint Did Not Match !!!";
#endif
					TripConCheck=false;//set TripConCheck to false if the first constraint didn't match
				}//END: if the nucleotide in the third position matches the constraint
			}//END: if the constraint is on xCoord.

			//if the constraint is on yCoord side
			else if(NucCompare(tripletArray[currConNum][currConPos-1][0][2],ct->numseq[yCoords[currConNum][currConPos-1]],0)){
#if defined(TRIPLET_VERBOSE_MODE)
				cerr << "    Triplet Constraint found on yCoord:\n    ";
#endif
				//if the nucleotide matches the constraint
				if((NucCompare(tripletArray[currConNum][currConPos-1][0][3],ct->numseq[yCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][0][0])) &&
				   (NucCompare(tripletArray[currConNum][currConPos-1][0][1],ct->numseq[yCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][0][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n      Match at (y+1)[" << yCoords[currConNum][currConPos] << "]=" << ct->numseq[yCoords[currConNum][currConPos]] << " and (y-1)[" << yCoords[currConNum][currConPos-2] << "]=" << ct->numseq[yCoords[currConNum][currConPos-2]] << endl;
					//QUESTION: temporary lines
					cerr << " foo " << tripletArray[currConNum][currConPos-1][0][0] << " " << yCoords[currConNum][currConPos] << "-" << yCoords[currConNum][currConPos-2] << endl;
#endif
					//if there is a second constraint for that match
					if(tripletArray[currConNum][currConPos-1][1][0]!=0){
						//if the second constraint doesnt' match
						if(!(NucCompare(tripletArray[currConNum][currConPos-1][1][3],ct->numseq[yCoords[currConNum][currConPos]],tripletArray[currConNum][currConPos-1][1][0])) &&
						   !(NucCompare(tripletArray[currConNum][currConPos-1][1][1],ct->numseq[yCoords[currConNum][currConPos-2]],tripletArray[currConNum][currConPos-1][1][0]))){
#if defined(TRIPLET_VERBOSE_MODE)
							cerr << "\n             !!! Second Constraint Did Not Match !!!";
#endif
							TripConCheck=false;//set TripConCheck to FALSE if the second constraint didn't match
						}//END: if the second constraint doesn't match
					}//END: if there is a second constraint
				}
				else{
#if defined(TRIPLET_VERBOSE_MODE)
					cerr << "\n                     !!! First Constraint Did Not Match !!!";
#endif
					TripConCheck=false;//set TripConCheck to false if the first constraint didn't match
				}//END: if the nucleotide in the third position matches the constraint
			}//END: if the constraint is on xCoord.

#if defined(TRIPLET_VERBOSE_MODE)
			if(TripConCheck==true)cerr << "\n**********MATCH*********\n";
			else cerr << "\n___________NO MATCH__________\n";
#endif
			
		}//END: if there is a triplet constraint for the previous position
	}//END: if the nucleotides are not separated by loops or bulges

	return TripConCheck;
	/*Return TRUE if either: 1. There were no triplet constraint for either the current or previous positions
	  2. There was a bulge or a loop in the middle of the triplet constraint
	  3. There were constraints and they all matched
	
	  Return FALSE if:       1. The constraint(s) didn't match
	*/
}

bool NucCompare(char ConNuc,int SeqNuc, char ConSign){//Given a constraint and a nucleotide, and a constraint sigh, the program returns 'true' if there is a match and 'false' if not
	
	/* Structure class nomenclature:
       A = 1 (R)
       C = 2 (Y)
       G = 3 (R)
       U = 4 (Y)
	   NAPSS Triplet constraint nomenclature
       Y = Pyrimidines
       R = Purines
	*/
#if defined(TRIPLET_VERBOSE_MODE)
	cerr << "  Comparing Nucleotides: Constraint=" << ConNuc << " Sequence=" << SeqNuc << " Constraint_Sign=" << ConSign ;
#endif

	//If the input for the constraint sign was 0, this lets the program know that it is comparing nucleotides, and not constraints
	//This operation is performed when the program tries to determine where is (xCoords or yCoords) 'A' or 'G' nucleotide of the triplet constraint located
	if(ConSign==0){
		if(ConNuc=='A' && SeqNuc==1 || ConNuc=='G' && SeqNuc==3) return true;
		else return false;
	}
	//If the comparison sigh is a '+', this means that it is a positive constraint (Can be ONLY this, and nothing else).
	else if(ConSign=='+'){
		if(ConNuc=='Y' && (SeqNuc==2 || SeqNuc==4))	return true;//Check if the nucleotide is a pyrimidine
		else if(ConNuc=='R' && (SeqNuc==1 || SeqNuc==3)) return true;//Check if the nucleotide is a purine
		else if(ConNuc=='X') return true;//Check if there is an 'X' in the triplet constraint. This means that there is no data.
		else return false;
	}
	//If the comparison sigh is a '-', this means that it is a negotive constraint (Can be anything BUT this)
	else if(ConSign=='-'){
		if(ConNuc=='Y' && (SeqNuc!=2 || SeqNuc!=4)) return true;//Check if the nucleotide is NOT a pyrimidine
		else if(ConNuc=='R' && (SeqNuc!=1 || SeqNuc!=3)) return true;//Check if the nucleotide is NOT a purine
		else if(ConNuc=='X') return true;//Check if there is an 'X' in the triplet constraint. This means that there is no data.
		else return false;
	}

	//RMW: This was added to napss-update-v4 and may change results. Try with/without this correction if there are discrepancies. 
	return true;
}

// Check for helical extension possibilities
void HelicalExtension(structure* ct, short** convertedDotplot, int* helixExtend, short** dgArray){
	for (int i = 1; i <= ct->GetNumberofStructures(); i++) {
		bool extensionAdded = true;
		while (extensionAdded) {
			extensionAdded = false;
			// Initialize the temporary array
			for (int iter2 = 0; iter2 <= ct->GetSequenceLength(); iter2++) {
				helixExtend[iter2] = 0;
			}
			
			int iter2 = 1;
			while (iter2 <= ct->GetSequenceLength()) {
				if (ct->GetPair(iter2,i) < iter2) {
					iter2++;
					continue;
				}
				else {
					// Check to see if the nucleotides downstream can pair with each other
					if (iter2 > 1 && ct->GetPair(iter2,i) < ct->GetSequenceLength() &&
					    ct->GetPair(iter2-1,i) == 0 && ct->GetPair(ct->GetPair(iter2,i)+1,i) == 0) { 
						// Nucleotides are currently unpaired 
						if (convertedDotplot[iter2-1][ct->GetPair(iter2,i)+1] !=0) { 
							// They could form a valid, stable pair
							if (helixExtend[iter2-1] == 0 && helixExtend[ct->GetPair(iter2,i)+1] == 0) {
								// There are no other possible extensions to these nucleotides yet - store this one
								helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
								helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
							}
							else {
								// One or both nucleotides already have a possible extension - if both do, 
								// make no changes because one new pair shouldn't break up two other possible pairs
								if (helixExtend[iter2-1] != 0 && helixExtend[ct->GetPair(iter2,i)+1] == 0) {
									// iter2 already has a possible helical extension to it - determine if the latest
									// pairing is more favorable in the original dotplot
									if (dgArray[iter2-1][helixExtend[iter2-1]] < 
									    dgArray[iter2-1][ct->GetPair(iter2,i)+1]) {
										helixExtend[helixExtend[iter2-1]] = 0;
										helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
										helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
									}
								}
								if (helixExtend[ct->GetPair(iter2,i)+1] != 0 && helixExtend[iter2-1] == 0) {
									// the nucleotide opposite iter2 already has a possible helical extension to it - 
									// determine if the latest pairing is more favorable in the original dotplot
									if (dgArray[ct->GetPair(iter2,i)+1][helixExtend[ct->GetPair(iter2,i)+1]] < 
									    dgArray[ct->GetPair(iter2,i)+1][iter2-1]) {
										helixExtend[helixExtend[ct->GetPair(iter2,i)+1]] = 0;
										helixExtend[ct->GetPair(iter2,i)+1] = iter2-1;
										helixExtend[iter2-1] = ct->GetPair(iter2,i)+1;
									}
								}
							}
						}
					}

					// Check to see if the nucleotides upstream can pair with each other
					if (ct->GetPair(iter2+1,i) == 0 && ct->GetPair(ct->GetPair(iter2,i)-1,i) == 0) { 
						// Nucleotides are currently unpaired 
						if (convertedDotplot[iter2+1][ct->GetPair(iter2,i)-1] !=0) { 
							// They could form a valid, stable pair
							if (helixExtend[iter2+1] == 0 && helixExtend[ct->GetPair(iter2,i)-1] == 0) {
								// There are no other possible extensions to these nucleotides yet - store this one
								helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
								helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
							}
							else {
								// One or both nucleotides already have a possible extension - if both do, 
								// make no changes because one new pair shouldn't break up two other possible pairs
								if (helixExtend[iter2+1] != 0 && helixExtend[ct->GetPair(iter2,i)-1] == 0) {
									// iter2 already has a possible helical extension to it - determine if the latest
									// pairing is more favorable in the original dotplot
									if (dgArray[iter2+1][helixExtend[iter2+1]] < 
									    dgArray[iter2+1][ct->GetPair(iter2,i)-1]) {
										helixExtend[helixExtend[iter2+1]] = 0;
										helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
										helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
									}
								}
								if (helixExtend[ct->GetPair(iter2,i)-1] != 0 && helixExtend[iter2+1] == 0) {
									// the nucleotide opposite iter2 already has a possible helical extension to it - 
									// determine if the latest pairing is more favorable in the original dotplot
									if (dgArray[ct->GetPair(iter2,i)-1][helixExtend[ct->GetPair(iter2,i)-1]] < 
									    dgArray[ct->GetPair(iter2,i)-1][iter2+1]) {
										helixExtend[helixExtend[ct->GetPair(iter2,i)-1]] = 0;
										helixExtend[ct->GetPair(iter2,i)-1] = iter2+1;
										helixExtend[iter2+1] = ct->GetPair(iter2,i)-1;
									}
								}
							}
						}
					}
					iter2++;
					continue;
				}
			}

			// Now loop through the temporary array and add any remaining potential pairs to the final structure
			for (int iter2 = 1; iter2 <= ct->GetSequenceLength(); iter2++) {
				if (helixExtend[iter2] != 0) {
					// Adding a pair - flip the switch
					extensionAdded = true;
					ct->SetPair(iter2,helixExtend[iter2],i);
				}
			}
		}
	}
}

string napss_error_message(int errorcode) {
	switch(errorcode) {
		case NAPSS_BAD_CONSTRAINT_FILE:
			return "The NMR Constraints file was not found or could not be read.\n";
		case NAPSS_ERR_CONSTRAINT_LEN: 
			return "Constraints cannot be less than two basepairs in length!\n";
		case NAPSS_ERR_NO_MATCHES:
			return "NO MATCHES FOUND\n";
		case NAPSS_BAD_RNA_COPY:
			return "The original RNA object could not be copied.\n";
		default: 
			return sfmt("Unknown error code: %i\n", errorcode);
	}
}