#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <ctime> 
//Usually the number of sequences in one Multilign object is small; thus,
//the pseudo-random number generated in random_shuffle function from STL::algorithm
//is good enough.
//#include <ctime>
//#include "../src/random.h" //random number generator
#include "Dynalign_object.h"
#include "Multilign_object.h"
#include "../src/outputconstraints.h" //required for readconstraints function
#include "../src/defines.h"

//Utility function.
bool file_exists(const char* const fileName);

//default constructor
Multilign_object::Multilign_object():dsvFiles(NULL), aliFiles(NULL), progress(NULL), dynobj(NULL), maxPairs(0), maxDsv(0), iterations(0), SHAPESlope(0), SHAPEIntercept(0),thermo(false) {
}

//constructor
Multilign_object::Multilign_object(const vector<vector<string> > &inputlist, const bool isrna, ProgressHandler *Progress): 
	inputList(inputlist), dsvFiles(NULL), aliFiles(NULL), progress(Progress), 
	dynobj(NULL), thermo(isrna,isrna?DT_RNA:DT_DNA),
	maxDsv(1), iterations(2), SHAPESlope(1.8), SHAPEIntercept(-0.6)  {
	ErrorCode = thermo.ReadThermodynamic();
	maxPairs = AverageLength();
}

Multilign_object::Multilign_object(const bool Multifind ,const string &outputmultifind, const vector<string> &ctfiles, ProgressHandler *Progress, const bool isrna):  
	dsvFiles(NULL), aliFiles(NULL), progress(Progress), 
	dynobj(NULL), thermo(isrna,isrna?DT_RNA:DT_DNA), 
	maxDsv(1), iterations(2), SHAPESlope(2.6), SHAPEIntercept(-0.8), output_multifind(outputmultifind),ct_files(ctfiles) {
		ErrorCode = thermo.ReadThermodynamic();
		maxPairs = AverageLength();
}
//Destructor
Multilign_object::~Multilign_object(){
    if(dsvFiles != NULL) {
        for (int i = 0; i < iterations; ++i){
            delete  [] dsvFiles[i];
        }
        delete [] dsvFiles;
    }

    if (aliFiles != NULL) {
        for (int i = 0; i < iterations; ++i){
            delete [] aliFiles[i];
        }
        delete [] aliFiles;
    }
}

// This function check the legality of the input filenames and
// prepare the parameters for the multilign calculations
// This function should be called before multilign calcultions at least once and
// whenever something related to seq/ct changes, e.g. SetIndexSeq, AddOneInput, RemoveOneInput, Randomize, etc.
int Multilign_object::PrepInput() {
  if(inputList.size() < 2) {return 5002;}
    // check the files exist.
    for(vvs_it rowIt = inputList.begin(); rowIt != inputList.end(); ++rowIt) {
        vector<string>::iterator it = rowIt->begin();
        // it points to seq filename currently
        if(!ifstream(it->c_str())) return 5012;
        // it points to ct filename currently
        if((++it)->empty()) return 5013;
        // it points to constraint filename currently
        if(!(++it)->empty())
            if(!ifstream(it->c_str())) return 5001;
        // it points to SHAPE filename curently
        if(!(++it)->empty())
            if(!ifstream(it->c_str())) return 5011;
    }

    if(0!=(ErrorCode=PairSeq1())) return ErrorCode;
    return 0;
}

int Multilign_object::PrepMultifindInput(){
  if( (input_alignment.size()<2) || ( (ct_files.size()!=input_alignment.size()) && (ct_files.size()!=0) ) ){return 5002;}
   for(vector<string>::iterator it=input_alignment.begin();it!=input_alignment.end();++it){
     if(it->find_first_not_of("atcguATCGU-")!=string::npos)return 5019;
     if(it->find_first_of("atcguATCGU")==string::npos)return 5020;
   }
  if(ct_files.size()!=0)for(vector<string>::iterator it=ct_files.begin();it!=ct_files.end();++it)if(it->empty())return 5013;
  if(0!=(ErrorCode=PairMultifindSeq1())) return ErrorCode;
  return 0;
}
  
void Multilign_object::ToHead(vvs_it first, vvs_it middle) {
    vvs_it previous = middle, tmp;
    while(middle != first){
        //cout << *(middle->begin()) << "  " << *(previous->begin()) << endl;
        swap(*middle--, *--previous);
        //*tmp = *middle; *middle = *previous; *previous = *tmp;
        //cout << *(middle->begin()) << "  " << *(previous->begin()) << endl;
    }
}


int Multilign_object::SetIndexSeq(size_t indexSeq) {
    // the first seq is the index seq; no need to make any change
    if(!--indexSeq) return 0;
    // the indexSeq is out of range
    if(indexSeq < 0 || indexSeq >= inputList.size()) return 5005;
    // move the indicated index seq to the beginning of the vectors by rotating
    ToHead(inputList.begin(), inputList.begin() + indexSeq);
    return 0;
}

int Multilign_object::SetIndexSeq(const string seqname) {
    if(seqname.empty()) return 5017;
    vvs_it it;
    for(it = inputList.begin(); it != inputList.end(); ++it)
        if(*(it->begin()) == seqname){
            SetIndexSeq(it - inputList.begin() + 1);
            break;
        }
    if(it==inputList.end())
        return 5018;
    return 0;
}



string Multilign_object::GetIndexSeq() const {
    if (inputList.empty()) return string();
    return *(inputList.begin()->begin());
}


//// initialize data member seq_ct
//int Multilign_object::MatchSeqCt(){
  //seq_ct = new pstr [Seq_List.size()];
  //// if running out of memory, seq_ct is set to NULL
  //if(!seq_ct) return 6000;
  //// initialize the seq_ct dynamic array.
  //for(vstr_index i = 0; i != Seq_List.size(); ++i){
    //*(seq_ct + i) = make_pair(Seq_List[i], Ct_List[i]);
    ////cout << (seq_ct+j)->second << " " << j <<" 44l" << endl;
  //}
  //return 0;
//}


// initialize data member seqPair
 int Multilign_object::PairSeq1(){
   if (inputList.size() < 2) {return 5002;}
   vs_index i = 0, j = i;
   seqPair.clear();
   while(++j != inputList.size()){
     seqPair.push_back(make_pair(i, j));
   }
   
   return 0;
}
 
 int Multilign_object::PairMultifindSeq1(){
   if(input_alignment.size()<2)return 5002;
   vs_index i=0,j=i;
   seqPair.clear();
   while(++j != input_alignment.size()){
     seqPair.push_back(make_pair(i, j));
   }
   
   return 0;
}
   


//return error messages based on the code from GetErrorCode and other error codes
string Multilign_object::GetErrorMessage(const int error) const {
    if (error==0) return "No Error.\n";
    else if (error==5001) return "The input constraint file can't be opened.\n";
    else if (error==5002) {cout<<"really?";return "You input less than 2 seq filenames.  Multilign requires 3 or more filenames\n";}
    else if (error==5003) return "One pairwise alignment can't be read for creating multiple alignment.\n";
    else if (error==5004) return "Sequence name is not found for removing.\n";
    else if (error==5005) return "The number indicating the index seq is out of range.\n";
    else if (error==5006) return "The value of MaxPairs is illegally less than -1.\n";
    else if (error==5007) return "The value of iterations is illegally less than 1.\n";
    else if (error==5008) return "The value of maxdsvchange is illegally less than zero.\n";
    else if (error==5009) return "The value of maxdsvchange is illegally larger than 99.\n";
    else if (error==5010) return "You are adding an empty string as seq filename or ct filename.\n";
    else if (error==5011) return "The input SHAPE file can't be opened.\n";
    else if (error==5012) return "At least one input seq file cannot be opened.\n";
    else if (error==5013) return "At least one output ct filename is not specified.\n";
    else if (error==5014) return "At least one intermediate dsv file cannot be deleted.\n";
    else if (error==5015) return "At least one intermediate aout file cannot be deleted.\n";
    else if (error==5016) return "An empty string is provided as multiple alignment file name.\n";
    else if (error==5017) return "An empty string is set as the index sequence.\n";
    else if (error==5018) return "The sequence file name to be set as index is not found.\n";
    else if (error==5019) return "The sequence contains abnormal symbols.\n";
    else if (error==5020) return "The sequence has no nucleotides.\n";
    else if (error==6000) return "Ran out of memory.\n";
    else if (error<100) {
        //This error range includes errors that derive from underlying RNA class instances
        //This can occur because the public pointer for RNA classes are accessed from the Dynalign_object instance
        if (dynobj == NULL) return RNA::GetErrorMessage(error);
        else return dynobj->GetRNA1()->GetErrorMessage(error);
    }

    else if (error < 5000) {
        if (dynobj == NULL) return "Error occured in Dynalign_object class.\n";
        else return dynobj->GetErrorMessage(error);
    }
    else return "Unknown Error\n";
}

void Multilign_object::ResetError() {
    this->ErrorCode = 0;
}


int Multilign_object::NameMultifindDsvFiles() {
  stringstream ss("");
  dsvFiles = new string * [iterations];
  for ( int j =0; j< iterations; ++j){
    dsvFiles[j] = new string [seqPair.size()];
    for ( size_t i = 0; i < seqPair.size(); ++i){
      dsvFiles[j][i]=output_multifind;
      dsvFiles[j][i]+="_";
      ss<<j<<"_"<<seqPair[i].first<<"_"<<seqPair[i].second<<".dsv";
      dsvFiles[j][i]+=ss.str();
      ss.str("");
      }
    }
  return 0;
}

int Multilign_object::NameDsvFiles() {
    dsvFiles = new string * [iterations];
    //if(dsvFiles) return 6000;
    stringstream ss;
    for ( int j =0; j< iterations; ++j){
      dsvFiles[j] = new string [seqPair.size()];
      //if(dsvFiles[j]) return 6000;
      for ( size_t i = 0; i < seqPair.size(); ++i){
        string one, two;
    //cout << seq_pair[i].first << " " << seq_pair[i].second << endl;
        size_t index = string::npos;

        // get the path of the ct filename for the index seq
        string path = inputList[seqPair[i].first][1];
        index = path.find_last_of("\\/");
        while(index!=string::npos && path.substr(index-1, 2)=="\\/")
            index=path.find_last_of("\\/", index);
        if (index!=string::npos)
            path = path.substr(0, index + 1);
        else
            path = "";


        // the first seq in the pair
        one = inputList[seqPair[i].first][0];
        if ( ( index = one.rfind(".seq") ) == (one.length() - 4) )
            // remove ".seq" suffix
            one=one.substr(0, index);

        // the second seq in the pair
        two = inputList[seqPair[i].second][0];
        if ( ( index = two.rfind(".seq") ) == (two.length() -4) )
            // remove ".seq" suffix
            two = two.substr(0, index);

        index = one.find_last_of("\\/");
        while(index!=string::npos && one.substr(index-1, 2)=="\\/")
            index=one.find_last_of("\\/", index);
        if (index!=string::npos)
            dsvFiles[j][i] = one.substr(index+1);
    else dsvFiles[j][i] = one;

        index = two.find_last_of("\\/");
        while(index!= string::npos && two.substr(index-1, 2)=="\\/")
            index = two.find_last_of("\\/", index);
        if(index!=string::npos)
            dsvFiles[j][i] = dsvFiles[j][i] + "_" + two.substr(index+1);
    else dsvFiles[j][i] = dsvFiles[j][i] + "_" + two;

        ss.str(""); ss << i+1;
        dsvFiles[j][i] = ss.str() + "_" + dsvFiles[j][i] ;
        ss.str(""); ss << j+1;
        dsvFiles[j][i] = path + ss.str() + "." + dsvFiles[j][i]+".dsv";
    //cout << dsvFiles[j][i] << endl;
      }
    }

    return 0;
}

int Multilign_object::NameMultifindAliFiles() {
  stringstream ss("");
  aliFiles = new string * [iterations];
  for ( int j =0; j< iterations; ++j){
    aliFiles[j] = new string [seqPair.size()];
    for ( size_t i = 0; i < seqPair.size(); ++i){
      aliFiles[j][i]=output_multifind;
      aliFiles[j][i]+="_";
      ss<<j<<"_"<<seqPair[i].first<<"_"<<seqPair[i].second<<".ali";
      aliFiles[j][i]+=ss.str();
      ss.str("");
      }
    }
  return 0;
}

int Multilign_object::NameAliFiles() {
    aliFiles = new string * [iterations];
    //if(aliFiles) return 6000;
    stringstream ss;
    for ( int j =0; j< iterations; ++j){
      aliFiles[j] = new string [seqPair.size()];
      for ( size_t i = 0; i < seqPair.size(); ++i){
        string one, two;
        size_t index = string::npos;
        // get the path of the ct filename for the index seq
        string path = inputList[seqPair[i].first][1];
        index = path.find_last_of("\\/");
        while(index!=string::npos && path.substr(index-1, 2)=="\\/")
            index=path.find_last_of("\\/", index);
        if (index!=string::npos)
            path = path.substr(0, index + 1);
        else
            path = "";

        one = inputList[seqPair[i].first][0];
        if ( ( index = one.rfind(".seq") ) == (one.length() - 4) )
          one=one.substr(0, index);

        two = inputList[seqPair[i].second][0];
        if ( ( index = two.rfind(".seq") ) == (two.length() -4) )
          two = two.substr(0, index);

        index = one.find_last_of("\\/");
        while(index!=string::npos && one.substr(index-1, 2)=="\\/")
          index=one.find_last_of("\\/", index);
        if(index!=string::npos)
          aliFiles[j][i] = one.substr(index+1);
    else aliFiles[j][i] = one;

        index = two.find_last_of("\\/");
        while(index!= string::npos && two.substr(index-1, 2)=="\\/")
          index = two.find_last_of("\\/", index);
        if(index!=string::npos)
          aliFiles[j][i] = aliFiles[j][i] + "_" + two.substr(index+1);
    else aliFiles[j][i] = aliFiles[j][i] + "_" + two;

        ss.str(""); ss << i+1;
        aliFiles[j][i] = ss.str() + "_" + aliFiles[j][i] ;
        ss.str(""); ss << j+1;
        aliFiles[j][i] = path + ss.str() + "." + aliFiles[j][i]+".aout";
    //cout << aliFiles[j][i] << endl;
      }
    }

    return 0;
}


int Multilign_object::CountBP(const int i, const int j, const double percent) const {
    // templating the dsv file
    Dynalign_object *tmp = new Dynalign_object(dsvFiles[j][i].c_str());
    tmp->GetRNA1()->SetTemperature(thermo.GetTemperature());
    // set a cutoff to count the number of base pairs with lowest free energy below this cutoff
    int bp_num = 0;
    int cutoff = int(tmp->GetLowestEnergy() * percent);
    // length is the length of index sequences
#ifndef MULTIFIND
    int length = RNA(inputList[seqPair[i].first][0].c_str(), FILE_SEQ, &thermo).GetSequenceLength();
#else
    int length = input_sequences[seqPair[i].first].size();
#endif
    // count the number of base pairs with lowest free energy below this cutoff
    for (int i=1; i<=length; ++i){
        for (int j=i; j<=length; ++j){
            if (tmp->GetBestPairEnergy(1, i, j) < cutoff )
        ++bp_num;
        }
    }
    //cout << bp_num;
    //// now bp_num are the total base pair number needed to be reduced.
    //bp_num -= maxPairs;
    //if (seqPair.size() > 1)
        //// now the bp_num is the number gradually disallowed for each dynalign calculation in 1st cycle
        //bp_num /= (seqPair.size()-1);
    //cout << "bp_num: " << bp_num << endl;
    delete tmp;
    return bp_num;
}


// the core function doing progressive dynalign calculations and templating
int Multilign_object::ProgressiveMultilign(
        const short int numProcessors,
        const bool Dsv, const bool Ali,
        const short int maxtrace,
        const short int bpwin, const short int awin,
        const short int percent,
        const short int imaxseparation,
        const float gap,
#ifdef DYNALIGN_II
		const float slope,
		const float intercept,
		const int max_elongation,
#else
		const bool singleinsert,
#endif
        const short int singlefold_subopt_percent,
        const bool local){

    // prepare the object to be ready for the real calculation.
#ifndef MULTIFIND
  if(0!=(ErrorCode=PrepInput())) return ErrorCode;
#else
  if(0!=(ErrorCode=PrepMultifindInput())) return ErrorCode;
#endif
    PartialProgress pp(progress); // This PartialProgress object allows us to pass fractional dynalign updates through to the main progress object.
    // percent of progress
    int iterationProgressPercent = 100 / (seqPair.size() * iterations);
    int dsvProgressPercent = (int)(iterationProgressPercent * 0.2); // percent advanced for a dsv template
    int dynCalcProgressPercent = iterationProgressPercent - dsvProgressPercent; // percent advanced for a dynalign calculation
	
	if (!thermo.VerifyThermodynamic()) return (ErrorCode=5); // read datatables if not already read.

    //   dGIndex.resize(input_alignment.size()-1);
    //  energies.resize(input_alignment.size()-1);
    // if Dsv is true, .dsv files will be stored as
    // j.i_<seq1_name>_<seq2_name>.aout,
    // where j is the cycle number and i is which number of calculation in the cycle
#ifndef MULTIFIND 
    if (Dsv) NameDsvFiles();
#else
    if(Dsv) NameMultifindDsvFiles();
#endif
    // if Ali is true, .aout files will be stored as
    // aliFiles[j][i] = j.i_<seq1_name>_<seq2_name>.aout,
    // where j is the cycle number and i is which number of calculation in the cycle.
#ifndef MULTIFIND
    if (Ali) NameAliFiles();
#else
    if (Ali) NameMultifindAliFiles();
#endif
    int stepBP = 0, totalBP = 0;
    string tmpfilename;
    int struct_num;
    RNA *rna;
    for (int j = 0; j < iterations; ++j){
        //cout << "Multilign Iteration: " << (j+1) << endl;
        for ( size_t i = 0; i < seqPair.size(); ++i){
            //cout << "\nPair " << i+1 << " in cycle " << j+1 << ':' << endl;
            //cout << '\t' << inputList[seqPair[i].first][0] << "<==>" << inputList[seqPair[i].second][0] << endl;
#ifndef MULTIFIND      
            dynobj = new Dynalign_object(inputList[seqPair[i].first][0].c_str(), FILE_SEQ,
                                    inputList[seqPair[i].second][0].c_str(), FILE_SEQ, &thermo);
#else
	    //   cerr<<input_sequences[seqPair[i].first]<<"\n";
	    //   cerr<<input_sequences[seqPair[i].second]<<"\n";
		    
	    dynobj = new Dynalign_object(input_sequences[seqPair[i].first].c_str(),input_sequences[seqPair[i].second].c_str());
		dynobj->GetRNA1()->SetTemperature(thermo.GetTemperature());
	  
#endif	  
            

            // read constraint file for the first seq if it exists
#ifndef MULTIFIND
            if(!inputList[seqPair[i].first][2].empty()){
                //cout << "\tConstraint 1: " << inputList[seqPair[i].first][2] << endl;
                if(0!=(ErrorCode=dynobj->GetRNA1()->ReadConstraints(inputList[seqPair[i].first][2].c_str()))) return ErrorCode;
            }

            // read constraint file for the second seq if it exists
            if(!inputList[seqPair[i].second][2].empty()){
                //cout << "\tConstraint 2: " << inputList[seqPair[i].second][2] << endl;
                if(0!=(ErrorCode=dynobj->GetRNA2()->ReadConstraints(inputList[seqPair[i].second][2].c_str()))) return ErrorCode;
            }

            structure *ct;
            ct = dynobj->GetRNA1()->GetStructure();
            // read SHAPE file for the first seq if it exists
            if(!inputList[seqPair[i].first][3].empty()){
                //cout << "\tSHAPE 1: " << inputList[seqPair[i].first][3] << endl;
                ct->SHAPEslope = SHAPESlope * conversionfactor;
                ct->SHAPEintercept = SHAPEIntercept * conversionfactor;
                ct->ReadSHAPE(inputList[seqPair[i].first][3].c_str());
                // DEBUG_SHAPE: ct->WriteSHAPE(sfmt("multi-shape1_i%i_j%i.txt", i,j));
            }

            ct = dynobj->GetRNA2()->GetStructure();
            // read SHAPE file for the second seq if it exists
            if(!inputList[seqPair[i].second][3].empty()){
                //cout << "\tSHAPE 2: " << inputList[seqPair[i].second][3] << endl;
                ct->SHAPEslope = SHAPESlope * conversionfactor;
                ct->SHAPEintercept = SHAPEIntercept * conversionfactor;
                ct->ReadSHAPE(inputList[seqPair[i].second][3].c_str());
                // DEBUG_SHAPE: ct->WriteSHAPE(sfmt("multi-shape2_i%i_j%i.txt", i,j));
            }
#else
#endif
            pp.setNextStep(dsvProgressPercent); // specify that the next step (Templatefromdsv) comprises a small amount (dsvProgressPercent) of the total progress.
            // doing dsv templating.
            if (i!=0) tmpfilename = dsvFiles[j][i-1];
            // i == 0 && j != 0
            else if (j!=0) tmpfilename = dsvFiles[j-1][seqPair.size()-1];
            // the first Dynalign is not templated
        if (!(i==0 && j==0) ){
                //cout << "\tdsv template file: " << tmpfilename << endl;
                if(0!=(ErrorCode=dynobj->Templatefromdsv(tmpfilename.c_str(), maxDsv))) return ErrorCode;
            }

             // add the work from the previous step and specify that the next step (Dynalign) comprises a percentage (dynCalcProgressPercent) of the total progress.
            pp.advanceToNextStep(dynCalcProgressPercent);
            dynobj->GetRNA1()->SetProgress(pp); // causes Dynalign to use the partial-progress object (as a proxy for this->progress)

            // set up the number of base pairs to template
            if (i != 0 && j != 0 && j < iterations - 1) {
                totalBP =- stepBP;
            }
            else totalBP = maxPairs;
	    // cerr<<dsvFiles[j][i]<<"\n";
	    //	    cerr<<"i "<<i<<" j "<<j<<"\n";

            #ifdef DYNALIGN_II
            //cout << "Running DYNALIGN_II" << endl;
			if (0!=(ErrorCode = dynobj->Dynalign(maxtrace, bpwin, awin,
                                               percent, imaxseparation,
                                               slope, intercept, gap, max_elongation,
                                               dsvFiles[j][i].c_str(),
                                               false, // Multilign has to set the optimal_only to false
                                               singlefold_subopt_percent,
                                               local, numProcessors,
                                               totalBP)) ){
			#else
            //cout << "Running DYNALIGN" << endl;
			if (0!=(ErrorCode = dynobj->Dynalign(maxtrace, bpwin, awin,
                                               percent, imaxseparation,
                                               gap, singleinsert,
                                               dsvFiles[j][i].c_str(),
                                               false, // Multilign has to set the optimal_only to false
                                               singlefold_subopt_percent,
                                               local, numProcessors,
                                               totalBP)) ){
			#endif

                //cout << "ERROR(cycle " << ++j << ' ' << (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): " << dynobj->GetErrorMessage(ErrorCode);
                //cout << "ERROR " << ErrorCode << endl;
                return ErrorCode;
            }
	        //cerr<<"we are out!\n";
            pp.stepComplete(); // update progress to indicate that the previous step (Dynalign) is complete.

            if (i==0 && j==0) {
                totalBP = CountBP();
                if (iterations==1)
                    stepBP = (totalBP - maxPairs) / (iterations * seqPair.size() );
                else
                    stepBP = (totalBP - maxPairs) / ( (iterations-1) * seqPair.size() );
            }

            dynobj->WriteAlignment( aliFiles[j][i].c_str() );

            // last iteration
            if (j == iterations-1){
                //cout << "Last iteration" << endl;
                if (i==seqPair.size()-1){
                    rna = dynobj->GetRNA1();
                    struct_num = rna->GetStructureNumber();
                    // calculate free energyies and write into ct variables
                    for(int n = 1; n <= struct_num; ++n) {
                        rna->CalculateFreeEnergy(n);
                    }
                    //cout << "295l   " << struct_num << endl;
#ifndef MULTIFIND
                    // cout << "WriteCt iteration i=" << i << " "  << flush;
                    // cout << "first=" << seqPair[i].first << " "   << flush;
                    // cout << "inputList=" << inputList[seqPair[i].first][1] << " "   << endl;
                    if(0!=(ErrorCode=dynobj->GetRNA1()->WriteCt( inputList[seqPair[i].first][1].c_str() ))) {
            		      return ErrorCode;
            		    }
#else
		    if(ct_files.size())if(0!=(ErrorCode=dynobj->GetRNA1()->WriteCt( ct_files[seqPair[i].first].c_str() ))) {

                        return ErrorCode;
                    }
#endif
                }
                ///// This is for debug/test.//////////
                ///// output the ct files of the index sequence in the ith pair of dynalign calculations.
                //else {
                    //rna = dynobj->GetRNA1();
                    //struct_num = rna->GetStructureNumber();
                    //// calculate free energyies and write into ct variables
                    //for(int n = 1; n <= struct_num; ++n) {
                        //rna->CalculateFreeEnergy(n);
                    //}
                    //if(0!=(ErrorCode=dynobj->GetRNA1()->WriteCt( inputList[seqPair[i].first][1] + '1').c_str() ))) return ErrorCode;
                //}
                ///// This is for debug/test./////////
		dGIndex.push_back(dynobj->GetRNA1()->CalculateFreeEnergy(1));
                rna = dynobj->GetRNA2();
                struct_num = rna->GetStructureNumber();
                // calculate free energyies and write into ct variables
                for(int n = 1; n <= struct_num; ++n) {
		  if (n==1)
    	            energies.push_back(rna->CalculateFreeEnergy(n));
		  else
                    rna->CalculateFreeEnergy(n);
                }
#ifndef MULTIFIND	      
                if(0!=(ErrorCode=dynobj->GetRNA2()->WriteCt( inputList[seqPair[i].second][1].c_str() ))) return ErrorCode;
#else 
		//	cerr<<ct_files.size()<<"\n";
		if(ct_files.size())
		if(0!=(ErrorCode=dynobj->GetRNA2()->WriteCt( ct_files[seqPair[i].second].c_str() ) )) return ErrorCode;
#endif	       
            }

            delete dynobj;

        }
    }

    //if(0!=(ErrorCode=WriteAlignment(allali))) return ErrorCode;
    pp.setWorkComplete(100); // set the progress to 100%

    return 0;
}

int Multilign_object::WriteAlignment(const string allali) const {
    // build multiple alignment based on pairwise alignments.
    // currently the written ali files in disk are read in to build multiple alignment.
    // if I know how to use the align array before they are written, it'll be better.
    if(allali.empty()) return 5016;
    int last_cycle = iterations - 1;
    vector<string> alignments;
    vector<string>::iterator beg, end;
    string sequence1, sequence2;
    ifstream ali_in;
    for (size_t i = 0; i < seqPair.size(); ++i) {
        ali_in.open(aliFiles[last_cycle][i].c_str());
        if (!ali_in){
            return 5003;
            //GetErrorMessage(ErrorCode);
            //exit(ErrorCode);
        }
        // first line is a comment
        getline(ali_in, sequence1);
        // get the alignments for sequence 1 and sequence 2
	vector<string> temp_pair_alignment;
	getline(ali_in, sequence1);
	//	cerr<<sequence1;
	temp_pair_alignment.push_back(sequence1);
	getline(ali_in, sequence2);
	//	cerr<<sequence2;
	temp_pair_alignment.push_back(sequence2);
	pair_alignments.push_back(temp_pair_alignment);
        if (alignments.empty()){
            alignments.push_back(sequence1);
            alignments.push_back(sequence2);
        }
        else {
            // Use a very simple alignment model:
            // whenever there is a gap needed, a gap will be inserted.
            for (size_t index = 0; index < sequence1.size() || index < alignments[0].size(); ++index){
                if (alignments[0].size() <= index)
                    for(beg = alignments.begin(); beg != alignments.end(); ++beg)
                        *beg += '-';
                else if (sequence1.size() <= index) {
                    sequence1 += '-';
                    sequence2 += '-';
                }
                else {
                    if (alignments[0][index] != sequence1[index]){
                        if (alignments[0][index] == '-'){
                            sequence1.insert(index, "-");
                            sequence2.insert(index, "-");
                        }
                        else {
                            for(beg = alignments.begin(); beg != alignments.end(); ++beg){
                                beg->insert(index, "-");
                            }
                        }
                    }
                }
            }
            alignments.push_back(sequence2);
        }
        ali_in.close();
    }

    ofstream out_ali(allali.c_str());
    //  int i=0;
#ifndef MULTIFIND      
    vvs_cit indexIt;
    for(beg = alignments.begin(), indexIt = inputList.begin();
            beg != alignments.end(); ++beg, ++indexIt){
     
      out_ali << *(indexIt->begin()) << ' '<<*beg;
      //	if (i==0)
  
      //	  for(int j=0; j < inputList.size() - 1; ++j)
#else
     for(beg = alignments.begin();
	     beg != alignments.end(); ++beg){
       //      cerr<<*beg;
	 out_ali << *beg ;
	 //	 if (i==0)
	 //	   for(int j=0; j < input_alignment.size() - 1; ++j)
#endif	   
	 //    out_ali <<' '<<dGIndex[j] << ' ';
	 //	else
	 //  out_ali << energies[i-1];
	 	out_ali << endl;
	//	++i;
    }

    return 0;

}


// To-Do: another style of Multilign calculation
int Multilign_object::MultiTempMultilign(){

    return 0;

}


int Multilign_object::SetMaxPairs(const int maxpairs){
    if(maxpairs < -1) return 5006;
    else if(maxpairs == -1) maxPairs = AverageLength();
    else maxPairs = maxpairs;
    return 0;
}


int Multilign_object::GetMaxPairs() const {
    return maxPairs;
}


int Multilign_object::AverageLength() const {
    int nt=0;
    vector<vector<string> >::const_iterator it;
    if(inputList.size()==0) return 0;// if there is no input seq at all.
    for (it = inputList.begin(); it != inputList.end(); ++it){
		RNA seq(it->begin()->c_str(), FILE_SEQ, &thermo);
		nt = nt + seq.GetSequenceLength();
    }
    return nt/inputList.size(); // the average length of all the input sequences
}


int Multilign_object::SetIterations(const int it){
    if(it < 1) return 5007;
    iterations = it;
    return 0;
}


int Multilign_object::GetIterations() const {
    return iterations;
}


int Multilign_object::SetMaxDsv(const float maxdsvchange){
    if(maxdsvchange <= 0) return 5008;
    if(maxdsvchange > 99) return 5009;
    maxDsv = maxdsvchange;
    return 0;
}


float Multilign_object::GetMaxDsv() const {
    return maxDsv;
}


int Multilign_object::GetSequenceNumber() const{
    return inputList.size();
}

void Multilign_object::SetTemperature(const double temp){
    thermo.SetTemperature(temp);
}

double Multilign_object::GetTemperature() const {
    return thermo.GetTemperature();
}


int Multilign_object::AddOneInput(const string seq, const string ct, const string constraint, const string shape) {
    if(seq.empty() || ct.empty()) return 5010;
    vector<string> tmp;
    tmp.push_back(seq); tmp.push_back(ct);
    tmp.push_back(constraint); tmp.push_back(shape);
    inputList.push_back(tmp);
    return 0;
}


int Multilign_object::RemoveOneInput(const string seq) {
    vvs_it it;
    bool found = false;
    for(it = inputList.begin(); it != inputList.end(); ++it) {
        if (*(it->begin()) == seq) {
            it = --inputList.erase(it);
            found = true;
        }
    }
    if (!found) return 5004;
    return 0;
}

void Multilign_object::SetSHAPESlope(const double slope) {
    SHAPESlope = slope;
}

double Multilign_object::GetSHAPESlope() const {
    return SHAPESlope;
}

void Multilign_object::SetSHAPEIntercept(const double intercept) {
    SHAPEIntercept = intercept;
}

double Multilign_object::GetSHAPEIntercept() const {
    return SHAPEIntercept;
}

//void Multilign_object::Randomize() {
    ////using DHM's random number generator
    //randomnumber rand;
    //rand.seed(time(NULL));
    //string tmpSeq, tmpCt;
    //size_t Seq_Num = Seq_List.size();
    //for (size_t i = 1; i < Seq_Num; ++i){
        //size_t r = i + size_t ( rand.roll() * Seq_Num ) % ( Seq_Num - i );
        //// randomize the order except the first seq/ct, i.e. the main/index seq
        //tmpSeq = Seq_List[i];
        //Seq_List[i] = Seq_List[r];
        //Seq_List[r] = tmpSeq;

        //tmpCt = Ct_List[i];
        //Ct_List[i] = Ct_List[r];
        //Ct_List[r] = tmpCt;
    //}
    //if(0!=(ErrorCode=MatchSeqCt())) return;
    //if(0!=(ErrorCode=PairSeq())) return;
//}


void Multilign_object::Randomize() {
	// randomly shuffle the order except the first 1, i.e. the index seq.
	std::srand(std::time(0));
	random_shuffle( ++inputList.begin(), inputList.end());
}

// Output the filenames to stdout in the order of progressive dynalign calculations.
// The first one are the index seq.
void Multilign_object::GetInputFilenames() {
    if(0!=(ErrorCode=PrepInput()))
        cout << GetErrorMessage(ErrorCode);

    cout << "Set Seq\tCt\tConstraint\tSHAPE:\n";
    for(vvs_cit rowIt = inputList.begin(); rowIt != inputList.end(); ++rowIt) {
        cout << "    ";
        for(vector<string>::const_iterator colIt = rowIt->begin();
                colIt != rowIt->end(); ++colIt) {
            cout << *colIt << " ";
        }
        cout << endl;
    }
}


void Multilign_object::GetPairs() {
    if(0!=(ErrorCode=PrepInput()))
        cout << GetErrorMessage(ErrorCode);

    cout << "Sequences are paired:\n";
    vector<pair<vs_index, vs_index> >::const_iterator it;
    for(it = seqPair.begin(); it != seqPair.end(); ++it){
        cout << inputList[it->first][0] << " <==> " << inputList[it->second][0] << endl;
    }
}

int Multilign_object::CleanupIntermediateFiles() const {
    if (dsvFiles != NULL){
        for ( int j =0; j< iterations; ++j)
            for ( size_t i = 0; i < seqPair.size(); ++i){
                if (remove(dsvFiles[j][i].c_str()) && file_exists(dsvFiles[j][i].c_str())) {
                    //cerr << "Cannot delete " << dsvFiles[j][i].c_str() << "  -- ErrNo: " << errno << endl;
                    //cerr << "File Exists: " << file_exists(dsvFiles[j][i].c_str()) << endl;
                    //return 5014;
                }
            }
    }
    if (aliFiles != NULL){
        for ( int j =0; j< iterations; ++j)
            for ( size_t i = 0; i < seqPair.size(); ++i){
                if (remove(aliFiles[j][i].c_str()))
                    return 5015;
            }
    }

    return 0;
}

// REMOVED: The nucleic acid type should not be changed after creation, because the thermo table has already been loaded.
//void Multilign_object::SetNucType(const bool isrna) {
//    isRNA = isrna;
//}
//
//bool Multilign_object::GetNucType() const {
//    return isRNA;
//}


void Multilign_object::SetProgress(ProgressHandler *Progress){
    progress = Progress;
}

void Multilign_object::StopProgress(){
    progress = NULL;
}

ProgressHandler* Multilign_object::GetProgress() const {
    return progress;
}

bool file_exists(const char* const fileName) {
    if (FILE *file = fopen(fileName, "r")) {
        fclose(file);
        return true;
    }
    return false;
}