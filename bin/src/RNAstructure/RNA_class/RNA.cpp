
//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications
// #define ENABLE_DEBUG_LOGS


#include "../src/debug_logging.h"

#include "RNA.h"

#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/DynProgArray.h"
#include "../src/forceclass.h"
#include "../src/pfunction.h"
#include "../src/boltzmann.h"
#include "../src/outputconstraints.h"

#ifndef PARTITION_CORE // exclude these from minimal partition-only builds (e.g. partition-rosetta)
    #include "../src/random.h"
    #include "../src/dotarray.h"
    #include "../src/stackclass.h"
    #include "../src/stackstruct.h"
    #include "../src/algorithm.h"
    #include "../src/alltrace.h"
    #include "../src/stochastic.h"
    #include "../src/MaxExpect.h"
    #include "../src/probknot.h"
    #include "../src/draw.h"
    #include "RsampleData.h"
#endif

#ifdef _CUDA_CALC_
extern "C"{
#include "../partition-smp/param.h"
#include "../partition-smp/prna.h"
#include "../partition-smp/util.h"
#include "../partition-smp/base.h"
}
#endif//_CUDA_CALC

#include <iostream>
#include <iomanip>

const float epsilon = 1e-6; // a small number for a tolerance in comparing floats

// Common constuctor logic:
//   1. set member fields to defaults.
//   2. If alphabetName is not NULL, call ReadThermodynamic to read the datatables
//   3. If sequenceOrFileName is not NULL:
//        - If type is SEQUENCE_STRING, load the sequence contained in sequenceOrFileName directly.
//        - Otherewise load the file specified by sequenceOrFileName.
void RNA::init(const char * sequenceOrFileName, const RNAInputType fileType, const bool allowUnknownBases /* default: false */, const bool skipThermoTables /* default: false */) {
    //set error status to zero
    ErrorCode=0;
    lastErrorDetails = "";

    //allocate ct
    ct = new structure();

    //Indicate that the partition function calculation has not been performed.
    partitionfunctionallocated = false;

    //Indicate that the energy data is not read.
    energyallocated = false;

    //Drawing coordinates have not been determined.
    drawallocated = false;

    //Do not report progress by default:
    progress = NULL;

    // Do not call ReadThermodynamic if any of the following are true:
    //  - The table is already loaded (which would also be the case if Thermodynamics was copied.)
    //  - No alphabet has been specified. (i.e. with the RNA(IsRNA) constructor)
    //  - We are loading a sav file (PFS or FSV/SAV)
    bool readThermo = !(IsAlphabetRead() || GetAlphabetName().empty() || 
                          fileType==FILE_PFS || fileType==FILE_SAV);
    if (readThermo) {
        this->skipThermoTables = skipThermoTables; // set this property in the Thermodynamics base class. If true, only the alphabet will be loaded..Not the energy files.
        //Now read the thermodynamic parameters immediately, so that the alphabet is known
        // Note that this uses the directory, alphabetName, and temperature already set in the inherited Thermodynamics object.
        if ((ErrorCode=ReadThermodynamic())!=0) return; //The thermodynamic parameter files were not found, so report an error.
        // Thermodynamics->data should now be loaded
        data->allowUnknownBases = allowUnknownBases;
    }
    //Now also add the pointer to the data_tables to the structure class:
    if (data!=NULL)ct->SetThermodynamicDataTable(data);

    if (sequenceOrFileName == NULL) return; // do NOT load sequence information if sequenceOrFileName is NULL

    // Load the sequence or file
    if (fileType == SEQUENCE_STRING) { // type==SEQUENCE_STRING indicates that sequenceOrFileName represents a sequence and should be loaded directly.
        // SetSequence returns 0 if successful. Otherwise it returns an error code (and may also set lastErrorDetails)
        ErrorCode = ct->SetSequence(sequenceOrFileName);
    } else {
        // sequenceOrFileName represents the path to a file
        ErrorCode = FileReader(sequenceOrFileName, fileType);
    }
}

//constructor where user provides a string with the sequence
RNA::RNA(const char sequence[], const bool IsRNA):Thermodynamics(IsRNA, IsRNA?DT_RNA:DT_DNA) {
    init(sequence, SEQUENCE_STRING);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

//RNA::RNA(const char sequence[], const char* const alphabetName, const bool allowUnknownBases) {
//  init(sequence, SEQUENCE_STRING, alphabetName, allowUnknownBases);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
//}

//  This constructor is deprecated. Instead use RNA(const char filename[], const RNAFileType type, const bool IsRNA)
//  For type: 1=>CT_FILEtype, 2 => FILE_SEQ, 3 => FILE_PFS, 4 => FILE_SAV, and 5 => FILE_DBN
RNA::RNA(const char filename[], const RNAInputType type, const bool IsRNA ):Thermodynamics(IsRNA, IsRNA?DT_RNA:DT_DNA) {
    init(filename, (RNAInputType)type);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

RNA::RNA(const char filepathOrSequence[], const RNAInputType type, const char* const alphabetName, const bool allowUnknownBases, const bool skipThermoTables) 
    : Thermodynamics(isAlphabetRNA(alphabetName), alphabetName) {
    init(filepathOrSequence, type, allowUnknownBases, skipThermoTables);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

// Copy Thermodynamics info and datatable from another Thermodynamics instance.
RNA::RNA(const char filepathOrSequence[], const RNAInputType type, const Thermodynamics *copyThermo) : Thermodynamics(*copyThermo) {
    init(filepathOrSequence, type);  // loads the sequence, but does not call ReadThermo if Thermodynamics is already loaded.
}

// Default constructor.
RNA::RNA(const bool IsRNA):Thermodynamics(IsRNA, NULL) {
    // Allocate the underlying structure class and nothing more.
    
    // The Thermodynamics alphabet name was set to NULL, so the thermodynamic parameter tables will NOT be loaded in init.
    // So the programmer is required to call ReadThermodynamic and set the sequence and structural information explicitly.
    // Setting the sequence to NULL in init indicates that no sequence should be loaded.
    init(NULL, SEQUENCE_STRING); 
}

// Override Thermodynamics::CopyThermo so we can copy the datatable to our CT.
void RNA::CopyThermo(Thermodynamics& copy) {
    Thermodynamics::CopyThermo(copy);
    ct->SetThermodynamicDataTable(copy.GetDatatable());
}

//Return the value of ErrorCode
int RNA::GetErrorCode() const {
    return ErrorCode;
}

//Return a c string that describes errors from GetErrorCode and other errors.
const char* RNA::GetErrorMessage(const int error) {
    switch(error) {
        case 0: return "No Error.\n";
        /* Standard File Errors */
        case 1: return "Input file not found.\n";
        case 2: return "Error opening file.\n"; // This should be returned to indicate a file IO error, NOT an error with the content of the file, which should be #28 or #29.
        /* Misc Errors */
        case 3: return "Structure number out of range.\n";
        case 4: return "Nucleotide number out of range.\n";
        case 5: return "Error reading thermodynamic parameters.\n";
        case 6: return "This would form a pseudoknot and is not allowed.\n";
        case 7: return "This pair is non-canonical and is therefore not allowed.\n";
        case 8: return "Too many restraints specified.\n";
        case 9: return "This nucleotide already under a conflicting constraint.\n";
        case 10: return "There are no structures to write to file.\n";
        case 11: return "Nucleotide is not a U.\n";
        case 12: return "Maximum pairing distance is too short.\n";
        case 13: return "Error reading constraint file.\n";
        case 14: return "A traceback error occurred.\n";
        case 15: return "No partition function data is available.\n";
        case 16: return "Wrong save file version used or file format not recognized.\n";
        case 17: return "This function cannot be performed unless a save file (.sav) was correctly loaded by the RNA constructor.\n";
        case 18: return "This threshold is too low to generate valid secondary structures.\n";
        case 19: return "The structure coordinates have not been determined, use DetermineDrawingCoordinates() to calculate the coordinates.\n";
        case 20: return "No sequence has been read.\n";
        case 21: return "Probabilities summed to greater than 1 in stochastic traceback.\n";
        case 22: return "Programming error.  Incorrect file type passed to constructor.\n";
        case 23: return "There are no structures present.\n";
        case 24: return "Too few iterations.  There must be at least one iteration.\n";
        case 25: return "Index is not a multiple of 10.\n";
        case 26: return "k, the equilibrium constant, needs to be greater than or equal to 0.\n";
        case 27: return "Lyngso O(N^3) internal loop search is not compatible with a parallel calculation.\n";  
        case 28: return "Error reading sequence.\n";
        case 29: return "Invalid file format.\n";
        case 30: return "Programming error: The thermodynamic parameters have not been read.\n";
        case 31: return "Length mismatch between sequence and annotation file.\n"; // Tried to load probability annotations from a PFS file obtained from folding a shorter sequence.
        case 32: return "Array size mismatch.\n"; // Caller passed in an array that is too small to hold requested information.
        case 33: return "Error opening pseudoknot penalty constants file.\n"; // Caller passed in an array that is too small to hold requested information.
        case 34: return "Error opening output file for writing.\n"; 
		case 35: return "Error writing output file.\n"; 
		case 36: return "Pairs must have probability greater than zero.  Therefore, the probknot threshold must be >= 0.";
        /* SHAPE, Experimental Pair Bonus, etc Restraint data */
        case 201: return "Restraint File Not Found (SHAPE or other experimental data).\n";
        case 202: return "Could not open or read Restraint file (SHAPE or other experimental data).\n";
        case 203: return "Wrong number of restraints in file (SHAPE or other experimental data).\n";
        case 204: return "Invalid Nucleotide Number in Restraint file (SHAPE or other experimental data).\n";
        case 215: return "This function is incompatible with other restraints (SHAPE or other experimental data).\n"; // for Rsample.
        /* General Errors */
        case 99: return "The calculation was canceled.\n";
        default: return "Unknown Error\n";
    }
}

const string RNA::GetErrorDetails() const {
    return lastErrorDetails.empty() ? ct->GetErrorDetails() : lastErrorDetails;
}
void RNA::SetErrorDetails(const string& details) {
    lastErrorDetails = details;
}

void RNA::SetSequenceLabel(const string& label) {
    GetStructure()->SetSequenceLabel(label);
}

//Return a string that describes errors from GetErrorCode and other errors.
//This uses GetErrorMessage to acrually get the errors.
string RNA::GetErrorMessageString(const int error) const {
    return GetErrorMessage(error);
}

//! If there was an error, this returns the error message, along with error details (if any).
//! if there was no error (i.e. GetErrorCode() returns 0 and GetErrorDetails() returns NULL), this function returns an empty string ("");
string RNA::GetFullErrorMessage() const {
    int code = GetErrorCode();
    string message(code==0?"":GetErrorMessage(code));
    string details = GetErrorDetails();

    // If message and details are both non-empty, 
    // combine them to "<message>: <details>"
    // This requires triming whitepace and the dot/period (.) from the end of message.
    if (!message.empty() && !details.empty()) {
        std::size_t last = message.find_last_not_of("\r\n\t .");
        if (last != string::npos)
            message.resize(last+1);
        message.append(": ");
    }
    message.append(details);
    // GetErrorMessage always ends with \n, so this should also, for consistency.
    if (!message.empty() && message[message.length()-1]!='\n')
        message+='\n';
    return message;
}

void RNA::ResetError() {
    this->ErrorCode = 0;
    this->lastErrorDetails="";
}


// Ensure that at a minumum number of structures have been created.
void RNA::EnsureStructureCapcacity(const int minimumStructures) {
    if (minimumStructures>ct->GetNumberofStructures()) {
        //Add one structure for each position between those available and the one being specified:
        for (int index=ct->GetNumberofStructures()+1;index<=minimumStructures;++index) ct->AddStructure();
    }
}

//User specifies a base pair between i and j in structure # structurenumber.
int RNA::SpecifyPair(const int i, const int j, const int structurenumber) {
    


    //start with error checking:
    if (i<0||i>ct->GetSequenceLength()||j<0||j>ct->GetSequenceLength()) return 4;
    else if (structurenumber<1) return 3;

    //also keep track of the maximum number of structures for which pairs have been specified
    EnsureStructureCapcacity(structurenumber);

    //now register the pair:
    ct->SetPair(i,j,structurenumber);


    return 0;

}
//Break a pair that i is involved in
int RNA::RemoveBasePair(const int i, const int structurenumber) {

    //start with error checking:
    if (i<0||i>ct->GetSequenceLength()) return 4;
    else if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 3;

    //Call the function for this in the underlying structure class.
    ct->RemovePair(i, structurenumber);


    //return that there was no error
    return 0;

}

double RNA::GetVprimeQ(const int i, const int j){
    cout << "Vprime\t" << v->dg[j][i+GetSequenceLength()] << endl;
    cout << "Q\t" << w5[GetSequenceLength()] << endl;
    cout << "Vprime/Q\t" << DIV(v->f(j,i+GetSequenceLength()),PROD(w5[GetSequenceLength()], PFSCALE(1, pfdata->scaling, 2))) << endl;
    
    return (double) DIV(v->f(j,i+GetSequenceLength()),PROD(w5[GetSequenceLength()], PFSCALE(1, pfdata->scaling, 2)));
}

double RNA::GetW(const int i, const int j){
    return (double) w->f(i,j);
}

//remove all pairs in structure # structurenumber.
//Also, roll back the number of specified structures if this is the last specified structure.
int RNA::RemovePairs(const int structurenumber, bool removeIfLastStructure) {
    
    //do some error checking
    if (structurenumber>ct->GetNumberofStructures()||structurenumber<1) return 3; //Structure number out of range


    //decrement the number of structures, if appropriate, i.e. this is the last structure
    if (removeIfLastStructure && structurenumber==ct->GetNumberofStructures()) {
        ct->RemoveLastStructure();
        return 0;
    }

    //otherwise, clean the selected structure of pairs:
    ct->CleanStructure(structurenumber);
    return 0;

}

#ifndef PARTITION_CORE //exclude from minimal partition-only builds (e.g. partition-rosetta)
//Calculate and return the folding free energy change for structure number structurenumber.
double RNA::CalculateFreeEnergy(const int structurenumber, const bool UseSimpleMBLoopRules) {
    //Do some simple error checking
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 0.0;

    if (!VerifyThermodynamic()) {
        ErrorCode = 5;//record an error
        return 0.0;//return 0.0 if a problem occurs
    }
    efn2(data,ct,structurenumber,UseSimpleMBLoopRules);

    //conversion factor is set in defines.h.  Free energies are multiplied by this factor internally so that integer math can be used.
    return (((double) ct->GetEnergy(structurenumber)/conversionfactor));
}

double RNA::ExteriorLoopCorrection(const int structurenumber, const bool UseSimpleMBLoopRules, int min_index, int max_index) {
    //Do some simple error checking
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 0.0;

    if (!VerifyThermodynamic()) {
        ErrorCode = 5;//record an error
        return 0.0;//return 0.0 if a problem occurs
    }
    //conversion factor is set in defines.h.  Free energies are multiplied by this factor internally so that integer math can be used.
    return (((double) ergexteriordiff(data, ct, structurenumber, UseSimpleMBLoopRules, min_index, max_index)));
}

//Write the details on the energy caclulation for all structures.
int RNA::WriteThermodynamicDetails(const char filename[], const bool UseSimpleMBLoopRules) {

    if (!VerifyThermodynamic()) return 5; //return non-zero if a problem occurs
    efn2(data,ct,0,UseSimpleMBLoopRules,filename);
    return 0;

}


#ifndef DYNALIGN_II
//Predict the secondary structure by free energy minimization.
//Also generate subooptimal solutions using a heuristic.
int RNA::FoldSingleStrand(const float percent, const int maximumstructures, const int window, const char savefile[], const int maxinternalloopsize, bool mfeonly, bool simple_iloops, bool disablecoax, bool allow_isolated) {
    char *savefilename;
    int percenti;
    int tracebackstatus;

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength() == 0) return 20;

    if (!VerifyThermodynamic()) return 5;

    //savefile will be passed to the function dynamic for structure prediction.
    //Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
    if (savefile == NULL || !strcmp(savefile, "")) savefilename = NULL;
    else 
    {
        savefilename = new char[((int) strlen(savefile)) + 1];
        strcpy(savefilename,savefile);
    }

    //right now, dynamic requires an integer specification of percent energy change.
    //FoldSingleStrand takes this is a float to provide the opportunity to reform this in the future.
    //For now, cast percent as an integer to pass to dynamic.
    percenti = (int) percent;

    //Predict the secondary structures.
    tracebackstatus = dynamic(ct, data, maximumstructures, percenti, window, progress, false, savefilename, maxinternalloopsize, mfeonly, simple_iloops, disablecoax,allow_isolated);

    //Clean up the memory use.
    if (savefilename != NULL) delete[] savefilename;
    if (progress && progress->canceled()) return 99; // This indicates that the folding operation was canceled by the user.

    if(tracebackstatus != 0) return 14;//This indicates a traceback error.
    else return 0;
}
#else
#endif

// Predict the lowest free energy secondary structure and generate all suboptimal structures.
int RNA::GenerateAllSuboptimalStructures(const float percent, const double deltaG) {

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    if (!VerifyThermodynamic()) return 5; //The thermodynamic data tables have not been read 

    //Call the alltrace function to do the work:
    alltrace(ct,data, ((short) percent), ((short) (deltaG*conversionfactor)),progress,NULL);

    return 0;


}

// Predict the structure with maximum expected accuracy and suboptimal structures.
int RNA::MaximizeExpectedAccuracy(const double maxPercent, const int maxStructures, const int window, const double gamma) {

    //first trap some possible errors
    if (!partitionfunctionallocated) {
        //There is no partition function data available.
        return 15;
    }


    //Past error trapping
    MaxExpectFill(ct, v, w5, pfdata, lfce, mod, fce, maxPercent, maxStructures, window, gamma, progress);

    return (progress&&progress->canceled()) ? 99 : 0;//no error return functionality right now

}
#endif // PARTITION_CORE

// This function predicts structures composed of probable base pairs.
int RNA::PredictProbablePairs(const float probability) {
    int i,j,count;
    char thresh[8];
    string label;//A string for making ct file labels

    //first trap some possible errors
    if (probability > epsilon && probability < 0.500-epsilon) {
        //The threshold is too low to be valie and not low enough that it will be considered zero, a the default
        return 18;
    }

    if (!partitionfunctionallocated) {
        //There is no partition function data available.
        return 15;
    }


    //Past error trapping

    
    
    

    
    if (probability>epsilon) {
        //The user specified a threshold, so use that and generate one structure
        

        //Get one clean structure
        if (ct->GetNumberofStructures()>0) {
            ct->CleanStructure(1);
            for (i=ct->GetNumberofStructures();i>1;--i) {
                ct->RemoveLastStructure();
            }
        }
        else {

            ct->AddStructure();
        }



        for (i=1;i<ct->GetSequenceLength();i++) {
            for (j=i+1;j<=ct->GetSequenceLength();j++) {

                if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce) > probability) {
                    //This pair exceeded the threshold, so add it to the list
                    ct->SetPair(i,j);
                    

                }

            }
        }

        //put the threshold in the ctlabel, which will appear in the header of a ct file
        sprintf(thresh,"%1.5f",probability);
        
        
        //Insert a label at the beging of ctlabel to provide the threshold.  This will be written in an output ct file.
        
        label = " >";
        label+=thresh;
        label+=" pairing probability; ";
        label+=ct->GetCtLabel(1);

        ct->SetCtLabel(label,1);
        
    }
    else {
        //The default threshold was specified, so create 8 structures, with thresholds of >=0.99, >=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.

        //Set up 8 clean structures:
        
        //Get 8 clean structures
        if (ct->GetNumberofStructures()>8) {
            for (i=ct->GetNumberofStructures();i>8;--i) {
                ct->RemoveLastStructure();
            }
            for (i=1;i<=8;++i) ct->CleanStructure(i);
        }
        else {
            for (i=1;i<=ct->GetNumberofStructures();++i) ct->CleanStructure(i);
            for (i=ct->GetNumberofStructures();i<8;++i) ct->AddStructure();
        }
        

        //loop over the structures and thresholds
        for (count=1;count<=8;count++) {

            for (i=1;i<ct->GetSequenceLength();i++) {
                for (j=i+1;j<=ct->GetSequenceLength();j++) {

                    if (count==1) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.99) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                            
                        }
                    }
                    else if (count==2) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.97) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==3) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.95) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==4) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.90) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==5) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.80) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==6) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.70) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==7) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.60) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                    else if (count==8) {
                        if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>.50) {
                            
                            //set this pair because it meets the threshold
                            ct->SetPair(i,j,count);
                        }
                    }
                }
            }

        }

        //add labels that would appear in a ct file header
        
        label = " >=97% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,2);
        
        label = " >=95% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,3);

        label = " >=90% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,4);

        label = " >=80% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,5);

        label = " >=70% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,6);

        label = " >=60% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,7);

        label = " >50% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,8);

        label = " >=99% probable pairs ";
        label+= ct->GetCtLabel(1);
        ct->SetCtLabel(label,1);

    }

    return 0;

}
#ifndef PARTITION_CORE //exclude from minimal partition-only builds (e.g. partition-rosetta)
int RNA::Rsample(const vector<double> &experimentalRestraints,  RsampleData &refdata, const int randomSeed, const char savefile[], const double cparam, const double offset, const int numsamples) {
        
const vector<double> &exp_rest = experimentalRestraints;
//vector<double> &ref_uu = refdata.react_uu;
//vector<double> &ref_pm = refdata.react_pmid;
//vector<double> &ref_pe = refdata.react_pend;

int errorCode, i, j;
rand64 rnd(randomSeed); // rand64 is a RNG defined in ../src/random.h  rnd(vector v) will return a random element from v

PartialProgress stepProgress(progress); // this is used so that partitionfunction will move the progress from 0 to 40% instead of 0 to 100% because there are subsequent steps.
ProgressHandler *oldProgress = progress; progress = &stepProgress; // store the current ProgressHandler and use stepProgress as a proxy for calls to PartitionFunction.
stepProgress.setNextStep(40); // specify that 40% of the work is to occur in the next step (PartitionFunction).

// Other restraints (e.g. added by ReadSHAPE or ReadOffset) are not compatible with Rsample.
// (Or maybe they are compatible...If so, remove this guard. This is here to prevent accidental loading of SHAPE data before calling Rsample)
if (ct->shaped) return 215;

// do 1st partition function
errorCode = PartitionFunction();
if (errorCode != 0)  return errorCode;

// indicate that the previous work is done and the next step (Stochastic) comprises another 20% of the work.
stepProgress.advanceToNextStep(20);

progress = NULL; // disable progress output during stochastic (otherwise the output is jibberish in SMP mode)

// do stochastic sampling
errorCode = Stochastic(numsamples, rnd.nextInt());
if (errorCode != 0) return errorCode;

progress = &stepProgress; // restore progress handler after stochastic
// assign shape reactivity to each nucleotide based on wether it is paired/unpaired
// by sampling from paired-mid/paired-end/unpaired shape reactivity distributions
  int seq_length;
  seq_length = GetSequenceLength();
  vector<double> shape_predicted(seq_length+1,0.0); //starts from 1 and element [0] is not used

  for ( i=1; i<= numsamples; i++) {
    // 1st nt
    if ( GetPair(1,i) ==0 ) { //u
      shape_predicted[1]+=rnd(refdata.react_uu);
    }
    else { //p-end
      shape_predicted[1]+=rnd(refdata.react_pend);
    }
    // nt's 2 to seq_length-1
    for(int j=2; j<=seq_length-1;j++) {
      if (GetPair(j,i) == 0 ) { //u
        shape_predicted[j]+=rnd(refdata.react_uu);
      }
      else if ( GetPair(j-1,i) != 0 && GetPair(j+1,i) != 0 ) { //p-mid
        shape_predicted[j]+=rnd(refdata.react_pmid);
      }
      else { //p-end
        shape_predicted[j]+=rnd(refdata.react_pend);
      }
    }
    //last nt
    if (  GetPair(seq_length,i) == 0 ) { //u
      shape_predicted[seq_length]+=rnd(refdata.react_uu);
    }
    else { //p-end
      shape_predicted[seq_length]+=rnd(refdata.react_pend);
    }
  }
   for (i=1;i<=seq_length;i++) {
    shape_predicted[i]=shape_predicted[i]/numsamples;
  }
// convert all experimental reactivity > -500 and less than -react_offset to -react_offset+0.01
// also convert original array to be [1] based
vector<double> ExpReactNorm(seq_length+1,0.0);
for (i=1;i<=seq_length;i++) {
    if (exp_rest[i-1] < -1.0*offset && exp_rest[i-1] > -500.0) ExpReactNorm[i] = -1.0*offset +0.01;
    else
       ExpReactNorm[i] = exp_rest[i-1];
  }

// find which nt have measured reactivity i.e are larger than -500.0
vector<int> shape_nucleotide(1, 0);
for (i=1;i<=seq_length;i++) {
    if (exp_rest[i-1] > -500.0) shape_nucleotide.push_back(i);
}

// convert all calculated reactivities that are < -react offset to -react_offset+0.01 
for (i=1;i<=seq_length;i++) {
    if (shape_predicted[i] < -1.0*offset )  shape_predicted[i] = -1.0*offset +0.01;
  }

  // calculate log_term
  vector<double> log_term(seq_length+1,0.0);
  for (i=1;i<shape_nucleotide.size();i++) {
    log_term[shape_nucleotide[i]] = log( (ExpReactNorm[shape_nucleotide[i]] + offset) / (shape_predicted[shape_nucleotide[i]]+offset) );
  }

  // make p which is pseudo-dg
  vector<double> p( seq_length+1, 0.0); 
 
  for (j=1;j<shape_nucleotide.size();j++) {
    p[shape_nucleotide[j]] = cparam * log_term[shape_nucleotide[j]];
  }
  // initialize SHAPE array and set all of its elements to 0
  ct->AllocateSHAPE();

  // assign elements of SHAPE array to be p
  for (j=1;j<=seq_length;j++) {
    ct->SHAPE[j] = p[j]*conversionfactor;
    ct->SHAPE[seq_length + j] = p[j]*conversionfactor;  
  }

    // Remove sampled structures
    ct->RemoveAllStructures();

    // indicate that the previous work is done and the next step (PartitionFunction) comprises another 40% of the work.
    stepProgress.advanceToNextStep(40);

   // call partition function again, now with pseudo-dg calculated using rsample algorithm
	errorCode = PartitionFunction(savefile);
    if (errorCode !=0 ) return errorCode;

    stepProgress.stepComplete(); // indicate that the last step is done (to put progress to 100%)
    progress = oldProgress; // restore previous progress handler.
    // remove the SHAPE data that was added above (and converted to equillibrium constants by PartitionFunction)
    ct->DeleteSHAPE();

    return 0;
}
#endif //PARTITION_CORE

//Calculate the partition function for the current sequence.
int RNA::PartitionFunction(const char savefile[], double temperature, bool disablecoax, bool restoreSHAPE, bool allowisolated, int maxinter) {
    int i,j;
    char *savefilename;
    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;
    
    if (!VerifyThermodynamic()) return 5; //The thermodynamic data tables have not been read 


    //savefile will be passed to the function dynamic for structure prediction.
    //Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
    if (is_blank(savefile)) 
	    savefilename=NULL;
    else {
        savefilename=new char[((int) strlen(savefile))+1];
        strcpy(savefilename,savefile);
    }

    if (partitionfunctionallocated) {
        delete v;
        delete w;
        delete wmb;
        delete wl;
        delete wlc;
        delete wmbl;
        delete wcoax;
        delete fce;
        delete[] lfce;
        delete[] mod;
        delete[] w3;
        delete[] w5;
        delete pfdata;
    }
    //Allocate the memory needed (only if this is the first call to pfunction):
    //indicate that the memory has been allocated so that the destructor will delete it.
    partitionfunctionallocated = true;

    //allocate space for the v and w arrays:
    w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    fce = new forceclass(ct->GetSequenceLength());
    lfce = new bool [2*ct->GetSequenceLength()+1];
    mod = new bool [2*ct->GetSequenceLength()+1];

    for (i=0;i<=2*ct->GetSequenceLength();i++) {
        lfce[i] = false;
        mod[i] = false;
    }

    for (i=0;i<ct->GetNumberofModified();i++) {

        if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
            mod[ct->GetModified(i)]=true;
            mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
        }
    }

    w5 = new PFPRECISION [ct->GetSequenceLength()+1];
    w3 = new PFPRECISION [ct->GetSequenceLength()+2];

    if (ct->intermolecular) {
        //take advantage of templating to prevent intramolecular base pairs

        ct->allocatetem();//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
        for (i=1;i<ct->inter[0];i++) {
            for (j=i+1;j<=ct->inter[2];j++) {

                //Set intermolecular pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
                ct->tem[j][i]=false;

            }
        }
        for (i=ct->inter[2]+1;i<ct->GetSequenceLength();i++) {
            for (j=i+1;j<=ct->GetSequenceLength();j++) {

                //Another set of intermolecular pairs to forbid.
                ct->tem[j][i]=false;

            }
        }
    }

    //Initialize the partition function datatable:
        //Ignore the setting of parameter temperature if it is less than zero.
        //Generally, this parameter should be left at the default.
    
    pfdata=new pfdatatable(data,TO_PFP(scalingdefinition),(temperature<0?GetTemperature():temperature));

    //This code converts the SHAPE array of data to equilibrium constants, which is
    //required for the partition function. If the restoreSHAPE parameter is true 
    // (which it is by default), then a backup copy of the SHAPE data is made and is
    // restored after partition, so that it can be used for Fold, etc.
    double *SHAPE_backup = NULL;
    if (ct->shaped) {
      if (restoreSHAPE) SHAPE_backup = ct->CopySHAPE(false/*includeSHAPEss*/);
      for (i=1;i<=2*ct->GetSequenceLength();i++) ct->SHAPE[i]=(double) boltzman( ct->SHAPE[i], pfdata->pftemp);
    }

    // Converts Experimental Pair Bonus array to equilibrium constants
    if (ct->experimentalPairBonusExists){
      for (i = 1; i <= 2*ct->GetSequenceLength(); i++){
        for(j = i; j <= 2*ct->GetSequenceLength(); j++){
          double avg = 0.5*(ct->EX[i][j] + ct->EX[j][i]);
          ct->EX[i][j] = boltzman(avg, pfdata->pftemp);
          ct->EX[j][i] = ct->EX[i][j]; //symmetrize as well...
        }
      }
    }

    //The next section handles the case where base pairs are not
            //not allowed to form between nucs more distant
            //than ct->GetPairingDistanceLimit()
    if (ct->DistanceLimited()) {
        //This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
        if (!ct->templated) ct->allocatetem();
        for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
            for (i=1;i<j;i++) {
                //Set distant pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
                if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
            }
        }
    }

#ifndef _CUDA_CALC_
    //cout << "Performing CPU Partition Function Code" << endl;
    //default behavior: calculate the partition function on the CPU
    calculatepfunction(ct,pfdata,progress,savefilename,false,&Q,w,v,wmb,wl,wlc,wmbl,wcoax,fce,w5,w3,mod,lfce,disablecoax,allowisolated, maxinter);
#else //ifdef _CUDA_CALC_
    //if cuda flag is set, calculate on GPU
    //this requires compilation with nvcc
    //it will work for pair probabilities but it will break stochastic traceback

//read nearest neighbor parameters at desired temperature
    const char *path = getDataPath();
    if (!path)
        die("%s: need to set environment variable $DATAPATH", "Could not find nearest neighbor parameters");
    else
        cout << endl << "Read in data tables" << endl;
    struct param par;
    //param *par = new param();
    //cout << "Size (cpp) : " << sizeof(*par) << endl;
    //cout << "Size (cpp int22) : " << sizeof(par->int22) << endl;
    

    param_read_from_text(path, &par, 0,0);
//    param_read_alternative_temperature(path, &par, pfdata->pftemp, !isrna);
//copy the sequence, correcting for 1-indexing
    
    char* buf = new char[ct->GetSequenceLength()+1];
    for(int i=0;i<GetSequenceLength();i++){
        buf[i] = ct->nucs[i+1];
    }
    buf[ct->GetSequenceLength()] = '\0'; //make it null-terminated
//calculate the partition function on the GPU
    int *bcp;
    bcp = ct->generate_constraint_matrix();

    prna_t p = prna_new(buf, &par, 1, bcp);
//transfer over the partition function from the GPU data structure
//because the GPU partition function is calculated in log space, we do not directly transfer the partition function
//instead, we put the pair probs, calculated in log space in V, 1 everywhere else
//so RNA::GetPairProbability(i,j) will return V(i,j)*1/1, which will be the pairing probability

	//define and build the inc array, to tracka allowed pairs:
	bool **inc = new bool* [data->alphabet.size()];
	for (int letters = 0; letters < data->alphabet.size(); ++letters) {
		inc[letters] = new bool[data->alphabet.size()];
		for (int letterpairs = 0; letterpairs < data->alphabet.size(); ++letterpairs) {
			inc[letters][letterpairs] = data->pairing[letters][letterpairs];

		}
	}

    for(int i=1;i<=ct->GetSequenceLength();i++){
        for(int j=i+1;j<=ct->GetSequenceLength();j++){
            //cout << "i:\t" << i << "\tj:\t" << j <<endl;
            // Get V(i,j)

			//make sure the nucleotides can pair, otherwise another conversion is needed:
			bool calculatev = true;
			if (!inc[ct->numseq[i]][ct->numseq[j]]) {
				//These are two nucleotides that cannot form a canonical pair
				calculatev = false;
			}
			if (calculatev) {
				int localbefore = 0;
				int localafter = 0;
				if ((i > 1 && j < ct->GetSequenceLength())) {

					localbefore = inc[ct->numseq[i - 1]][ct->numseq[j + 1]];

				}
                else localbefore = 0;

				//after = 0 if a stacked pair cannot form 3' to locali
				if ((j - i) > minloop + 2) {
					localafter = inc[ct->numseq[i + 1]][ct->numseq[j - 1]];

				}
				else localafter = 0;

				//if there are no stackable pairs to locali.localj then don't allow a pair locali,localj
				if ((localbefore == 0) || (localafter == 0)) {
					//v->f(locali,localj)= 0;
					calculatev = false;
				}

			}

			if (calculatev) {
				v->f(i, j) = XLOG_TO_PF(-1 * get_v_array(p, i - 1, j - 1)); //V(i,j) = P(i,j)

				//Get V'(i,j)
				//In the partition-cuda code, V'(i,j) is stored in V(j,i)
				v->f(j, i + ct->GetSequenceLength()) = XLOG_TO_PF(-1 * get_v_array(p, j - 1, i - 1)); //V'(i,j) = 1
			}
			else {
				v->f(i, j) = ZERO;
				v->f(j, i + ct->GetSequenceLength()) = ZERO;
			}
        }
        //cout << "i:\t" << i << endl;
        w5[i] = PFSCALE(XLOG_TO_PF(-1*get_w5_array(p, i-1)), pfdata->scaling, -2); //have to divide by scaling*2
        w3[i] = PFSCALE(XLOG_TO_PF(-1*get_w3_array(p, i-1)), pfdata->scaling, -2); //because RNA::calculateprobability expects scaling values
    }

	//clean up the inc array, which tracked what can pair
	for (int letters = 0; letters < data->alphabet.size(); ++letters) delete[] inc[letters];
	delete[] inc;

    delete[] buf;
    prna_delete(p);
    delete[] bcp;
    //delete &par;
#endif
	if (savefilename!=NULL) {
        if (!progress||!progress->canceled()) writepfsave(savefilename,ct,w5,w3,v,w,wmb,wl,wlc,wmbl,wcoax,fce,mod,lfce,pfdata);

        //clean up some memory use:
        delete[] savefilename;
    }
    if (SHAPE_backup!=NULL) { ct->LoadSHAPE(SHAPE_backup, false/*includeSHAPEss*/); delete[] SHAPE_backup; }
    return progress&&progress->canceled() ? 99 : 0;
}

#ifndef PARTITION_CORE //exclude from minimal partition-only builds (e.g. partition-rosetta)
//Predict maximum expected accuracy structures that contain pseudoknots from either a sequence or a partition function save file.
int RNA::ProbKnot(int iterations, int MinHelixLength, double threshold) {

    //first trap some possible errors
    if (!partitionfunctionallocated) {
        //There is no partition function data available.
        return 15;
    }

    if (iterations < 1) {
        //there can't be fewer than one iteration
        return 24;

    }

	if (threshold < 0) {
		//There can't be pairs with less than zero probablity
		return 36; 
	}

    //Past error trapping
    //Call the ProbKnot Program:
    return ProbKnotAssemble(v, w5, ct, pfdata, lfce, mod, pfdata->scaling, fce, iterations, MinHelixLength, threshold );


}

//Predict maximum expected accuracy structures that contain pseudoknots from a file containing ensemble of structures.
int RNA::ProbKnotFromSample(int iterations, int MinHelixLength, double threshold) {

    if (iterations < 1) {
        //there can't be fewer than one iteration
        return 24;

    }

	if (threshold < 0) {
		//There can't be pairs with less than zero probablity
		return 36; 
	}

    //Past error trapping
    //Call the ProbKnot Program:
    return ProbKnotAssemble( ct, iterations, MinHelixLength, threshold );

}


//Refold a sequence using data from a save file.
int RNA::ReFoldSingleStrand(const float percent, const int maximumstructures, const int window) {

    if (!energyallocated) {
        //A .sav file was not read by the constructor.  Therefore, this function cannot be performed.
        return 17;

    }

    //Now do the refolding.
    return traceback(ct, data, ev, ew, ewmb, ew2, ewmb2, ew3, ew5, fce, lfce, vmin, maximumstructures, (int) percent, window,mod);

}


//Sample structures from the Boltzman ensemble.
int RNA::Stochastic(const int structures, const int seed) {

    if (!partitionfunctionallocated) {
        //There is no partition function data available.
        return 15;
    }

	// Remove existing structures (in case stochastic or another operation was 
	// previously performed).
	ct->RemoveAllStructures();

    //Past error trapping, call the stochastic traceback function
    return stochastictraceback(w,wmb,wmbl,wcoax,wl,wlc,v,
        fce, w3,w5,pfdata->scaling, lfce, mod, pfdata, structures,
        ct, seed, progress);


}
#endif

//Force a nucleotide to be double stranded (base paired).
//Return an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceDoubleStranded(const int i) {
    int index;

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //Check that nucleotide is valid
    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

    //Check for conflicting restraints (anything specifying i to be unpaired):
    for (index=0;index<ct->GetNumberofSingles();index++) {

        if(ct->GetSingle(index)==i) return 9;
    }

    
    ct->AddDouble(i);
    
    return 0;

    


}

//Function to specify a nucleotide, i, that is a U in a GU pair
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint, 11 = nucleotide not U).
int RNA::ForceFMNCleavage(const int i) {
    int index;

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //Check that nucleotide is valid
    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

    //Check to make sure the nucleotide is U.
    if (ct->numseq[i]!=4) return 11;

    //Check for conflicting restraints (anything specifying i to be unpaired):
    for (index=0;index<ct->GetNumberofSingles();index++) {

        if(ct->GetSingle(index)==i) return 9;
    }

    //Check for a nucleotide already forced to be in a pair that is not a GU pair.
    for (index=0;index<ct->GetNumberofPairs();index++) {

        if (i==ct->GetPair5(index)&&ct->numseq[ct->GetPair3(index)]!=3) return 9;
        else if (i==ct->GetPair3(index)&&ct->numseq[ct->GetPair5(index)]!=3) return 9;

    }


    ct->AddGUPair(i);
    
    return 0;

    


}

//Specify the maximum distance allowed between paired nucleotides in subsequent structure prediction.
//return An integer that indicates an error code (0 = no error, 12 = too long or too short).
int RNA::ForceMaximumPairingDistance(const int distance) {

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    if (distance < minloop+1) return 12;
    else {
        
        ct->SetPairingDistance(distance);
        return 0;
    }


}


//Indicate a nucleotide that is accessible to chemical modification.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified).
int RNA::ForceModification(const int i) {

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //Check that nucleotide is valid
    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range


    


    
    //Go ahead and record the constraint.
    ct->AddModified(i);
    
    return 0;
    


}

//Force a base pair between nucleotides i and j.
//Returns an error code: (0 = no error, 4 = nucleotide out of range, 6 = pseudoknot formation, 7 = non-canonical pair, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForcePair(const int i, const int j) {
    bool allowedpairs[6][6]={{false,false,false,false,false,false},{false,false,false,false,true,false},{false,false,false,true,false,false},{false,false,true,false,true,false},
    {false,true,false,true,false,false},{false,false,false,false,false,false}};
    int index;
    int locali,localj;

    //First perform the error checking:

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //Note: In structure, forced pairs run between index of 1 and a maximum of maxforce-1.

    //if (ct->npair==(maxforce-1)) return 8;//This means there are too many pair constraints.

    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range
    if (j<1||j>ct->GetSequenceLength()) return 4;//j is out of range

    if (!allowedpairs[ct->numseq[i]][ct->numseq[j]]) return 7;//non-canonical pair

    //sort indexes from 5' to 3':
    if (i>j) {
        locali=j;
        localj=i;
    }
    else {
        locali=i;
        localj=j;
    }

    //check for pseudoknots with any other forced pair or the same nucleotide forced into two pairs:
    for (index=0;index<ct->GetNumberofPairs();index++) {
        if (locali<ct->GetPair5(index)&&ct->GetPair5(index)<localj&&localj<ct->GetPair3(index)) return 6;//a pseudoknot

        if (locali==ct->GetPair5(index)||locali==ct->GetPair3(index)||localj==ct->GetPair5(index)||localj==ct->GetPair3(index)) return 9;//i or j is in another forced pair

    }

    //now check for other conflicting restraints:
    for (index=0;index<ct->GetNumberofForbiddenPairs();index++) {
        if(ct->GetForbiddenPair5(index)==locali && ct->GetForbiddenPair3(index)==localj ) return 9;//The pair was forbidden.

    }
    for (index=0;index<ct->GetNumberofSingles();index++) {

        if(ct->GetSingle(index)==locali||ct->GetSingle(index)==localj)  return 9;//i or j was previously forced single-stranded.
    }

    //Now register the restraint because the error checking was clear or errors.
    ct->AddPair(locali,localj);

    return 0;

}

//Prohibit a pair between two nucleotides in subsequent structure prediction.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = nucleotide in conflicting restraint).
int RNA::ForceProhibitPair(const int i, const int j) {
    int index,locali,localj;

    //First perform the error checking:

    //Note: In structure, forced pairs run between index of 0 and a maximum of maxforce-1.

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;


    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range
    if (j<1||j>ct->GetSequenceLength()) return 4;//j is out of range


    //sort indexes from 5' to 3':
    if (i>j) {
        locali = j;
        localj = i;
    }
    else {
        locali=i;
        localj=j;
    }

    //check to make sure this pair hasn't been forced:
    for (index=0;index<ct->GetNumberofPairs();index++) {

        if (locali==ct->GetPair5(index)&&localj==ct->GetPair3(index)) return 9;//i or j is in a forced pair

    }


    //Now register the restraint because the error checking was clear or errors.
    ct->AddForbiddenPair(locali,localj);

    return 0;

}

//Force a nucleotide to be single stranded in subsequent structure prediction.
//An integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceSingleStranded(const int i) {
    int index;

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //Check that nucleotide is valid
    if (i<1||i>ct->GetSequenceLength()) return 4;//i is out of range

    //Check for conflicting constraints; anything forcing a nucleotide to be paired.

    for (index=0;index<ct->GetNumberofPairs();index++) {//check all the forced pairs
        if (i==ct->GetPair5(index)||i==ct->GetPair3(index)) return 9;//i is in a forced pair
    }
    for (index=0;index<ct->GetNumberofDoubles();index++) {//check all the force doubles
        if(ct->GetDouble(index)==i) return 9;
    }
    for (index=0;index<ct->GetNumberofGU();index++) {//check all the force FMN
        if(ct->GetGUpair(index)==i) return 9;
    }


    //Register the constraint:
    ct->AddSingle(i);

    return 0;

}

//Return a nucleotide that is forced double stranded.
int RNA::GetForcedDoubleStranded(const int constraintnumber) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofDoubles()) return 0;

    //Now return the constraint.
    return ct->GetDouble(constraintnumber);

}

//Return a nucleotide that is accessible to FMN cleavage.
int RNA::GetForcedFMNCleavage(const int constraintnumber) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofGU()) return 0;

    //Now return the constraint.
    return ct->GetGUpair(constraintnumber);

}

//Return a nucleotide that is accessible to modification.
int RNA::GetForcedModification(const int constraintnumber) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofModified()) return 0;

    //Now return the constraint.
    return ct->GetModified(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.

}

//Return a nucleotide in a forced pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedPair(const int constraintnumber, const bool fiveprime) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofPairs()) return 0;

    //Now return the constraint.
    if (fiveprime) return ct->GetPair5(constraintnumber);
    else return ct->GetPair3(constraintnumber);

}

//Return a nucleotide in a prohibited pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedProhibitedPair(const int constraintnumber, const bool fiveprime) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofForbiddenPairs()) return 0;

    //Now return the constraint.
    if (fiveprime) return ct->GetForbiddenPair5(constraintnumber);
    else return ct->GetForbiddenPair3(constraintnumber);

}

//Return a nucleotide that is forced single stranded.
int RNA::GetForcedSingleStranded(const int constraintnumber) {

    //First make sure the constraintnumber is valid.
    if (constraintnumber<0||constraintnumber>=ct->GetNumberofSingles()) return 0;

    //Now return the constraint.
    return ct->GetSingle(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.


}


//Return the maximum pairing distance.
//Return an integer that indicates the maximum distance allowed between paired nucleotides, where -1 indicates that the maximum distance is not set.
int RNA::GetMaximumPairingDistance() {

    if (ct->DistanceLimited()) return ct->GetPairingDistanceLimit();
    else return -1;

}


//Return the number of nucletides forced to be paired.
int RNA::GetNumberOfForcedDoubleStranded() {

    return ct->GetNumberofDoubles();

}

// Add an experimental bonus to a pair of nucleotides
//void RNA::SetExperimentalBonus(const int i, const int j, const double bonus){

//  ct->EX[i][j] = bonus;

//}


//!Return the number of nucleotides accessible to FMN cleavage.
int RNA::GetNumberOfForcedFMNCleavages() {

    return ct->GetNumberofGU();

}

//!Return the number of nucleotides accessible to chemical modification.
int RNA::GetNumberOfForcedModifications() {

    return ct->GetNumberofModified();

}

//!Return the number of forced base pairs.
int RNA::GetNumberOfForcedPairs() {

    return ct->GetNumberofPairs();

}

//!Return the number of prohibited base pairs.
int RNA::GetNumberOfForcedProhibitedPairs() {

    return ct->GetNumberofForbiddenPairs();

}

//!Return the number of nucleotides that are not allowed to pair.
int RNA::GetNumberOfForcedSingleStranded() {

    return ct->GetNumberofSingles();

}

//Read a set of folding constraints to disk in a plain text file.
//filename is a c string that is the file name to be read.
//Returns an integer that indicates an error code (0 = no error, 1 = file not found, 13 = error reading constraint file).
int RNA::ReadConstraints(const char filename[]) {
    FILE *check;

    //check that the file exists.
    if ((check = fopen(filename, "r"))== NULL) {
        //the file is not found
        fclose(check);
        return 1;
    }
    fclose(check);

    //Now read the constraints
    if (readconstraints(filename, ct)) return 0;
    else return 13;
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
//slope is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
//intercept is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double slope, const double intercept, RestraintType modifier, const bool IsPseudoEnergy) {
    // ct->ReadSHAPE will verify that the SHAPE input file exists
    int code;
    if (IsPseudoEnergy) {
        //This is the pseudo energy version
        ct->SHAPEslope=slope*conversionfactor;//register the slope in tenths of kcal/mol
        ct->SHAPEintercept=intercept*conversionfactor;//register the intercept in tenths of a kcal/mol
        code = ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
    }
    else
        code = ct->ReadSHAPE(filename,(float) slope,(float) intercept);//call ReadSHAPE() with parameters to parse thresholds
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code;
}

//Read DMS data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains DMS data.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadDMS(const char filename[], const bool bynt) {
    
    int code;
    if (bynt) {
      code = ct->ReadSHAPE(filename, RESTRAINT_DMSNT); //call ReadSHAPE() to read the file and determine pseudo energies
    } else {
      code = ct->ReadSHAPE(filename, RESTRAINT_DMS); //call ReadSHAPE() to read the file and determine pseudo energies
    }
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code;    
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions - overloaded version for including single-stranded SHAPE.
//filename is a c string that indicates a file that contains SHAPE data.
//dsSlope is the double-stranded slope.
//dsIntercept is the double-stranded intercept.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
//ssSlope is the single-stranded slope.
//ssIntercept in the single-stranded intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double dsSlope, const double dsIntercept, const double ssSlope, const double ssIntercept, RestraintType modifier) {
    // ct->ReadSHAPE will verify that the SHAPE input file exists

    ct->SHAPEslope=dsSlope*conversionfactor;//register the slope in tenths of kcal/mol
    ct->SHAPEintercept=dsIntercept*conversionfactor;//register the intercept in tenths of a kcal/mol
    ct->SHAPEslope_ss=ssSlope*conversionfactor;//register the slope in tenths of kcal/mol
    ct->SHAPEintercept_ss=ssIntercept*conversionfactor;//register the intercept in tenths of a kcal/mol
    int code = ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the slope.
//parameter2 is the intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadExperimentalPairBonus(const char filename[], double const experimentalOffset, double const experimentalScaling ) {
    // ct->ReadExperimentalPairBonus will verify that the input file exists
    int code = ct->ReadExperimentalPairBonus(filename, experimentalOffset, experimentalScaling );
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}

//Read Double Strand Offset
int RNA::ReadDSO(const char filename[]) {
    // ct->ReadOffset will verify that the SHAPE input file exists
    int code = ct->ReadOffset(NULL,filename);
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}

//Read Single Strand Offset
int RNA::ReadSSO(const char filename[]) {
    // ct->ReadOffset will verify that the SHAPE input file exists
    int code = ct->ReadOffset(filename,NULL);
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}


//Remove all previously defined constraints.
void RNA::RemoveConstraints() {
    ct->RemoveConstraints();
    ct->min_gu=0;
    ct->min_g_or_u=0;
    ct->nneighbors=0;
    ct->nregion=0;
    ct->nmicroarray=0;
}

void RNA::SetConstraints(vector<int> ss) { //HERE TONY
   for (int i=0; i<ss.size(); i++) ct->AddSingle(ss[i]);
}


//Add extrinsic restraints for partition function calculations.
int RNA::SetExtrinsic(int i, int j, double k) {
    int locali, localj;

    //First do the error checking:
    //check the indices
    if (i<1||i>ct->GetSequenceLength()||j<1||j>ct->GetSequenceLength()) return 4;

    //make sure the equilibrium constant is not less than zero.
    if (k<0) return 26;

    //Now past error checking:

    //sort indexes from 5' to 3':
    if (i>j) {
        locali = j;
        localj = i;
    }
    else {
        locali=i;
        localj=j;
    }



    if (ct->constant==NULL) {
        //allocate the space needed in structure to store the constants
        ct->allocateconstant();
    }


    ct->constant[localj][locali] = TO_XLOG(k);

    return 0;
}

//Write the set of folding constraints to disk.
int RNA::WriteConstraints(const char filename[]) {

    outputconstraints(filename,ct);
    return 0;

}

//Specify a comment for inclusion in subsequently written .ct files.
int RNA::AddComment(const char comment[], const int structurenumber) {
    string label;

    //start with error checking:
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) return 3;

    //now register the comment (at the end of existing string, but before the newline there by default):
    label = ct->GetCtLabel(structurenumber);
    
    //Remove the existing newline character, if it exists:
    if (label.length()>0) if (label[label.length()-1]=='\n') label.erase(label.length()-1);

    //Add the comment
    label+=comment;

    //Add back the newline
    label+="\n";

    //add a newline at the end of the comment -- required when ct is written
    ct->SetCtLabel(label,structurenumber);

    return 0;

}

//Write a ct file of the structures
int RNA::WriteCt(const char filename[], bool append, CTCommentProvider &commentProvider) const {
    if (ct->GetNumberofStructures()>0)
        return ct->ctout(filename,append,commentProvider);
    return 10; //an error code
}

//Write a dot-bracket file of the structures.
int RNA::WriteDotBracket(const char filename[], const int structurenumber, 
						 const DotBracketFormat format, 
						 CTCommentProvider &commentProvider) const {
    if (ct->GetNumberofStructures()>0)
        return ct->writedotbracket(filename, structurenumber, format, commentProvider);
    return 10; //an error code
}

#ifndef PARTITION_CORE //exclude from minimal partition-only builds (e.g. partition-rosetta)
//Break any pseudoknots that might be in a structure.
int RNA::BreakPseudoknot(const bool minimum_energy, const int structurenumber, const bool useFastMethod) {
    if (useFastMethod && !minimum_energy) {
        // use the fast method of breaking pseudoknots.
        // It doesn't require loading any energy tables and is completely alphabet and energy agnostic.
        if (structurenumber > 0)
        ct->BreakPseudoknots(structurenumber);
        else
        for(int i = 1; i <= GetStructureNumber(); i++)
        ct->BreakPseudoknots(i);
    } else {
    int i,j,structures;
    structure *tempct;


    double **bpProbArray; //contains the raw bp probabilities for each bp
    double *bpSSProbArray; //contains the raw single strand probability for a base
    double **vwArray;  //contains v and w recursion values
    double **vwPArray; //the v' and w' recursion values

    double *w3Array=0;//w3Array[i] is the maximum score from nucletides i to ct->GetSequenceLength()
    double *w5Array=0;//w5Array[i] is the maximum score from nucleotides 1 to i
                            

    int Length;
    int start,stop;

    // WARNING sumPij is set to 0, but is never accessed.
    double sumPij; //holds the sum of probabilities for base pairs based on a specific i over js


    //Make sure there are structures...
    if (ct->GetNumberofStructures()<=0) {
        return 23;
    }

    //If a specific structure is requested, i.e. structurenumber!=0, make sure the structure is valid.
    if (structurenumber!=0) {
        if (structurenumber<0||structurenumber>ct->GetNumberofStructures()) return 3;

    }


    // verify that the datatable has been read (otherwise, read it now).
    // This is required no matter which branch is used below. (The datatable's 'pairing' member is used by MEAFill.)
    if (!VerifyThermodynamic()) return 5;


    //allocate tempct:
    tempct = new structure(2);
    tempct->allocate(ct->GetSequenceLength());
    tempct->SetSequenceLabel("temp\n");
    tempct->SetThermodynamicDataTable(data);
    if (minimum_energy) {
        //mnimum_energy==true, so minimize free energy in the pseudoknot free structure
        for (i=1;i<=ct->GetSequenceLength();i++) {
            tempct->numseq[i]=ct->numseq[i];
        }
        //allocate a template of allowed pairs
        tempct->allocatetem();

        //Break pairs for each structure (first check if only one structure should be treated)

        if (structurenumber!=0) {
            start = structurenumber;
            stop = structurenumber;

        }
        else {
            start = 1;
            stop = ct->GetNumberofStructures();
        }
        for (structures=start;structures<=stop;structures++) {

            //initialize all base pairs as unallowed:
            for (i=0;i<=ct->GetSequenceLength();i++) {
                for (j=i+1;j<=ct->GetSequenceLength();j++) {
                    tempct->tem[j][i] = false;
                }
           }

            //now allow the pairs that are in the loaded ct:
            for (i=1;i<=ct->GetSequenceLength();i++) {
                if (ct->GetPair(i,structures)>i) {
                    tempct->tem[ct->GetPair(i,structures)][i] = true;
                }
            }

            if (tempct->GetNumberofStructures() > 0) {
                //strip the structure from the ct at this point
                for (int index=tempct->GetNumberofStructures();index>0;--index) tempct->RemoveLastStructure();

            }

            //Predict the secondary structures.
            alltrace(tempct,data,0,0,progress,NULL,false);//
            //dynamic(tempct, data, 1, 0, 0, progress, false, NULL, 30);



            //copy the pairs back to ct
            ct->CleanStructure(structures);
            for (i=1;i<=ct->GetSequenceLength();i++) {
                if (tempct->GetPair(i)>i) ct->SetPair(i,tempct->GetPair(i),structures);

            }
            //copy the energy back to ct (this energy is probably inaccurate and should be updated by calling efn2(data, ct, structures, false); before saving.
            ct->SetEnergy(structures,tempct->GetEnergy(1));
        }
    }//end of minimum_energy == true
    else {
        //This is minimum_energy == false, so maximize pairs

        for (i=1;i<=ct->GetSequenceLength();i++) {
            tempct->numseq[i]=ct->numseq[i];

        }


        //loop over all structures unless the programmer requested a specific structure
        if (structurenumber!=0) {
            start = structurenumber;
            stop = structurenumber;

        }
        else {
            start = 1;
            stop = ct->GetNumberofStructures();
        }
        for (structures=start;structures<=stop;structures++) {



             
            //strip the structure from the ct at this point
            if (tempct->GetNumberofStructures()>0) {
                for(int index=tempct->GetNumberofStructures();index>0;index--) tempct->RemoveLastStructure();
            }



            //allocate main arrays and initialize the Arrays to defaults
            bpProbArray = new double *[ct->GetSequenceLength()+1];
            bpSSProbArray = new double [ct->GetSequenceLength()+1];
            vwArray = new double *[ct->GetSequenceLength()+1];
            vwPArray = new double *[ct->GetSequenceLength()+1];

            // WARNING: sumPij is never used!!! the following line is ineffective.
            sumPij = 0;

            for (i=0;i<=ct->GetSequenceLength();i++) {
                bpProbArray[i] = new double [ct->GetSequenceLength()+1];
                vwArray[i] = new double [ct->GetSequenceLength()+1];
                vwPArray[i] = new double [ct->GetSequenceLength()+1];

                //tempct->basepr[1][i] = 0;
                bpSSProbArray[i] = 0;

                for (j=0;j<=ct->GetSequenceLength();j++) {
                    bpProbArray[i][j]=0;
                    vwArray[i][j]=-0;
                    vwPArray[i][j]=-0;
                }
            }




            tempct->nucs[0] = ' ';

            // Recursion rules investigate for vwArray
            //    1)  if the base pair (BP) probability value is 0, skip that pair
            //    2)  hairpin turns at 5 BP
            //    3)  stack/internal/bulge pairing at 7 BPs
            //    4)  multibranching at 12 BPs (2 hairpins and a stack)
            // Because of the rules for hairpin, the probabilities
            //    can be taken from the bpProbArray directly until j-i > 5

            // Calculate the single stranded probabilities for each base
            // Pi = 1 - (for all j, sum(Pij)
            // fill in w for the diagonal for the Pi,i
            for (i=1; i<=ct->GetSequenceLength(); i++)
            {


                if (ct->GetPair(i,structures)!=0) bpSSProbArray[i] = 0.0;
                else bpSSProbArray[i] = 1.0;



                vwArray[i][i] = bpSSProbArray[i];
            } // end loop over each base pair


            //Calculate the base pair probabilities to start...
            for (Length = 2; Length <=ct->GetSequenceLength(); Length++)
            {

                //begin populating v and w along the diagonal starting with the
                //   shortest hairpin length
                for (i=1, j=i+Length-1; j<=ct->GetSequenceLength(); i++, j++)
                {
                    if (ct->GetPair(i,structures)==j) {
                        bpProbArray[j][i]=1.0;
                    }
                    else {
                        bpProbArray[j][i]=-1.0;
                    }

                }
            }


            //Call the MEAFill routine.
            //Note the false at the end "allows" non-canonical pairs.  This is required so that
            //non-canonical pairs aren't spuriosly broken
            MEAFill(tempct, bpProbArray, bpSSProbArray, vwArray, vwPArray, w5Array, w3Array, &data->pairing, 1.0, 0, progress, false);

            // start traceback
            trace(tempct, vwArray, vwPArray, bpProbArray, 1.0, 0, 1, 0);






            // Deallocate memory for the MaxExpect calculation
            //Arrays with functionality in the fill step
            for (i=0;i<=ct->GetSequenceLength();i++) {
                delete[] bpProbArray[i];
            }
            delete[] bpProbArray;

            delete[] bpSSProbArray;

            for (i=0; i<=ct->GetSequenceLength(); i++) {
                delete[] vwArray[i];
                delete[] vwPArray[i];
            }
            delete[] vwArray;
            delete[] vwPArray;




            //copy the pairs back to ct

            //clean the structure of pairs first
            ct->CleanStructure(structures);

            //loop over all bases to find pairs
            for (i=1;i<=ct->GetSequenceLength();i++) {
                if(tempct->GetPair(i,1)>i) ct->SetPair(i,tempct->GetPair(i,1),structures);

            }
            //copy the energy back to ct (this energy is probably inaccurate and should be updated by calling efn2(data, ct, structures, false); before saving.
            ct->SetEnergy(structures,tempct->GetEnergy(1));
        }
    }
    delete tempct;
    } // from else clause of:  if(useFastMethod && !minimum_energy)
    return 0;
}
#endif //PARTITION_CORE

// Report if there are any pseudoknots in a structure.
bool RNA::ContainsPseudoknot(const int structurenumber) {
    //make sure structurenumber is a valid structure
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
        ErrorCode=3;
        return false;
    }
    else
        return ct->HasPseudoknots(structurenumber);
}


//Get the ensemble folding free energy change as determined by the partition function.
double RNA::GetEnsembleEnergy() {

    //check to see if partition function data is present.
    if (!partitionfunctionallocated) {
        ErrorCode = 15;
        return 0.0;
    }

    //past error trapping, so set no error found
    ErrorCode = 0;

    //calculate the ensemble folding free energy change, -RT ln (Q).
    //One needs to also account for the fact that a scaling is applied per nucleotide.
    return -RKC*pfdata->pftemp*(PF_LOG_PFPRECISION(w5[ct->GetSequenceLength()])-ct->GetSequenceLength()*PF_LOG_PFPRECISION(pfdata->scaling));
    // return -RKC*GetTemperature()*log(exp(LOGPF(w5[ct->GetSequenceLength()])-ct->GetSequenceLength()*LOGPF(pfdata->scaling)));
    // return -RKC*GetTemperature()*log(exp(LOGPF(w5[ct->GetSequenceLength()])) / exp(ct->GetSequenceLength()*LOGPF(pfdata->scaling)));
    // return -RKC*GetTemperature()*log(exp(LOGPF(w5[ct->GetSequenceLength()])) / exp(LOGPF(pfdata->scaling))*ct->GetSequenceLength());
    // return -RKC*GetTemperature()*log(exp(LOGPF(w5[ct->GetSequenceLength()])) / pow(exp(LOGPF(pfdata->scaling)),ct->GetSequenceLength()));
    // return -RKC*GetTemperature()*log(TO_LINEAR(w5[ct->GetSequenceLength()])) / pow(TO_LINEAR(pfdata->scaling),ct->GetSequenceLength());
    // return -RKC*GetTemperature()*log(TO_LINEAR(w5[ct->GetSequenceLength()])) / pow(TO_LINEAR(pfdata->scaling),ct->GetSequenceLength());

    // linear: exp(LOGPF(X) = exp(log(X)) = X
    // log_double: exp(LOGPF(X) = exp(log(X_ld)) = (double)X_ld
    // LOG_CALC: exp(LOGPF(X) = exp(X)
    
}

// Gets the ensemble defect.
double RNA::GetEnsembleDefect(const int structurenumber, const int start, const int end, const char filename[]) {
	ofstream output;

	//if filename is specified, open the file and write a one-line header
	if (strcmp(filename,"")!=0) {
		output.open(filename,ios_base::app);
		output << "Structure " << structurenumber << "\n";
		output << "nucleotide\tdefect\n";
	}
	
	//check to see if partition function data is present.
	if (!partitionfunctionallocated)
		PartitionFunction();

	int pos_start = start;
	int pos_end = end;

	//The default is to do a global calculaton (i.e. start == 0 and end == 0)
	if (pos_start == 0) pos_start = 1;
	if (pos_end == 0) pos_end = GetSequenceLength();

	double defect=0.0;//initialize the 'defect'
	ct->BreakPseudoknots(structurenumber); // break pseudoknots in the structure (because they are not represented in the pair probabilities generated by partition)
	
	//Calculate Ensemble Defect (ED)
	for(int i=pos_start;i<=pos_end;++i){//for every nucleotide in the structure
		double nucdefect = 0.0;

		//Calculate ED for unpaired nucleotides
		if (GetPair(i, structurenumber) == 0) {
			//Sum all probabilities of nucleotide 'i' to form pairs with all other nucleotides - 'j'
			for (int j = 1; j <= GetSequenceLength(); ++j) {
				if (i < j) {
					
					nucdefect += GetPairProbability(i, j);//add the probabilities to the 'defect'
				}
				else if (j < i)
					nucdefect += GetPairProbability(j, i);
			}
		}
		else {  //GetPair(i, structurenumber) !=0
			if (GetPair(i, structurenumber) < i) {
				nucdefect += (1.0 - GetPairProbability(GetPair(i, structurenumber), i));
			}
			else {//GetPair(i, structurenumber) > i
				nucdefect += (1.0 - GetPairProbability(i,GetPair(i, structurenumber)));
			}
		}
		defect += nucdefect;

		//if filename is specified, write the nucdefect to the file:
		if (strcmp(filename, "") != 0) {
			output<<i<<"\t"<<nucdefect<<"\n";
		}

	}//END: Calculate Ensemble Defect

	//if filename is specified, close the file
	if (strcmp(filename, "") != 0) {
		output.close();
	}

	return defect;
}


//Get the folding free energy change for a predicted structure.
double RNA::GetFreeEnergy(const int structurenumber) {
    //make sure structurenumber is a valid structure
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
        ErrorCode=3;
        return 0.0;
    }
    else {
        //error trapping complete
        return ct->GetEnergy(structurenumber)/(double)conversionfactor;
    }
}

// Get the nucleotide to which the specified nucleotide is paired.
int RNA::GetPair(const int i, const int structurenumber) {

    //check to make sure i is a valid nucleotide
    if (i<1||i>ct->GetSequenceLength()) {
        ErrorCode=4;
        return 0;
    }
    //make sure there is structure information for this RNA
    else if (ct->GetNumberofStructures()==0){
        ErrorCode = 23;
        return 0;
    }
    //now make sure structurenumber is a valid structure
    else if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
        ErrorCode=3;
        return 0;
    }
    else {
        //error trapping complete
        return ct->GetPair(i,structurenumber);
    }
}

//Extract the lowest free energy for a structure that contains the i-j pair using data from a save file (.sav).
double RNA::GetPairEnergy(const int i, const int j) {
    int locali, localj;

    //check whether the a save file was read (by the constructor)
    if (!energyallocated) {
        ErrorCode = 17;
        return 0.0;
    }


    if (i<1||i>ct->GetSequenceLength()) {
        //i is out of range
        ErrorCode = 4;
        return 0.0;
    }
    if (j<1||j>ct->GetSequenceLength()) {
        //j is out of range
        ErrorCode = 4;
        return 0.0;
    }


    //sort indexes from 5' to 3':
    if (i>j) {
        locali = j;
        localj = i;
    }
    else {
        locali=i;
        localj=j;
    }

    //No error;
    ErrorCode = 0;

    //calculate and return the energy
    return ((((double) (ev->f(locali,localj)+ev->f(localj,locali+ct->GetSequenceLength())))/conversionfactor));



}

//Get the total number of specified structures
int RNA::GetStructureNumber() const {

        return ct->GetNumberofStructures();

}

//return the base pairing probability between nucleotides i and j
double RNA::GetPairProbability(const int i, const int j) {

    //check to see if partition function data is present.
    if (!partitionfunctionallocated) {
        ErrorCode = 15;
        return 0.0;
    }

    //check that the nucleotides are in the correct range
    if (i<1||j>ct->GetSequenceLength()||j<0||j>ct->GetSequenceLength()) {
        ErrorCode = 4;
        return 0.0;
    }

    //past error trapping, so set no error found
    ErrorCode = 0;

    //calculate the base pair probability
    return (double) calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce);


}

// Fills an array with all basepair probabilities
//  Function requires that the partition function data be present either because PartitionFunction() 
// has been called or the constructor that reads a partition function save was used.  
// This function generates internal error codes that can be accessed by GetErrorCode(): 0 = no error, nonzero = error.
// The errorcode can be resolved to a c string using GetErrorMessage.
// param arr: A pointer to an array of doubles to receive the probabilities. The array should be allocated to hold (N-1)*N/2 elments, where N is the sequence length.
//        The array will be filled with probabilities in upper-triangular format:  e.g. for N=5 the array would be: 
//       P(1:2) P(1:3) P(1:4) P(1:5)   P(2:3) P(2:4) P(2:5)   P(3:4) P(3:5)   P(4:5)   (i.e. 5*4/2 = 10 elements)
//       And the position of P(I, J)  where I < J is at 
//                    (2*N-I)*(I-1)/2+J-I-1 (if I and J are 1-based) or 
//                    (2*N-I-1)*I/2+J-I-1   (if I and J are 0-based)
// param size: The size of the array of doubles. If this is less than the required size, the array will not be modified and the function will return the required array size.
// return: The total size required to hold the array of probabilities. This is equal to (N-1)*N/2 elments, where N is the sequence length. 
//        The function returns a negative number if an error occurred. Call GetErrorMessage(-returnValue) to get a textual description of the error.
int RNA::GetPairProbabilities(double* arr, const int size){
    int len = ct->GetSequenceLength();
    int required = (len-1)*len/2;
    if (size < required)
        return required;

    //check to see if partition function data is present.
    if (!partitionfunctionallocated)
        return -15; // return negative error code.
    
    int pos = 0;
    for(int i=0; i < len; i++)
        for(int j=i+1; j < len; j++)
            arr[pos++] = calculateprobability(i+1,j+1,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce);

    return required;
}

//Provide the comment from the ct file as a string.
std::string RNA::GetCommentString(const int structurenumber) {
    //Add some code for backwards compatibility:
    //In the past, all labels were associated with structures.  Now, labels can come with the sequence.
    //If there are no structures, return the sequence label:

    if (ct->GetNumberofStructures()==0||structurenumber==-1) 
        return ct->GetSequenceLabel();


    //start with error checking:
    if (structurenumber<1||structurenumber>ct->GetNumberofStructures()){
        //The request is for a structure that is out of range
        ErrorCode = 3;
        return "";
    }

    return ct->GetCtLabel(structurenumber);
}

#ifndef PARTITION_CORE //exclude from minimal partition-only builds (e.g. partition-rosetta)
//Determine the coordinates for drawing a secondary structure.
int RNA::DetermineDrawingCoordinates(const int height, const int width, const int structurenumber) {

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;

    //First check for errors in the specification of the structure number:
    if (structurenumber<0||structurenumber>ct->GetNumberofStructures()) return 3;

    if (!drawallocated) {
        //If the memory has not been allocated for coordinates, go ahead and allocate it now
        structurecoordinates = new coordinates(ct->GetSequenceLength());
        drawallocated = true;
    }

    //now perform the calculation for coordinate determination:
    place(structurenumber, ct, structurecoordinates, height, width);

    //no errors, return 0
    return 0;


}

//Get the X coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideXCoordinate(const int i) {

    if (!drawallocated) {
        //The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
        ErrorCode = 19;
        return 0;
    }

    if (i<0||i>ct->GetSequenceLength()) {
        //The nucleotide is invalid, indicate an error.
        ErrorCode = 4;
        return 0;
    }

    //Fetch the coordinate from the coordinates structure, structurecoordinates.
    return structurecoordinates->x[i];


}

//Get the Y coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideYCoordinate(const int i) {

    if (!drawallocated) {
        //The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
        ErrorCode = 19;
        return 0;
    }

    if (i<0||i>ct->GetSequenceLength()) {
        //The nucleotide is invalid, indicate an error.
        ErrorCode = 4;
        return 0;
    }

    //Fetch the coordinate from the coordinates structure, structurecoordinates.
    return structurecoordinates->y[i];


}

// Get the X coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelXCoordinate(const int i) {

    if (!drawallocated) {
        //The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
        ErrorCode = 19;
        return 0;
    }

    if (i<0||i>ct->GetSequenceLength()) {
        //The nucleotide is invalid, indicate an error.
        ErrorCode = 4;
        return 0;
    }

    if (i%10!=0) {
        //The nucleotide index is not a multiple of 10, return an error
        ErrorCode = 25;
        return 0;

    }

    //Fetch the coordinate from the coordinates structure, structurecoordinates.
    return structurecoordinates->num[i/10][0];


}

// Get the Y coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelYCoordinate(const int i) {

    if (!drawallocated) {
        //The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
        ErrorCode = 19;
        return 0;
    }

    if (i<0||i>ct->GetSequenceLength()) {
        //The nucleotide is invalid, indicate an error.
        ErrorCode = 4;
        return 0;
    }

    if (i%10!=0) {
        //The nucleotide index is not a multiple of 10, return an error
        ErrorCode = 25;
        return 0;

    }

    //Fetch the coordinate from the coordinates structure, structurecoordinates.
    return structurecoordinates->num[i/10][1];

}
#endif // PARTITION_CORE

//Get the identity of nucleotide i.
char RNA::GetNucleotide(const int i) {

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) {
        ErrorCode = 20;
        return '-';
    }

    //Check that nucleotide is valid
    if (i<1||i>ct->GetSequenceLength()) {
        ErrorCode = 4;//i is out of range
        return '-';
    }
    return ct->nucs[i];
}

//Get the total length of the sequence
int RNA::GetSequenceLength() const {
    return ct->GetSequenceLength();
}

const char* RNA::GetSequence() const {
    return ct->GetSequence();
}

std::string RNA::GetSequence(size_t start, size_t length) const {
    if (start<1) start=1;
    if (start>GetSequenceLength()) return "";
    if (length==string::npos)
        length=GetSequenceLength()-start;
    length = std::min(length, GetSequenceLength()-start);
    return string(ct->nucs+start, length);
}

//Return the type of backbone (true = RNA, false = DNA).
bool RNA::GetBackboneType() const {
    return isrna;
}


//Access the underlying structure class.
structure *RNA::GetStructure() {
    return ct;
}


//Provide a ProgressHandler for following calculation progress.
//A ProgressHandler class has a public function void update(int percent) that indicates the progress of a long calculation.
void RNA::SetProgress(ProgressHandler& Progress) {
    progress = &Progress;
    return;
}


//Provide a means to stop using a ProgressHandler.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the ProgressHandler is deleted, so that this class does not make reference to it.
void RNA::StopProgress() {
    progress=NULL;
    return;
}

ProgressHandler* RNA::GetProgress() {
        return progress;
}


RNA::~RNA() {





    if (partitionfunctionallocated) {
        //The partition function calculation was performed, so the memory allocated for the partition function needs to be deleted.
        delete[] lfce;
        delete[] mod;
        delete[] w5;
        delete[] w3;
        delete v;
        delete w;
        delete wmb;
        delete wl;
		delete wlc;
        delete wmbl;
        delete wcoax;
        delete fce;
        delete pfdata;

    }

    if (energyallocated) {
        //A folding save file was opened, so clean up the memory use.

        delete[] lfce;
        delete[] mod;

        delete[] ew5;
        delete[] ew3;

        if (ct->intermolecular) {
            delete ew2;
            delete ewmb2;
        }

        delete ev;
        delete ew;
        delete ewmb;
        delete fce;
    }

    #ifndef PARTITION_CORE // exclude these from minimal partition-only builds (e.g. partition-rosetta)
    if (drawallocated) {
        //The drawing coordiantes have been determined, so the memory needs to be cleaned up.
        delete structurecoordinates;

    }
    #endif //PARTITION_CORE

    delete ct;//delete the structure
}

// string getFileNameWithoutExt(const char*)


//This is a protected function for handling file input.
int RNA::FileReader(const char filename[], const RNAInputType type) {

    if (!isStdIoFile(filename) && !fileExists(filename)) {
        SetErrorDetails(sfmt("The path '%s' is invalid or does not exist.", filename));
        return 1; // file not found.
    }

    if (type==FILE_CT||type==FILE_SEQ||type==FILE_DBN)
        if (!IsAlphabetRead()) return 30;

    // RMW 2015-03-12: try/catch to prevent any unexpected errors from crashing the program. At least show an error message.
    // (Previous to this, passing the wrong type of file to openct caused an unhandled memory allocation exception.)
    try {
        //open the file based on type:
        switch(type) {
            case FILE_MSA: // ct file   
                ct->openseqx(filename);
                return ct->open_alignment(filename); 
            case FILE_CT: // ct file
                return ct->openct(filename);  // openct returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            case FILE_DBN: // dot-bracket file
                return ct->opendbn(filename); // opendbn returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            case FILE_SEQ:
                //type indicates a .seq file
                return ct->openseqx(filename);
            case FILE_PFS: { //partition function save file
                LOG_INFO("Reading PFS file")
                short vers;
                //allocate the ct file by reading the save file to get the sequence length:
                std::ifstream sav(filename,std::ios::binary);

                read(&sav,&(vers));//read the version of the save file

                if (vers!=pfsaveversion) {
                    //Wrong version!
                    sav.close();
                    return 16;
                }

                int sequencelength;
                //read the length of the sequence
                read(&sav,&(sequencelength));
                sav.close();

                LOG_INFO("Allocate Arrays")
                //allocate everything
                ct->allocate(sequencelength);

                w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                fce = new forceclass(ct->GetSequenceLength());

                w5 = new PFPRECISION [ct->GetSequenceLength()+1];
                w3 = new PFPRECISION [ct->GetSequenceLength()+2];

                lfce = new bool [2*ct->GetSequenceLength()+1];
                mod = new bool [2*ct->GetSequenceLength()+1];

                pfdata = new pfdatatable();
				data = new datatable();

                //indicate that the memory has been allocated so that the destructor will delete it.
                partitionfunctionallocated = true;
                LOG_INFO("Reading save data")

                //load all the data from the pfsavefile:
                readpfsave(filename, ct, w5, w3,v, w, wmb,wl, wlc,wmbl, wcoax, fce,&pfdata->scaling,mod,lfce,pfdata,data);
                return 0;
            }  // FILE_PFS
            #ifndef PARTITION_CORE // exclude these from minimal partition-only builds (e.g. partition-rosetta)
            case FILE_SAV: { // folding save file.
                short vers; 

                //peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
                std::ifstream sav(filename,std::ios::binary);

                //Read the save file version information.
                read(&sav,&vers);

                //check the version
                if (vers!=safiversion) {
                    //Wrong version!
                    sav.close();
                    return 16;
                }

                int sequencelength;
                read(&sav,&(sequencelength));
                read(&sav,&(ct->intermolecular));
                sav.close();

                //indicate that everything is allocated and needs to be deleted in the destructor
                energyallocated = true;

                //allocate everything
                ct->allocate(sequencelength);

                ew = new DynProgArray<integersize>(ct->GetSequenceLength());
                ev = new DynProgArray<integersize>(ct->GetSequenceLength());
                ewmb = new DynProgArray<integersize>(ct->GetSequenceLength());
                fce = new forceclass(ct->GetSequenceLength());


                lfce = new bool [2*ct->GetSequenceLength()+1];
                mod = new bool [2*ct->GetSequenceLength()+1];

                ew5 = new integersize [ct->GetSequenceLength()+1];
                ew3 = new integersize [ct->GetSequenceLength()+2];

                if (ct->intermolecular) {
                    ew2 = new DynProgArray<integersize>(ct->GetSequenceLength());
                    ewmb2 = new DynProgArray<integersize>(ct->GetSequenceLength());

                    for (int i=0;i<3;i++) read(&sav,&(ct->inter[i]));

                }
                else {
                    ew2 = NULL;
                    ewmb2 = NULL;
                }

                //indicate that the thermodynamic parameters are read and available (and need to be deleted).
                data = new datatable();

                //now read the file.
                readsav(filename, ct, ew2, ewmb2, ew5, ew3, lfce, mod, data,
                         ev, ew, ewmb, fce, &vmin);

                //set error status to zero
                 return 0;
            }  //FILE_SAV
            #endif //PARTITION_CORE
            default:
                return 22; // error - invalid file type
        } // SWITCH type
    } catch (std::exception* ex) {
        SetErrorDetails(ex->what());
        return 2;
    }
}


//The following should not be included for compilations for Windows:
#ifndef _WINDOWS_GUI

//A global function for error reporting
void errmsg(int err,int erri) {

if (err==30) {
    std::cerr << "End Reached at traceback #"<<erri<<"\n";
   return;
}
if (err==100) {
    std::cerr << "error # "<<erri;
   return;
}
switch (err) {
    case 1:
    std::cerr << "Could not allocate enough memory";
      break;
   case 2:
    std::cerr << "Too many possible base pairs";
      break;
   case 3:
    std::cerr << "Too many helixes in multibranch loop";
   case 4:
    std::cerr << "Too many structures in CT file";
   default:
    std::cerr << "Unknown error";
}
//std::cin >> err;
return;

}

#endif

#ifdef COUNTING
bool RNA::WriteDataCounters(string count_File, std::vector<double> counts) {
	if (!VerifyThermodynamic()) {
		return 0;
	}

    ofstream ofile(count_File.c_str());
    if (!ofile.good()) return false;
    
    for (int i=0; i<counts.size(); i++){
        ofile << setprecision(4) << fixed << counts[i] << endl;
    }

    return ofile.good();
}

std::vector<double> RNA::GetDataCounters(){
    return data->get_data_counts();
}

bool RNA::ResetDataCounters(){
    if (!VerifyThermodynamic()) {
		return 0;
	}

    data->clear_parameter_usage_data();

    return 1;
}
double RNA::CalculateUncertainty(){
    return data->calculate_uncertainty();
}

#endif //COUNTING
