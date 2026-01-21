//The TwoRNA class code.

#include "TwoRNA.h"



//!Constructor - user provides a sequences as c strings.
TwoRNA::TwoRNA(const char sequence1[], const char sequence2[], bool isRNA) {
	lastErrorDetails="";

	//Use sequence1 in the constructor for rna1.
	rna1=new RNA(sequence1, isRNA);

	//Use sequence2 in the constructor for rna2.
	rna2=new RNA(sequence2, SEQUENCE_STRING, rna1); // use the RNA(sequence,fileType,*Thermo) constructor to use the datatable from rna1 (so that it isn't re-read from disk)

	ErrorCodeTwo=0;
	if (rna1->GetErrorCode()!=0) ErrorCodeTwo=1000;
	if (rna2->GetErrorCode()!=0) ErrorCodeTwo+=2000;
}


//Constructor
//Open filename1 for the first sequence and filename2 for the second sequence.
//These are type1 for the first seequence and type2 for the second sequence.
TwoRNA::TwoRNA(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, bool isRNA) {
	lastErrorDetails="";

	//Use the rna2 constructor for the second sequence.
	rna1 = new RNA(filename1,type1,isRNA);

	//Use the rna2 constructor for the second sequence.
	rna2 = new RNA(filename2,type2,rna1); // use the RNA(sequence,fileType,*Thermo) constructor to use the datatable from rna1 (so that it isn't re-read from disk)

	ErrorCodeTwo=0;
	if (rna1->GetErrorCode()!=0) ErrorCodeTwo=1000;
	if (rna2->GetErrorCode()!=0) ErrorCodeTwo+=2000;
}

//Constructor
//Open filename1 for the first sequence and filename2 for the second sequence.
//These are type1 for the first seequence and type2 for the second sequence.
TwoRNA::TwoRNA(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, const Thermodynamics *copyThermo) {
	lastErrorDetails="";

	//Use the rna2 constructor for the second sequence.
	rna1 = new RNA(filename1,type1,copyThermo);

	//Use the rna2 constructor for the second sequence.
	rna2 = new RNA(filename2,type2,copyThermo);

	ErrorCodeTwo=0;
	if (rna1->GetErrorCode()!=0) ErrorCodeTwo=1000;
	if (rna2->GetErrorCode()!=0) ErrorCodeTwo+=2000;
}


//!Constructor
//!Default constructor that requires no parameters.
TwoRNA::TwoRNA() {
	ErrorCodeTwo=0;
	lastErrorDetails="";

	//Use the rna1 constructor for the second sequence.
	rna1 = new RNA();

	//Use the rna2 constructor for the second sequence.
	rna2 = new RNA();
}

//Set the current temperature in K for calculations.
int TwoRNA::SetTemperature(double temperature) {


	//Temperature only needs to be set in rna1 because the thermodynamic parameters from rna1 are used.
	return rna1->SetTemperature(temperature);


}

//Get the current folding temperature in K.
double TwoRNA::GetTemperature() const {

	//Temperature only needs to be checked in rna1 because the thermodynamic parameters from rna1 are used.
	return rna1->GetTemperature();

}


//Return the value of ErrorCode
int TwoRNA::GetErrorCode() const {
	return ErrorCodeTwo;
}

//Reset the error codes to 0
void TwoRNA::ResetError() {
    rna1->ResetError();
    rna2->ResetError();
	lastErrorDetails="";
	ErrorCodeTwo=0;
}

//Return error messages based on code from GetErrorCode and other error codes.		

//		0 = no error
//		1000 = Error associated with sequence 1 or with procedure, get message from sequence 1
//		2000 = Error associated with sequence 2, get message from sequence 2
//		3000 = Errors with each sequence, get messages from each


const char* TwoRNA::GetErrorMessage(const int error) {
	

	if (error==0) return "No Error.\n";
	else if (error==1000) {
			//This is a message that refers to sequence 1 (file i/o).
			strcpy(compoundmessage,"Error in Sequence 1: ");
			strcat(compoundmessage,rna1->GetErrorMessage(error-1000));
	}
	
	
	else if (error==2000) {
		//2 indicates that there is an error from the second sequence that needs to be accessed from rna2:
		strcpy(compoundmessage,"Error in Sequence 2: ");
		strcat(compoundmessage,rna2->GetErrorMessage(error-2000));

	}
	else if (error==3000) {
		//3 indicates errors in each sequence
		strcpy(compoundmessage,"Error in Sequence 1: ");
		strcat(compoundmessage,rna1->GetErrorMessage(rna1->GetErrorCode()));
		strcat(compoundmessage,"Error in Sequence 2: ");
		strcat(compoundmessage,rna2->GetErrorMessage(rna2->GetErrorCode()));

	}
	else {

		strcpy(compoundmessage,"Unknown Error Occurred\n");
	}
	
	//This is only reached in the state when there is an error.
	return compoundmessage;
}

//Return error messages based on code from GetErrorCode and other error codes.		
std::string TwoRNA::GetErrorMessageString(const int error) {
	std::string temp;
	temp = GetErrorMessage(error);
	return temp;
}

const string TwoRNA::GetErrorDetails() const {
	if (!lastErrorDetails.empty()) return lastErrorDetails;
	return  !rna1->GetErrorDetails().empty() ? rna1->GetErrorDetails() : rna2->GetErrorDetails();
}
void TwoRNA::SetErrorDetails(const string& details) {
	lastErrorDetails = details;
}

//!return A pointer to the underlying structure class for sequence 1.
RNA *TwoRNA::GetRNA1() {
	return rna1;
}


//!return A pointer to the underlying structure class for sequence 2.
RNA *TwoRNA::GetRNA2() {
	return rna2;
}


//Destructor
TwoRNA::~TwoRNA() {

	
	delete rna1;
	delete rna2;


}
