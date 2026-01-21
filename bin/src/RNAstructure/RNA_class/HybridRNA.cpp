//HybridRNA source code.


#include "HybridRNA.h"
#include "../src/bimol.h"


//Constructors:
//Pass the constructor parameters directly to the parent TwoRNA class.

HybridRNA::HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, Thermodynamics *copyThermo) {
	
	sequences = new TwoRNA(filename1, type1, filename2, type2, copyThermo);
	//By default, allow intramolecular pairs
	forbidunimolecular = false;
}

HybridRNA::HybridRNA(const char sequence1[], const char sequence2[], const bool IsRNA)
	: RNA(NULL/*sequence*/, IsRNA) { // this constructor loads the thermodynamic datatable (whereas RNA(bool) does not. No Sequence is loaded.
	// Load the two sub-sequences
	commonconstructor(sequence1, SEQUENCE_STRING, sequence2, SEQUENCE_STRING);
}

HybridRNA::HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, const bool IsRNA) 
  : RNA(NULL/*sequence*/, IsRNA) { // this constructor loads the thermodynamic datatable (whereas RNA(bool) does not. No Sequence is loaded.
	// Load the two sub-sequences
	commonconstructor(filename1, type1, filename2, type2);
}

//This constructor read the thermodynamic parameters associated with library
HybridRNA::HybridRNA(const char filename1[], const int type1, const char filename2[], const int type2, const char *alphabet)
  : RNA(NULL/*sequence*/,0/*file-type*/,alphabet) { // this constructor loads the thermodynamic datatable (whereas RNA(bool) does not. No Sequence is loaded.
	// Load the two sub-sequences
	commonconstructor(filename1, type1, filename2, type2);
}

//perform common constructor tasks
void HybridRNA::commonconstructor(const char filename1[], const int type1, const char filename2[], const int type2) {
	// Read both sequences and use THIS as the Thermodynamics reference (i.e. use the datatable that was already loaded in the constructor)
	sequences = new TwoRNA(filename1, type1, filename2, type2, this);
	//By default, allow intramolecular pairs
	forbidunimolecular = false;	
}

//Predict the bimolecular structure, no intramolecular pairs, and consider accessibility
int HybridRNA::AccessFold(const double gamma, const float percent, const int maximumstructures, const int window, const int maxinternalloopsize) {

	//check to make sure that a sequence has been read for each sequence
	if (GetRNA1()->GetStructure()->GetSequenceLength()==0) return 20;
	if (GetRNA2()->GetStructure()->GetSequenceLength()==0) return 20;

	//now read the thermodynamic parameters, using the inheritance of the base class.
	if (!VerifyThermodynamic()) return 5;//record an error		

	accessfold(GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), GetStructure(), maxinternalloopsize, maximumstructures, ((int)percent), window, 
		GetDatatable(), gamma, GetBackboneType(), GetTemperature());

	//also record the fact that the RNADuplex structure is intermolecular, in case the user manipulates the underlying ct
	GetStructure()->intermolecular = true;
	GetStructure()->inter[0]=GetRNA1()->GetStructure()->GetSequenceLength()+1;
	GetStructure()->inter[1]=GetRNA1()->GetStructure()->GetSequenceLength()+2;
	GetStructure()->inter[2]=GetRNA1()->GetStructure()->GetSequenceLength()+3;

	return 0;
}
		

//Predict the bimolecular secondary structure
int HybridRNA::FoldBimolecular(const float percent, const int maximumstructures, const int window, const char savefile[], const int maxinternalloopsize) {
	int i,j;
	
	if (!VerifyThermodynamic()) return 5;

	//First do some error trapping:
	
	//check to make sure that a sequence has been read for each sequence
	if (GetRNA1()->GetStructure()->GetSequenceLength()==0) return 20;
	if (GetRNA2()->GetStructure()->GetSequenceLength()==0) return 20;

	
	//Now populate the sequence data from each sequence to the common structure in the RNA base class 
	SetupBimolecular();



   if (forbidunimolecular) {
		//forbid unimolecular pairs
	   GetStructure()->allocatetem();
	   for (i=1;i<GetRNA1()->GetStructure()->GetSequenceLength();i++) {
		   for (j=i+1;j<=GetRNA1()->GetStructure()->GetSequenceLength();j++) {
				GetStructure()->tem[j][i]=false;
		   }

	   }
	   for (i=GetRNA1()->GetStructure()->GetSequenceLength()+3;i<GetStructure()->GetSequenceLength();i++) {
		   for (j=i+1;j<=GetStructure()->GetSequenceLength();j++) {
				GetStructure()->tem[j][i]=false;
		   }

	   }
   }

   //Now that the sequence has been set up, the base class, single sequence method can be used
   return RNA::FoldSingleStrand(percent,maximumstructures,window,savefile,maxinternalloopsize);

    


}

//Perform a simple bimolecular structure prediction with no intramolecular pairs
int HybridRNA::FoldDuplex(const float percent, const int maximumstructures, const int window, const int maxinternalloopsize) {

	
	//check to make sure that a sequence has been read for each sequence
	if (GetRNA1()->GetStructure()->GetSequenceLength()==0) return 20;
	if (GetRNA2()->GetStructure()->GetSequenceLength()==0) return 20;

	//now read the thermodynamic parameters, using the inheritance of the base class.
	if (!VerifyThermodynamic()) return 5;//record an error		

	bimol(GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), GetStructure(), maxinternalloopsize, maximumstructures, ((int)percent), window, GetDatatable());

	//also record the fact that the RNADuplex structure is intermolecular, in case the user manipulates the underlying ct
	GetStructure()->intermolecular = true;
	GetStructure()->inter[0]=GetRNA1()->GetStructure()->GetSequenceLength()+1;
	GetStructure()->inter[1]=GetRNA1()->GetStructure()->GetSequenceLength()+2;
	GetStructure()->inter[2]=GetRNA1()->GetStructure()->GetSequenceLength()+3;
	
	return 0;

}


//Perform a bimolecular partition function calculations.
int HybridRNA::PartitionFunctionBimolecular(const char savefile[]) {


	//First do some error trapping:
	
	//check to make sure that a sequence has been read for each sequence
	if (GetRNA1()->GetStructure()->GetSequenceLength()==0) return 20;
	if (GetRNA2()->GetStructure()->GetSequenceLength()==0) return 20;

	
	//Now populate the sequence data from each sequence to the common structure in the RNA base class 
	SetupBimolecular();

	//Now that the composite sequence is set up, the base class function can be used:
	return PartitionFunction(savefile);

}

// Provide a pointer to the underlying RNA class for sequence 1.
RNA *HybridRNA::GetRNA1() {

	return sequences->GetRNA1();

}



// Provide a pointer to the underlying RNA class for sequence 2.
RNA *HybridRNA::GetRNA2(){

	return sequences->GetRNA2();

}




//Return the value of ErrorCode
int HybridRNA::GetErrorCode() {
	return sequences->GetErrorCode()+ErrorCode;
}

//Return error messages based on code from GetErrorCode and other error codes.		
 const char* HybridRNA::GetErrorMessage(const int error) {
	

	if (error==0) return "No Error.\n";
	else if (error>=1000) {
		//This error was generated by an TwoRNA base class event.
		return sequences->GetErrorMessage(error);
	}
	else if (error<100) {
		//This error was generated by the RNA base class or a base class message applies:
		return RNA::GetErrorMessage(error);

	}
	else return "Unknown Error\n";
	
}

// Get whether intramolecular pairs are allowed.		
bool HybridRNA::GetForbidIntramolecular() {
	return forbidunimolecular;
}


// Set whether intramolecular pairs are allowed.
void HybridRNA::SetForbidIntramolecular(const bool forbid) {

	forbidunimolecular = forbid;
}


//Provide a TProgressDialog for following calculation progress.
void HybridRNA::SetProgress(ProgressHandler& Progress) {

	progress = &Progress;

	return;
}


//Provide a means to stop using a TProgressDialog.
void HybridRNA::StopProgress() {

		
	

	progress=NULL;
	return;

}


//Destructor
HybridRNA::~HybridRNA() {
	
	//Delete the instance of TwoRNA:
	delete sequences;

}


//This function makes a composite sequence in the RNA::ct with each of the sequences in TwoRNA
void HybridRNA::SetupBimolecular() {
	int i,j;

	RNA *r1 = sequences->GetRNA1(), *r2 = sequences->GetRNA2();
	structure *ct = GetStructure();

	//Copy the data table information here:
	// If this datatable has NOT been read, copy it (if possible) or read it.

	// if (GetEnergyRead())
	// 	cout << "Hybrid Energy already read. Temperature: " << GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "Hybrid Energy NOT read." << endl;

	// if (r1-GetEnergyRead())
	// 	cout << "RNA1 Energy already read. Temperature: " << r1->GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "RNA1 Energy NOT read." << endl;

	// if (r2-GetEnergyRead())
	// 	cout << "RNA2 Energy already read. Temperature: " << r2->GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "RNA2 Energy NOT read." << endl;

	if (!GetEnergyRead()) {
		if (r1->GetEnergyRead())
			CopyThermo(*r1);
		else if (r2->GetEnergyRead())
			CopyThermo(*r2);
	} else {
		// Set the datatables in RNA1 and RNA2 to match this, in case SetTemperature etc has been called on THIS HybridRNA object.
		if (r1->GetDatatable()!=data)
			r1->CopyThermo(*this);
		if (r2->GetDatatable()!=data)
			r2->CopyThermo(*this);
	}

	// if (GetEnergyRead())
	// 	cout << "AFTER: Hybrid Energy read. Temperature: " << GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "AFTER: Hybrid Energy NOT read." << endl;

	// if (r1-GetEnergyRead())
	// 	cout << "AFTER: RNA1 Energy read. Temperature: " << r1->GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "AFTER: RNA1 Energy NOT read." << endl;

	// if (r2-GetEnergyRead())
	// 	cout << "alreadyRNA2 Energy read. Temperature: " << r2->GetDatatable()->GetLoadedTemperature() << endl;
	// else
	// 	cout << "alreadyRNA2 Energy NOT read." << endl;

	if (GetDatatable()!=GetStructure()->GetThermodynamicDataTable())
		cout << "Datatable mismatch!" << endl;

	if (r1->GetDatatable()!=GetDatatable()||r1->GetStructure()->GetThermodynamicDataTable()!=GetDatatable())
		cout << "RNA1 Datatable mismatch!" << endl;

	if (r2->GetDatatable()!=GetDatatable()||r2->GetStructure()->GetThermodynamicDataTable()!=GetDatatable())
		cout << "RNA2 Datatable mismatch!" << endl;


	//First check if this has been done:
	if (ct->allocated) return;

	//copy the labels over
	string label;
	label = r1->GetStructure()->GetSequenceLabel();

	//Make sure there is no newline at the end of this line.
	if (label[label.size()-1]=='\n') {
		label.erase(label.size()-1,1);		
	}

	label+="_";
	label+=r2->GetStructure()->GetSequenceLabel();


	ct->SetSequenceLabel(label);
	

	//Sequence length is total of each sequence plus a three nuc linker
	
	ct->allocate(r1->GetStructure()->GetSequenceLength()+r2->GetStructure()->GetSequenceLength()+3);

	for (i=1;i<=r1->GetStructure()->GetSequenceLength();i++) {
		ct->numseq[i] = r1->GetStructure()->numseq[i];
		ct->nucs[i] = r1->GetStructure()->nucs[i];
		ct->hnumber[i] = r1->GetStructure()->hnumber[i];

	}
	
	for (i=1;i<=r2->GetStructure()->GetSequenceLength();i++) {
		ct->numseq[i+r1->GetStructure()->GetSequenceLength()+3] = r2->GetStructure()->numseq[i];
		ct->nucs[i+r1->GetStructure()->GetSequenceLength()+3] = r2->GetStructure()->nucs[i];
		ct->hnumber[i+r1->GetStructure()->GetSequenceLength()+3] = r2->GetStructure()->hnumber[i];

	} 	
      
   

     //add intermolecular linker
  for (i=r1->GetStructure()->GetSequenceLength()+1; i<(r1->GetStructure()->GetSequenceLength()+4); ++i) {

	//ToDo: Could guard this code to make sure the linker exists in this dataset.

		//use the first linker entry (likely the only linker entry)
		ct->numseq[i]=data->basetonum(data->linker[0]);
		ct->nucs[i] = data->linker[0];
		ct->hnumber[i] = 0;
  }

 


   ct->inter[0] = r1->GetStructure()->GetSequenceLength()+1;
   ct->inter[1] = r1->GetStructure()->GetSequenceLength()+2;
   ct->inter[2] = r1->GetStructure()->GetSequenceLength()+3;

   ct->intermolecular = true;


   //Also copy information about nucleotides that must be single stranded
		//(These were entered as lowercase by the user.)

   for (i=0;i<r1->GetStructure()->GetNumberofSingles();++i) {
		ct->AddSingle(r1->GetStructure()->GetSingle(i));
   }
   for (i=0;i<r2->GetStructure()->GetNumberofSingles();++i) {
		ct->AddSingle(r2->GetStructure()->GetSingle(i)+r1->GetStructure()->GetSequenceLength()+3);
   }

}

