//Dynalign-object source code.


#include "Dynalign_object.h"
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"
//#include "../src/phmm/aln_env_utils.h"
//#include "../src/align_services.h"
//#include "../src/hmm_arrays.h"

#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h" 
#else
	#include "../src/platform.h"
#endif //_WINDOWS


#include <cstdio>
#include <iostream>


//Constructors:

//Default constructor:
Dynalign_object::Dynalign_object()
	:TwoRNA() {

}

//Pass the constructor parameters directly to the parent TwoRNA class.
Dynalign_object::Dynalign_object(const char sequence1[], const char sequence2[], const bool IsRNA)
	:TwoRNA(sequence1, sequence2, IsRNA) {

	//Do some initilaizations
	CommonConstructor();

}

Dynalign_object::Dynalign_object(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, const bool IsRNA)
	: TwoRNA(filename1, type1, filename2, type2, IsRNA) {

	//Do some initilaizations
	CommonConstructor();


}

Dynalign_object::Dynalign_object(const char filename1[], const RNAInputType type1, const char filename2[], const RNAInputType type2, const Thermodynamics *copyThermo)
	: TwoRNA(filename1, type1, filename2, type2, copyThermo) {

	//Do some initilaizations
	CommonConstructor();


}

//Constructor for reading a dynalign save file.
Dynalign_object::Dynalign_object(const char filename[]): TwoRNA() {
	
	//Open the save file, make the output file
	short i;
	//integersize en1;
	
	short maxsep;
	
	//integersize crit;
	
	
	dynalignheap heap;

	//Do some initializations
	CommonConstructor();

	//Make sure the filename exists:
	if ((fopen(filename, "r"))== NULL) {
		//the file is not found
		ErrorCodeTwo=106;
		return;	
	}

	//open the save file to peek at the sizes needed to allocate arrays:
	std::ifstream sav(filename, std::ios::binary);

	//Check the save file version
	int versioncheck;
	read(&sav, &versioncheck);
	//Make sure the right version is being used:
	if (versioncheck != dynalignsavefileversion) {
		//the file is not found
		ErrorCodeTwo = 111;
		return;
	}
	
	data = new datatable;


#ifdef DYNALIGN_II
        int max_elongation;
        short islope,iintercept;
#else
        bool singleinsert;
#endif
        bool local;
	bool **allowed_alignments;
		//indicate that a save file was read so that the destructor can do the cleanup
	savefileread=true;
	




	read(&sav, &modificationflag);

	int length1,length2;


	//start with structure information
	read(&sav,&(length1));
	read(&sav,&(length2));
	read(&sav,&maxsep);
	sav.close();

#ifdef DYNALIGN_II
 

  int single1_vmin=0;
  int single2_vmin=0;
  integersize *single1_w5,*single1_w3,*single2_w5,*single2_w3;
  bool *single1_lfce,*single2_lfce,*single1_mod,*single2_mod;
  // DynProgArray<integersize> *single1_w2,*single1_wmb2,*single2_w2,*single2_wmb2;
  
  //allocate everything
  DynProgArray<integersize> single1_w(length1),single2_w(length2);
  DynProgArray<integersize> single1_v(length1),single2_v(length2);
  DynProgArray<integersize> single1_wmb(length1),single2_wmb(length2);
  DynProgArray<integersize> single1_we(length1,0),single2_we(length2,0);
  forceclass single1_fce(length1),single2_fce(length2);


  single1_lfce = new bool [2*length1+1];
  single2_lfce = new bool [2*length2+1];
  single1_mod = new bool [2*length1+1];
  single2_mod = new bool [2*length2+1];
  single1_w5 = new integersize [length1+1];
  single1_w3 = new integersize [length1+2];
  single2_w5 = new integersize [length2+1];
  single2_w3 = new integersize [length2+2];
  //the following are used for chemical modfication cases
#else
#endif

	


	
	if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[length1+1];
	  for (i=0;i<=length1;++i) allowed_alignments[i]=new bool [length2+1];

	}
	else allowed_alignments=NULL;
	//fill the low and highend arrays:
	//allocate space:
	lowend = new short [2*length1];
	highend = new short [2*length1];


	if (modificationflag==1) {
		vmod = new dynalignarray();
	}
	else vmod=NULL;

	v = new varray();
	w = new dynalignarray();

	

	w3 = new wendarray();
	w5 = new wendarray();

#ifdef DYNALIGN_II
        opendynalignsavefile(filename,GetRNA1()->GetStructure(), GetRNA2()->GetStructure(),  v, w, vmod, w3, w5, &single1_w,&single2_w,&single1_v,&single2_v,&single1_wmb,&single2_wmb,&single1_we,&single2_we,
		       single1_lfce,single2_lfce,single1_mod,single2_mod,&single1_fce,&single2_fce,&single1_vmin,&single2_vmin,data, &max_elongation,  &maxsep, &islope, &iintercept, &gap, &lowest,&local,allowed_alignments,lowend,highend);
#else
	opendynalignsavefile(filename, GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), v, w, vmod, w3, w5, data, &singleinsert, &maxsep, &gap, &lowest, &local, allowed_alignments,lowend,highend);
#endif	



	if (maxsep<0) {
	  
		for (i=0;i<=length1;++i) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
	}

	
	
	//No errors encounted
	ErrorCodeTwo=0;
	return;


}

//Constructor
//This constructor is used to perform Dynalign refolding.
Dynalign_object::Dynalign_object(const char* filename, const short maxtrace, const short bpwin, const short awin, const short percent): TwoRNA() {
	int i;

	CommonConstructor();

	//Make sure the filename exists:
	if ((fopen(filename, "r"))== NULL) {
		//the file is not found
		ErrorCodeTwo=106;
		return;
		
	}



	//open the save file to peek at the sizes needed to allocate arrays:
	std::ifstream sav(filename,std::ios::binary);


	read(&sav, &modificationflag);

	//start with structure information
	int length;
	read(&sav,&length);
	
	sav.close();

	//allocate space for the alignment
	align = new short *[maxtrace];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<maxtrace;i++)  align[i] = new short [length+1];
#ifdef DYNALIGN_II
	refolddynalign(filename,GetRNA1()->GetStructure(),GetRNA2()->GetStructure(),align,maxtrace,bpwin,awin,percent, false, NULL);
#else
	refolddynalign(filename,GetRNA1()->GetStructure(),GetRNA2()->GetStructure(),align,maxtrace,bpwin,awin,percent);
#endif

	//No errors encountered.
	ErrorCodeTwo=0;
	return;

}


//Run a dynalign calculation.
//Return 0 if there is no error, return 107 if the dat files aren't found, 
#ifdef DYNALIGN_II
int Dynalign_object::Dynalign(const short int maxtrace, 
                              const short int bpwin, const short int awin, const short int percent, const short int imaxseparation, 
							  const float slope, const float intercept, const float gap, const int max_elongation,  const char savefile[], 
							  const bool optimalonly, const short int singlefold_subopt_percent, const bool local, 
	const short int numProcessors,
      const int maxpairs)
#else
int Dynalign_object::Dynalign(const short int maxtrace, 
	                          const short int bpwin, const short int awin, const short int percent, const short int imaxseparation, 
							  const float gap, const bool singleinsert, const char savefile[], 
							  const bool optimalonly, const short int singlefold_subopt_percent, const bool local, 
	const short int numProcessors,
      const int maxpairs)
#endif

 {
	bool constraints= false;
	bool **allowed_alignments;//An array for storing those nucleotide pairs that can align or not.
	int errormessage;
	int i;

	structure *ct;//have a structure pointer in case the user does template from ct

	//Read the thermodynamic parameters, if necessary, and store in RNA1:
	errormessage = GetRNA1()->VerifyThermodynamic()?0:5;

	if (errormessage!=0) return 110;

	//decide if there are folding or alignment constraints:

	else if (GetRNA1()->GetStructure()->GetNumberofPairs()>0) constraints=true;
	else if (GetRNA2()->GetStructure()->GetNumberofPairs()>0) constraints=true;
	else if (GetRNA1()->GetStructure()->GetNumberofForbiddenPairs()>0) constraints=true;
	else if (GetRNA2()->GetStructure()->GetNumberofForbiddenPairs()>0) constraints=true;
	else if (GetRNA1()->GetStructure()->GetNumberofSingles()>0) constraints=true;
	else if (GetRNA2()->GetStructure()->GetNumberofSingles()>0) constraints=true;
	else if (GetRNA1()->GetStructure()->GetNumberofModified()>0) constraints=true;
	else if (GetRNA2()->GetStructure()->GetNumberofModified()>0) constraints=true;
	else if (GetRNA1()->GetStructure()->GetNumberofGU()>0) constraints=true;
	else if (GetRNA2()->GetStructure()->GetNumberofGU()>0) constraints=true;
	//convert the gap parameter to 10*kcal/mol
	short igapincrease = (int)(gap * 10.0);
#ifdef DYNALIGN_II
        short islope = (int)(slope * 10.0);
        short iintercept = (int)(intercept * 10.0);
#else
#endif
	//do the dynalign calculation
  	
	//This section folds thie individual sequences to find insignificant pairs to be ingnored by Dynalign.
	GetRNA1()->GetStructure()->allocatetem();
	GetRNA2()->GetStructure()->allocatetem();

	//The default behavior is to fold the sequences to determine pair hat should be not allowed in dynalign
	
	if(dsv_templated) {
		if(templatefromdsv(GetRNA1()->GetStructure(), templatefilename, MAXDSV, maxpairs))
          return 109;
	}
	else if (ct_templated) {
		ct = new structure();

		ct->openct(templatefilename);

		templatefromct(ct);

		delete ct;
	}
	else templatefromfold(GetRNA1()->GetStructure(), GetRNA1()->GetDatatable(), singlefold_subopt_percent);
	templatefromfold( GetRNA2()->GetStructure(), GetRNA1()->GetDatatable(), singlefold_subopt_percent );

	//This next section determined the allowed nucleotide alignments if the HMM forward-backward is used:
	if (imaxseparation < 0) {
		//allocate space in allowed_alignments
		allowed_alignments = new bool *[GetRNA1()->GetStructure()->GetSequenceLength()+1];
		for (i=0;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			allowed_alignments[i] = new bool [GetRNA2()->GetStructure()->GetSequenceLength()+1];	
	
		}

		// Needed for having nucleotide sequences as c strings.
		GetRNA1()->GetStructure()->nucs[GetRNA1()->GetStructure()->GetSequenceLength() + 1] = 0;
		GetRNA2()->GetStructure()->nucs[GetRNA2()->GetStructure()->GetSequenceLength() + 1] = 0;

		// Sequence1 nucleotides.
	calculate_coinc_probs_env(GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), allowed_alignments, forcealign);	
	}
	else allowed_alignments = NULL;

	//allocate space for the alignment
	align = new short *[maxtrace];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<maxtrace;i++)  align[i] = new short [GetRNA1()->GetStructure()->GetSequenceLength()+1];
	//also store the maxtrace for deleting later
	Maxtrace = maxtrace;



	//The heart of the calculation is here:

#ifdef DYNALIGN_II
	errormessage = dynalign(GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), align, imaxseparation, islope, iintercept, igapincrease, GetRNA1()->GetDatatable(),
                                 maxtrace, bpwin, awin, percent, forcealign, max_elongation, allowed_alignments, GetRNA1()->GetProgress(),
           savefile, optimalonly, local,
           /*force =*/ constraints, numProcessors); 
#else	
	errormessage = dynalign(GetRNA1()->GetStructure(), GetRNA2()->GetStructure(), align, imaxseparation, igapincrease, GetRNA1()->GetDatatable(),
           singleinsert, maxtrace, bpwin, awin, percent, forcealign, allowed_alignments, GetRNA1()->GetProgress(),
           savefile, optimalonly, local,
           /*force =*/ constraints, numProcessors); 
#endif
	

	if (imaxseparation < 0) {
		//delete space in allowed_alignments
		
		for (i=0;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			delete[] allowed_alignments[i];	
	
		}

		delete[] allowed_alignments;

	}
	return errormessage;

}




// Write the alignment to disk.	
void Dynalign_object::WriteAlignment(const char filename[]) {

	//Write the alignment
	alignout(align,filename,GetRNA1()->GetStructure(),GetRNA2()->GetStructure());
	return;

}



//Force an alignment during a Dynalignc calculation).
//return An integer that indicates an error code (0 = no error, 100 = nucleotide i out of range, 101 = nucleotide k out of range). 
int Dynalign_object::ForceAlignment(const int i, const int k) {

	//check the indexes
	if (i<1||i>GetRNA1()->GetSequenceLength()) return 100;
	else if (k<1||k>GetRNA2()->GetSequenceLength()) return 101;



	//OK Error checking came up clean:

	//Allocate the space for forcing alignments if not already done.
	if (forcealign==NULL) AllocateForceAlign();

	forcealign[0][i]=k;
	forcealign[1][k]=i;


	return 0;

}

//Get an alignment constraint.
int Dynalign_object::GetForcedAlignment(const int i, const int seq) {

	//check the indexes
	if (!(seq==1||seq==2)) return 0;
	if (seq==1) if (i<1||i>GetRNA1()->GetSequenceLength()) return 0;
	if (seq==2) if (i<1||i>GetRNA2()->GetSequenceLength()) return 0;

	//Make sure some constraints have been defined.
	if (forcealign==NULL) return 0;

	//OK, passed error checks.
	return forcealign[seq-1][i];


}



//Read alignment constraints from disk.
int Dynalign_object::ReadAlignmentConstraints(const char filename[]) {

	//Make sure the file exists or return an error
	if ((fopen(filename, "r"))== NULL) {
		//the file is not found
		return 102;
		
	}

	


	//Allocate space in forcealign, if needed
	if (forcealign==NULL) AllocateForceAlign();

	
	
	readalignmentconstraints(filename, forcealign, GetRNA1()->GetStructure(),  GetRNA2()->GetStructure()); 

	return 0;

	//Note that space is reserved in the error list for an error reading the file (103), but the back end does not yet support this.

}

//Read a ct file to determine what pairs will be allowed for sequence 1 in a subsequent dynalign calculation.
//return An integer that indicates an error code (0=no error, 104=file not found, 105=template is already specified)
int Dynalign_object::Templatefromct(const char ctfilename[]) {

	//Make sure the file exists or return an error
	if ((fopen(ctfilename, "r"))== NULL) {
		//the file is not found
		return 104;
		
	}
	
	if (templatefilename!=NULL) return 105;

	storetemplatefilename(ctfilename);

	ct_templated = true;

	return 0;
}


		
//Read a dynalign save file to determine what pairs will be allowed for sequence 1 in a subsequent dynalign calculation.
//return An integer that indicates an error code (0=no error, 106=file not found, 105=template is already specified)
int Dynalign_object::Templatefromdsv(const char dsvfilename[], const float maxdsvchange) {

	//Make sure the file exists or return an error
	if ((fopen(dsvfilename, "r"))== NULL) {
		//the file is not found
		return 106;
		
	}

	if (templatefilename!=NULL) return 105;

	storetemplatefilename(dsvfilename);

	dsv_templated = true;

	MAXDSV = maxdsvchange;

	return 0;
}


//Report the best energy for pair i-j from sequence #sequence.
double Dynalign_object::GetBestPairEnergy(const int sequence,const int a, const int b) {
	int i,j,k,l;


	if (!savefileread) {
		//There is an error if the savefile was not read by the constructor.

		ErrorCodeTwo=107;
		return 0.0;

	}

	//check the indexes
	if (sequence==1) {
		if (a<1||a>GetRNA1()->GetSequenceLength()) {
			ErrorCodeTwo=108;
			return 0.0;
		}
		else if (b<1||b>GetRNA1()->GetSequenceLength()) {
			ErrorCodeTwo=108;
			return 0.0;
		}
	}
	else {//sequence==2
		if (a<1||a>GetRNA2()->GetSequenceLength()) {
			ErrorCodeTwo=108;
			return 0.0;
		}
		else if (b<1||b>GetRNA2()->GetSequenceLength()) {
			ErrorCodeTwo=108;
			return 0.0;
		}

	}


	//Past error checking
	ErrorCodeTwo=0;


	if (array==NULL) {
		//This is the first call of this function.  The array needs to be filled.

		array = new double **[2];
		array[0] = new double *[GetRNA1()->GetStructure()->GetSequenceLength()+1];
		for (i=0;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			array[0][i] = new double [i+1];
		}
		for (i=0;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) {
				array[0][i][j]=INFINITE_ENERGY;
			}
		}
		array[1] = new double *[GetRNA2()->GetStructure()->GetSequenceLength()+1];
		for (i=0;i<=GetRNA2()->GetStructure()->GetSequenceLength();i++) {
			array[1][i] = new double [i+1];
		}
		for (i=0;i<=GetRNA2()->GetStructure()->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) {
				array[1][i][j]=INFINITE_ENERGY;
			}
		}


		////
		for (i=1;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			for (j=i+minloop;j<=GetRNA1()->GetStructure()->GetSequenceLength();j++) {
				for (k=max(lowend[i],1);k<=(min(GetRNA2()->GetStructure()->GetSequenceLength(),highend[i]));k++) {
					for (l=max(lowend[j],k);l<=(min(GetRNA2()->GetStructure()->GetSequenceLength(),highend[j]));l++) {
						
						//use a and b when refering to the energy arrays
						//a = k-i+maxsep;
						//b = l-j+maxsep;	

						//if ((a+ct2->numofbases-ct1->numofbases>=0)&&(a+ct2->numofbases-ct1->numofbases<=2*maxsep)) {

							
					
						array[0][j][i]=min(array[0][j][i],((double) (v->f(i,j,k,l)+v->f(j,i+GetRNA1()->GetStructure()->GetSequenceLength(),l,k+GetRNA2()->GetStructure()->GetSequenceLength())))/conversionfactor);
						
						//allow single BP inserts
						if (i>1&&j<GetRNA1()->GetStructure()->GetSequenceLength()&&k>1&&l<GetRNA1()->GetStructure()->GetSequenceLength()&&k+1>=lowend[i+1]&&k+1<=highend[i+1]&&k<=highend[i-1]&&k>=lowend[i-1]
							&&l>=lowend[j+1]&&l<=highend[j+1]&&l-1<=highend[j-1]&&l-1>=lowend[j-1]) {
							array[0][j][i]=min(array[0][j][i],((double) (v->f(i+1,j-1,k+1,l-1)+v->f(j+1,(i-1)+GetRNA1()->GetStructure()->GetSequenceLength(),l,k+GetRNA2()->GetStructure()->GetSequenceLength())+
								erg1(i-1,j+1,i,j,GetRNA1()->GetStructure(),data)+erg1(i,j,i+1,j-1,GetRNA1()->GetStructure(),data)+erg1(k,l,k+1,l-1,GetRNA2()->GetStructure(),data)+2*gap))/conversionfactor);	
						}
						



					
						array[1][l][k]=min(array[1][l][k],((double) (v->f(i,j,k,l)+v->f(j,i+GetRNA1()->GetStructure()->GetSequenceLength(),l,k+GetRNA2()->GetStructure()->GetSequenceLength())))/conversionfactor);

						//allow single BP inserts
						if (i>1&&j<GetRNA1()->GetStructure()->GetSequenceLength()&&k>1&&l<GetRNA2()->GetStructure()->GetSequenceLength()&&k+1>=lowend[i+1]&&k+1<=highend[i+1]&&k-1<=highend[i]&&k-1>=lowend[i]
							&&l+1>=lowend[j]&&l+1<=highend[j]&&l-1<=highend[j-1]&&l-1>=lowend[j-1]) {
							array[1][l][k]=min(array[1][l][k],(double) ((v->f(i+1,j-1,k+1,l-1)+v->f(j,i+GetRNA1()->GetStructure()->GetSequenceLength(),l+1,k-1+GetRNA2()->GetStructure()->GetSequenceLength())+
								erg1(k-1,l+1,k,l,GetRNA2()->GetStructure(),data)+erg1(k,l,k+1,l-1,GetRNA2()->GetStructure(),data)+erg1(i,j,i+1,j-1,GetRNA1()->GetStructure(),data)+2*gap))/conversionfactor);	
						}

						

					}
				}
			}

		}


	}

	return array[sequence-1][b][a];




}


//Return the lowest free energy change from a previous Dynalign calculation
double Dynalign_object::GetLowestEnergy() {


	if (!savefileread) {
		//There is an error if the savefile was not read by the constructor.

		ErrorCodeTwo=107;
		return 0.0;

	}

	//convert from an integer representation to a double with conversionfactor (defined in defines.h)
	return (((double) lowest)/conversionfactor);

}

//Return error messages based on code from GetErrorCode and other error codes.		

//		0 = no error
//		100-999 = Error associated with Dynalign, to be handled here.
//		>=1000 = Errors for underlying sequence, get message from TwoRNA base class.
const char* Dynalign_object::GetErrorMessage(const int error) {

	if (error > 999) return TwoRNA::GetErrorMessage(error);
	else if (error == 100) return "Nucleotide from sequence 1 is out of range.\n";
	else if (error == 101) return "Nucleotide from sequence 2 is out of range.\n";
	else if (error == 102) return "Alignment constraint file not found.\n";
	else if (error == 103) return "Error reading alignment constraint file.\n";
	else if (error == 104) return "CT file not found.\n";
	else if (error == 105) return "A template has already been specified; only one is allowed.\n";
	else if (error == 106) return "DSV file not found.\n";
	else if (error == 107) return "Data not available to calculate energy.\n";
	else if (error == 108) return "Nucleotide out of range.\n";
	else if (error == 109) return "Value of maxpairs is too large to be achievable.\n";
	else if (error == 110) return "Error reading thermodynamic parameters.\nPlease set environment variable DATAPATH to the location of the thermodynamic parameters.\n";
	else if (error == 111) return "DSV file is the wrong version.";
	else if (error==0) return "No Error.\n";
	else return "Unknown Error.\n";

}

//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void Dynalign_object::SetProgress(ProgressHandler& Progress) {

	//use the TProgressDialog pointer from RNA1
	GetRNA1()->SetProgress(Progress);

	return;

}


//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void Dynalign_object::StopProgress() {

	//use the TProgressDialog pointer from RNA1
	GetRNA1()->StopProgress();
	return;

}

//Destructor
Dynalign_object::~Dynalign_object() {
	int i;


	//clean up memory allocation
	if (align!=NULL) {
		for (i=0;i<Maxtrace;i++)  delete[] align[i];
		delete[] align;
	}

	//clean up forcedalignments if necessary:
	if (forcealign != NULL) {
		delete[] forcealign[0];
		delete[] forcealign[1];
		delete[] forcealign;
	}

	//clean up filename if needed:
	if (templatefilename!=NULL) delete[] templatefilename;

	if (savefileread) {
		//If the savefile was read, do the emory cleanup
		if (modificationflag) delete vmod;
		delete v;
		delete w;
		delete w3;
		delete w5;
		delete[] lowend;
		delete[] highend;
		delete data;
	}

	//delete array if need be
	if (array!=NULL) {

		
		for (i=0;i<=GetRNA1()->GetStructure()->GetSequenceLength();i++) {
			delete[] array[0][i];
		}
		delete[] array[0];
		
		
		
		for (i=0;i<=GetRNA2()->GetStructure()->GetSequenceLength();i++) {
			delete[] array[1][i];
		}
		delete[] array[1];
		delete[] array;

	}
}


void Dynalign_object::CommonConstructor() {

	//By default, the alignment information is not allocated
	align=NULL;


	//By default, there are no alignment constraints
	forcealign = NULL;

	//By default, dsv templating is off
	dsv_templated = false; 

	//By default, ct templatimg is off
	ct_templated = false;

	//By default, no filename is needed
	templatefilename = NULL;

	//By default, a dynalign save file was not read
	savefileread=false;

	//By default, pair energies are not needed
	array = NULL;


}

//Allocate space in the forcealign array for storing alignment constraints.
void Dynalign_object::AllocateForceAlign() {
	int index;
	
	forcealign=new short *[2];

		
	forcealign[0]=new short [GetRNA1()->GetStructure()->GetSequenceLength()+1];
	forcealign[1]=new short [GetRNA2()->GetStructure()->GetSequenceLength()+1];
	for (index=1;index<=GetRNA1()->GetStructure()->GetSequenceLength();index++) {
		forcealign[0][index]=0;
    }
	for (index=1;index<=GetRNA1()->GetStructure()->GetSequenceLength();index++) {
		forcealign[1][index]=0;
    }

}

void Dynalign_object::storetemplatefilename(const char *name) {

	templatefilename = new char [strlen(name)+1];
	strcpy(templatefilename,name);
}
