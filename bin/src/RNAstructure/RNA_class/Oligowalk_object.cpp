#include "Oligowalk_object.h"
#include "RNA.h"

#include "../src/OligoScreenCalc.h"
#include "../src/intermolecular.h"
#include "../src/thermo.h"

#include <cstdlib>


//Constructor - user provides a sequence as a c string.
Oligowalk_object::Oligowalk_object(const char sequence[])
	:RNA(sequence) {
		CommonConstructor();
}

//Constructor - user provides a filename for existing file as a c string.
Oligowalk_object::Oligowalk_object(const char filename[], const int type)
	:RNA(filename, type) {

		CommonConstructor();


}

//Default Constructor - user provides nothing.
Oligowalk_object::Oligowalk_object(const bool IsRNA)
	:RNA(IsRNA) {

		CommonConstructor();

}


//Perform an OligoWalk calculation
//usesub = 0 -> no suboptimal; 3 -> heuristic
int Oligowalk_object::Oligowalk(const int oligo_length, const bool isDNA, const int option, const double oligo_concentration, const int usesub, const int start, const int stop   ) {

	
	datatable *enthalpy;
	thermo *helixstack;
	int test=-1;
	string stackf,datapath = getDataPath();
	int error;
	datatable *dnadata;
	//use a character array to turn off SHAPE:
	char shapefile='\0';

	int i,j,k,l;

	// if it's a CT file, it has no substructures, so add one.
	if (GetStructure()->GetNumberofStructures()==0)
		GetStructure()->AddStructure();

	Thermodynamics *ddata;
	rddata *hybriddata,*enthalpyhybrid;

	//Make sure this is the first (and only allowed call of OligoWalk)
	if (table!=NULL) return 101;


	length = oligo_length;//save oligo_length for use in the destructor and for error checking later

	//Ensure that the thermodynamic data tables have been read.
	if (!VerifyThermodynamic()) return 5;//return non-zero if a problem occurs
    
	//Read the enthalpy data from disk using the underlying thermodynamics class GetEnthalpyData function:
	enthalpy = GetEnthalpyTable();

	//make sure that the thermodynamics parameters could be read from disk
	if (enthalpy==NULL) {
		return 5;//5 is the error code that indicates the thermodynamic parameters were not found
	}


	//Although prefiltering isn't used in this version of OligoWalk, it needs to be initialized and passed to function.
	//Note that the true indicates that empirical scores would be used a
	prefilter = new siPREFILTER(*GetDatatable(),*enthalpy,0,true,GetStructure()->GetSequenceLength() - oligo_length + 2,isDNA);


	//Allocate the tables needed for storing the results
	table = new int*[GetStructure()->GetSequenceLength() - oligo_length + 2];

	for (i = 0; i < GetStructure()->GetSequenceLength() - oligo_length + 2; i++) {
		//DHM commented out these lines for now.  They need to be restored later.
		//if (siRNA) table[i] = new int[7];
   		//else table[i] = new int[6];
		table[i] = new int[6];
	}

	//allocate memory of number of suboptimal structures
	numofsubstructures= new int*[GetStructure()->GetSequenceLength() - oligo_length +2];
	for (i = 0; i < GetStructure()->GetSequenceLength() - oligo_length + 2; i++)	{
		numofsubstructures[i]= new int [2];
		numofsubstructures[i][0]=0;
		numofsubstructures[i][1]=0;
	}

	//Allocate helixstack:
	helixstack = new thermo(datapath);

	//Now read the DNA parameters if the oligos are DNA:
	if (isDNA) {
		ddata = new Thermodynamics(false);//allocate space for DNA parameters
		
		//set the temperature for the DNA parameters 
		error = ddata->SetTemperature(GetTemperature());
		
		//Check fo an error 
		if (error!=0) {
			delete ddata;
			delete prefilter;
			delete helixstack;
			return error;
		}

		//Read the dna thermodynamic parameters
		error = ddata->ReadThermodynamic();
		
		//Check fo an error 
		if (error!=0) {
			delete ddata;
			delete prefilter;
			delete helixstack;
			return error;
		}

		//Read the hybrid data as well
		hybriddata = new rddata;

		stackf=datapath+"/stackdr.dat";

		//Check for errors
		if (readrd (hybriddata,stackf)==0) {
      		delete ddata;
			delete prefilter;
			delete helixstack;
			return 5;
		}

		if (GetTemperature()<310||GetTemperature()>311) {
		
			//The temperature is simgificantly different from 37 dgrees C, so read and use the enthalpy data.

			stackf=datapath+"/stackdr.dh";
			enthalpyhybrid = new rddata;

			//Check for errors
			if (readrd (enthalpyhybrid,stackf)==0) {
				delete ddata;
				delete prefilter;
				delete enthalpyhybrid;
				delete helixstack;
				return 5;
			}

			for (i=0;i<5;i++) {
				for (j=0;j<5;j++) {
					for (k=0;k<5;k++) {
						for (l=0;l<5;l++) {
							hybriddata->stack[i][j][k][l]=Tscale(GetTemperature(),hybriddata->stack[i][j][k][l],enthalpyhybrid->stack[i][j][k][l]);
						}
					}
				}
			}
			hybriddata->init=Tscale(GetTemperature(),hybriddata->init,enthalpyhybrid->init);
			delete enthalpyhybrid;

		}

		helixstack->DH = datapath+"/stackdr.dh"; 
		helixstack->DS = datapath+"/stackdr.ds";
		helixstack->HELIX = datapath+"/helixdr.dat";

		dnadata = ddata->GetDatatable();
	}
	else {
		dnadata = NULL;
	}

	if (helixstack->read()==0) {
		//This means an error occurred reading the helixstack parameters	
		if (isDNA) delete ddata;
		delete prefilter;
		delete helixstack;
		
		return 5;
	}

	//For now, siRNA is off.
	//if (siRNA) {
	//	//mask will store whether an oligo meets siRNA design criteria
	//	mask = new bool [ct.numofbases - length + 2];

	//}
	//else mask = NULL;

	//note that foldsize is temporarilly set to zero so that there is no folding size limit
	//note that distance is set to zero so that there is no maximum pairping distance
	//note that test is set to -1 so there is no testing
	//note that write is set to FALSE to turn off writing
	olig(isDNA, option, GetStructure(), oligo_length, oligo_concentration, usesub,
		start, stop, 0/*foldsize*/, 0/*distance*/,
		table, numofsubstructures, &shapefile, &test, false,
		*GetDatatable(), *dnadata, helixstack, hybriddata,
		prefilter,GetProgress());

	//siRNA is off right now
	//if (oligoobject->siRNA) {
	//	filterbysirna(&oligoobject->ct,oligoobject->table,oligoobject->length,&oligoobject->data,oligoobject->mask,oligoobject->asuf,oligoobject->tofe,oligoobject->fnnfe);
	//}

	//note that foldsize is temporarilly set to zero here.
	//report(oligoobject->outputfile, &oligoobject->ct, oligoobject->table, oligoobject->numofsubstructures, 
	//	oligoobject->length, oligoobject->isdna, oligoobject->c, oligoobject->usesub,oligoobject->start,oligoobject->stop,prefilter,0,
	//	oligoobject->mask,oligoobject->asuf,oligoobject->tofe,oligoobject->fnnfe);

	//Clean up if this is DNA
	if (isDNA) {
		delete ddata;
	}
	delete helixstack;
	return 0;
}

// Get the breaking target DG for a given nucleotide.
double Oligowalk_object::GetBreakTargetDG(const int index) {

	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}

	//note that the free eneergy changes are stored internally as integers, where they are multiplied by conversionfactor
	return (((double) table[index][2])/conversionfactor);



}

//Get the duplex DG for a given nucleotide.
double Oligowalk_object::GetDuplexDG(const int index) {
	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}

	//note that the free eneergy changes are stored internally as integers, where they are multiplied by conversionfactor
	return (((double) table[index][1])/conversionfactor);

}

// Get the bimolecular oligo-oligo DG for a given nucleotide.

double Oligowalk_object::GetOligoOligoDG(const int index) {
	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}

	//note that the free eneergy changes are stored internally as integers, where they are multiplied by conversionfactor
	return (((double) table[index][4])/conversionfactor);

}

//Get the oligo-self DG for a given nucleotide.
double Oligowalk_object::GetOligoSelfDG(const int index) {
	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}
	//note that the free eneergy changes are stored internally as integers, where they are multiplied by conversionfactor
	return (((double) table[index][3])/conversionfactor);

}

//Get the overall DG for a given nucleotide.
double Oligowalk_object::GetOverallDG(const int index) {
	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}

	//note that the free eneergy changes are stored internally as integers, where they are multiplied by conversionfactor
	return (((double) table[index][0])/conversionfactor);

}

// Get the Tm for a given nucleotide.
double Oligowalk_object::GetTm(const int index) {
	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		ErrorCode = 100;
		return 0.0;

	}
	else if (index<1||index>GetStructure()->GetSequenceLength() - length + 1) {
		//This means an invalid index is being used
		ErrorCode = 3;
		return 0.0;
	}

	//note that the Tm is stored as an integer, multiplied by 10.
	return (((double) table[index][5])/10);

}


//Write a report to disk using the OligoWalk data.
int Oligowalk_object::WriteReport(const char outputfilename[], const int oligo_length, const bool isDNA, const int option, const double oligo_concentration, const int usesub, const int start, const int stop) {

	//Perform error checking:
	if (table==NULL) {
		//This means that OligoWalk() was not called successfully.
		return 100;

	}

	//The last three 0's indicate this is not siRNA design
	ofstream out(outputfilename);
	report(out, isDNA, option, GetStructure(), oligo_length, oligo_concentration, usesub,
		start,stop,0/*foldsize*/,0/*distance*/,
		table, numofsubstructures, NULL/*shapefile*/, prefilter,
		NULL,0,0,0,
		false/*html*/,true/*header*/,true/*body*/);

	//No errors:
	return 0;

}


//Perform an OligoScreen calculation.
int Oligowalk_object::OligoScreen(const char infilename[], const char outfilename[]) {
	rddata *hybriddata,*enthalpyhybrid;
	char stackf[maxfil];
	int i,j,k,l;


	//Check that the inputfilename exists:
	FILE *check;

	//check that the file exists.
	if ((check = fopen(infilename, "r"))== NULL) {
		//the file is not found
		//fclose(check);
		return 1;		
	}
	fclose(check);

	//Ensure that the thermodynamic data tables have been read.
	if (!VerifyThermodynamic()) return 5;//return non-zero if a problem occurs

	//Now read the RNA-DNA hybrid parameters, if necessary:
	if (!isrna) {
		//This is DNA oligos

		//Get the information from $DATAPATH, if available
		strcpy(stackf,getDataPath()); // getDataPath will return DATAPATH if it exists or "." otherwise.
		strcat(stackf,"/stackdr.dat");

		//check that the datafile exists:
		if ((check = fopen(stackf, "r"))
			== NULL) {
			return 5;
		}

		hybriddata = new rddata;
		readrd(hybriddata,stackf);
		
		if (GetTemperature()<310||GetTemperature()>311) {
			//The temperature has changed.
			//Read the enthalpy data into a rddata.
			
			//Get the information from $DATAPATH, if available
			strcpy(stackf,getDataPath());
			strcat(stackf,"/stackdr.dh");

			//check that the datafile exists:
			if ((check = fopen(stackf, "r"))
				== NULL) {
				delete hybriddata;
				return 5;
			}
			
			enthalpyhybrid = new rddata;
			readrd(enthalpyhybrid,stackf);
      		

			for (i=0;i<5;i++) {
				for (j=0;j<5;j++) {
					for (k=0;k<5;k++) {
						for (l=0;l<5;l++) {
							hybriddata->stack[i][j][k][l]=Tscale(GetTemperature(),hybriddata->stack[i][j][k][l],enthalpyhybrid->stack[i][j][k][l]);
						}
					}
				}
			}
			delete enthalpyhybrid;
		}

	}
	else {
		//This is RNA oligos
		hybriddata= NULL;
	}

	//call the backend function
	OligoScreenCalc(infilename, outfilename, data, hybriddata);


	//cleanup:
	if (!isrna) {
		delete hybriddata;

	}

	return 0;

}


//Return a c string that describes errors from GetErrorCode and other errors. 
const char* Oligowalk_object::GetErrorMessage(const int error) {

	if (error==100) return "No OligoWalk data present.  Perform an OligoWalk calculation first using a call to OligoWalk()\n";
	if (error==101) return "OligoWalk has been performed.  Only one OligoWalk calculation can be performed.\n";

	//Return the base class's message if not recognized as an Oligowalk_object-specific error.
	else return RNA::GetErrorMessage(error);


}


//Destructor
Oligowalk_object::~Oligowalk_object() {

	int i;

	if (table!=NULL) {
		//Table and numofsubstructures were allocated, so delete them
		for (i = 0; i < GetStructure()->GetSequenceLength() - length + 2; ++i) {
			delete[] table[i];
			delete[] numofsubstructures[i];
		}
		delete[] table;
		delete[] numofsubstructures;

		delete prefilter;
	}

}

//Common code for all constructors
void Oligowalk_object::CommonConstructor() {
	
	table = NULL;
	numofsubstructures = NULL;


}


