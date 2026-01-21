#include <cstdlib>
#include <cstring>
#include "thermodynamics.h"

// Constructor:
Thermodynamics::Thermodynamics(const bool isRNA, const char* const alphabetName, const double temperature) {
	isrna=isRNA; 		  //store the backbone type (this is superceded by alphabetName, but some code may still use it.)
	enthalpy=NULL; 		  //set the enthalpy parameters to an unread status
	data=NULL; 
	nominal_temperature=temperature;  //set the folding temperature
	nominal_alphabetName=as_str(alphabetName); // the name of the alphabet opened by a call to opendat  (This is stored in case the datatable needs to be reopened after a call to SetTemperature)
	skipThermoTables=false;
	copied=false; 		// whether the energy data has been copied from another Thermodynamics class instance (and therefore should NOT be deleted by THIS class).
}

// copy Constructor
Thermodynamics::Thermodynamics(const Thermodynamics& copy) {
	// copyThermo calls ClearEnergies and ClearEnthalpies, so `data` and `enthalpies` must be set to NULL before calling it.
	data = enthalpy = NULL; copied = false;
	CopyThermo(copy);
}

// This replaces ReadThermodynamic.
// IMPORTANT: Thermodynamics::CopyThermodynamic should not be used, except possibly in constructors of derived classes.
//    Use the Thermodyanamic(Thermodynamic* copy) constructor instead.
//    Reason: Since the extended alphabet, the datatable is loaded in the RNA class constructor, 
//    so copying the datatable must also happen in the constructor to avoid reading a 
//    (possibly different) datatable first.
void Thermodynamics::CopyThermo(const Thermodynamics& copy) {
	ClearEnergies();
	ClearEnthalpies();
	isrna=copy.isrna; //store the backbone type (this is superceded by alphabetName, but some code may still use it.)
	enthalpy=copy.enthalpy; //set the enthalpy parameters 
	data=copy.data;
	nominal_temperature=copy.GetTemperature();	//set the folding temperature
	nominal_alphabetName=copy.GetAlphabetName(); // the name of the alphabet opened by a call to opendat  (This is stored in case the datatable needs to be reopened after a call to SetTemperature)
	skipThermoTables=copy.skipThermoTables;
	copied=true; // 'true' indicates the datatable was copied. i.e. this does NOT own the datatable (and therefore should NOT delete it).
}

// Returns true if the two temperatures are within +/- TOLERANCE of each other. (e.g. +/- 0.01)
inline bool equal_within_tolerance(double temp1, double temp2) {
	return abs(temp1-temp2)<TOLERANCE;
}

//Set the folding temperature:
//Return an error code pertaining to reading the thermodynamic parameters:
int Thermodynamics::SetTemperature(double temperature) {
	// If the thermodynamic parameter files were already read at some point, re-read them from disk.
	// (Skip reading if the temperature is the same as the already-loaded values)
	// Note that if the data tables have not yet been loaded, we do NOT need to load them, because
	// ReadThermodynamic will be called by the client code later. In constrast, if the client has already
	// called ReadThermodynamic, we SHOULD re-read the data-tables because the client will expect that 
	// SetTemperature modifies the existing tables.
	if (GetEnergyRead()&&!equal_within_tolerance(data->GetLoadedTemperature(),temperature))
		return ReloadDataTables(temperature); // reload the data (from the original directory and using the original alphabetName)
	else
		// Either (A) The data tables have not yet been loaded. or 
		//        (B) The temperature of the data tables is already the same as the specified temperature.
		// In either case just set the nominal folding temperature.
		// If (A) then the nominal temperature will be applied the next time ReadThermodynamic is called.
		nominal_temperature = temperature;
	return 0;
}

//! Re-loads the free energy data tables from original
//!  data directory and using the original alphabet name.
//! If new_temperature is less than 0 it is ignored. Otherwise, the data tables will be read in and
//! scaled to the specified temperature.
int Thermodynamics::ReloadDataTables(const double new_temperature) {
	if (data==NULL) return 30;
	// make a copy of these datatable members because datatable may modify them in opendat.
	string dir(data->data_directory), alphabet(data->alphabetName);
	return ReadThermodynamic(dir.c_str(), alphabet.c_str(), new_temperature<0 ? data->temperature : new_temperature );
}

// deletes the thermodynamic data tables.
void Thermodynamics::ClearEnergies() {
	if (data!=NULL && !copied) delete data;
	data = NULL;
	copied = false;
}
void Thermodynamics::ClearEnthalpies() {
	if (enthalpy!=NULL) delete enthalpy;
	enthalpy = NULL;
}

//return the current folding temperature
// (see comments regarding nominal_temperature to understand why it might differ from the result of this function)
double Thermodynamics::GetTemperature() const {
	// If the datatables have already been loaded, return the TRUE temperature 
	// (i.e. the one that represents the currently loaded thermodynamic parameters.)
	// Otherwise, return the nominal temperature (which *will* be applied when
	// the datatables are read in at later time (by a call to ReadThermodynamic);
	return GetEnergyRead()?data->GetLoadedTemperature() : nominal_temperature;
}

//return the current folding temperature
// (see comments regarding nominal_alphabetName to understand why it might differ from the result of this function)
string Thermodynamics::GetAlphabetName() const {
	// If the datatables have already been loaded, return the TRUE alphabetName 
	// (i.e. the one that represents the currently loaded thermodynamic parameters.)
	// Otherwise, return the nominal alphabetName (which *will* be applied when
	// the datatables are read in at later time (by a call to ReadThermodynamic);
	return IsAlphabetRead()?data->GetAlphabetName() : nominal_alphabetName;
}

//Destructor:
Thermodynamics::~Thermodynamics() {

	//If the thermodynamic parameter files were read at some point, delete them now:
	ClearEnergies();

	//If the enthalpy parameters were read from disk, they must be deleted now:
	ClearEnthalpies();

}


datatable *Thermodynamics::GetDatatable() {
	return data;
}



bool Thermodynamics::VerifyThermodynamic() { 
	// Open the energy datatable if it hasn't been opened already (or if it was closed by SetTemperature etc)
	skipThermoTables = false; // force loading the full tables (not just the alphabet)
	return GetEnergyRead() || ReadThermodynamic()==0; // use the current alphabetName and directory.
}

//read the thermodynamic parameters from disk at location $DATAPATH or pwd
int Thermodynamics::ReadThermodynamic(const char *const directory, const char *alphabetName, const double set_temperature) {
	// Only allocate the datatable if energyread is false, meaning that no parameters are loaded
	//	This is important because the user might alter the temperature with SetTemperature(), triggering a re-read of the parameters.
	if (data==NULL) data = new datatable();

	if (!is_blank(alphabetName)) nominal_alphabetName = alphabetName;

	if (nominal_alphabetName.empty())
		// The client has not specified an alphabetname. Use the default.
		nominal_alphabetName=isrna?DT_RNA:DT_DNA; // use "rna" or "dna" depending on isrna

	if (set_temperature>=0) nominal_temperature = set_temperature; //update the nominal temperature. it will be applied below.

	// cout << sfmt("ReadThermo T=%f  dir=%s  alpha=%s", nominal_temperature, directory==NULL?"":directory, alphabetName==NULL?"":alphabetName)  << endl;

	//open the data files -- must reside in pwd or $DATAPATH. opendat will auto-detect the data-path if directory is NULL
	//open the thermodynamic data tables
	int error = 0;
	if (data->opendat(directory, nominal_alphabetName.c_str(), false, skipThermoTables)==0) { 
		error=5; // 5=error reading thermo tables
	} else {
		//now check to see if the temperature is other than 310.15 K:
		if (!equal_within_tolerance(nominal_temperature,TEMP_37C)) { 
			// causes the data table to read in the *.dh enthalpy tables 
			// and rescale the thermodynamic parameters for the specified temperature
			error=data->ScaleToTemperature(nominal_temperature); // returns an error code if the enthalpy tables can't be loaded. the code corresponds to RNA::GetErrorMessage, so we can return it directly
		}
	}
	if (error!=0) ClearEnergies(); // delete energy tables because they were not loaded correctly (or at the desired temperature).
	return error;
}


// This function is used to provide an enthalpy table.
datatable *Thermodynamics::GetEnthalpyTable(const char* alphabetName) {
	//start by determining if the parameters need to be read or whether they have already been read:
	if (enthalpy==NULL) {
		string alpha=as_str(alphabetName);
		
		//allocate a table to storte enthalpy parameters
		enthalpy = new datatable();

		// If the user did not specify an alphabet, see if one is currently loaded (i.e. in alphabetName). If so, use that. Otherwise use isrna to decide whether to use the RNA or DNA alphabet.
		if (alpha.empty())
			alpha = GetAlphabetName(); // get the currently loaded alphabet name (if data is loaded) or nominal_alphabetName otherwise.
		
		// if alphabetName is still empty, set it to "rna" or "dna" (this would be the case if data is NOT currently loaded and nominal_alphabetName=NULL)
		if (alpha.empty())
			alpha=nominal_alphabetName=isrna?DT_RNA:DT_DNA; // use "rna" or "dna" depending on isrna

		//open the enthlpy parameters and check for errors
		if (enthalpy->opendat(data->data_directory.c_str(), alpha.c_str(), true)==0)
			ClearEnthalpies(); //an error has occured. (enthalpy is set to NULL and returned).
	}
	// Either the parameters have been read or an error has occured in which case enthalpy is NULL. Return it either way.
	return enthalpy;
}


bool Thermodynamics::IsAlphabetRead() const {
	return data!=NULL && data->loadedAlphabet;
}

//Return whether this instance of Thermodynamics has the paremters populated (either from disk or from another Thermodynamics class).
bool Thermodynamics::GetEnergyRead() const {
	return data!=NULL && data->loadedTables; // Note: loadedAlphabet is always true if loadedTables is true.
}

