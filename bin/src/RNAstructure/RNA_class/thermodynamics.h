//Class Thermodynamics -- Wrapper for RNAstructure functions that open and store thermodynamic parameters
//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(THERMODYNAMICS_H)
#define THERMODYNAMICS_H


//Import the definition of struct datatable.
#include "../src/rna_library.h"


//TOLERANCE is the maximum deviation from 310.15 K before which the enthalpy parameters are read from disk to adjust the free energyy changes from 310.15 K.
#define TOLERANCE 0.01 


//! Thermodynamics Class.
/*!
	The RNA class provides an encapsulation of the functions and struct for reading and storing thermodynamic parameters.
	This includes methods for changing folding temperatures
	This class is intended for use in inheritance for classes that provide functionality. 
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.

class Thermodynamics {
	
	public:
		//**********************************************
		//Constructors:
		//**********************************************
		Thermodynamics(const bool isRNA=true, const char* const alphabetName=NULL, const double temperature=310.15);
		Thermodynamics(const Thermodynamics& copyThermo);  // copy constructor. Use this in place of Thermodynamics::CopyThermodynamic

		//**********************************************
		//Functions that pertain to folding temperature:
		//**********************************************

		//!Set the temperature of folding in K.

		//!This function allows the user to specify folding temperatures other than 310.15 K (37 degrees C).
		//!This changes folding free energy changes that would be returned for existing structures and would alter the set of structures predicted.
		//!When this function is called, the thermodynamic parameter files are immediately read from disk.  These include both enthalpy
		//!parameters (.dh files) and free energy changes at 310.15 (.dat files).  The files must either be at a location indicated by the $DATAPATH
		//!environment variable or in pwd.  Changing the temperature only alters subsequent calculations.  For example, if a structure prediction
		//!method has been called, the set of predicted structures are not changed at the time of a call to SetTemperature.  Likewise, SetTemperature
		//!must be called before calling a structure prediction method if a temperature other than the 310.15 K default is desired.
		//!The function returns an error code where 0 is no error and non-zero errors can be parsed by by GetErrorMessage() or GetErrorMessageString() in an inheriting class.
		//\param temperature is double that indicates the folding temperature in K.
		//\return An integer that indicates an error code for reading the thermodynamic parameters.  0 = no error.  5 = files not found.
		int SetTemperature(double temperature);

		//!Get the current folding temperature in K.
		//!\return A double that indicates the folding temperature in K.
		double GetTemperature() const;

		//! Get the name of the extended alphabet for which thermodynamic parameters should be loaded.
		string GetAlphabetName() const;

		//**********************************************
		//Functions that read parameters from disk:
		//**********************************************

		//! Function to read the thermodynamic parameters.

		//! This function depends on temp, the current temperature, to determine in the folding free energies need to be set to other than those read
		//!	in files for 310.15 K.
		//!	Return of zero => no error and a return of non-zero indicates error.
		//! Public functions that need the thermodynamic parameters call this function automatically. 
		//! By default, the path to the thermodynamic paramaters is fetched from the $DATAPATH environment variable.
		//! If a specific path is needed, $DATAPATH is overridden by specifying the pathname explicitly here as a parameter.
		//! \return An int that indicates whether an error occured.
		//! \param directory is a pointer to a cstring that indicates the directory that contains the thermodynamnic parameter files.  By default, this is NULL and the environment variable $DATAPATH is consulted to get this path. 
		//! \param alphabet is a pointer to a cstring that indicates the name of the nucleic acid alphabet (e.g. "rna", "dna" etc);
		//! \param temperature specify a custom temperature. If this parameter is less than zero it is ignored. Otherwise the free energy parameters are recalculated to apply to the specified temperature.
		int ReadThermodynamic(const char *directory=NULL, const char *alphabet=NULL, const double temperature=-1.0);

		//! Re-loads the free energy data tables from the original directory and 
		//!    using the original alphabet name.
		//! \param temperature specify a custom temperature. If this parameter is less than zero it is ignored. 
		//!        Otherwise the free energy parameters are recalculated to apply to the specified temperature.
		//! 	   scaled to the specified temperature.
		int ReloadDataTables(const double new_temperature=-1.0);

		//! Force the datatables to be read if they haven't already. Return true if the tables were already loaded or if the attempt to (re)open them succeded. 
		bool VerifyThermodynamic(); 

		//**********************************************
		//Functions that provide accessibility of the underlying tables:
		//**********************************************


		//! This function is used during inheritance to provide access to the free energy change parameters.
		//! This function generates no error codes.  (Error checking was done for this during construction).
		//!\return A pointer to datatable with free energy change parameters.
		datatable *GetDatatable();

		

		//! This function is used to provide an enthalpy table.
		//! This function will return a NULL pointer if there is an error reading the tables from disk.
		//! It is important that programs check the status of the pointer before using it, i.e. make sure it is not NULL.
		//! \param alphabet Specify which alphabet to use. If this is NULL (the default) the currently loaded alphabet, if there is one, will be used. Otherwise either default RNA or DNA alphabet will be used, depending on the value of isrna.
		//! \return A pointer to datatable with the enthalpy change parameters.
		datatable *GetEnthalpyTable(const char* alphabet=NULL);

		//! Clear the currently loaded energy datatable and release its resources.
		void ClearEnergies();
		//! Clear the currently loaded enthalpy datatable and release its resources.
		void ClearEnthalpies();

		//!Return whether this instance of Thermodynamics has the paremters populated (either from disk or from another Thermodynamics class).
		
		//!\return A bool that indicates whether the parameters are populated (true = yes).
		bool GetEnergyRead() const;
		//!\return A bool that indicates whether the alphabet was populated.
		bool IsAlphabetRead() const;

		//**********************************************
		//Destructor:
		//**********************************************
		~Thermodynamics();	
	
	
		//Keep track of whether the sequence is RNA or DNA.
		bool isrna;

	protected:

		//!Copy thermodynamic parameters from an instance of Thermodynamics class.

		//!This is generally not needed because functions automatically populate
		//!the parameters from disk.  It is helpful, however, when a large number
		//!of calculations with be performed because the parameters can then
		//!be read from disk only once.
		//!Note that the source Thermodynamics class must have been initialized with the "correct" ISRNA value and correct temperature.
		//!Return 0 if no error and non-zero errors can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!\param thermo is a pointer to Thermodynamics class.  That must have already called the ReadThermodynamics() function.  
		//! 
		//! IMPORTANT: CopyThermo should not be used, except possibly in constructors of derived classes.
		//!    Use the Thermodyanamic(Thermodynamic* copy) constructor instead.
		//!    Reason: Since the extended alphabet, the datatable is loaded in the RNA class constructor, 
		//!    so copying the datatable must also happen in the constructor to avoid reading a 
		//!    (possibly different) datatable first.
		//! This replaces Thermodynamics::CopyThermodynamic
		virtual void CopyThermo(const Thermodynamics& copy);

		////Get the set of file names for reading parameters.
		////The c strings will all receive file names that start with the path name, datapath.
		////isRNA indicates whether the DNA or RNA parameters are being read (true = RNA, false = DNA).
		////isEnthlpy indicates whether the enthalpy parameters are being read (true = enthalpy, false = free energy at 37 degrees C).
		
		//Keep track of whether the parameter files were read.
		// bool energyread;

		//Class to store thermodynamic parameters.
		datatable *data;

		//Class to store enthalpy parameters.
		datatable *enthalpy;

		bool copied; // whether the datatable was copied from another Thermodynamics class. (in which case, this class should NOT delete it. The owner should.)

		//The desired folding temperature in K (310.15 by default)
		// Note that the ACTUAL data table may have a DIFFERENT temperature than this if the 
		// data is shared among one or more another Thermodynamics objects and another 
		// one has modified the data.
		// Therfore, if the data tables have been loaded, the "TRUE" temperature should be obtained from this->GetTemperature() (which internally queries the data object)
		double nominal_temperature;

		// The name of the desired extended alphabet to open.
		// Note that the ACTUAL data table may have a DIFFERENT alphabet than this if the 
		// data is shared among one or more another Thermodynamics objects and another 
		// one has loaded the data with a different alphabet.
		// Therfore, if the data tables have been loaded, the "TRUE" alphabet should be obtained from this->GetAlphabetName() (which internally queries the data object)
		string nominal_alphabetName;


		bool skipThermoTables; // if true, only the alphabet will be loaded. The energy tables are skipped. This should ONLY be set to true if the program makes no use of energies (e.g. ct2dot etc)
};

#endif
