// RNAstructureBackendCalculator class.
//
// (c) 2010 Mathews Lab, University of Rochester Medical Center.
// Written by Jessica S. Reuter.
//
// The RNAstructureBackendCalculator class handles all calculations and data
// access for the Java GUI interface.

// #include <sys/syscall.h>
// #include <unistd.h>
#include "RNAstructureBackendCalculator.h"
// DEBUG TIMING: #include "Windows.h"



////////////////////////////////////////////////////
// non-member utility functions.
////////////////////////////////////////////////////
const char* const SAVFILE_EXT = "sav";

// If ext is null, SAVFILE_EXT will be used instead.
string chageFileExt(const string ctFile, const char* ext=NULL) { 
	int found = ctFile.find_last_of( "." );
	if (ext==NULL) ext=SAVFILE_EXT;
	int len=strlen(ext);
	if (found == string::npos) 
		return ctFile + "." + ext;
	return ctFile.substr(0, found + 1) + ext;
}

/*
 * Constructor and Destructor.
 */

//////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
RNAstructureBackendCalculator::RNAstructureBackendCalculator() {
	// ErrorCode = 0;
	// errorDetails="";

	// Initialize default values for all non-forced constraints.
	constraintsHolder.maxLoop = 30;
	constraintsHolder.maxPair = -1;
	constraintsHolder.shapeEnergy = false;
	constraintsHolder.shapeFile = "";
	constraintsHolder.shapeParam1 = 0.0;
	constraintsHolder.shapeParam2 = 0.0;
	constraintsHolder.temperature = 310.15;

	// Initialize the Dynalign object and its error checker to null pointers.
	dynalign = NULL;
	dynalignChecker = NULL;

	// Initialize the HybridRNA object and its error checker to null pointers.
	hybrid = NULL;
	hybridChecker = NULL;

	// RMW - removed ProgressMonitor, as it was not necessary: monitor = 0; // Initialize the progress monitor to a null pointer.
	progress = new SimpleProgressHandler(); // Initialize the progress dialog. (set the outputStream to NULL so that progress isn't written to std out.)
	ownsProgress = true;

	// Set the Multilign and TurboFold structs to be inactive.
	turboFold.isActive = false;
	multilign.isActive = false;

	// Initialize the oligo object and its error checker to null pointers.
	oligo = NULL;
	oligoChecker = NULL;

	// Initialize the single strand RNA object and its error checker to null
	// pointers.
	rna = NULL;
	rnaChecker = NULL;

	// Initialize the sequence handler.
	for( int i = 1; i <= 3; i++ ) { sequenceHandler.push_back( "" ); }
}

//////////////////////////////////////////////////////////////////////////////
// Destructor.
//////////////////////////////////////////////////////////////////////////////
RNAstructureBackendCalculator::~RNAstructureBackendCalculator() {

	// If the Dynalign object is not null, delete it.
	// If the Dynalign error checker is not null, delete it.
	if( dynalign != NULL ) { delete dynalign; }
	if( dynalignChecker != NULL ) { delete dynalignChecker; }

	// If the HybridRNA object is not null, delete it.
	// If the HybridRNA error checker is not null, delete it.
	if( hybrid != NULL ) { delete hybrid; }
	if( hybridChecker != NULL ) { delete hybridChecker; }

	// If the progress monitor is not null, delete it.
	// RMW - removed ProgressMonitor : if( monitor != 0 ) { delete monitor; }
	if( progress != NULL && ownsProgress ) { delete progress; }

	// If the oligo object is not null, delete it.
	// If the oligo error checker is not null, delete it.
	if( oligo != NULL ) { delete oligo; }
	if( oligoChecker != NULL ) { delete oligoChecker; }

	// If the single strand object is not null, delete it.
	// If the single strand error checker is not null, delete it.
	if( rna != NULL ) { delete rna; }
	if( rnaChecker != NULL ) { delete rnaChecker; }
}

// Creates a new progress dialog, deleting any existing one.
void RNAstructureBackendCalculator::resetProgress() {
	//delete progress;
	//progress = new TProgressDialog();
	progress->reset();
}

/*
 * AllSub Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for AllSub calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildAllSubDataStructure(string file, bool isRNA ) { 
	return loadSingleRNA(file, isRNA);
}

//////////////////////////////////////////////////////////////////////////////
// Get the absolute energy difference.
//////////////////////////////////////////////////////////////////////////////
double RNAstructureBackendCalculator::getSuboptimalAbsoluteDiff() {

	// Get the sequence length.
	int length = rna->GetSequenceLength();

	// Return an absolute energy difference based on the length.
	if( length > 1200 ) { return 0.25; }
	else if( length > 800 ) { return 0.5; }
	else if( length > 500 ) { return 0.75; }
	else if( length > 300 ) { return 1; }
	else if( length > 120 ) { return 1.5; }
	else if( length > 50 ) { return 3; }
	else { return 10; }
}

//////////////////////////////////////////////////////////////////////////////
// Get the maximum percent energy difference.
//////////////////////////////////////////////////////////////////////////////
float RNAstructureBackendCalculator::getSuboptimalPercentDiff() {

	// Get the sequence length.
	int length = rna->GetSequenceLength();

	// Return a percent energy difference based on the length.
	if( length > 1200 ) { return 5; }
	else if( length > 800 ) { return 8; }
	else if( length > 500 ) { return 10; }
	else if( length > 300 ) { return 15; }
	else if( length > 120 ) { return 20; }
	else if( length > 50 ) { return 25; }
	else { return 50; }
}

//////////////////////////////////////////////////////////////////////////////
// Run AllSub Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runAllSub(
	string ctFile, float percent, double absolute ) {

	// Set the calculation temperature.
	rna->SetTemperature( constraintsHolder.temperature );

	// Add SHAPE data, if necessary. If an error occurred, return it.
	if( constraintsHolder.shapeFile != "" ) {
		int shapeErr = rna->ReadSHAPE(
			constraintsHolder.shapeFile.c_str(),
			constraintsHolder.shapeParam1, constraintsHolder.shapeParam2, RESTRAINT_SHAPE,
			constraintsHolder.shapeEnergy );
		if( shapeErr != 0 ) { return rnaChecker->returnError( shapeErr ); }
	}

	// Create a variable to handle error codes.
	int error = 0;

	// Create the progress monitor.
	// RMW - removed ProgressMonitor, as it was not necessary:// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	rna->SetProgress( *progress );
	// DEBUG TIMING: DWORD startTime = GetTickCount();
	// Run suboptimal structure generation.
	error = rna->GenerateAllSuboptimalStructures( percent, absolute );

	// Delete the progress monitor.
	rna->StopProgress();
	//delete dialog;

	// TODO: Add a TProgressDialog to rna->WriteCt because it can take a HUGE amount of time for AllSub! The progress monitor could be updated every 10 structures or so.
	// DEBUG TIMING: cerr << "before WriteCt" << GetTickCount()-startTime << endl;
	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = rna->WriteCt( ctFile.c_str() ); }

	// DEBUG TIMING: cerr << "done with WriteCt" << GetTickCount()-startTime << endl;;

	// Return the error string.
	return rnaChecker->returnError( error );
}



/*
 * Bifold Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for bifold calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildBifoldDataStructure(string file1, string file2, bool isRNA ) {
	return loadDoubleRNA(file1, file2, isRNA, FILE_SEQ);
}

//////////////////////////////////////////////////////////////////////////////
// Run Bifold Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runBifold(
	string ctFile, float percent, int maxStructures, int windowSize,
	bool saveFile, bool intraForbidden ) {

	// Set the calculation temperature.
	hybrid->SetTemperature( constraintsHolder.temperature );

	// Create a variable to handle error codes.
	int error = 0;

	// Determine the save file name.
	string saveFileName = saveFile ? generateSavFile(ctFile) : "";

	// Set whether intramolecular pairs are allowed or not.
	hybrid->SetForbidIntramolecular( intraForbidden );

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	hybrid->SetProgress( *progress );

	// Run bifold.
	error = hybrid->FoldBimolecular(
		percent, maxStructures, windowSize, saveFileName.c_str(),
		constraintsHolder.maxLoop );

	// Delete the progress monitor.
	hybrid->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = hybrid->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return hybridChecker->returnError( error );
}



/*
 * Bipartition Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for bipartition calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildBipartitionDataStructure(string file1, string file2, bool isRNA ) {
	return loadDoubleRNA(file1, file2, isRNA, FILE_SEQ);
}

//////////////////////////////////////////////////////////////////////////////
// Run Bipartition Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runBipartition( string pfsFile ) {

	// Set the calculation temperature.
	hybrid->SetTemperature( constraintsHolder.temperature );

	// Create the progress monitor.
	resetProgress(); // RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	// TProgressDialog* progress = new TProgressDialog( *(monitor) );
	hybrid->SetProgress( *progress );

	// Run the bipartition calculation.
	int mainError = hybrid->PartitionFunctionBimolecular( pfsFile.c_str() );

	// Delete the progress monitor.
	hybrid->StopProgress();
	//delete progress;

	// Return any error.
	return hybridChecker->returnError( mainError );
}


//////////////////////////////////////////////////////////////////////////////
// Dynalign Module
//////////////////////////////////////////////////////////////////////////////
// Create the Dynalign_object
string RNAstructureBackendCalculator::buildDynalignDataStructure(
	string file1, string file2, bool isRNA ) {
	if( dynalign != 0 ) { delete dynalign; }
	if( dynalignChecker != 0 ) { delete dynalignChecker; }
	dynalign = new Dynalign_object( file1.c_str(), 2, file2.c_str(), 2, isRNA );
	dynalignChecker = new ErrorChecker<Dynalign_object>( dynalign );
	return dynalignChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Clear Dynalign alignment constraints.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::clearDynalignAlignmentConstraints() {
	dynalignAlignments.clear();
}

//////////////////////////////////////////////////////////////////////////////
// Get the Dynalign alignment constraints as a string.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getDynalignAlignmentConstraints() {
	stringstream stream( stringstream::in | stringstream::out );
	stream << "<html>Alignment Constraints<br/><br/>";

	map<int, int>::iterator it = dynalignAlignments.begin();
	while( it != dynalignAlignments.end() ) {
		stream << "(1) " << it->first << " -- (2) " << it->second << "<br>";
		it++;
	}

	return stream.str();
}

//////////////////////////////////////////////////////////////////////////////
// Get the Dynalign alignment window size.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getDynalignAlignmentWindowSize() {

	// Get the sequence length.
	int length = dynalign->GetRNA1()->GetSequenceLength();

	// Return a window size based on the length.
	if( length > 1200 ) { return 20; }
	else if( length > 800 ) { return 15; }
	else if( length > 500 ) { return 11; }
	else if( length > 300 ) { return 7; }
	else if( length > 120 ) { return 5; }
	else if( length > 50 ) { return 3; }
	else { return 2; }
}

//////////////////////////////////////////////////////////////////////////////
// Get the Dynalign structure window size.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getDynalignStructureWindowSize() {

	// Get the sequence length.
	int length = dynalign->GetRNA1()->GetSequenceLength();

	// Return a window size based on the length.
	if( length > 500 ) { return 3; }
	else if( length > 300 ) { return 2; }
	else if( length > 50 ) { return 1; }
	else { return 0; }
}

//////////////////////////////////////////////////////////////////////////////
// Read a Dynalign alignment constraints file.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::readDynalignAlignmentConstraintsFile(
	string file ) {

	string result = "";
	ifstream fileStream( file.c_str() );
	string line;
	stringstream lineStream( stringstream::in | stringstream::out );

	int len1 = dynalign->GetRNA1()->GetSequenceLength();
	int len2 = dynalign->GetRNA2()->GetSequenceLength();

	while( fileStream.good() ) {
		getline( fileStream, line );
		lineStream.str( "" );
		lineStream << line;

		int first, second;
		lineStream >> first;
		lineStream >> second;

		bool nuc1OK =
			( first >= 1 ) &&
			( first <= len1 );
		bool nuc2OK =
			( second >= 1 ) &&
			( second <= len2 );

		if( !( nuc1OK && nuc2OK ) ) {
			result = "Error reading Dynalign alignment constraints file.";
		} else {
			dynalignAlignments[first] = second;
		}
	}
	fileStream.close();

	return result;
}

//////////////////////////////////////////////////////////////////////////////
// Run Dynalign calculations.
//////////////////////////////////////////////////////////////////////////////
#ifdef DYNALIGN_II
string RNAstructureBackendCalculator::runDynalign(
	string ctFile1, string ctFile2, string saveFile, string alignFile,
	double percent, int structures, int windowStr, int windowAli, float gap,
	float slope, float intercept, int max_elongation ) {
#else
 string RNAstructureBackendCalculator::runDynalign(
        string ctFile1, string ctFile2, string saveFile, string alignFile,
        double percent, int structures, int windowStr, int windowAli, float gap,
        bool isInsert ) {
#endif
	// Set the temperature for the Dynalign calculation.
	dynalign->SetTemperature( constraintsHolder.temperature );

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
    resetProgress();  //RMW Removed ProgressMonitor:  TProgressDialog* progress = new TProgressDialog( *(monitor) );
	dynalign->SetProgress( *progress );
	
	// Run the Dynalign calculation.
#ifdef DYNALIGN_II
	int mainError = dynalign->Dynalign(structures, windowStr, windowAli, 
		percent, (short)-99, slope, intercept, gap, max_elongation,
		saveFile.c_str() );
#else
	int mainError = dynalign->Dynalign(structures, windowStr, windowAli, 
        percent, (short)-99, gap, isInsert,
		saveFile.c_str() );
#endif

	// Delete the progress monitor.
	dynalign->StopProgress();
	//  delete progress;

	// If an error occurred in the main calculation, return it.
	if( mainError != 0 ) { return dynalignChecker->returnError( mainError ); }

	// Write the CT output files. If an error occurs with either of these
	// writings, return it.
	int ctErr = dynalign->GetRNA1()->WriteCt( ctFile1.c_str() );
	if( ctErr != 0 ) { return dynalignChecker->returnError( ctErr ); }
	ctErr = dynalign->GetRNA2()->WriteCt( ctFile2.c_str() );
	if( ctErr != 0 ) { return dynalignChecker->returnError( ctErr ); }

	// Write the alignment file.
	dynalign->WriteAlignment( alignFile.c_str() );

	// Return the empty string.
	return "";
}

//////////////////////////////////////////////////////////////////////////////
// Set a Dynalign alignment constraint.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setDynalignAlignmentConstraint(
	int nuc1, int nuc2 ) {

	// Check to make sure the first strand nucleotide is in range; if not,
	// return an error.
	int length1 = dynalign->GetRNA1()->GetSequenceLength();
	if( ( nuc1 < 1 ) || ( nuc1 > length1 ) ) {
		return "Nucleotide in sequence 1 out of range.";
	}

	// Check to make sure the second strand nucleotide is in range; if not,
	// return an error.
	int length2 = dynalign->GetRNA2()->GetSequenceLength();
	if( ( nuc2 < 1 ) || ( nuc2 > length2 ) ) {
		return "Nucleotide in sequence 2 out of range.";
	}

	// Set the nucleotides aligned and return the empty string.
	dynalignAlignments[nuc1] = nuc2;
	return "";
}

//////////////////////////////////////////////////////////////////////////////
// Write a Dynalign alignment constraints file.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::writeDynalignAlignmentConstraintsFile(
	string file ) {

	ofstream out( file.c_str() );
	map<int, int>::iterator it = dynalignAlignments.begin();
	while( it != dynalignAlignments.end() ) {
		out << it->first << "\t" << it->second << endl;
		it++;
	}
	out << "-1\t-1" << endl;
	out.close();
}



/*
 * Efn2 Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for efn2 calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildEfn2DataStructure(string file, bool isRNA ) {
	return loadSingleRNA(file, isRNA, FILE_CT);
}

//////////////////////////////////////////////////////////////////////////////
// Run Efn2 Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runEfn2(
	string outFile, bool writeDetails ) {

	// Create a variable to handle errors.
	string error = "";

	// Set the calculation temperature.
	rna->SetTemperature( constraintsHolder.temperature );

	// Calculate the structure type. If the type isn't radial, set the
	// structure type result as an error and return it.
	string typeResult = getStructureType();
	if( typeResult == "Circular" ) {
		stringstream errorStream( stringstream::in | stringstream::out );
		errorStream << "Structures contain pseudoknots, so Efn2 cannot run.";
		error = errorStream.str();
	} else if( typeResult != "Radial" ) { error = typeResult; }
	if( error != "" ) { return error; }

	// If a thermodynamic details file should be written, write that file.
	if( writeDetails ) {

		// Write the thermodynamic details file, and return completion result.
		int code = rna->WriteThermodynamicDetails( outFile.c_str() );
		return rnaChecker->returnError( code );
	}

	// If a simple list file should be written, write that file.
	else {

		// Create a vector to hold energies.
		vector<double> energies;

		// For each structure, calculate energy.
		int structures = rna->GetStructureNumber();
		for( int i = 1; i <= structures; i++ ) {

			// Calculate the free energy and check for errors.
			energies.push_back( rna->CalculateFreeEnergy( i ) );
			if( ( error = rnaChecker->returnError() ) != "" ) {
				i += structures;
			}
		}

		// If all energies were calculated correctly, continue.
		if( error == "" ) {

			// Open a stream to the energy list file.
			ofstream out( outFile.c_str() );

			// For each structure, write the energy to the file.
			for( int i = 1; i <= structures; i++ ) {
				out << "Structure: " << i << "   Energy = " << fixed
				    << setprecision( 1 ) << energies.at( i - 1 ) << endl;
			}

			// Close the energy list file.
			out.close();
		}
	}

	// Return the error string.
	return error;
}



/*
 * Fold Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for Fold calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildFoldDataStructure(string file, bool isRNA ) { return loadSingleRNA(file, isRNA); }

//////////////////////////////////////////////////////////////////////////////
// Get the folding window size.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getFoldWindowSize() {

	// Get the sequence length.
	int length = rna->GetSequenceLength();

	// Return a window size based on the length.
	if( length > 1200 ) { return 20; }
	else if( length > 800 ) { return 15; }
	else if( length > 500 ) { return 11; }
	else if( length > 300 ) { return 7; }
	else if( length > 120 ) { return 5; }
	else if( length > 50 ) { return 3; }
	else { return 2; }
}

//////////////////////////////////////////////////////////////////////////////
// Run Fold Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runFold(string ctFile, int percent, int structures, int window, bool saveFile) {

	// Set the calculation temperature.
	rna->SetTemperature( constraintsHolder.temperature );

	// Set the maximum pairing distance, if applicable.
	if( constraintsHolder.maxPair != -1 ) {
		rna->ForceMaximumPairingDistance( constraintsHolder.maxPair );
	}

	// Determine the save file name.
	string saveFileName = saveFile ? generateSavFile(ctFile) : "";

	// Add SHAPE data, if necessary. If an error occurred, return it.
	if( constraintsHolder.shapeFile != "" ) {
		int shapeErr = rna->ReadSHAPE(
			constraintsHolder.shapeFile.c_str(),
			constraintsHolder.shapeParam1, constraintsHolder.shapeParam2, RESTRAINT_SHAPE,
			constraintsHolder.shapeEnergy );
		if( shapeErr != 0 ) { return rnaChecker->returnError( shapeErr ); }
	}

	// Create a variable to handle error codes.
	int error = 0;

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	rna->SetProgress( *progress );

	// Run single strand folding.
	error = rna->FoldSingleStrand(
		percent, structures, window, saveFileName.c_str(),
		constraintsHolder.maxLoop );

	// Delete the progress monitor.
	rna->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = rna->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return rnaChecker->returnError( error );
}



/*
 * MaxExpect Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for MaxExpect calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildMaxExpectDataStructure(string file ) {
	return loadSingleRNA(file, true, FILE_PFS);
}

//////////////////////////////////////////////////////////////////////////////
// Run MaxExpect Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runMaxExpect(
	string ctFile, double score, int structures, int window, double gamma ) {

	// Create a variable to handle error codes.
	int error = 0;

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	rna->SetProgress( *progress );

	// Run maximum expected accuracy.
	error = rna->MaximizeExpectedAccuracy( score, structures, window, gamma );

	// Delete the progress monitor.
	rna->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = rna->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return rnaChecker->returnError( error );
}



/*
 * Multilign Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Set the Multilign module active.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::activateMultilign() {

	multilign.isActive = true;
}

//////////////////////////////////////////////////////////////////////////////
// Add a tuple to the Multilign file array.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::addMultilignTuple(
	string seqFile, string ctFile ) {

	vector<string> row;
	row.push_back( seqFile );
	row.push_back( ctFile );
	row.push_back( "" );
	row.push_back( "" );
	multilign.files.push_back( row );
}

//////////////////////////////////////////////////////////////////////////////
// Delete a tubple from the Multilign file array.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::deleteMultilignTuple(
	unsigned int index ) {

	multilign.files.erase( multilign.files.begin() + ( index - 1 ) );
}

//////////////////////////////////////////////////////////////////////////////
// Get the CT structure file name for a particular Multilign sequence.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getMultilignCT( int index ) {

	return multilign.files[index-1][1];
}

//////////////////////////////////////////////////////////////////////////////
// Get the max pairs in the underlying Multilign object.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getMultilignMaxPairs() {

	// If less than two sequences have been input, return 0.
	if( multilign.files.size() < 2 ) { return 0; }

	// Otherwise, return the number of pairs calculated.
	else {

		// Create the Multilign object and check it for errors
		// If an error occurred, return 0.
		Multilign_object* object = new Multilign_object( multilign.files );
		ErrorChecker<Multilign_object>* checker =
			new ErrorChecker<Multilign_object>( object );
		if( checker->returnError() != "" ) { return 0; }

		// Set and get the max pairs.
		object->SetMaxPairs( -1 );
		int pairs = object->GetMaxPairs();

		// Delete the Multilign object and error checker.
		delete object;
		delete checker;

		// Return the max pairs.
		return pairs;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Get the Multilign sequence set as a formatted string.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getMultilignSequenceSetData() {

	// Initialize the sequence set stream.
	stringstream seqSet( stringstream::in | stringstream::out );

	// Put each sequence set pair in the stream.
	int size = multilign.files.size();
	for( int i = 1; i <= size; i++ ) {
		seqSet << i << ".\t" << multilign.files[i-1][0] << endl
		       << "\t" << multilign.files[i-1][1];
		if( i != size ) { seqSet << endl; }
	}

	// Return the sequence set string.
	return seqSet.str();
}

//////////////////////////////////////////////////////////////////////////////
// Get the number of sequences used in Multilign.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getNumMultilignSequences() {

	return multilign.files.size();
}

//////////////////////////////////////////////////////////////////////////////
// Run Multilign calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runMultilign(
	int percent, int maxStructures, int bpWindow, int alignWindow, double gap,
	bool insert, double maxDsvChange, int maxPairs, int cycles,
	string alignFile, bool saveFiles, bool isRNA ) {

	// Create a variable to handle error codes.
	string errorString = "";

	// Create the Multilign object and check it for errors.
	Multilign_object* object = new Multilign_object( multilign.files, isRNA );
	ErrorChecker<Multilign_object>* checker =
		new ErrorChecker<Multilign_object>( object );
	errorString = checker->returnError();
	if( errorString != "" ) { return errorString; }

	// If no error occurred, set the maximum dsv change, the maxpairs, the
	// number of cycles, and the temperature.
	if( errorString == "" ) {
		object->SetMaxDsv( maxDsvChange );
		object->SetMaxPairs( maxPairs );
		object->SetIterations( cycles );
		object->SetTemperature( constraintsHolder.temperature );
	}

	// If no error occurred, run the main calculation.
	if( errorString == "" ) {
		// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
		resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
		object->SetProgress( progress );

		// Run a Multilign calculation.
		int mainError = object->ProgressiveMultilign(
			1, true, true, maxStructures, bpWindow, alignWindow, percent,
			(short)-99, gap, insert, (short)30, false );

		// Delete the progress monitor.
		object->StopProgress();
		// delete dialog;

		// Check the main calculation for errors.
		errorString = checker->returnError( mainError );
	}

	// If no error occurred, write the alignment file.
	if( errorString == "" ) {
		int alignError = object->WriteAlignment( alignFile );
		errorString = checker->returnError( alignError );
	}

	// If intermediate files should not be saved, clean up the intermediate
	// files generated by Multilign.
	if( !saveFiles ) { object->CleanupIntermediateFiles(); }

	// Delete the Multilign object and error checker.
	delete object;
	delete checker;

	// Return the calculation result.
	return errorString;
}



/*
 * OligoScreen Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for OligoScreen calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildOligoScreenDataStructure() {
	if( oligo != 0 ) { delete oligo; }
	if( oligoChecker != 0 ) { delete oligoChecker; }
	oligo = new Oligowalk_object();
	oligoChecker = new ErrorChecker<Oligowalk_object>( oligo );
	return oligoChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Run OligoScreen Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runOligoScreen(
	string in, string out, bool isRNA ) {

	// Create a variable that handles errors.
	string error = "";

	// Set the oligomer chemistry and the calculation temperature.
	oligo->isrna = isRNA;
	oligo->SetTemperature( constraintsHolder.temperature );

	// Run OligoScreen and check for errors.
	int oligoError = oligo->OligoScreen( in.c_str(), out.c_str() );
	error = oligoChecker->returnError( oligoError );

	// Return the error string.
	return error;
}



/*
 * OligoWalk Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for OligoWalk calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildOligoWalkDataStructure(
	string file, bool isRNA ) {
	if( oligo != 0 ) { delete oligo; }
	if( oligoChecker != 0 ) { delete oligoChecker; }
	
	string ext = toLower(getFileExt(file));
	RNAInputType type = (ext=="seq"||ext=="fasta"||ext=="fa") ? FILE_SEQ : FILE_CT;
	oligo = new Oligowalk_object( file.c_str(), type);
	oligo->isrna = isRNA;
	oligoChecker = new ErrorChecker<Oligowalk_object>( oligo );
	return oligoChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Get whether an oligo can be folded bimolecular.
//////////////////////////////////////////////////////////////////////////////
bool RNAstructureBackendCalculator::canFoldOligoOligo( int index ) {

	return ( oligo->GetOligoOligoDG( index ) != 0.0 );
}

//////////////////////////////////////////////////////////////////////////////
// Get whether an oligo can be folded unimolecular.
//////////////////////////////////////////////////////////////////////////////
bool RNAstructureBackendCalculator::canFoldOligoSelf( int index ) {

	return ( oligo->GetOligoSelfDG( index ) != 0.0 );
}

//////////////////////////////////////////////////////////////////////////////
// Determine the maximum number of oligos allowed.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::determineOligoMaximum( int length ) {

	return getOligoTargetLength() - ( length - 1 );
}

//////////////////////////////////////////////////////////////////////////////
// Fold an oligo.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::foldOligo(
	string sequence, int index, bool bimolecular, bool isRNA, string file ) {

	// Flip the sequence.
	string flip( sequence );
	reverse( flip.begin(), flip.end() );

	// If the fold should be unimolecular, use an RNA object to do the fold.
	if( bimolecular == false ) {
		RNA* seq = new RNA( flip.c_str(), isRNA );
		seq->FoldSingleStrand( 100, 20, 0 );
		seq->WriteCt( file.c_str() );
		delete seq;
	}

	// Otherwise, use a HybridRNA object to do the fold.
	else {
		HybridRNA* seq = new HybridRNA( flip.c_str(), flip.c_str(), isRNA );
		seq->FoldBimolecular( 100, 20, 0 );
		seq->WriteCt( file.c_str() );
		delete seq;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Get the data for all oligos.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getAllOligoData( int height ) {

	// Create the string stream that holds all oligo data.
	// Also, get the bounds of the target region.
	stringstream data( stringstream::in | stringstream::out );
	int start = getGraphRegionBegin();
	int end = getGraphRegionEnd();

	// Add the data for the overall/duplex overlaid bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {

			// Get the next overall value and set it as the minimum or maximum,
			// if applicable.
			double nextOverall = oligo->GetOverallDG( i );
			if( nextOverall < minimum ) { minimum = nextOverall; }
			if( nextOverall > maximum ) { maximum = nextOverall; }

			// Get the next duplex value and set it as the minimum or maximum,
			// if applicable.
			double nextDuplex = oligo->GetDuplexDG( i );
			if( nextDuplex < minimum ) { minimum = nextDuplex; }
			if( nextDuplex > maximum ) { minimum = nextDuplex; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the overall and duplex values.
			double overall = oligo->GetOverallDG( i );
			double duplex = oligo->GetDuplexDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( overall > 0 ) { upper = ( overall / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( overall < 0 ) { middle = ( overall / minimum ) * minSpan; }
			else { middle = ( duplex / minimum ) * minSpan; }

			// Calculate the lower bar height.
			double lower = 0;
			if( overall < 0 ) { lower = ( duplex / minimum ) * minSpan; }

			// If the middle and lower bars are stacked, adjust the lower bar
			// length to allow them to stack properly.
			if( ( middle > 0 ) && ( lower > 0 ) ) { lower -= middle; }

			// Add the bar heights to the data stream.
			data << upper << "," << middle << "," << lower << " ";
		}

		// Close off the section.
		data << ";";
	}

	// Add the data for the overall bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {
			double next = oligo->GetOverallDG( i );
			if( next < minimum ) { minimum = next; }
			if( next > maximum ) { maximum = next; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the next value.
			double value = oligo->GetOverallDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( value > 0 ) { upper = ( value / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( value < 0 ) { middle = ( value / minimum ) * minSpan; }

			// Add the bar heights to the data stream.
			// The lowest bar is always set to 0 because it's never used.
			data << upper << "," << middle << ",0 ";
		}

		// Close off the section.
		data << ";";
	}

	// Add the data for the duplex bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {
			double next = oligo->GetDuplexDG( i );
			if( next < minimum ) { minimum = next; }
			if( next > maximum ) { maximum = next; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the next value.
			double value = oligo->GetDuplexDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( value > 0 ) { upper = ( value / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( value < 0 ) { middle = ( value / minimum ) * minSpan; }

			// Add the bar heights to the data stream.
			// The lowest bar is always set to 0 because it's never used.
			data << upper << "," << middle << ",0 ";
		}

		// Close off the section.
		data << ";";
	}

	// Add the data for the broken target bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {
			double next = oligo->GetBreakTargetDG( i );
			if( next < minimum ) { minimum = next; }
			if( next > maximum ) { maximum = next; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the next value.
			double value = oligo->GetBreakTargetDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( value > 0 ) { upper = ( value / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( value < 0 ) { middle = ( value / minimum ) * minSpan; }

			// Add the bar heights to the data stream.
			// The lowest bar is always set to 0 because it's never used.
			data << upper << "," << middle << ",0 ";
		}

		// Close off the section.
		data << ";";
	}

	// Add the data for the unimolecular bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {
			double next = oligo->GetOligoSelfDG( i );
			if( next < minimum ) { minimum = next; }
			if( next > maximum ) { maximum = next; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the next value.
			double value = oligo->GetOligoSelfDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( value > 0 ) { upper = ( value / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( value < 0 ) { middle = ( value / minimum ) * minSpan; }

			// Add the bar heights to the data stream.
			// The lowest bar is always set to 0 because it's never used.
			data << upper << "," << middle << ",0 ";
		}

		// Close off the section.
		data << ";";
	}

	// Add the data for the bimolecular bars.
	{

		// Determine the minimum and maximum values these bars will reach.
		// Then add those values to the data stream.
		double minimum = numeric_limits<double>::infinity();
		double maximum = minimum * -1;
		for( int i = start; i <= end; i++ ) {
			double next = oligo->GetOligoOligoDG( i );
			if( next < minimum ) { minimum = next; }
			if( next > maximum ) { maximum = next; }
		}
		if( maximum < 0 ) { maximum = 0; }
		data << minimum << " " << maximum << " ";

		// Calculate the range covered by the bars.
		double totalSpan = maximum - minimum;
		double minSpan = height * ( abs( minimum ) / totalSpan );
		double maxSpan = height * ( abs( maximum ) / totalSpan );

		// Determine the sizes of the bars.
		for( int i = start; i <= end; i++ ) {

			// Get the next value.
			double value = oligo->GetOligoOligoDG( i );

			// Calculate the upper bar height.
			double upper = 0;
			if( value > 0 ) { upper = ( value / maximum ) * maxSpan; }

			// Calculate the middle bar height.
			double middle = 0;
			if( value < 0 ) { middle = ( value / minimum ) * minSpan; }

			// Add the bar heights to the data stream.
			// The lowest bar is always set to 0 because it's never used.
			data << upper << "," << middle << ",0 ";
		}

		// Close off the section.
		data << ";";
	}

	// Return the data string.
	return data.str();
}

//////////////////////////////////////////////////////////////////////////////
// Get the oligo for a particular index.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getDisplayedOligo( int index ) {

	// Initialize the oligo string and determine the oligo chemistry.
	string oligoString = "";
	bool isRNA = ( oligoDataHolder.oligoType == "RNA Oligo" );

	// Add any whitespace necessary before the oligo.
	for( int i = 1; i < index; i++ ) { oligoString += " "; }

	// Build the oligo as a complement to the target with the proper length.
	oligoString += "3";
	for( int i = 1; i <= oligoDataHolder.oligoLength; i++ ) {
		char nuc = oligo->GetNucleotide( index + ( i - 1 ) );
		if( ( nuc == 'A' ) && isRNA ) { oligoString += "U"; }
		else if( ( nuc == 'A' ) && !isRNA ) { oligoString += "T"; }
		else if( nuc == 'C' ) { oligoString += "G"; }
		else if( nuc == 'G' ) { oligoString += "C"; }
		else if( nuc == 'T' ) { oligoString += "A"; }
		else if( nuc == 'U' ) { oligoString += "A"; }
		else if( nuc == 'X' ) { oligoString += "X"; }
	}
	oligoString += "5";

	// Return the oligo string.
	return oligoString;
}

//////////////////////////////////////////////////////////////////////////////
// Get the index where the graph region starts.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getGraphRegionBegin() {

	return oligoDataHolder.graphBegin;
}

//////////////////////////////////////////////////////////////////////////////
// Get the index where the graph region ends.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getGraphRegionEnd() {

	return oligoDataHolder.graphEnd;
}

//////////////////////////////////////////////////////////////////////////////
// Get the most stable oligo.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getMostStableOligo() {

	return oligoDataHolder.mostStable;
}

//////////////////////////////////////////////////////////////////////////////
// Get all data for a particular oligo.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getOligoLabelData( int index ) {

	// Create the string stream that holds oligo data.
	stringstream stream( stringstream::in | stringstream::out );

	// Build the data for a particular oligo in the stream.
	stream << "Oligo Number: " << index << ";"
	       << "<html>Overall" << oligoDataHolder.deltaG
	       << oligo->GetOverallDG( index ) << ";"
	       << "<html>Break Target" << oligoDataHolder.deltaG
	       << oligo->GetBreakTargetDG( index ) << ";"
	       << oligoDataHolder.concentration << ";"
	       << "<html>Duplex" << oligoDataHolder.deltaG
	       << oligo->GetDuplexDG( index ) << ";"
	       << "<html>Oligo-Self" << oligoDataHolder.deltaG
	       << oligo->GetOligoSelfDG( index ) << ";"
	       << oligoDataHolder.oligoType << ";"
	       << "<html>T<sub>m</sub> = " << oligo->GetTm( index ) << ";"
	       << "<html>Oligo-Oligo" << oligoDataHolder.deltaG
	       << oligo->GetOligoOligoDG( index ) << ";";

	// Return the data string.
	return stream.str();
}

//////////////////////////////////////////////////////////////////////////////
// Get the oligo target length.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getOligoTargetLength() {

	return oligo->GetSequenceLength();
}

//////////////////////////////////////////////////////////////////////////////
// Get the oligo target sequence.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getOligoTargetSequence() {

	string sequence = "<html>&nbsp;";
	for( int i = 1; i <= getOligoTargetLength(); i++ ) {
		char nuc = oligo->GetNucleotide( i );
		if( oligo->GetPair( i ) == 0 ) { sequence += nuc; }
		else {
			sequence += "<span style=\"color:red\">";
			sequence += nuc;
			sequence += "</span>";
		}
	}
	sequence += "&nbsp;";
	return sequence;
}

//////////////////////////////////////////////////////////////////////////////
// Run OligoWalk calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runOligoWalk(
	string report, int mode, string chemistry, bool suboptimal, int length,
	int amount, string unit, int start, int stop ) {

	// Do any necessary variable conversions.
	bool isDNA = ( chemistry == "DNA" );
	double concentration =
		( unit == "mM" ) ? (double)amount * 0.001 :
		( unit == "uM" ) ? (double)amount * 0.000001 :
		( unit == "nM" ) ? (double)amount * 0.000000001 :
		( unit == "pM" ) ? (double)amount * 0.000000000001 :
		0.0;
	int suboptimalOption = ( suboptimal ) ? 3 : 0;

	// Create a variable that monitors errors.
	int error = 0;

	// Set the calculation temperature.
	oligo->SetTemperature( constraintsHolder.temperature );

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	oligo->SetProgress( *progress );

	// Run OligoWalk.
	error = oligo->Oligowalk(
		length, isDNA, mode, concentration, suboptimalOption, start, stop );

	// Delete the progress monitor.
	oligo->StopProgress();
	// delete dialog;

	// If there was an error running OligoWalk, return the error.
	if( error != 0 ) { return oligoChecker->returnError( error ); }

	// Write the report file, and if an error occurred, return.
	error = oligo->WriteReport(
		report.c_str(), length, isDNA, mode, concentration,
		suboptimalOption, start, stop );

	// If there was an error writing the report, return the error.
	if( error != 0 ) { return oligoChecker->returnError( error ); }

	// Create the string that holds the delta G symbol and save it in the
	// oligo data handler.
	stringstream temperature( stringstream::in | stringstream::out );
	temperature << fixed << setprecision( 2 )
	            << ( constraintsHolder.temperature - 273.15 );
	oligoDataHolder.deltaG = " \u0394G\u00B0<sub>" + temperature.str() + "</sub> = ";

	// Save the oligo concentration in the oligo data handler.
	stringstream concStream( stringstream::in | stringstream::out );
	concStream << amount << " " << unit;
	oligoDataHolder.concentration = concStream.str();

	// Save the oligo type and oligo length in the oligo data handler.
	oligoDataHolder.oligoType = chemistry + " Oligo";
	oligoDataHolder.oligoLength = length;

	// Save the binding region bounds in the oligo data handler.
	oligoDataHolder.graphBegin = start;
	oligoDataHolder.graphEnd = stop;

	// Determine the most stable oligo and save it in the oligo data handler.
	double minValue = numeric_limits<double>::infinity();
	for( int i = start; i <= stop; i++ ) {
		double value = oligo->GetOverallDG( i );
		if( value < minValue ) {
			minValue = value;
			oligoDataHolder.mostStable = i;
		}
	}

	// Return an empty string to show that the calculation finished correctly.
	return "";
}



/*
 * Partition Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for partition calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildPartitionDataStructure(string file, bool isRNA ) {
	return loadSingleRNA(file, isRNA);
}

//////////////////////////////////////////////////////////////////////////////
// Run Partition Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runPartition( string pfsFile ) {

	// Set the calculation temperature.
	rna->SetTemperature( constraintsHolder.temperature );

	// Set the maximum pairing distance, if applicable.
	if( constraintsHolder.maxPair != -1 ) {
		rna->ForceMaximumPairingDistance( constraintsHolder.maxPair );
	}

	// Add SHAPE data, if necessary. If an error occurred, return it.
	if( constraintsHolder.shapeFile != "" ) {
		int shapeErr = rna->ReadSHAPE(
			constraintsHolder.shapeFile.c_str(),
			constraintsHolder.shapeParam1, constraintsHolder.shapeParam2,RESTRAINT_SHAPE,
			constraintsHolder.shapeEnergy );
		if( shapeErr != 0 ) { return rnaChecker->returnError( shapeErr ); }
	}

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	rna->SetProgress( *progress );

	// Do the partition function calculation.
	int partError = rna->PartitionFunction( pfsFile.c_str() );

	// Delete the progress monitor.
	rna->StopProgress();
	// delete dialog;

	// Return any error that occurred.
	return rnaChecker->returnError( partError );
}



/*
 * ProbKnot Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for ProbKnot calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildProbKnotDataStructure(string file ) {
	return loadSingleRNA(file, true, FILE_PFS);
}

//////////////////////////////////////////////////////////////////////////////
// Run ProbKnot Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runPseudoknotPrediction(
	string ctFile, int iterations, int helix ) {

	// Create a variable that handles errors.
	string error = "";

	// Run pseudoknot predictioon and check for errors.
	int pseudoError = rna->ProbKnot( iterations, helix );
	error = rnaChecker->returnError( pseudoError );

	// If no error occurred, write a CT file of pseudoknotted structures, then
	// check for errors after writing.
	if( error == "" ) {
		int writeError = rna->WriteCt( ctFile.c_str() );
		error = rnaChecker->returnError( writeError );
	}

	// Return the error string.
	return error;
}

/*
 * Refold Dynalign Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for Dynalign refolding calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildRefoldDynalignDataStructure(
	string file ) {
	if( dynalign != 0 ) { delete dynalign; }
	if( dynalignChecker != 0 ) { delete dynalignChecker; }
	dynalign = new Dynalign_object( file.c_str() );
	dynalignChecker = new ErrorChecker<Dynalign_object>( dynalign );
	return dynalignChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Run Dynalign Refolding Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runDynalignRefold(
	string ctFile1, string ctFile2, string alignFile, int percent,
	int structures, int windowStr, int windowAli ) {

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	dynalign->SetProgress( *progress );

	// Run the Dynalign refolding calculation.
	int mainError =
		dynalign->Dynalign( structures, windowStr, windowAli, percent );

	// Delete the progress monitor.
	dynalign->StopProgress();
	// delete dialog;

	// Return an error if the Dynalign calculation didn't finish correctly.
	if( mainError != 0 ) { return dynalignChecker->returnError( mainError ); }

	// Write the output files, and return an error if one occurs.
	int ctError1 = dynalign->GetRNA1()->WriteCt( ctFile1.c_str() );
	if( ctError1 != 0 ) { return dynalignChecker->returnError( ctError1 ); }
	int ctError2 = dynalign->GetRNA2()->WriteCt( ctFile2.c_str() );
	if( ctError2 != 0 ) { return dynalignChecker->returnError( ctError2 ); }
	dynalign->WriteAlignment( alignFile.c_str() );

	// Return the empty string if no error occurred.
	return "";
}



/*
 * Refold Single Structure Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for single structure refolding calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildRefoldSingleDataStructure(string file ) {
	return loadSingleRNA(file, true, FILE_SAV);
}

//////////////////////////////////////////////////////////////////////////////
// Get the refolding window size.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getRefoldWindowSize() {

	// Get the sequence length.
	int length = rna->GetSequenceLength();

	// Return a window size based on the length.
	if( length > 1200 ) { return 20; }
	else if( length > 800 ) { return 15; }
	else if( length > 500 ) { return 11; }
	else if( length > 300 ) { return 7; }
	else if( length > 120 ) { return 5; }
	else if( length > 50 ) { return 3; }
	else { return 2; }
}

//////////////////////////////////////////////////////////////////////////////
// Run Refold Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runRefold(
	string ctFile, int percent, int structures, int window ) {

	// Create a variable to handle error codes.
	int error = 0;

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	rna->SetProgress( *progress );

	// Run refolding.
	error = rna->ReFoldSingleStrand( percent, structures, window );

	// Delete the progress monitor.
	rna->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = rna->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return rnaChecker->returnError( error );
}



/*
 * Remove Pseudoknots Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for pseudoknot removal.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildRemovePseudoknotsDataStructure(
	string file, bool isRNA ) {
	return loadSingleRNA(file, isRNA, FILE_CT);
}

//////////////////////////////////////////////////////////////////////////////
// Run Pseudoknot Removal Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runPseudoknotRemoval(
	string ctFile, bool minimize ) {

	// Set the calculation temperature.
	rna->SetTemperature( constraintsHolder.temperature );

	// Create a variable that handles errors.
	string error = "";

	// Run pseudoknot removal and check for errors.
	int removeError = rna->BreakPseudoknot( minimize );
	error = rnaChecker->returnError( removeError );

	// If no error occurred, write a CT file of pseudoknot-free structures,
	// then check for errors after writing.
	if( error == "" ) {
		int writeError = rna->WriteCt( ctFile.c_str() );
		error = rnaChecker->returnError( writeError );
	}

	// Return the error string.
	return error;
}



/*
 * Sequence Display Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Get the sequence comments.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getSequenceComment() {

	return sequenceHandler[1];
}

//////////////////////////////////////////////////////////////////////////////
// Get the sequence data.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getSequenceData() {

	return sequenceHandler[2];
}

//////////////////////////////////////////////////////////////////////////////
// Get the sequence title.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getSequenceTitle() {

	return sequenceHandler[0];
}

//////////////////////////////////////////////////////////////////////////////
// Read the sequence data, depending on the type.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::readSequenceData( string file ) {

	// Open the file and initialize the line variable.
	ifstream in( file.c_str() );
	string line = "";

	// Read the first line of the file, then use that line to determine what
	// type of file it is.
	getline( in, line );
	char firstChar = line[0];
	string type =
		( firstChar == '>' ) ? "fasta" :
		( firstChar == 'L' ) ? "genbank" :
		( firstChar == ';' ) ? "sequence" :
		"text";

	// If the type is a FASTA file, read the FASTA file.
	if( type == "fasta" ) {

		// Clip off the header line character, then set the header in the
		// sequence data object.
		line = line.erase( 0, 1 );
		sequenceHandler[0] = line;

		// For all the other lines in the file, put them in the data vector.
		while( !in.eof() ) {
			getline( in, line );
			sequenceHandler[2] += line;
		}
	}

	// If the type is a Genbank file, read the Genbank file.
	else if( type == "genbank" ) {

		// Until the "DEFINITION" line of the file, skip all lines, and when
		// the "DEFINITION" line is found, save it as the title. Then, until
		// the "ACCESSION" line is reached, append all lines to the title.
		// Once this is done, save the constructed title string.
		while( line.substr( 0, 10 ) != "DEFINITION" ) { getline( in, line ); }
		string title = line;
		while( line.substr( 0, 9 ) != "ACCESSION" ) {
			title += line;
			getline( in, line );
		}
		sequenceHandler[0] = title;

		// Until the "VERSION" line of the file, skip all lines, and when the
		// "VERSION" line is found, save it as the comment.
		while( line.substr( 0, 7 ) != "VERSION" ) { getline( in, line ); }
		sequenceHandler[1] = line.substr( 7 );

		// Until the "ORIGIN" line of the file, skip all lines, and then save
		// all the lines directly after that as the sequence data.
		while( line.substr( 0, 6 ) != "ORIGIN" ) { getline( in, line ); }
		getline( in, line );
		while( line.find_first_of( "/" ) == string::npos ) {
			int noSpace = line.find_first_not_of( " " );
			sequenceHandler[2] += line.substr( noSpace );
		}
	}

	// If the type is a sequence file, read the sequence file.
	else if( type == "sequence" ) {

		// Read in all the comment lines.
		while( line.substr( 0, 1 ) == ";" ) {
			sequenceHandler[1] += line.substr( 1 );
			getline( in, line );
		}

		// Save the sequence title.
		sequenceHandler[0] = line;
		getline( in, line );
		// Read in the sequence.
		while( !in.eof() ) {
			
			if (line.length()>0) {
				char last = line.at( line.length() - 1 );
				if( last == '1' ) { line = line.substr( 0, line.length() - 1 ); }
			}
			sequenceHandler[2] += ( line + "\n" );
			getline( in, line );
		}
	}

	// If the type is a text file, read the text file.
	else {

		// Simply read the file as sequence data.
		sequenceHandler[2] += ( line + "\n" );
		while( !in.eof() ) {
			getline( in, line );
			sequenceHandler[2] += ( line + "\n" );
		}

	}

	// Close the file.
	in.close();
}

//////////////////////////////////////////////////////////////////////////////
// Set the sequence comments.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSequenceComment( string comments ) {

	sequenceHandler[1] = comments;
}

//////////////////////////////////////////////////////////////////////////////
// Set the sequence data.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSequenceData( string data ) {

	sequenceHandler[2] = data;
}

//////////////////////////////////////////////////////////////////////////////
// Set the sequence title.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSequenceTitle( string title ) {

	sequenceHandler[0] = title;
}

//////////////////////////////////////////////////////////////////////////////
// Write a FASTA file.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::writeFastaFile( string file ) {

	// Open the file stream.
	ofstream out( file.c_str() );

	// Write the title and description together as the FASTA header line,
	// then write the sequence into the file.
	out << ">" + sequenceHandler[0] + " " + sequenceHandler[1] << endl
	    << sequenceHandler[2] << endl;

	// Close the file and return the result string.
	out.close();
	return "FASTA file written successfully.";
}

//////////////////////////////////////////////////////////////////////////////
// Write a sequence file.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::writeSequenceFile( string file ) {

	// Open the file stream.
	ofstream out( file.c_str() );

	// Replace all newlines ("\n") in the comments with newlines followed by
	// comment line markers ("\n; ").
	string comments = "; " + sequenceHandler[1];
	for( unsigned int i = 1; i <= comments.size(); i++ ) {
		if( comments[i-1] == '\n' ) {
			comments.replace( i - 1, 1, "\n; " );
			i += 3;
		}
	}

	// Write the comments into the file, followed by the title line, the
	// sequence itself, and the file ending line.
	out << comments << "\n";

	if (sequenceHandler[0]!="") {
	    out << sequenceHandler[0] << "\n";
	}
	else {
		//If the title is an empty string, write a space in its place.  This fixes the problem that a title of some sort is required in a .seq.
		out << " " << "\n";
	}
	out  << sequenceHandler[2] << endl
	    << "1" << "\n";

	// Close the file and return the result string.
	out.close();
	return "Sequence file written successfully.";
}



/*
 * Stochastic Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for stochastic calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildStochasticDataStructure(
	string file ) {
	return loadSingleRNA(file, true, FILE_PFS);
}

//////////////////////////////////////////////////////////////////////////////
// Run Stochastic Sampling Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runStochastic(
	string outFile, int ensemble, int seed ) {

	// Create a variable that handles errors.
	string error = "";

	// Run stochastic sampling and check for errors.
	int stochasticError = rna->Stochastic( ensemble, seed );
	error = rnaChecker->returnError( stochasticError );

	// If no error occurred, write a CT file of stochastic sampled structures,
	// then check for errors after writing.
	if( error == "" ) {
		int writeError = rna->WriteCt( outFile.c_str() );
		error = rnaChecker->returnError( writeError );
	}

	// Return the error string.
	return error;
}



/*
 * TurboFold Module.
 */

//////////////////////////////////////////////////////////////////////////////
// Set the TurboFold module active.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::activateTurboFold() {

	turboFold.isActive = true;
}

//////////////////////////////////////////////////////////////////////////////
// Add a tuple to the TurboFold file array.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::addTurboFoldTuple(
	string seqFile, string ctFile, string pfsFile ) {

	// Create the partition function save file name.
	
	if (pfsFile.empty()) 
		pfsFile = ctFile.substr( 0, ctFile.find_last_of( "." ) ) + ".pfs";

	// Add the proper data.
	vector<string> row;
	row.push_back( seqFile );
	row.push_back( ctFile );
	row.push_back( pfsFile );
	turboFold.files.push_back( row );
}

//////////////////////////////////////////////////////////////////////////////
// Delete a tuple from the TurboFold file array.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::deleteTurboFoldTuple(
	unsigned int index ) {

	turboFold.files.erase( turboFold.files.begin() + ( index - 1 ) );
}

//////////////////////////////////////////////////////////////////////////////
// Get the number of sequences used in TurboFold.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getNumTurboFoldSequences() {

	return turboFold.files.size();
}

//////////////////////////////////////////////////////////////////////////////
// Get the CT file associated with a particular TurboFold sequence.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getTurboFoldCT( int index ) {

	return turboFold.files[index-1][1];
}

//////////////////////////////////////////////////////////////////////////////
// Get the partition function save file associated with a particular TurboFold
// sequence.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getTurboFoldSaveFile( int index ) {

	return turboFold.files[index-1][2];
}

//////////////////////////////////////////////////////////////////////////////
// Get the TurboFold sequence set as a formatted string.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getTurboFoldSequenceSetData() {

	// Initialize the sequence set stream.
	stringstream seqSet( stringstream::in | stringstream::out );

	// Put each sequence set pair in the stream.
	int size = turboFold.files.size();
	for( int i = 1; i <= size; i++ ) {
		seqSet << i << ".\t" << turboFold.files[i-1][0] << endl
		       << "\t" << turboFold.files[i-1][1];
		if( i != size ) { seqSet << endl; }
	}

	// Return the sequence set string.
	return seqSet.str();
}

//////////////////////////////////////////////////////////////////////////////
// Run TurboFold calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runTurboFold(
	double gammaT, int turboIter, string mode, double percent, int structures,
	int window, double gammaM, int pkIter, int helix, double cutoff, string outputAlnFile) {

	// Get the number of sequences to read.
	int size = turboFold.files.size();

	// Initialize the temporary FASTA file name and the error string.
	const char* tempFileName = tmpnam( NULL );
	string errorString = "";

	// Build the sequences vector and the pfs files vector.
	vector<string> sequences, pfsFiles;
	for( int i = 1; i <= size; i++ ) {
		sequences.push_back( turboFold.files[i-1][0] );
		pfsFiles.push_back( turboFold.files[i-1][2] );
	}

	// Create a TurboFold object and its error checker.
	TurboFold* turbo = new TurboFold( &sequences, &pfsFiles, outputAlnFile );
	ErrorChecker<TurboFold>* checker = new ErrorChecker<TurboFold>( turbo );
	errorString = checker->returnError();

	resetProgress();
	progress->update(12);
	turbo->SetProgress(*progress);
	// Set the temperature of the TurboFold calculation.
	turbo->SetTemperature( constraintsHolder.temperature );

	// If no error occurred in the constructor, run TurboFold folding.
	if( errorString == "" ) {
		int mainError = turbo->fold( gammaT, turboIter );
		errorString = checker->returnError( mainError );
	}

	// If no error occurred in the main calculation, run a specific mode
	// calculation.
	if( errorString == "" ) {
		for( int i = 1; i <= size; i++ ) {
			int modeError = 0;
			if( mode == "MEA" ) {
				modeError = turbo->MaximizeExpectedAccuracy(
					i, percent, structures, window, gammaM );
			} else if( mode == "ProbKnot" ) {
				modeError = turbo->ProbKnot( i, pkIter, helix );
			} else {
				modeError = turbo->PredictProbablePairs( i, cutoff );
			}

			errorString = checker->returnError( modeError );
			if( errorString != "" ) { i += size; }
		}
	}

	// If no error occurred in mode calculation, write the CT files.
	if( errorString == "" ) {
		for( int i = 1; i <= size; i++ ) {
			int writeError =
				turbo->WriteCt( i, turboFold.files[i-1][1].c_str() );
			errorString = checker->returnError( writeError );
			if( errorString != "" ) { i += size; }
		}
	}

	// Return the error string.
	return errorString;
}

//////////////////////////////////////////////////////////////////////////////
// Run TurboFold calculations in maximum expected accuracy mode.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runTurboFoldMaximumExpectedAccuracy(
	double turboGamma, int turboIterations, double percent, int structures,
	int window, double meaGamma, string outAlnFile ) {

	// Run TurboFold in the workhorse method in MEA mode.
	return runTurboFold(
		turboGamma, turboIterations, "MEA",
		percent, structures, window, meaGamma, -1, -1, -1, outAlnFile );
}

//////////////////////////////////////////////////////////////////////////////
// Run TurboFold calculations in pseudoknot mode.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runTurboFoldPseudoknot(
	double turboGamma, int turboIterations, int pkIterations, int helix, string outAlnFile ) {

	// Run TurboFold in the workhorse method in ProbKnot mode.
	return runTurboFold(
		turboGamma, turboIterations, "ProbKnot",
		-1, -1, -1, -1, pkIterations, helix, -1, outAlnFile );
}

//////////////////////////////////////////////////////////////////////////////
// Run TurboFold calculations in Threshold mode.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runTurboFoldThreshold(
	double turboGamma, int turboIterations, double threshold, string outAlnFile ) {

	// Run TurboFold in the workhorse method in Threshold mode.
	return runTurboFold(
		turboGamma, turboIterations, "Threshold",
		-1, -1, -1, -1, -1, -1, threshold, outAlnFile );
}


//////////////////////////////////////////////////////////////////////////////
// Build an RNA data structure for modules that take a single sequence as input.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::loadSingleRNA(string file, const bool isRNA, const RNAInputType fileType) {
	if( rna != 0 ) { delete rna; }
	if( rnaChecker != 0 ) { delete rnaChecker; }
	rna = new RNA( file.c_str(), fileType, isRNA );
	rnaChecker = new ErrorChecker<RNA>( rna );
	return rnaChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Build a HybridRNA data structure for modules that take two sequences as input (e.g. bifold, bipartition, DuplexFold, AccessFold calculations).
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::loadDoubleRNA(
	string file1, string file2, const bool isRNA, const RNAInputType fileType ) {
	if( hybrid != 0 ) { delete hybrid; }
	if( hybridChecker != 0 ) { delete hybridChecker; }
	hybrid = new HybridRNA( file1.c_str(), fileType, file2.c_str(), fileType, isRNA );
	hybridChecker = new ErrorChecker<HybridRNA>( hybrid );
	return hybridChecker->returnError();
}

//////////////////////////////////////////////////////////////////////////////
// Build a data structure for AccessFold calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildAccessFoldDataStructure(string file1, string file2, bool isRNA ) {
	return loadDoubleRNA(file1, file2, isRNA);
}

//////////////////////////////////////////////////////////////////////////////
// Run AccessFold Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runAccessFold(string ctFile, 
		float percent=50, int maxStructures=20, double gamma=0.4,
		int windowSize=0, bool saveFile=false ) {
			
	// Set the calculation temperature.
	hybrid->SetTemperature( constraintsHolder.temperature );

	// Create a variable to handle error codes.
	int error = 0;

	// Determine the save file name.
	string saveFileName = saveFile ? generateSavFile(ctFile) : "";

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	hybrid->SetProgress( *progress );

	progress->update(20);

	// Run AccessFold.
	error = hybrid->AccessFold(gamma, percent, maxStructures, windowSize, constraintsHolder.maxLoop );

	// Delete the progress monitor.
	hybrid->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = hybrid->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return hybridChecker->returnError( error );
}


//////////////////////////////////////////////////////////////////////////////
// Build a data structure for DuplexFold calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::buildDuplexFoldDataStructure(string file1, string file2, bool isRNA ) {
	return loadDoubleRNA(file1, file2, isRNA);
}


//////////////////////////////////////////////////////////////////////////////
// Run DuplexFold Calculations.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::runDuplexFold(string ctFile, 
		float percent=50, int maxStructures=20,
		int windowSize=0, bool saveFile=false ) {
			
	// Set the calculation temperature.
	hybrid->SetTemperature( constraintsHolder.temperature );

	// Create a variable to handle error codes.
	int error = 0;

	// Determine the save file name.
	string saveFileName = saveFile ? generateSavFile(ctFile) : "";

	// Create the progress monitor.
	// RMW - removed ProgressMonitor:   monitor = new ProgressMonitor();
	resetProgress();  // RMW - removed ProgressMonitor:   TProgressDialog* dialog = new TProgressDialog( *(monitor) );
	hybrid->SetProgress( *progress );

	// Run DuplexFold.
	error = hybrid->FoldDuplex(percent, maxStructures, windowSize, constraintsHolder.maxLoop );

	// Delete the progress monitor.
	hybrid->StopProgress();
	// delete dialog;

	// If no error occurred, write a CT file, then check for errors after
	// writing.
	if( error == 0 ) { error = hybrid->WriteCt( ctFile.c_str() ); }

	// Return the error string.
	return hybridChecker->returnError( error );
}


/*
 * Constraints Mutators and Accessors.
 */

//////////////////////////////////////////////////////////////////////////////
// Remove all folding constraints.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::clearFoldingConstraints( int strand ) {

	// Remove constraints based on the current data structure.
	if( rna != 0 ) { rna->RemoveConstraints(); }
	else if( hybrid != 0 ) { hybrid->RemoveConstraints(); }
	else if( dynalign != 0 ) {
		if( strand == 1 ) {
			dynalign->GetRNA1()->RemoveConstraints();
		} else {
			dynalign->GetRNA2()->RemoveConstraints();
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Get the folding constraints.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getFoldingConstraints( int strand ) {

	// Get the strand whose constraints to show.
	RNA* data = 0;
	if( rna != 0 ) { data = rna; }
	else if( hybrid != 0 ) { data = hybrid; }
	else if( dynalign != 0 ) {
		if( strand == 1 ) { data = dynalign->GetRNA1(); }
		else { data = dynalign->GetRNA2(); }
	}
	if( data == 0 ) { return ""; }

	// Create a stringstream to hold the constraints data.
	stringstream stream( stringstream::in | stringstream::out );
	stream << "<html>";

	// Add data for forced base pairs.
	stream << "Forced Base Pairs <br>";
	int numType = data->GetNumberOfForcedPairs();
	if( numType < 1 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedPair( i - 1, true ) << " - "
			       << data->GetForcedPair( i - 1, false );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}
	stream << "<br><br>";

	// Add data for prohibited base pairs.
	stream << "Forced Prohibited Pairs <br>";
	numType = data->GetNumberOfForcedProhibitedPairs();
	if( numType == 0 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedProhibitedPair( i - 1, true ) << " - "
			       << data->GetForcedProhibitedPair( i - 1, false );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}
	stream << "<br><br>";

	// Add data for chemical modifications.
	stream << "Forced Modifications <br>";
	numType = data->GetNumberOfForcedModifications();
	if( numType == 0 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedModification( i - 1 );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}
	stream << "<br><br>";

	// Add data for FMN cleavages.
	stream << "Forced FMN Cleavages <br>";
	numType = data->GetNumberOfForcedFMNCleavages();
	if( numType == 0 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedFMNCleavage( i - 1 );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}
	stream << "<br><br>";

	// Add data for forced single stranded nucleotides.
	stream << "Forced Single Stranded <br>";
	numType = data->GetNumberOfForcedSingleStranded();
	if( numType == 0 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedSingleStranded( i - 1 );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}
	stream << "<br><br>";

	// Add data for forced double stranded nucleotides.
	stream << "Forced Double Stranded <br>";
	numType = data->GetNumberOfForcedDoubleStranded();
	if( numType == 0 ) { stream << "None"; }
	else {
		for( int i = 1; i <= numType; i++ ) {
			stream << data->GetForcedDoubleStranded( i - 1 );
			if( i != numType ) {
				if( ( i % 5 ) != 0 ) { stream << ", "; }
				else { stream << "<br>"; }
			}
		}
	}

	// Return the stringstream data.
	return stream.str();
}

//////////////////////////////////////////////////////////////////////////////
// Get the maximum constrainable index.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getMaxConstraintIndex( int strand ) {

	RNA* data = 0;
	if( rna != 0 ) { data = rna; }
	else if( hybrid != 0 ) { data = hybrid; }
	else if( dynalign != 0 ) {
		if( strand == 1 ) { data = dynalign->GetRNA1(); }
		else { data = dynalign->GetRNA2(); }
	}
	if( data == 0 ) { return 0; }
	else { return data->GetSequenceLength(); }
}

//////////////////////////////////////////////////////////////////////////////
// Get the maximum bulge/internal loop size.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getMaxLoop() {

	// Get the maximum loop size.
	return constraintsHolder.maxLoop;
}

//////////////////////////////////////////////////////////////////////////////
// Get the maximum pairing distance allowed between nucleotides.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getMaxPair() {

	// Get the maximum distance.
	return constraintsHolder.maxPair;
}

//////////////////////////////////////////////////////////////////////////////
// Get the calculation temperature.
//////////////////////////////////////////////////////////////////////////////
double RNAstructureBackendCalculator::getTemperature() {

	// Get the calculation temperature.
	return constraintsHolder.temperature;
}

//////////////////////////////////////////////////////////////////////////////
// Read in a folding constraints file.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::readFoldingConstraintsFile(
	string file, int strand ) {

	// Read a constraints file based on the current data structure.
	if( rna != 0 ) {
		int code = rna->ReadConstraints( file.c_str() );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->ReadConstraints( file.c_str() );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->ReadConstraints( file.c_str() ) :
			dynalign->GetRNA1()->ReadConstraints( file.c_str() );
		return dynalignChecker->returnError( code );
	}

	// If reading wasn't possible, return an error.
	return "Error reading folding constraints file.";

}

//////////////////////////////////////////////////////////////////////////////
// Set a cleaved nucleotide.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setCleavedNucleotide(
	int nuc, int strand ) {

	// Set a cleaved nucleotide based on the current data structure.
	if( rna != 0 ) {
		int code = rna->ForceFMNCleavage( nuc );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->ForceFMNCleavage( nuc );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->ForceFMNCleavage( nuc ) :
			dynalign->GetRNA1()->ForceFMNCleavage( nuc );
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This nucleotide cannot currently be cleaved.";
}

//////////////////////////////////////////////////////////////////////////////
// Set a double stranded nucleotide.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setDoubleStrandedNucleotide(
	int nuc, int strand ) {

	// Set a double stranded nucleotide based on the current data structure.
	if( rna != 0 ) {
		int code = rna->ForceDoubleStranded( nuc );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->ForceDoubleStranded( nuc );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->ForceDoubleStranded( nuc ) :
			dynalign->GetRNA1()->ForceDoubleStranded( nuc );
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This nucleotide cannot currently be forced double stranded.";
}

//////////////////////////////////////////////////////////////////////////////
// Set a forced helix.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setForcedHelix(
	int nuc1, int nuc2, int length, int strand ) {

	// Set a helix based on the current data structure.
	if( rna != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			code = rna->ForcePair( nuc1 + increment, nuc2 - increment );
			if( code != 0 ) { break; }
		}
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			code = hybrid->ForcePair( nuc1 + increment, nuc2 - increment );
			if( code != 0 ) { break; }
		}
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			int pair1 = nuc1 + increment;
			int pair2 = nuc2 + increment;
			code = ( strand != 2 ) ?
				dynalign->GetRNA1()->ForcePair( pair1, pair2 ) :
				dynalign->GetRNA2()->ForcePair( pair1, pair2 );
			if( code != 0 ) { break; }
		}
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This helix cannot currently be forced.";
}

//////////////////////////////////////////////////////////////////////////////
// Set the maximum bulge/internal loop size.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setMaxLoop( int newLoop ) {

	// Set the maximum loop size.
	constraintsHolder.maxLoop = newLoop;
}

//////////////////////////////////////////////////////////////////////////////
// Set the maximum pairing distance allowed between nucleotides.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setMaxPair( int newPair ) {

	// Set the maximum distance.
	constraintsHolder.maxPair = newPair;
}

//////////////////////////////////////////////////////////////////////////////
// Set a modified nucleotide.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setModifiedNucleotide(
	int nuc, int strand ) {

	// Set a modified nucleotide based on the current data structure.
	if( rna != 0 ) {
		int code = rna->ForceModification( nuc );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->ForceModification( nuc );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->ForceModification( nuc ) :
			dynalign->GetRNA1()->ForceModification( nuc );
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This nucleotide cannot currently be modified.";
}

//////////////////////////////////////////////////////////////////////////////
// Set a forced helix.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setProhibitedHelix(
	int nuc1, int nuc2, int length, int strand ) {

	// Set a prohibited helix based on the current data structure.
	if( rna != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			int pair1 = nuc1 + increment;
			int pair2 = nuc2 + increment;
			code = rna->ForceProhibitPair( pair1, pair2 );
			if( code != 0 ) { break; }
		}
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			int pair1 = nuc1 + increment;
			int pair2 = nuc2 + increment;
			code = hybrid->ForceProhibitPair( pair1, pair2 );
			if( code != 0 ) { break; }
		}
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = 0;
		for( int i = 1; i <= length; i++ ) {
			int increment = i - 1;
			int pair1 = nuc1 + increment;
			int pair2 = nuc2 + increment;
			code = ( strand != 2 ) ?
				dynalign->GetRNA1()->ForceProhibitPair( pair1, pair2 ) :
				dynalign->GetRNA2()->ForceProhibitPair( pair1, pair2 );
			if( code != 0 ) { break; }
		}
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This helix cannot currently be prohibited.";
}

//////////////////////////////////////////////////////////////////////////////
// Set the SHAPE constraints file.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSHAPEFile( string file ) {

	// Set the SHAPE file.
	constraintsHolder.shapeFile = file;
}

//////////////////////////////////////////////////////////////////////////////
// Set the first SHAPE constraints parameter.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSHAPEParam1( double value ) {

	// Set the SHAPE parameter.
	constraintsHolder.shapeParam1 = value;
}

//////////////////////////////////////////////////////////////////////////////
// Set the second SHAPE constraints parameter.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSHAPEParam2( double value ) {

	// Set the SHAPE parameter.
	constraintsHolder.shapeParam2 = value;
}

//////////////////////////////////////////////////////////////////////////////
// Set the type of SHAPE constraints used.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setSHAPEType( bool isEnergy ) {

	// Set the SHAPE constraints type.
	constraintsHolder.shapeEnergy = isEnergy;
}

//////////////////////////////////////////////////////////////////////////////
// Set a single stranded nucleotide.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::setSingleStrandedNucleotide(
	int nuc, int strand ) {

	// Set a single stranded nucleotide based on the current data structure.
	if( rna != 0 ) {
		int code = rna->ForceSingleStranded( nuc );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->ForceSingleStranded( nuc );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->ForceSingleStranded( nuc ) :
			dynalign->GetRNA1()->ForceSingleStranded( nuc );
		return dynalignChecker->returnError( code );
	}

	// If setting wasn't possible, return an error.
	return "This nucleotide cannot currently be forced single stranded.";
}

//////////////////////////////////////////////////////////////////////////////
// Set the calculation temperature.
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::setTemperature( double newTemp ) {

	// Set the calculation temperature.
	constraintsHolder.temperature = newTemp;
}

//////////////////////////////////////////////////////////////////////////////
// Write in a folding constraints file.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::writeFoldingConstraintsFile(
	string file, int strand ) {

	// Write a constraints file based on the current data structure.
	if( rna != 0 ) {
		int code = rna->WriteConstraints( file.c_str() );
		return rnaChecker->returnError( code );
	} else if( hybrid != 0 ) {
		int code = hybrid->WriteConstraints( file.c_str() );
		return hybridChecker->returnError( code );
	} else if( dynalign != 0 ) {
		int code = ( strand == 2 ) ?
			dynalign->GetRNA2()->WriteConstraints( file.c_str() ) :
			dynalign->GetRNA1()->WriteConstraints( file.c_str() );
		return dynalignChecker->returnError( code );
	}

	// If writing wasn't possible, return an error.
	return "Error writing folding constraints file.";

}



/*
 * Utility Methods.
 */

//////////////////////////////////////////////////////////////////////////////
// Get the percent progress for this calculation.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::getProgressNumber() {
	// Return the percent progress.
	return progress==NULL ? 0 : progress->progress();
}

//////////////////////////////////////////////////////////////////////////////
// Cancel this calculation (when implemented in native method).
//////////////////////////////////////////////////////////////////////////////
void RNAstructureBackendCalculator::cancelOperation() {
	if (progress!=NULL)
		progress->cancel();
}

bool RNAstructureBackendCalculator::wasCanceled() {
	return progress!=0&&progress->canceled();
}


// Typically functions that can generate a saveFile (e.g. Fold)
// infer the name of the saveFile from the output file name.
// however, the user can override this by calling setSaveFile to 
// set the name explicitly before calling runFold etc.
void RNAstructureBackendCalculator::setSaveFile(string path) {
	saveFilePath = path;
}

// Returns the string set by setSaveFile or the empty string if 
// it has not been set.
string RNAstructureBackendCalculator::getSaveFile() {
	return saveFilePath;

}

// If saveFilePath is NOT empty (i.e. the user has set saveFilePath explicitly
// via setSaveFile) then saveFilePath is returned as-is.
// Otherwise, this function generates a new saveFilePath derived from 
// outputFile by changing the file extension to ".sav". If outputFile
// has no extension, then ".sav" is appended.
string RNAstructureBackendCalculator::generateSavFile(const string outputFile, const char* const savExtension) { 
	return saveFilePath.empty() ? 
		chageFileExt(outputFile, savExtension) :
		saveFilePath;
}

//////////////////////////////////////////////////////////////////////////////
// Get the type of structure being handled.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getStructureType() {

	// Get the sequence length and the number of structures.
	int length = rna->GetSequenceLength();
	int numStructures = rna->GetStructureNumber();

	// Go through the length of each structure and check for pairs.
	// If pairs are found, set the pairs flag.
	bool hasPairs = false;
	for( int i = 1; i <= numStructures; i++ ) {
		for( int j = 1; j <= length; j++ ) {
			if( rna->GetPair( j, i ) > j ) {
				hasPairs = true;
				i = numStructures + 1;
				j = length + 1;
			}
		}
	}

	// If no pairs were found, return an error.
	if( !hasPairs ) { return "This structure contains no pairs."; }

	// Check each individual structure for pseudoknots.
	for( int i = 1; i <= numStructures; i++ ) {

		// Check if the structure contains a pseudoknot.
		bool pseudoknot = rna->ContainsPseudoknot( i );

		// If an error occurred trying to find a pseudoknot, return it.
		string error = rnaChecker->returnError();
		if( error != "" ) { return error; }

		// If the structure contains a pseudoknot, return "Circular".
		if( pseudoknot ) { return "Circular"; }
	}

	// If the method got here, the structures have pairs but no pseudoknots,
	// so return "Radial".
	return "Radial";
}

// string showPID() {
// 	char buffer [255];
// 	sprintf( buffer, "pid: %d  ppid: %d  tid: %ld", 
// 		getpid(), getppid(), syscall(SYS_gettid));
// 	return string(buffer);
// }


//////////////////////////////////////////////////////////////////////////////
// Set a system environment variable for the current process, by passing a string in the form "NAME=VALUE".
// Returns 0 if successful or non-zero otherwise.
//////////////////////////////////////////////////////////////////////////////
int RNAstructureBackendCalculator::setEnvVar(const string& envstr) {
	// A copy of the c_str must be made because putenv on linux 
	//  and Mac is not putenv(const char*), it is putenv(char*)
	//  I.e. the OS does NOT make a copy of the char*, but instead stores the
	//  pointer to it. This means that modifying the *str will modify the 
	//  environment.  A corrollary is that the *str must NOT be freed, or 
	//  it will leave the environment with a pointer to invalid memory.
	// see: http://www.greenend.org.uk/rjk/tech/putenv.html	
	char *str = copy_cstr(envstr.c_str()); 
	if (str == NULL) return 1; //error
	return putenv(str);
}

//////////////////////////////////////////////////////////////////////////////
// Get a system environment variable for the current process, 
// by passing a string in the name of the variable.
// Returns an empty string if not found.
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getEnvVar(const string& envvar) {
	return as_str(getenv(envvar.c_str())); //creates a COPY of the char*
}

//////////////////////////////////////////////////////////////////////////////
// Get the official build version of the RNAstructure library 
// (defined in src/version.h), e.g. "5.8"
//////////////////////////////////////////////////////////////////////////////
string RNAstructureBackendCalculator::getVersion() {
	return string(RNASTR_BUILD_VERSION); //creates a COPY of the char*
}

int RNAstructureBackendCalculator::GetErrorCode() {
	//return ErrorCode;
	if (rna!=NULL&&rna->GetErrorCode()) return rna->GetErrorCode();
	if (dynalign!=NULL&&dynalign->GetErrorCode()) return dynalign->GetErrorCode();
	if (hybrid!=NULL&&hybrid->GetErrorCode()) return hybrid->GetErrorCode();
	if (oligo!=NULL&&oligo->GetErrorCode()) return oligo->GetErrorCode();
	// if (rnaChecker!=NULL&&rnaChecker->isErrorStatus(false)) return rnaChecker->isErrorStatus(false);
	// if (dynalignChecker!=NULL&&dynalignChecker->isErrorStatus(false)) return dynalignChecker->isErrorStatus(false);
	// if (hybridChecker!=NULL&&hybridChecker->isErrorStatus(false)) return hybridChecker->isErrorStatus(false);
	// if (oligoChecker!=NULL&&oligoChecker->isErrorStatus(false)) return oligoChecker->isErrorStatus(false);
}
string RNAstructureBackendCalculator::GetFullErrorMessage() {
	if (rna!=NULL&&rna->GetErrorCode()) return rnaChecker->returnError();
	if (dynalign!=NULL&&dynalign->GetErrorCode()) return dynalignChecker->returnError();
	if (hybrid!=NULL&&hybrid->GetErrorCode()) return hybridChecker->returnError();
	if (oligo!=NULL&&oligo->GetErrorCode()) return oligoChecker->returnError();
}

// Note: The ProgressHandler that would be used below would be different from that created by RNA.i due to the 
// separate invocations of SWIG.
// ProgressHandler& RNAstructureBackendCalculator::GetProgress() { return *progress; }
// void RNAstructureBackendCalculator::SetProgress(ProgressHandler &new_progress, const bool deleteWhenFinished) {
// 	if (progress!=NULL&&ownsProgress) delete progress;
// 	progress=&new_progress;
// 	ownsProgress=deleteWhenFinished;
// }
// void RNAstructureBackendCalculator::SetError(const int code, const string &details) {
// 	if (ErrorCode==0) ErrorCode = code;
// 	if (!details.empty()) {
// 		if (errorDetails.empty())
// 			errorDetails = details;
// 		else
// 			errorDetails = errorDetails + "\n" + details;
// 	}
// }