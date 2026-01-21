/*
 * A program that compares two structures in a circular layout.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "CircleCompare_Interface.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
CircleCompare_Interface::CircleCompare_Interface() {

	// Initialize the calculation type description.
	calcType = "Circular structure comparison";

	// Initialize the optional annotation file names.
	probabilityFile = "";
	SHAPEFile = "";

	// Set most boolean flags to false.
	alternative = false;
	exact = false;
	filenames = false;
	levorotatory = false;
	probabilityForPredicted = false;
	isSVG = false;

	// Set the encircling of nucleotides to true.
	encircle = true;

	// The specific structure in the predicted CT file that should be compared to the accepted CT.
	// Set this to -1 as the default, which compares all structures individually.
	number = -1;

	// Set the number of structures in the predicted strand to a default of 1.
	// No matter what, if the program is run correctly, there will be at least one predicted structure compared.
	predictedCount = 1;

	// Initialize the text annotation flag to false.
	textAnnotation = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool CircleCompare_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "CircleCompare" );
	parser->addParameterDescription( "predicted ct", "The name of a file containing CT data for the predicted structure." );
	parser->addParameterDescription( "accepted ct", "The name of a file containing CT data for the accepted structure." );
	parser->addParameterDescription( "output file", "The name of an image file to which output will be written. Usually, this is a Postscript image file, although the user can specify that it be an SVG image file instead. To specify SVG images, the \"--svg\" flag must be specified in conjunction with a structure number flag." );

	// Add the alternative option.
	vector<string> alternativeOptions;
	alternativeOptions.push_back( "-a" );
	alternativeOptions.push_back( "-A" );
	alternativeOptions.push_back( "--alternative" );
	parser->addOptionFlagsNoParameters( alternativeOptions, "Specify that an alternative color scheme should be used. Default is not to use the alternative color scheme." );

	// Add the exact option.
	vector<string> exactOptions;
	exactOptions.push_back( "-e" );
	exactOptions.push_back( "-E" );
	exactOptions.push_back( "--exact" );
	parser->addOptionFlagsNoParameters( exactOptions, "Specify exact comparison when structure comparison is scored. Default is to allow flexible pairings." );

	// Add the file option.
	vector<string> fileOptions;
	fileOptions.push_back( "-f" );
	fileOptions.push_back( "-F" );
	fileOptions.push_back( "--file" );
	parser->addOptionFlagsNoParameters( fileOptions, "Specify that structure file names should be shown in addition to their descriptions. Default is not to show structure file names." );

	// Add the levorotatory option.
	vector<string> levorotatoryOptions;
	levorotatoryOptions.push_back( "-l" );
	levorotatoryOptions.push_back( "-L" );
	levorotatoryOptions.push_back( "--levorotatory" );
	parser->addOptionFlagsNoParameters( levorotatoryOptions, "Specify that the drawn structure is rendered counterclockwise. Default is to render drawn structures clockwise." );

	// Add the number option.
	vector<string> numberOptions;
	numberOptions.push_back( "-n" );
	numberOptions.push_back( "-N" );
	numberOptions.push_back( "--number" );
	parser->addOptionFlagsWithParameters( numberOptions, "Specify the index of a particular structure in the predicted CT to be compared with the accepted CT, one-indexed. Default is -1, which signifies all structures output to one file." );

	// Add the probability option for the predicted structure.
	vector<string> probabilityOptions;
	probabilityOptions.push_back( "-p" );
	probabilityOptions.push_back( "-P" );
	probabilityOptions.push_back( "--probability" );
	parser->addOptionFlagsWithParameters( probabilityOptions, "Specify the name of the file from which base pairing probability data will be read for annotation. This file should describe pairing data for the predicted structure, not the accepted structure. Default is no probability annotation file used." );

	// Add the probability option for the accepted structure.
	vector<string> probability2Options;
	probability2Options.push_back( "-p2" );
	probability2Options.push_back( "-P2" );
	probability2Options.push_back( "--probability2" );
	parser->addOptionFlagsWithParameters( probability2Options, "Specify the name of the file from which base pairing probability data will be read for annotation. This file should describe pairing data for the accepted structure, not the predicted structure. Default is no probability annotation file used." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-s" );
	shapeOptions.push_back( "-S" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify the name of the file from which SHAPE data will be read for annotation. Default is no SHAPE annotation file used." );

	// Add the SVG option.
	vector<string> svgOptions;
	svgOptions.push_back( "--svg" );
	parser->addOptionFlagsNoParameters( svgOptions, "Specify that the output file should be an SVG image file, rather than a Postscript image file. Note that only one SVG image can be written into a particular file, so the structure number flag must also be specified when writing an SVG document." );

	// Add the text annotation option.
	vector<string> textOptions;
	textOptions.push_back( "-t" );
	textOptions.push_back( "-T" );
	textOptions.push_back( "--text" );
	parser->addOptionFlagsWithParameters( textOptions, "Specify the name of the text file from which base pairing probability data will be read for annotation. This file should describe pairing data for the predicted structure, not the accepted structure. Default is no probability annotation file used." );

	// Add the uncircled option.
	vector<string> uncircledOptions;
	uncircledOptions.push_back( "-u" );
	uncircledOptions.push_back( "-U" );
	uncircledOptions.push_back( "--uncircled" );
	parser->addOptionFlagsNoParameters( uncircledOptions, "Specify that no circles should surround nucleotides when drawing. Default is to surround nucleotides with circles." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		predicted = parser->getParameter( 1 );
		accepted = parser->getParameter( 2 );
		output = parser->getParameter( 3 );
	}

	// Get the alternative option.
	if( !parser->isError() ) { alternative = parser->contains( alternativeOptions ); }

	// Get the exact option.
	if( !parser->isError() ) { exact = parser->contains( exactOptions ); }

	// Get the file option.
	if( !parser->isError() ) { filenames = parser->contains( fileOptions ); }

	// Get the levorotatory option.
	if( !parser->isError() ) { levorotatory = parser->contains( levorotatoryOptions ); }

	// Get the number option.
	if( !parser->isError() ) {
		parser->setOptionInteger( numberOptions, number );
		bool badNumber =
			( number != -1 ) &&
			( number < 0 );
		if( badNumber ) { parser->setError( "structure number" ); }
	}

	// Get the probability annotation file option.
	// This can be read either as a partition function save file or a text file.
	if( !parser->isError() ) {
		if( parser->contains( probabilityOptions ) ) {
			probabilityFile = parser->getOptionString( probabilityOptions, true );
			probabilityForPredicted = true;
		}

		if( parser->contains( probability2Options ) ) {
			probabilityFile = parser->getOptionString( probability2Options, true );
			probabilityForPredicted = false;
		}

		if( parser->contains( textOptions ) ) {
			probabilityFile = parser->getOptionString( textOptions, true );
			textAnnotation = true;
			probabilityForPredicted = true;
		}
	}

	// Get the SHAPE annotation file option.
	if( !parser->isError() ) { SHAPEFile = parser->getOptionString( shapeOptions, true ); }

	// Get the SVG option.
	if( !parser->isError() ) { isSVG = parser->contains( svgOptions ); }

	// Get the uncircled option.
	if( !parser->isError() ) { encircle = !parser->contains( uncircledOptions ); }

	// If in SVG mode and no structure number is specified, show an error message.
	if( ( !parser->isError() ) && parser->contains( svgOptions ) && ( !parser->contains( numberOptions ) ) ) {
		cerr << "No structure number specified with SVG structure." << endl
		     << "Please specify the structure to draw as an SVG image file." << endl;
		parser->setError();
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void CircleCompare_Interface::run() {

	/**
	 * Create and validate the comparison object.
	 */

	// Print a message saying that validation has started.
	cout << "Validating structure files..." << flush;

	// Create a comparison object and validate it.
	StructureComparedImageHandler* compare = new StructureComparedImageHandler();
	string validResult = compare->validateFiles( predicted, accepted, number );
	int numPredictedStructures = compare->getNumPredictedStructures();

	// Delete the comparison object and print a message saying that validation has finished.
	// If an error occurred, show it and return.
	delete compare;
	if( validResult == "" ) { cout << "done." << endl; }
	else {
		cerr << endl << validResult << endl
		     << calcType << " complete with errors." << endl;
		return;
	}

	/**
	 * Assemble the structure comparisons if the comparison object is valid.
	 */

	// Print a message saying that structure comparison has started.
	cout << "Assembling circular structure comparison file..." << endl;

	// If the output file already exists, delete the file so a new file with that name can be generated.
	ifstream test( output.c_str() );
	bool exists = test.good();
	test.close();
	if( exists ) { remove( output.c_str() ); }

	// Determine which structures are going to be drawn.
	int start = ( number == -1 ) ? 1 : number;
	int end = ( number == -1 ) ? numPredictedStructures : number;

	// Read in the data for each individual structure and write them to the appropriate output file.
	for( int i = start; i <= end; i++ ) {

		// Print a message saying that a specific structure comparison has started.
		cout << "    Assembling comparison " << i << "..." << flush;

		// Create a new comparison object and set its defaults.
		StructureComparedImageHandler* compare = new StructureComparedImageHandler();
		compare->setColorScheme( alternative );
		compare->setNucleotidesCircled( encircle );

		// Read the pairs.
		string compareResult = "";
		if( ( probabilityForPredicted == true ) && ( probabilityFile != "" ) ) {
			compareResult = compare->readPredictedStructureProbability( predicted, i, probabilityFile, textAnnotation );
		} else if( SHAPEFile != "" ) {
			compareResult = compare->readPredictedStructureSHAPE( predicted, i, SHAPEFile );
		} else {
			compareResult = compare->readPredictedStructure( predicted, i );
		}

		// Read in the accepted structure.
		if( compareResult == "" ) {
			if( ( probabilityForPredicted == false ) && ( probabilityFile != "" ) ) {
				compareResult = compare->readAcceptedStructureProbability( accepted, probabilityFile, textAnnotation );
			} else {
				compareResult = compare->readAcceptedStructure( accepted );
			}
		}

		// Overlay the comparison and set its data.
		if( compareResult == "" ) {
			compare->overlayStructures();
			compare->addComparisonData( predicted, accepted, i, exact, isSVG, filenames );
		}

		// Flip the structure, if asked.
		if( ( compareResult == "" ) && ( levorotatory == true ) ) { compare->flipHorizontally(); }

		// Write the structure image, as either SVG or Postscript.
		if( compareResult == "" ) {
			if( isSVG ) { compare->writeSVG( output ); }
			else { compare->writePostscript( output, ( number == -1 ) ); }
		}

		// Delete the comparison object.
		delete compare;

		// If no error occurred, print a message saying that a specific structure comparison has finished.
		// If any errors occurred in this iteration, show the error and return.
		if( compareResult == "" ) { cout << "done." << endl; }
		else {
			cerr << endl << "Error comparing accepted structure with predicted structure " << i << ": " << compareResult << endl
			   	 << calcType << " complete with errors." << endl;
			return;
		}
	}

	/*
	 * Print confirmation of run finishing successfully, if it did.
	 */

	cout << calcType << " complete." << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	CircleCompare_Interface* runner = new CircleCompare_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
