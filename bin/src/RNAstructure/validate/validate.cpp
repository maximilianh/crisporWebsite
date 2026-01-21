/*
 * A program that verifies the format of several types of files used by the 
 * RNAstructure sofware package.
 *
 * (c) 2017-2018 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */

#include "validate.h"
#include <limits.h>

#define MAX_SEQUENCE_LABEL_LENGTH 500 // maximum characters in a sequence label. (Arbitrary limit to prevent malicious users of the web interface from uploading garbage.)

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
validate_Interface::validate_Interface() {

	// Initialize the calculation type description.
	fileType = "auto"; //automatically determine type from extension.
	knownType = RFT_UNKNOWN;
	file = "";
	maxLength = INT_MAX; // maximum sequence length;
	alphabetName = "rna";
	allowUnknownBases = false;
	errorMessage = "";
	show_info = false;
	quiet = false;

	// tracks count of sequences for which info is added (if show_info is true)
	info_seq_count=0;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool validate_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "validate" );
	parser->addParameterDescription( "file path", "The path to the file that should be validated." );

	// Add the type option.
	vector<string> typeOptions;
	typeOptions.push_back( "-t" );
	typeOptions.push_back( "--type" );
	parser->addOptionFlagsWithParameters(  typeOptions, "The type of file validation to perform. "
		"This can be numeric (1=CT, 2=SEQ/FASTA, 3=PFS, 4=SAV, 5=DotBracket, 6=CON, 7=SHAPE/CHEM, 8=OLIS)"
		" or text (one of: ct, seq|fasta, pfs, sav, dot|dbn|braket, con, shape|chem, lis|list|olis)"
		" or \"auto\" (the default) which uses the file extension to determine the type." );

	// Add the sequence length limit option
	vector<string> lengthOptions;
	lengthOptions.push_back( "-l" );
	lengthOptions.push_back( "--length" );
	parser->addOptionFlagsWithParameters(  lengthOptions, "Limit the length of the validated sequence to the specified number of bases." );

	// Add the alphabet option
	vector<string> alphabetOptions;
	alphabetOptions.push_back( "-a" );
	alphabetOptions.push_back( "--alphabet" );
	parser->addOptionFlagsWithParameters(  alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag." );

	// Add the unknownBases option
	vector<string> unknownBasesOptions;
	unknownBasesOptions.push_back( "-u" );
	unknownBasesOptions.push_back( "--unknown" );
	parser->addOptionFlagsNoParameters(  unknownBasesOptions, "Suppress errors due to unknown bases in sequences." );

	vector<string> quietOption;
	quietOption.push_back( "-q" );
	quietOption.push_back( "--quiet" );
	parser->addOptionFlagsNoParameters(  quietOption, "Suppress unnecessary output. The process exit code indicates the result of validation (as usual)." );

	vector<string> infoOption;
	infoOption.push_back( "-i" );
	infoOption.push_back( "--info" );
	parser->addOptionFlagsNoParameters(  infoOption, "Display information about the sequence. (This also implies --quiet.) Currently the only information shown is the sequence length of each validated sequence." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	string numberString;
	if( !parser->isError() ) {

		// get the name of the file to validate
		file = parser->getParameter( 1, true );

		// Get the type option.
		parser->setOptionString(typeOptions, fileType);

		// Get the length limit option.
		parser->setOptionInteger(lengthOptions, maxLength);

		parser->setOptionString(alphabetOptions, alphabetName);

		allowUnknownBases = parser->contains( unknownBasesOptions );

		show_info = parser->contains( infoOption );

		quiet = show_info || parser->contains( quietOption );
	}

	if( !parser->isError() ) {
		knownType = parseFileType(file, fileType);
		if (knownType == RFT_UNKNOWN) {
			if (fileType == "auto" || fileType=="AUTO")
				cerr << "File type could not be determined from file extension. Please specify the type explicitly with the " << typeOptions.back() << " option. (File: \"" << file << "\")." << endl;
			else
				cerr << "Invalid or unknown file type: \"" << fileType << "\"." << endl;
			parser->setError();
		}
	}


	// Delete the parser and return whether the parser encountered an error.
	bool hadError = parser->isError();
	delete parser;
	return !hadError;
}

// Make a string upper-case
void strToUpper(string& s) {
	for (int i = 0; i < s.size(); i++)
		 s[i] = toupper(s[i]);
}

RNAFileType validate_Interface::parseFileType(const string& file, const string& userType) {
	string type = userType; // get upper-case version of type.
	strToUpper(type);

	if (type=="AUTO") return fileTypeFromPath(file);

	// test for numeric values.
	if (type=="1") return RFT_CT;
	if (type=="2") return RFT_SEQ;
	if (type=="3") return RFT_PFS;
	if (type=="4") return RFT_SAV;
	if (type=="5") return RFT_DBN;
	if (type=="6") return RFT_CON;
	if (type=="7") return RFT_SHAPE;
	if (type=="8") return RFT_LIS;
	if (type=="9") return RFT_MULTISEQ;

	return fileTypeFromExt(userType);

	return RFT_UNKNOWN;
}

void validate_Interface::addInfo(const string& seq) {
	if (!show_info) return;
	info_seq_count++;
	cout << "Seq " << info_seq_count << "\tLength:" << seq.size() << endl;
}
void validate_Interface::addInfo(RNA& rna) {
	// In the future, show other info, such as energy etc.
	addInfo(rna.GetSequence());
}

bool validate_Interface::validateOligoList() {
	// first validate the file as a single text sequence. 
	// This will identify any invalid characters.
	RNA rna(file.c_str(), RFT_SEQ, alphabetName.c_str(), allowUnknownBases, SKIP_THERMO);
	if( rna.GetErrorCode() != 0 ) {
		errorMessage = rna.GetFullErrorMessage(); // already ends in "\n"
		return false;
	}
	
	// now validate each sequence line-by-line to verify sequence length.
	string line;
	ifstream input(file.c_str());
	long linenum = 0;
	while(getline(input, line)) {
		linenum++;
		line.erase( std::remove_if( line.begin(), line.end(), ::isspace ), line.end() ); // remove whitespace (even embedded within the line)
		if (line.length()==0)
			continue;
		if (line.length() > maxLength) {
			errorMessage = sfmt("Sequence exceeds the maximum allowed number of bases (%d). (line %ld).", maxLength, linenum);
			return false;
		}
		addInfo(line); // Add information about each sequences if the user requested it (show_info)
	}
	return true;
}

bool validate_Interface::validateFoldConstraintFile() {
	if (!fileExists(file.c_str())) {
		errorMessage = RNA::GetErrorMessage(1);
		return false;
	}
	ifstream input(file.c_str());
	if (!input) {
		errorMessage = RNA::GetErrorMessage(2);
		return false;
	}
	
	string line;
	string sectionName;
	int linenum = 0;
	int sectionNum = 0;
	bool sectionEnded = true;
	const char* sections[] = {"DS:","SS:","Mod:","Pairs:","FMN:","Forbids:"};
	stringstream s;
	while (!getline(input, line).fail()) {
		linenum++;
		if (line.find_first_not_of("\r\t ")==string::npos)
			continue; // blank line or just whitespace.
		
		//cout << "parsing line " << linenum << "\t[" << line << "]" << endl;

		s.clear();
		s.str(line); // assign line to stringstream.
		if (sectionEnded) {
			if (sectionNum > 5) {
				errorMessage = sfmt("Invalid data after the last constraint section ended (line %d).", linenum);
				return false;
			}

			s >> sectionName;
			if (sectionName != sections[sectionNum]) {
				errorMessage = sfmt("Missing constraint section %s on line %d.", sections[sectionNum], linenum);
				return false;
			}
			sectionEnded = false;
			continue;
		}
		int i1, i2;
		s >> i1;
		if (sectionNum == 3 || sectionNum == 5)
			s >> i2;
		
		if (s.fail()) {
			errorMessage = sfmt("Invalid data on line %d. Expected a numeric index.", linenum);
			return false;
		}
		if (i1 == -1) {
			sectionNum++;
			sectionEnded = true;
		}
	}
	if (sectionNum < 6) {
		errorMessage = sfmt("Missing constraint section %s.", sections[sectionNum]);
		return false;
	}
	return true;
}

// Validates a FASTA or SEQ file that contains multiple sequences (instead of the usual single sequence)
bool validate_Interface::validateMultipleSequenceFile() {
	if (!fileExists(file.c_str())) {
		errorMessage = RNA::GetErrorMessage(1);
		return false;
	}
	ifstream input(file.c_str());
	if (!input) {
		errorMessage = RNA::GetErrorMessage(2);
		return false;
	}
	
	datatable data;
	data.opendat(NULL, alphabetName.c_str(), false, SKIP_THERMO);
	data.allowUnknownBases = this->allowUnknownBases;

	  // handle both FASTA and SEQ format. 
      // FASTA file have the following format:
      // '>' <SEQUENCE_LABEL>                   // A greater-than symbol (>) starts each new sequence. The symbol is followed on the same line by the sequence label
      // <SEQUENCE...>                          // The sequence immediately follows the label-line. Technically it should not contain whitespace, but this is tolerated.
      // '>' <SEQUENCE_LABEL>                   // Subsequence sequences follow, each starting with '>'
      // <SEQUENCE...>
      // ...
      // SEQ file have the following format:
      // ';' <COMMENT>                          // SEQ files always start with at least one comment line, starting with ';'
      // [ ';' <COMMENT> ...  ]                 // multiple subsequent comment lines are allowed.
      // <SEQUENCE_LABEL>                       // The first line that doesn't start with ';' is the sequence label.
      // <SEQUENCE...>                          // The sequence starts on the line after the label.

	string line;
	int linenum = 0, count = 0;
	string sequence;
	bool started=false,inComment=false,isSEQ=false;

	int state;	
	while (!getline(input, line).fail()) {
		linenum++;
		bool verifyPrev=false, isLabel=false, startNext=false;
		if (line.find_first_not_of("\r\t ")==string::npos)
			continue; // blank line or just whitespace.
		
		if (line[0]=='>')  { // Check for next FASTA sequence
			verifyPrev=true; // verify the previous sequence if one was started.
			isLabel=true;
			startNext=true;
		} else if (line[0] == ';') { // Check for next SEQ sequence
			verifyPrev=true; // verify the previous sequence if one was started.
			inComment=true;  // This line and possibly subsequent lines are comments.
			isSEQ=true; // allows '1' to end the sequence.
		} else if (inComment) {
			// previous line(s) were SEQ comment(s) (starting with ';')
			// So this line must be the sequence label.
			isLabel = true;
			inComment = false;
			startNext = true;
		}

		if (isLabel) {
			// this is the beginning of a new FASTA sequence. The sequence label is on this line, so just ignore it, as long as it is a reasonable size.
			if (line.size() > MAX_SEQUENCE_LABEL_LENGTH) {
				errorMessage = sfmt("Error in sequence %d: The sequence label is too long (on line %d).", count+1, linenum); // use count+1 because count is not updated until the startNext routine below
				return false;
			}
		}

		if (verifyPrev & started) {
			// the previous sequence ended. validate it.
			if (sequence.size() > this->maxLength) {
				errorMessage = sfmt("Sequence exceeds the maximum allowed number of bases (Size: %d, Maximum: %d).", sequence.size(), maxLength);
				return false;
			}
			if (sequence.empty()) {
				errorMessage = sfmt("Sequence %d is empty.", count); // here count refers to the PREVIOUS sequence (count has not yet been incremented for the current one)
				return false;
			}
			addInfo(sequence); // Add information about each sequences if the user requested it (show_info)
			started = false;
		}

		if (!(isLabel|inComment)) {
			bool ended = false;
			for (int i = 0; i < line.size(); i++) {
				char c = line[i];
				if (::isspace(c)) continue;
				
				if (isSEQ && c=='1') {
					ended = true;
					continue;
				} 
				if (ended) {
					errorMessage = sfmt("Invalid data past end of SEQ file (i.e. past the terminating '1' character).  (Found '%c' in sequence %d on line %d column %d).", c, count, linenum, i+1);
					return 0;
				}

				int num = data.basetonum(c);
				if (num == -1) {
					errorMessage = sfmt("Invalid or unknown nucleobase '%c' (in sequence %d on line %d column %d).", c, count, linenum, i+1);
					return 0;
				}

				sequence.push_back(c);
				started=true; // just in case this is a plain-text sequence.
			}
		}
		
		if (startNext) {
			count++; // we are officially in the next sequence.
			started = true;
			sequence.clear();
		}
	}

	// reading is done. Verify state.
	// Note: comments are allowed at the end. As long as there is no label-line, it should be OK.
	if (!inComment && started) {  
		// the previous sequence ended. validate it.
		if (sequence.size() > this->maxLength) 
			errorMessage = sfmt("Sequence exceeds the maximum allowed number of bases (Size: %d, Maximum: %d).", sequence.size(), maxLength);
		else if (sequence.empty())
			errorMessage = sfmt("Sequence %d is empty.", count);
		addInfo(sequence);
	}
	return errorMessage.empty();
}

bool validate_Interface::validateChemFile() {
	if (!fileExists(file.c_str())) {
		errorMessage = RNA::GetErrorMessage(1);
		return false;
	}
	ifstream input(file.c_str());
	if (!input) {
		errorMessage = RNA::GetErrorMessage(2);
		return false;
	}
	
	string line;
	int linenum = 0;
	stringstream s;
	while (!getline(input, line).fail()) {
		linenum++;
		if (line.find_first_not_of("\r\t ")==string::npos)
			continue; // blank line or just whitespace.

		s.clear();
		s.str(line); // assign line to stringstream.
		int i; 	double val;
		s >> i; 
		if (s.fail()) {
			errorMessage = sfmt("Invalid data on line %d. Expected a numeric index.", linenum);
			return false;
		}
		s >> val;
		if (s.fail()) {
			errorMessage = sfmt("Invalid data on line %d. Expected a decimal value.", linenum);
			return false;
		}
	}
	return true;
}

bool validate_Interface::validateRNAFile() {
	// validate the file using RNA's built-in file reading.
	RNA rna(file.c_str(), knownType, alphabetName.c_str(), allowUnknownBases, SKIP_THERMO);
	if( rna.GetErrorCode() != 0 ) 
		errorMessage = rna.GetFullErrorMessage(); // already ends in "\n"
	else if (rna.GetSequenceLength() > maxLength)
		errorMessage = sfmt("Sequence exceeds the maximum allowed number of bases (Size: %d, Maximum: %d).", rna.GetSequenceLength(), maxLength);
	else {
		addInfo(rna);
		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool validate_Interface::run() {
	string errmsg = "";
	bool success = false;

	// Show a message saying that conversion has started.
	if (!quiet) {
		cout << "Validating file..." << flush;
		if (allowUnknownBases) cout << "(unknown bases ignored)...";
	}

	switch (knownType) {
		case RFT_LIS:
			success = validateOligoList();
			break;
		case RFT_SHAPE:
			success = validateChemFile();
			break;
		case RFT_CON:
			success = validateFoldConstraintFile();
			break;
		case RFT_MULTISEQ:
			success = validateMultipleSequenceFile();
			break;
		default:
			success = validateRNAFile();
			break;
	}

	if (!quiet) cout << "done." << endl;

	if (success) {
		if (!quiet) cout << "Validation successful." << endl;
	} else {
		if (errorMessage[errorMessage.length()-1]!='\n') errorMessage.append("\n");
		cerr << "Validation Error: " << errorMessage; // errmsg should end with \n
	}

	//delete rna;
	return success;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	validate_Interface* runner = new validate_Interface();
	bool success = false;
	if( runner->parse( argc, argv ) ) success = runner->run();
	delete runner;
	return success?0:1;
}
