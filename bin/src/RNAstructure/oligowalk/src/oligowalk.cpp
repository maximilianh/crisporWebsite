/*====================================================================================================

oligowalk.cpp : This is the file including main function
oligowalk_test: This is the program to explore the database finding correlation between inhibition efficacy and energy

OligoWalk calculate the binding affinity of oligomer with strucutured RNA
They are revised based on Mathews' code from RNAStructure

oligomer (such as siRNA) is L in lengther
optimal and sub-optimal structures of target can be calculated
partition calculation and stochastic smapling method can also be 
used to predict the breaking energy of target and oligo structure

intermolecular.cpp: handle the options of folding and binding
pclass.cpp: includes the important classed(Pclass for normal partion function, OligoPclass for refilling 
of constrained sequence and scan folding at different site				
Created:Nov. 2005						Modified: 
zhi_lu@urmc.rochester.edu

Jan 2018 - Significantly revised interface - Richard Watson
=======================================================================================================*/
#include <iostream>
#include <iomanip>
#include <cstring>
//#include <cstdarg>
// #include <cctype>
// #include <cstdlib>
#include <stdio.h>
#include "../../RNA_class/RNA.h"
#include "../../src/intermolecular.h"
#include "../../src/ParseCommandLine.h"

using namespace std;

// Open a sub-section of an RNA sequence.
// Returns 0 on error or 1 on success.
int openseq_frac (structure *ct, RNA *rna, int start, int end, bool include_structure) {
	if (end==0)	end = rna->GetSequenceLength(); // end-position=zero means use up to the actual end of the sequence.
	if (start > rna->GetSequenceLength()){
		cerr<<"Sequence is shorter than start-position: "<<start<<"\n";
	} else if (end > rna->GetSequenceLength()){
		cerr<<"Sequence is shorter than end-position: "<<end<<"\n";
	} else if (end < start){
		cerr<<"End-position "<<end<<" is less than start-position " <<start<<".\n";
	} else {
		// indices are valid.
		ct->RemoveConstraints();
		ct->SetThermodynamicDataTable(rna->GetDatatable());
		ct->SetSequence(string(rna->GetSequence()).substr(start-1,end-start+1));
		ct->SetSequenceLabel(rna->GetStructure()->GetSequenceLabel());

		if (include_structure) {
			//also copy over structure information
			for (int structures = 1; structures <= rna->GetStructureNumber(); ++structures) {
				ct->AddStructure();
				for (int position = start; position <= end; ++position) {
					if (rna->GetPair(position,structures) > 0) {
						if (rna->GetPair(position, structures) >= start && rna->GetPair(position, structures <= end)) {
							ct->SetPair(position - start + 1, rna->GetPair(position, structures) - start + 1,structures);
						}

					}

				}
			}

		}

		return 1; // true indicates success.
	}
	return 0; // failed. error message was already output above.
}

// Parse a string to a double, but include special unit-factor 
// abbreviations: mM uM nM pM
bool parseConcentration(const string& input, double& out);
// Parse a 'unit' descriptor which is either an integer (e.g. -6 for uM)
// or an abbreviation: M, mM, uM, nM, or pM
bool parseUnit(const string& input, int& out);

//void show_help(); // show command-line help message.
//void show_param(const char* const flag, const char*const paramName, const char* description);

int main(int argc, char *argv[]) {
	string seqfilename, reportfilename, shapefile;
	int i,j,
		scanStart=1, // aka 'M': start position of scanning on folded region of target (default: nuc 1)
		scanEnd=0, // aka 'N': start position of scanning on folded region of target (default: 0 - indicates the end of the sequence)
		// start and end are used to truncate the sequence.
		start=1, // start position of folding region of target (default: start at nuc 1)
		end=0,   // end position of folding region of target (default: 0 - the end of the sequence)
		foldsize=0, // see foldsize description below.
		unit = 0; // concentration unit factor e.g. 0=Molar, -6=micromolar
	int distance=0;
	int mode=1, // default mode: 1 - break local target structure to bind oligo 
		suboptimal=0,
		length=-1, // -1 indicates not entered (show error message)
		useprefilter=0; // 0 indicates prefilter not used.
	int **table,**numofsubstructures;
	int *TEST,TESTnum=-1; vector<int> TESTon;
	double conc = 1E-6;
	bool isdna=false;
	bool isstructure = false;
	bool scoreit=false;
	bool writesav=false,webserver=false,html=false,include_header=true;
	structure *ct;
	datatable *data,*dhdata,*ddata = NULL;
	rddata *hybriddata = NULL;
	thermo *helixstack;
	//TProgressDialog PD;
	siPREFILTER *prefilter;
	//define the datapath from shell
	const string datapath = getDataPath(DT_RNA);
	Thermodynamics dnaThermo(false);
	// These are additional selection criteria that could be used in the future.
	const bool *mask=NULL; const double asuf=0; const double tofe=0; const double fnnfe=0;
	TProgressDialog* progress=NULL;
	
	FixWindowsExponentFormat();
/*--------------------------------------------------------------------------------------------------
mode 1 - break local target structure to bind oligo 
mode 2 - refold target RNA after oligo binding
mode 3 - no target structure considered

suboptimal 0 - only consider optimal structure
suboptimal 1 - like choice 3,using suboptimal structures,but the whole set from alltrace() function prediction
suboptimal 2 - using partition function considering every possible structure  of target
			     suboptimal 2 can only used with mode 2
suboptimal 3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
suboptimal 4 - using stochastic sampling method to sample 1000 structures

useprefilter 1 - using criteria to prefill functional siRNA; (-- you may not want to type -test )
foldsize >0  - only folding a fragment with size=foldsize+binding length, 
			   which is centered on the siRNA binding region
			   when foldsize>1, only option 2 plus Usesub 2 is the available option
test 0         testing all sites without prefilter score
test 3 2 5 8   testing 3 sites 2,5 and 8 without prefilter score -- do not type -test if you want prefiter score
----------------------------------------------------------------------------------------------------*/
	int error = 0;

	// Create the command line parser
	ParseCommandLine* parser = new ParseCommandLine( "OligoWalk" );
	parser->addParameterDescription( "sequence file", "The name of a sequence file (SEQ, FASTA) or structure file (CT, DBN) containing the input sequence." );
	parser->addParameterDescription( "report file", "The name of a report file to which results will be written." );

	// Add the constraint file option.
	vector<string> concOptions = parser->addFlag(true, "-c -co --conc", "Molar concentration of oligos. "
		"E.g. \"1.5E-6\", \"1.5uM\", or \"0.0000015\".\n"
		"Unit Abbreviations: mM=10^-3 uM=10^-6, nM=10^-9, pM=10^-12\n"
		"This may be used in conjunction with the '--unit' flag: --conc 1.5 --unit -6 (1.5 micromolar)."
		);

	vector<string> modeOptions = parser->addFlag(true, "-m --mode", 
		"Determines how target structure is used:\n"
		"  1=Break Local Structure.\n"
		"  2=Refold target RNA after oligo binding.\n"
		"  3=No target structure considered.");

	vector<string> lengthOptions = parser->addFlag(true, "-l --length", "Length of oligomers for hybridization.");
	vector<string> startPosOptions = parser->addFlag(true, "-st --start", "Start position of folding region of target.");
	vector<string> endPosOptions = parser->addFlag(true, "-en --end", "End position of folding region of target.");
	vector<string> scanStartOptions = parser->addFlag(true, "-ss --from", "Start position of scanning on folded region of target.");
	vector<string> scanEndOptions = parser->addFlag(true, "-se --to", "End position of scanning on folded region of target.");
	vector<string> dnaOptions = parser->addFlag(false, "-d --dna",  "Indicate that the oligomers are DNA (as opposed to RNA).");
	vector<string> suboptimalOptions = parser->addFlag(true, "-s --suboptimal", "Determines suboptimal structure options:\n"
		"  0=Do not consider suboptimal structures.\n"
		"  1=Use AllSub to generate suboptimal structures.\n"
		"  2=Use Partition Function to generate all suboptimal structures.\n"
		"  3=Use a heuristic method for both oligo-free and oligo-bound RNA.\n"
		"  4=Use stochastic sampling to generate 1000 structures.");

	vector<string> unitOptions = parser->addFlag(true, "-u --unit", "Specifies a power-of-ten to multiply the concentration by.\n"
		"For example \"-co 3 -unit -6\" is equal to \"-co 3uM\" or \"-co 3E-6\"");

	vector<string> prefilterOptions = parser->addFlag(true, "-fi --filter", "Whether to use the siRNA prefilter to prefill functional siRNA.\n"
		"0=No Prefilter; 1=Use Prefilter");

	vector<string> scoreOptions = parser->addFlag(false, "-score", "Score the siRNA prefilter.");

	vector<string> structureOptions = parser->addFlag(false, "--structure", "Read a structure rather than a sequence.");

	vector<string> foldsizeOptions = parser->addFlag(true, "-fold", "Only fold a fragment with the specified size (plus the length of \n"
		"the oligomer), which is centered on the binding region.\n"
		"When FOLD_SIZE > 1, MODE (-m) and SUBOPTIMAL (-s) must both be 2.");

	vector<string> distanceOptions = parser->addFlag(true, "-dist", "Limit the maximum distance between nucleotides that can pair.");
	vector<string> shapefileOptions = parser->addFlag(true, "-sh --shape", "Specify a SHAPE data file.");
	vector<string> testOptions = parser->addFlag(true, "-test", "Perform self-tests. The parameter should be a list of space-separated test numbers, e.g.: -test '1 2 5'");
	vector<string> savOptions = parser->addFlag(false, "-write", "Write sav files to save time in test mode.");
	vector<string> webserverOptions = parser->addFlag(false, "-w --webserver", "Use special output for running OligoWalk as a webserver.\n"
		"This implies HTML=true and it sends the header (which lists the \n"
		"parameters that were used in the calulation) to STDOUT instead of \n"
		"including them in the report. It also turns off progress updates.");

	vector<string> noHeaderOptions = parser->addFlag(false, "--no-header", "Do not include a list of the parameters used in the \n"
		"calulation in the output report file.");

	vector<string> htmlOptions = parser->addFlag(false, "--html", "Write the report in HTML mode instead of plain text.");

	// Parse the command line.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( parser->isError() ) return 1;

	seqfilename = parser->getParameter( 1 );
	reportfilename = parser->getParameter( 2 );

	parser->setOptionInteger(modeOptions, mode);
	parser->setOptionInteger(lengthOptions, length);
	parser->setOptionInteger(startPosOptions, start);
	parser->setOptionInteger(endPosOptions, end);
	parser->setOptionInteger(scanStartOptions, scanStart);
	parser->setOptionInteger(scanEndOptions, scanEnd);
	parser->setOptionInteger(suboptimalOptions, suboptimal);
	parser->setOptionInteger(prefilterOptions, useprefilter);
	parser->setOptionInteger(foldsizeOptions, foldsize); //fold substructure size centered on binding
	parser->setOptionInteger(distanceOptions, distance); //Limit the maximum distance between nucleotides that can pair


	if (parser->contains(unitOptions))
		if (!parseUnit(parser->getOptionString(unitOptions, false), unit))
			parser->setErrorSpecialized("The specified 'unit' parameter is invalid.");

	if (parser->contains(concOptions))
		if (!parseConcentration(parser->getOptionString(concOptions, false), conc))
			parser->setErrorSpecialized("The specified 'concentration' parameter is invalid.");

	parser->setOptionString(shapefileOptions, shapefile);

	isdna = parser->contains(dnaOptions);
	scoreit = parser->contains(scoreOptions);
	isstructure = parser->contains(structureOptions);
	writesav = parser->contains(savOptions); //write sav files to save time in test mode
	html = parser->contains(htmlOptions); //write sav files to save time in test mode
	include_header = !parser->contains(noHeaderOptions);

	// Use special output for running OligoWalk as a webserver.
	//   This implies HTML=true and it sends the header information (i.e. output 
	//   listing the parameters that were used in the calulation) to STDOUT instead of
	//   including them in the report.
	//   It also turns off realtime progress updates.
	if (parser->contains(webserverOptions)) {
		webserver=true;
		html=true;
		include_header=false;
	}

	// Read in test options, e.g. -test '1 3 5'
	if (parser->contains(testOptions)) {
		stringstream ss(parser->getOptionString(testOptions, false));
		int num;
		while(!ss.eof()) {
			ss >> num;
			if (ss.fail()&&ss.eof()) 
				break; // no invalid data, but we reached the end after attempting the input (because only whitespace was left after the last input)
			
			if (ss.fail()) {
				cerr << "Invalid test number input: " << ss.str() << endl;
				return 1;
			}
			TESTon.push_back(num);
		}
	}

	
	if(parser->isError()) return 1;

	// if no parsing errors occurred above, perform additional input verification.
    if (seqfilename.empty())
    	parser->setErrorSpecialized("Please enter a valid sequence file name.");

    if (reportfilename.empty())
		parser->setErrorSpecialized("Please enter the output report file name.");

	if (length==-1)
		// -1 indicates not entered (show error message)
		parser->setErrorSpecialized("Please enter the length of hybridization (e.g. -l 19 )");

	if (suboptimal==2 && mode != 2)
		parser->setErrorSpecialized("Partition function (-s 2) can only used with refolding target (-m 2).");

	if (mode<1 || mode>3)
		parser->setErrorSpecialized("Invalid local-structure mode parameter (-m). Valid values are 1-3.");

	if (suboptimal<0 || suboptimal>4)
		parser->setErrorSpecialized("Invalid suboptimal parameter (-s). Valid values are 0-4.");

	if( parser->isError() ) return 1;

	//Check for compatability of parameters:
	//mode 1 requires an input structure and requires suboptimal == 0 or suboptimal == 3.
	if (mode == 1) {
		if (!(suboptimal == 0 || suboptimal == 3)) {
			cerr << "Mode 1 (--mode 1) is only compatible with --suboptimal 0 or --suboptimal 3." << endl;
			return 1;
		}
		if (!isstructure) {
			cerr << "Mode 1 requires reading a structure instead of a sequence, the --structure option." << endl;
			return 1;
		}

	}

	if (unit!=0)
		conc *= pow(10, unit);

	//open the sequence into ct 
	int flag;
	if (isstructure) {
		//open a structure
		flag = FILE_CT;	
	}
	else {
		//open a sequence
		flag = FILE_SEQ;
	}
	RNA rna(seqfilename.c_str(),flag, true); // always open as RNA. (DNA data table is passed separately)

	if (rna.GetErrorCode()!=0) {
		cerr << rna.GetFullErrorMessage();
		return 1;
	}
	ct=new structure;

	// copy fragment of structure 
	if (!openseq_frac(ct, &rna, start, end, isstructure)) {
		cerr << "Failed to open sequence file."<<endl;
		delete ct;
		return 1;
	}

	if (scanEnd==0) scanEnd=ct->GetSequenceLength();
	scanEnd=scanEnd-length+1; // User inputs scanEnd as the end of the region to scan. But since the oligo itself takes up `length` bases, we have to subtract it from scanEnd.

	if (TESTnum==-1)	scoreit=true;
	else				scoreit=false;
	TEST = new int[ct->GetSequenceLength()+1];
	for (i=1;i<=ct->GetSequenceLength();i++)	{
		if (TESTnum==0)		TEST[i]=1;
		else				TEST[i]=0;
		for (j=1;j<=TESTnum;j++) {
			if (i==TESTon[j])	TEST[i]=1;
		}
	}
	//allocate memory of energy tables
	table = new int*[ct->GetSequenceLength() - length + 2];
	for (i = 0; i < ct->GetSequenceLength() - length + 2; i++) {
			table[i] = new int[6];
	}
	//allocate memory of number of suboptimal structures
	numofsubstructures= new int*[ct->GetSequenceLength() - length +2];
	for (i = 0; i < ct->GetSequenceLength() - length + 2; i++)	{
		numofsubstructures[i]= new int [2];
		numofsubstructures[i][0]=0;
		numofsubstructures[i][1]=0;
	}
	
	// RMW 2017-02-27  - Prevent Segfault in intermolecular::olig where GetPair() was called before any structure was allocated.
	if (ct->GetNumberofStructures()==0)
		ct->AddStructure();
	
	helixstack = new thermo(datapath);

	// open the enthalpy files
	data=rna.GetDatatable();
	dhdata=rna.GetEnthalpyTable();
	if (dhdata==NULL) {
		cerr<<"Error reading enthalpy parameters."<<endl;
		return 1;
	}

	if (isdna) {
		// Read the DNA datatable
		if (dnaThermo.ReadThermodynamic()!=0) {
			cerr<<"Error reading DNA thermodynamic parameters."<<endl;
			return 1;
		}
		ddata = dnaThermo.GetDatatable();
		hybriddata = new rddata;
		string stackf = datapath+"/stackdr.dat";
		if (!readrd(hybriddata, stackf)) {
			cerr<<"Error reading hybrid stacking parameters." << endl
			    <<"Please verify the following file:\n\t"
			    <<stackf<<endl;
			return 1;
		}
		helixstack->DH = datapath+"/stackdr.dh"; 
		helixstack->DS = datapath+"/stackdr.ds";
		helixstack->HELIX = datapath+"/helixdr.dat";
	}

	if (!helixstack->read()) {
		cerr<<"Error reading thermo helix stacking parameters." << endl
		    <<"Please verify the following files:\n\t"
		    <<helixstack->DH<<"\n\t"
		    <<helixstack->DS<<"\n\t"
		    <<helixstack->HELIX<<endl;
		return 1;
	}

//define prefilter
	prefilter = new siPREFILTER(*data,*dhdata,useprefilter,scoreit,ct->GetSequenceLength() - length + 2,isdna);

// output header information		
	// Call report() here, but specify that only the header (which lists the parameters used in the calculation) 
	// should be output, not the table body.
	report(cout/*use stdout*/, 
		isdna, mode, ct, length, conc, suboptimal, scanStart, scanEnd, foldsize, distance, 
		table, numofsubstructures, shapefile.c_str(), prefilter,
		mask, asuf, tofe, fnnfe, 
		html, true/*header*/, false/*body*/);

	if (!webserver) {
		cout << "OligoWalk Calculation Started..."<<endl;
		// Create the progress monitor.
		progress = new TProgressDialog();
	}

// run oligowalk
	olig(isdna, mode, ct, length, conc, suboptimal, 
		scanStart, scanEnd, foldsize, distance, 
		table, numofsubstructures, shapefile.c_str(), TEST, writesav, 
		*data, *ddata, helixstack, hybriddata, 
		prefilter, progress);
// output report file
	
	ofstream reportfile(reportfilename.c_str());
	report(reportfile, 
		isdna, mode, ct, length, conc, suboptimal, scanStart, scanEnd, foldsize, distance, 
		table, numofsubstructures, shapefile.c_str(), prefilter,
		mask, asuf, tofe, fnnfe, 
		html, include_header, true);

//delete the dynamic parameters
	for (i = 0; i < ct->GetSequenceLength() - length + 2; i++) {
			delete[] table[i];
		delete[] numofsubstructures[i];
	}
	delete[] table;
	delete[] numofsubstructures;
	delete[] TEST;	
	
	if (isdna) {
		//delete ddata;
		delete hybriddata;
	}
	
	delete helixstack;
	delete ct;
	delete prefilter;
	
	if (progress!=NULL) delete progress;

	if (!webserver) {
		cout << endl << "OligoWalk Calculation Complete." << endl <<
		                "Wrote report file " << reportfilename << endl;
	}
	return 0;
}

// Show command-line help for a flag.
// This displays the flag followed by a parameter-name (which can be NULL
// if the flag does not accept a parameter) followed by a description
// of the flag.
// The description may contain multiple lines separated by 
// linefeed (\n). This character causes the text to be wrapped 
// as well as properly indented.
// void show_param(const char* const flag, const char*const paramName, const char* description) {
// 	const int FLAG_WIDTH=10; // maximum length of flag-names.
// 	const char*const val = is_blank(paramName) ? "" : paramName;
// 	// output flag name and (optional) parameter.
// 	cerr << setw(FLAG_WIDTH) << flag << ' ' << val << endl; //fprintf(stderr, "%10s %-10s\n", flag, val);
// 	char line[200]; const char *start=description, *end;
// 	do {
// 		end = strchr(start, '\n'); // search for linefeed
// 		strncpy(line, start, end==NULL?198:end-start);
// 		if (end!=NULL) line[end-start]='\0'; // since we didn't copy to the end, we have to terminate the string.
// 		cerr << setw(FLAG_WIDTH+1) << "" << line << endl; //fprintf(stderr, "%10s %s\n", "", line);
// 		start=end+1;
// 	} while(end);
//   	cerr << endl;
// }

// void show_help() {
// 	//The parameters in command line were wrong, explain that to user
// 	cerr <<"Usage: OligoWalk [-modifiers]"<<endl;
// 	cerr <<"Example:" << endl << "\tOligoWalk -type r -seq foo.seq -o foo.rep -m 1 -l 18 -s 1 -co 3uM" << endl;
// 	cerr <<"Modifiers:"<< endl;
// 	show_param("-h",NULL,"Show this help message.");
// 	show_param("-seq","INPUT_FILE","Input sequence file.");
// 	show_param("-o","REPORT_FILE","Output report file.");
// 	show_param("-m","MODE","Determines how target structure is used:\n"
// 		"  1=Break Local Structure.\n"
// 		"  2=Refold target RNA after oligo binding.\n"
// 		"  3=No target structure considered.");
// 	show_param("-st","START_POS","Start position of folding region of target.");
// 	show_param("-en","END_POS","End position of folding region of target.");
// 	show_param("-M","SCAN_START","Start position of scanning on folded region of target.");
// 	show_param("-N","SCAN_END","End position of scanning on folded region of target.");
// 	show_param("-l","LENGTH","Length of hybridization.");
// 	show_param("-type","OLIGO_TYPE", "Nucleic Acid Type of oligomers: d=DNA, r=RNA.");
// 	show_param("-s","SUBOPTIMAL","Determines suboptimal structure options:\n"
// 		"  0=Do not consider suboptimal structures.\n"
// 		"  1=Use AllSub to generate suboptimal structures.\n"
// 		"  2=Use Partition Function to generate all suboptimal structures.\n"
// 		"  3=Use a heuristic method for both oligo-free and oligo-bound RNA.\n"
// 		"  4=Use stochastic sampling to generate 1000 structures.");
// 	show_param("-co","CONCENTRATION","Molar concentration of oligo, \n"
// 		"e.g. \"1.5E-6\" or \"1.5uM\" or \"0.0000015\".\n"
// 		"Unit Abbreviations: mM=10^-3 uM=10^-6, nM=10^-9, pM=10^-12");
// 	show_param("-unit","POWER_OF_TEN","Specifies a power-of-ten to multiply the concentration by.\n"
// 		"For example \"-co 3 -unit -6\" is equal to \"-co 3uM\" or \"-co 3E-6\"");
// 	show_param("-fi","USE_PREFILTER","Whether to use the siRNA prefilter to prefill functional siRNA.\n"
// 		"0=No Prefilter; 1=Use Prefilter");
// 	show_param("-score",NULL,"Score the prefilter of siRNA.");
// 	show_param("-fold","FOLD_SIZE","Only fold a fragment with the specified size (plus the length of \n"
// 		"the oligomer), which is centered on the binding region.\n"
// 		"When FOLD_SIZE > 1, MODE (-m) and SUBOPTIMAL (-s) must both be 2.");
// 	show_param("-dist","DISTANCE","Limit the maximum distance between nucleotides that can pair.");
// 	show_param("-shape","SHAPE_FILE","Specify a SHAPE data file.");
// 	show_param("-test",NULL,"Perform self-tests.");
// 	show_param("-write",NULL,"Write sav files to save time in test mode.");
// 	show_param("-webserver",NULL,"Use special output for running OligoWalk as a webserver.\n"
// 		"This implies HTML=true and it sends the header (which lists the \n"
// 		"parameters that were used in the calulation) to STDOUT instead of \n"
// 		"including them in the report. It also turns off progress updates.");
// 	show_param("-no-header",NULL,"Do not include a list of the parameters used in the \n"
// 		"calulation in the output report file.");
// 	show_param("-html",NULL,"Write the report in HTML mode instead of plain text.");
// }

//trim whitespace from a string.
// void ctrim(char *out, size_t max_len, const char *str) {
//   const char *end; size_t out_size = 0;
//   while(isspace(*str)) str++; // Trim leading space
//   if(*str != 0) {
// 	  // Trim trailing space
// 	  end = str + strlen(str) - 1;
// 	  while(end > str && isspace(*end)) end--;
// 	  end++;
// 	  // Set output size to minimum of trimmed string length and buffer size minus 1
// 	  out_size = (end - str) < max_len-1 ? (end - str) : max_len-1;
// 	  memcpy(out, str, out_size); // Copy trimmed string
//   }
//   out[out_size] = 0; //add null terminator
// }

bool parseUnit(const string& input, int& out) {
	string unit(input);
	trim(unit);
	// check for abbreviations first.
	if (unit=="M")       out= 0;
	else if (unit=="mM") out=-3;
	else if (unit=="uM") out=-6;
	else if (unit=="nM") out=-9;
	else if (unit=="pM") out=-12;
	else if (!parseInt(unit, out))
		return false; // cerr <<  "Invalid 'unit' parameter: \"" << unit << "\"" << endl;
	return true; // one of the abbreviations above was successful.
}

// same as parseDbl, but include special unit-factor 
// abbreviations: mM uM nM pM
bool parseConcentration(const string& input, double& out) {
	double num;
	stringstream ss(input);
	ss >> num;
	if (ss.fail()) {
		cerr <<  "ERROR: cannot convert string to number: \"" << input << "\"" << endl;
		return false;
	}
	out = num;

	string unit; int iunit=0;
	ss >> unit;
	if (!trim(unit).empty()) {
		if (!parseUnit(unit, iunit)) {
			cerr <<  "Invalid 'unit' abbreviation in 'concentration' parameter: \"" << unit << "\"" << endl;
			return false;
		}
		cout << "unit: " << iunit << endl;
		if (iunit!=0)
			out *=pow(10, iunit);
	}
	return true;
}
