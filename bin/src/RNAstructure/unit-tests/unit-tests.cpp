// very basic unit testing of the RNAstructure library

#include "unit-tests.h"

using namespace std;

int total_passed = 0, total_failed = 0; // global counters of test results
void post_result(bool passed) { if (passed) total_passed++; else total_failed++; }

// Utility function to copy the basepr vector from a structure.
void getBasePairs(structure &ct, vector<int> &pairs, int structurenumber=1) {
	pairs.resize(ct.GetSequenceLength()+1);
	pairs[0]=0;
	for(unsigned int i=1; i < pairs.size(); i++)
		pairs[i]=ct.GetPair(i, structurenumber);
}

void test_CommonUtils() {
	// Test trim functions
	string s, s1("   Hello   \r\n"), s2("\n"), s3("\n \t \r  "), s4("0");
	cout << "---- TESTING escapeChars in " __FILE__ " ----" << endl;
	escapeChars(s="Hello   \r\n\t"); TEST_EQUAL(s,"Hello   \\r\\n\\t");
	escapeChars(s="Hello   \x0D\x0A\x09"); TEST_EQUAL(s,"Hello   \\r\\n\\t");
	escapeChars(s="Hello   \x01\x02\x11"); TEST_EQUAL(s,"Hello   \\x01\\x02\\x11");

	cout << "---- TESTING escapeChars in " __FILE__ " ----" << endl;
	s=createSafeFilename("    E.coli; 12345adeAg!    "); TEST_EQUAL(s,"E.coli; 12345adeAg!");
	s=createSafeFilename("    E.coli: 12345adeAg! ?/\\ : <> | ***   "); TEST_EQUAL(s,"E.coli_ 12345adeAg! ___ _ __ _ ___");
	s=createSafeFilename("\n\r\tE.coli: 12345adeAg_\n hi"); TEST_EQUAL(s,"E.coli_ 12345adeAg__ hi");
	s=createSafeFilename("\n\r\tE.coli: 12345adeAg_\n hi", ".ct", true); TEST_EQUAL(s,"E.coli__12345adeAg___hi.ct");

	cout << "---- TESTING trimLeft in " __FILE__ " ----" << endl;
	trimLeft(s=s1); TEST_EQUAL(s,"Hello   \r\n");
	trimLeft(s=s2); TEST_EQUAL(s,"");
	trimLeft(s=s3); TEST_EQUAL(s,"");
	trimLeft(s=s4); TEST_EQUAL(s,"0");
	cout << "---- TESTING trimRight in " __FILE__ " ----" << endl;
	trimRight(s=s1); TEST_EQUAL(s,"   Hello");
	trimRight(s=s2); TEST_EQUAL(s,"");
	trimRight(s=s3); TEST_EQUAL(s,"");
	trimRight(s=s4); TEST_EQUAL(s,"0");
	cout << "---- TESTING trim in " __FILE__ " ----" << endl;
	trim(s=s1); TEST_EQUAL(s,"Hello");
	trim(s=s2); TEST_EQUAL(s,"");
	trim(s=s3); TEST_EQUAL(s,"");
	trim(s=s4); TEST_EQUAL(s,"0");

	cout << "---- TESTING toLower and toUpper in " __FILE__ " ----" << endl;
	TEST_EQUAL(toLower(s="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"),"abcdefghijklmnopqrstuvwxyz0123456789");
	TEST_EQUAL(toUpper(s="abcdefghijklmnopqrstuvwxyz0123456789"),"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");
	// test string functions that accept a const string& and return a copy of it, with the desired modification.
	TEST_EQUAL(s=toLower((const string)"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"),"abcdefghijklmnopqrstuvwxyz0123456789");
	TEST_EQUAL(s=toUpper((const string)"abcdefghijklmnopqrstuvwxyz0123456789"),"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

	cout << "---- TESTING parseInt in " __FILE__ " ----" << endl;
	bool b; int i; double d;
	i=7;b=parseInt("12", i); TEST_EQUAL(b,true); TEST_EQUAL(i,12);
	i=7;b=parseInt("   12    ", i); TEST_EQUAL(b,true); TEST_EQUAL(i,12); // whitespace
	i=7;b=parseInt(" 12 blah", i); TEST_EQUAL(b,false); TEST_EQUAL(i,7); // fails due to extra text
	i=7;b=parseInt(" 12 blah", i, false); TEST_EQUAL(b,true); TEST_EQUAL(i,12); // ignore extra text after number
	i=7;b=parseInt("blah 12 blah", i, false); TEST_EQUAL(b,false); TEST_EQUAL(i,7); 
	i=7;b=parseInt("0xFF",i); TEST_EQUAL(b,true); TEST_EQUAL(i,255); // hex
	i=7;b=parseInt("0777",i); TEST_EQUAL(b,true); TEST_EQUAL(i,511); // octal

	cout << "---- TESTING parseDbl in " __FILE__ " ----" << endl;
	d=7;b=parseDbl("-12.75E-3", d); TEST_EQUAL(b,true); TEST_EQUAL(d,-12.75E-3);
	d=7;b=parseDbl("\t\n -12.75E-3 \t\n", d); TEST_EQUAL(b,true); TEST_EQUAL(d,-12.75E-3); // whitespace
	d=7;b=parseDbl(" -12.75E-3blah", d); TEST_EQUAL(b,false); TEST_EQUAL(d,7.0); // fails due to extra text
	d=7;b=parseDbl(" -12.75E-3blah", d, false); TEST_EQUAL(b,true); TEST_EQUAL(d,-12.75E-3); // ignore extra text after number
	d=7;b=parseDbl("blah -12.75E-3 blah", d, false); TEST_EQUAL(b,false); TEST_EQUAL(d,7.0);

	cout << "---- TESTING join in " __FILE__ " ----" << endl;
	vector<int> test;
	TEST_EQUAL(join(test),"");
	test.push_back(0);
	TEST_EQUAL(join(test),"0");
	test.push_back(1);
	TEST_EQUAL(join(test),"0,1");

	cout << "---- TESTING split in " __FILE__ " ----" << endl;
	vector<string> v1 = split("Hello there you ", " ");
	vector<string> vexp; vexp.push_back("Hello"); vexp.push_back("there"); vexp.push_back("you"); vexp.push_back("");
	TEST_EQUAL(v1,vexp);

	cout << "---- TESTING sfmt in " __FILE__ " ----" << endl;
	TEST_EQUAL(s=sfmt("Test %d is %.2f%% %s!", 25, 100.0, "nice"),"Test 25 is 100.00% nice!");
	#define LONG_MSG "A very long message that is bigger than the 256 character initial limit. This message is so big that it might just break the buffers of sfmt unless sfmt is working properly, in which case it shouldn't bother it too much. Ok let's just see if it works. It should, I guess. Let's keep our fingers crossed, shall we? Wonderful!"
	TEST_EQUAL(s=sfmt("%s is %d %s long", LONG_MSG, strlen(LONG_MSG), "characters"),LONG_MSG " is 324 characters long");

	cout << "---- TESTING getFileExt in " __FILE__ " ----" << endl;
	TEST_EQUAL(getFileExt("hello.txt"),"txt");
	TEST_EQUAL(getFileExt(s="file.dir/hello_txt"),"");
	TEST_EQUAL(getFileExt(s="file.dir/hello.txt.banana"),"banana");
	TEST_EQUAL(getFileExt(s="file"),"");
	TEST_EQUAL(getFileExt(s="file\\test.JPEG"),"JPEG");
	TEST_EQUAL(getFileExt(s="file\\test."),"");
}

void test_Thermodynamics() {
	Thermodynamics thermo;
	TEST_EQUAL(thermo.GetEnergyRead(), false);
	TEST_EQUAL(thermo.IsAlphabetRead(), false);
	TEST_EQUAL(thermo.GetDatatable(), (datatable*)NULL);

	TEST_EQUAL(thermo.ReadThermodynamic(), 0); // load datatables and verify non-zero exit code.
	TEST_EQUAL(thermo.GetEnergyRead(), true);
	TEST_EQUAL(thermo.IsAlphabetRead(), true);
	TEST_EQUAL(thermo.GetAlphabetName(), "rna");
	TEST_EQUAL(thermo.GetTemperature(), TEMP_37C);

	// test the copy constructor.
	Thermodynamics t2(thermo);
	TEST_EQUAL(t2.GetEnergyRead(), true);
	TEST_EQUAL(t2.IsAlphabetRead(), true);
	TEST_EQUAL(t2.GetAlphabetName(), "rna");
	TEST_EQUAL(t2.GetTemperature(), TEMP_37C);
	TEST_EQUAL((size_t)t2.GetDatatable(), (size_t)thermo.GetDatatable());

	TEST_EQUAL(thermo.SetTemperature(200), 0); // set temperature and verify non-zero exit code.
	TEST_EQUAL(thermo.GetTemperature(), 200.0);
	TEST_EQUAL(thermo.GetDatatable()->temperature, 200.0);
	
	// t2 should pick up the change in temperature
	TEST_EQUAL(t2.GetDatatable()->temperature, 200.0);
	TEST_EQUAL(t2.GetTemperature(), 200.0);

	TEST_EQUAL(t2.ReloadDataTables(), 0); // reload data tables and verify that the return code is 0 (success)
	TEST_EQUAL(t2.GetTemperature(), 200.0);
	TEST_EQUAL(t2.GetAlphabetName(), "rna");

	TEST_EQUAL(t2.ReadThermodynamic(NULL, "dna", TEMP_37C), 0); // load new data tables and verify that the return code is 0 (success)

	TEST_EQUAL(thermo.GetEnergyRead(), true);	
	TEST_EQUAL(thermo.IsAlphabetRead(), true);
	TEST_EQUAL(thermo.GetAlphabetName(), "dna");
	TEST_EQUAL(thermo.GetTemperature(), TEMP_37C);
}


void test_Structure() {
	cout << "---- TESTING structure class in " __FILE__ " ----" << endl;
	Thermodynamics thermo; thermo.ReadThermodynamic();
	structure &str = *(new structure(1));
	TEST_EQUAL(str.IsAlphabetLoaded(), false);
	TEST_EQUAL(str.IsThermoDataLoaded(), false);
	TEST_EQUAL(str.SetSequence("ATCGATCGAG"), 30);
	str.SetThermodynamicDataTable(thermo.GetDatatable());

	TEST_EQUAL(str.IsAlphabetLoaded(), true);
	TEST_EQUAL(str.IsThermoDataLoaded(), true);	
	TEST_EQUAL(str.SetSequence("ATCGATCGAG"), 0);
	TEST_EQUAL(str.GetSequence(), "ATCGATCGAG");

	TEST_EQUAL(str.GetSequenceLength(), 10);
	str.AddStructure();
	str.SetPair(1, 7);
	str.SetPair(2, 6);
	str.SetPair(3, 9);
	str.SetPair(4, 8);
	str.SetPair(5, 10);

	vector<int> pairs, optimal, crossing, levels;
	str.FindPseudoknots(1, &crossing, &optimal);

	getBasePairs(str, pairs);
	TEST_EQUAL(join(pairs),    "0,7,6,9,8,10,2,1,4,3,5");
	TEST_EQUAL(join(optimal),  "0,7,6,0,0,0,2,1,0,0,0");
	TEST_EQUAL(join(crossing), "0,0,0,9,8,10,0,0,4,3,5");
	
	str.GetPseudoknotRanks(levels);
	TEST_EQUAL(join(levels), "0,1,1,2,2,3,1,1,2,2,3");

	//str.ctout("unit-test.ct");
	const char*const DBN_FILE = "~unit-test.dbn";
	remove(DBN_FILE);
	TEST_EQUAL(str.writedotbracket(DBN_FILE), 0);

	structure &str2 = *(new structure(1));
	str2.SetThermodynamicDataTable(thermo.GetDatatable());
	int result;
	TEST_EQUAL(result=str2.opendbn(DBN_FILE), 0);
	if (result != 0)
		cerr << RNA::GetErrorMessage(result) << str2.GetErrorDetails();
	else {
		getBasePairs(str, pairs);
		TEST_EQUAL(join(pairs), "0,7,6,9,8,10,2,1,4,3,5");
	}
	std::remove(DBN_FILE);
	
	delete &str;
	delete &str2;
}

void test_Misc() {
	string s;
	eraseEnergyLabel(s="  ENERGY = -14.9  Banana"); TEST_EQUAL(s, "Banana");
	eraseEnergyLabel(s="  SCORE = 40  Frog", "SCORE"); TEST_EQUAL(s, "Frog");

	structure*str = new structure(2);
	str->AddStructure();
	str->AddStructure();
	str->SetEnergy(2,36);
	TEST_EQUAL(s=CTComments::None.getComment(str, 1), "");
	TEST_EQUAL(s=CTComments::None.getComment(str, 2), "");
	TEST_EQUAL(s=CTComments::Energy.getComment(str, 1), "");
	TEST_EQUAL(s=CTComments::Energy.getComment(str, 2), string("ENERGY = ") + (conversionfactor==10?"3.6":"0.36") );
	delete str;
}



// Run all tests.
int main(int argc,  char* argv[]) {
	cout << "Unit test results are shown below." << endl << flush;
	cout << "Tests are identified by line-number (in square brackets)." << endl << flush;

	// test basic unit-test features
	TEST_EQUAL("Hello","Hello");
	TEST_EQUAL(string("Hello"),"Hello");
	TEST_EQUAL("Hello",string("Hello"));
	TEST_EQUAL((int)true,1);
	TEST_EQUAL((double)1,1.0);

	test_CommonUtils();
	test_Thermodynamics(); 
	test_Structure();
	test_Misc();

	cout << "Done testing. Passed: " << total_passed << " Failed: " << total_failed << endl;
	return (total_failed == 0 && total_passed != 0) ? 0 : 1;
}