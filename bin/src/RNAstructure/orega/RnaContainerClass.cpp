#include "RnaContainerClass.h"
#include <iostream>
#include <sstream>

using namespace std;

RNA_Container::RNA_Container() {
	
	
	//set the pointer to RNA to NULL.
	rna = NULL;
	
}


RNA *RNA_Container::Return_RNA() {
	
	return rna;
}

int RNA_Container::Allocate(const char* filenameOrSequence, RNAInputType type, RNA* copyThermo) {
	rna = new RNA(filenameOrSequence, type, copyThermo);
	//Also specify a dummy pair to add the first structure
	rna->SpecifyPair(1,2);
	rna->RemoveBasePair(1);
	return rna->GetErrorCode();
}

int RNA_Container::Allocate(const char* filename, int type, bool IsRNA) {
	if (type==0)//.seq file
		rna = new RNA(filename,2,IsRNA);
	if (type==1)//save file
		rna = new RNA(filename,IsRNA);
	//Also specify a dummy pair to add the first structure
	rna->SpecifyPair(1,2);
	rna->RemoveBasePair(1);
	return rna->GetErrorCode();
}

RNA_Container::~RNA_Container() {
	
	//If rna was allocated, delete it.
	if (rna!=NULL) {
		
		delete rna;
	}
}




vector<string> seqsFromFile(string filename){
	vector<string> seqs;
	std::fstream input(filename.c_str());
	if(!input.is_open()){
		cerr <<"failed to open sequence file"<<std::endl;
		exit(1);
	}
	for(string line; getline(input, line);){
		if(!line.empty()){
			seqs.push_back(line);
		}
	}
	return seqs;
}

//take path to a file with mutated seqs from a previous run
//delimited by newlines
RNA_Container* fromSaveFile(string filename, RNA* copyThermo){
	vector<string> seqs = seqsFromFile(filename);
	RNA_Container* rna = new RNA_Container[2*seqs.size()];
	for(int i=0;i<seqs.size();i++){
		int error = rna[i].Allocate(seqs[i].c_str(),SEQUENCE_STRING,copyThermo);
		if (error!=0){ 
			cerr <<"failed to allocate rna classes from sequence file\n";
			exit(1);
		}
	}
	for (int i = seqs.size(); i < 2*seqs.size(); i++) {
		int error = rna[i].Allocate(seqs[i- seqs.size()].c_str(), SEQUENCE_STRING, copyThermo);
		if (error != 0) {
			cerr << "failed to allocate rna classes from sequence file\n";
			exit(1);
		}
	}
	return rna;
}

//take path to a file in .seq format
RNA_Container* fromSeqFile(const char* filename, int numberofseqs, RNA* copyThermo){
	RNA_Container* rna = new RNA_Container[numberofseqs*2];
	for (int i=0;i<numberofseqs*2;++i) {
		int error = rna[i].Allocate(filename,FILE_SEQ,copyThermo);
		if (error!=0) {
			//An error occurred:
			string errorname = "Error in Allocate ";
			errorname +=  rna[i].Return_RNA()->GetErrorMessage(error);
			cerr << errorname.c_str() << endl;
			exit(1);
		}
	}
	return rna;
}

string nucstostring(char* nucs, int n){
	vector<char> v;
	for(int i=0;i<n;i++)
		v.push_back(nucs[i+1]);
	return string(v.begin(),v.end());
}

void toFile(string outfilename,RNA_Container* r,int numberofseqs){
	std::ofstream f(outfilename.c_str());
	if(!f.is_open()){
		cerr <<"failed to open output file, sequence not written"<<std::endl;
		return;
	}
	for(int i=0;i<numberofseqs;i++){
		RNA* strand = r[i].Return_RNA();
		char* nucs = strand->GetStructure()->nucs;
		const int n = strand->GetSequenceLength();
		f<<nucstostring(nucs,n)<<std::endl;
	}
}
