#include "RNAFileType.h"
#include <iostream>

using namespace std;

// This function attempts to determine the type of RNAstructure file based on a filename.
// It simply determines the extension of the file and passes it to fileTypeFromExt
RNAFileType fileTypeFromPath(const string& filePath) {
	size_t slash = filePath.find_last_of("/\\"); // find last occurance of a slash
	if (slash == string::npos) slash = 0;
	size_t dot = filePath.rfind('.'); // find last dot
	if (dot == string::npos || dot < slash) return RFT_UNKNOWN; // either there is no dot or the last dot is in the directory name. Either way, the file has no extension.
	return fileTypeFromExt(filePath.substr(dot+1));
}

// This function attempts to determine the type of RNAstructure file based on either a file extension or 
//    a string indicator/abbreviation of the file type. (for example the 'validate' program allows
//    the user to specify the type of file valildation to perform by providing one of the abbreviations 
//    tested for in this function.
//    e.g.  validate <file> --type 'ct'   or   validate <file> --type 'mseq'
// It's important to note that a file type/intent cannot always be inferred from the actual file extension.
//   for example, a FASTA file could contain multiple sequences, so it could represent either a 'RFT_SEQ' (i.e. single-sequence file) 
//   or a 'RFT_MULTISEQ' (i.e. multiple-sequence file)
RNAFileType fileTypeFromExt(const string& extension) {
	string ext = extension; // create a copy we can convert to upper-case
	for (int i=0;i<extension.size();i++) ext[i]=toupper(extension[i]); // get upper-case ext
	
	if (ext=="CT") return RFT_CT;
	if (ext=="SEQ"||ext=="FASTA"||ext=="TSEQ") return RFT_SEQ;
	if (ext=="PFS") return RFT_PFS;
	if (ext=="SAV"|ext=="FSF") return RFT_SAV;
	if (ext=="DBN"||ext=="DOT"||ext=="BRACKET") return RFT_DBN;
	if (ext=="CON"||ext=="FCON") return RFT_CON;
	if (ext=="CHEM"||ext=="SHAPE") return RFT_SHAPE;
	if (ext=="OLIS"||ext=="LIS"||ext=="LIST") return RFT_LIS;
	// Re: MSEQ - this is not really an extension that is used. Instead multiple sequences are typically stored in a FASTA file. 
	// So the "MSEQ" indicator is only used is by code that expects multiple sequences or by a program like 'validate' 
	// in which the user can specify the type explicitly.
	if (ext=="MSEQ"||ext=="MFASTA"||ext=="MULTISEQ") return RFT_MULTISEQ; 
	return RFT_UNKNOWN;
}

// Returns a file extension (e.g. "ct", "fasta" etc) that would be 
// appropriate for the given RNAFileType.
const char* const extFromFileType(RNAFileType type) {
	switch(type) {
		case RFT_CT: return "ct";
		case RFT_SEQ: return "fasta";
		case RFT_PFS: return "pfs";
		case RFT_SAV: return "sav";
		case RFT_DBN: return "dbn";
		case RFT_CON: return "con";
		case RFT_SHAPE: return "shape";
		case RFT_LIS: return "lis";
		case RFT_MULTISEQ: return "mseq";
		default: return "unknown";
	}
}