// This file lists some common file types used in RNAstructure
// as well as some utility functions for determining the file type
// based on a file extension.

#if !defined(RNA_FILE_TYPE_H)
#define RNA_FILE_TYPE_H

#include <string>
#include <cstring>

//! An enumeration of some common file types used in RNAstructure.
//! Defining these helps keep numeric values in sync between c++ code 
//! and client code (e.g. Java, python etc) that access the RNAstructure 
//! class library via a native interface (e.g. with SWIG-generated glue code).
enum RNAFileType {
	RFT_UNKNOWN = -1, // Invalid or unknown value. Sometimes used to indicate auto-detect.
	//! A string containing the sequence -- NOT a file name.
	RFT_SEQUENCE=0,
	//! CT Structure file
	RFT_CT=1,
	//! Sequence file (FASTA, SEQ or plain-text sequence)
	RFT_SEQ=2,
	//! Partition function save file (*.pfs)
	RFT_PFS=3,
	//! Folding Save file (*.sav)
	RFT_SAV=4,
	//! Dot-Bracket Notation Structure file (*.dbn, *.dot, *.bracket)
	RFT_DBN=5,
	//! Folding constraints file (*.con)
	RFT_CON=6,
	//! Chemical modification constraints file  (e.g. SHAPE)
	RFT_SHAPE=7,
	//! Oligo list file (*.lis)
	RFT_LIS=8,
	//! File containing multiple sequences (e.g. used by Turbofold etc)
	RFT_MULTISEQ=9,
};

//! This function attempts to determine the type of RNAstructure file based on a filename.
//! It simply determines the extension of the file and passes it to fileTypeFromExt
RNAFileType fileTypeFromPath(const std::string& filePath);

//! This function attempts to determine the type of RNAstructure file based on either a file extension or 
//!    a string indicator/abbreviation of the file type. (for example the 'validate' program allows
//!    the user to specify the type of file valildation to perform by providing one of the abbreviations 
//!    tested for in this function.
//!    e.g.  validate <file> --type 'ct'   or   validate <file> --type 'mseq'
//! It's important to note that a file type/intent cannot always be inferred from the actual file extension.
//!   for example, a FASTA file could contain multiple sequences, so it could represent either a 'FILE_SEQ' (i.e. single-sequence file) 
//!   or a 'FILE_MULTISEQ' (i.e. multiple-sequence file)
RNAFileType fileTypeFromExt(const std::string& ext);

//! Returns a file extension (e.g. "ct", "fasta" etc) that would be 
//! appropriate for the given RNAFileType.
const char* const extFromFileType(RNAFileType type);

#endif // RNA_FILE_TYPE_H