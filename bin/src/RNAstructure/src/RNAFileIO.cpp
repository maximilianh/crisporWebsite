// This file contains common File I/O routines used in RNAstructure.

//Opens a file and attempts to determine whether it is a CT file or a SEQ file.
//Returns FILE_SEQ if the file is likely to be a sequence file, FILE_CT if it is
//a CT or DBN file.
//Returns 0 if the file cannot be opened.
RNAInputType guessInputFileType (const string& inputfile) {
	ifstream in(inputfile);
	if (!in.good()) return 0;
	// 1. Allow leading blank lines and ;-comments in ALL file types.
	// 2. When text is encountered, detect starting characters:
	//    A) '>' suggests either FASTA or DBN (but could be SEQ)
	//    B) A number suggests CT (but could be SEQ)
	//    C) A valid BASE character suggests a text-file (but could be SEQ)
	//    D) Note:: SEQ file can start with *anything* because it is the sequence name.
	//         But unlike the other formats, a SEQ file MUST have a preceding ;-comment.
	//         So if '>' or a BASE is encountered *before* any comment, it must NOT be a SEQ.
	//    A.1) Both FASTA and DBN have BASES on the first line.
	//    A.2) DBN has a second line that should start with a symbol. ",-.([{<"
	//    A.3) FASTA may or may not have a second line of bases.
}
