#if !defined(RnaContainerClass_H)
#define RnaContainerClass_H

#include "../RNA_class/RNA.h"
#include "../src/random.h"
#include "../src/PseudoParser.h"
#include <vector>
#include <iomanip>


using namespace std;

//This class provides an RNA class, and allows an array because there is a default constructor.
class RNA_Container {

public:

	//Default constructor
	RNA_Container();

    
	//Return a pointer to the start of the underlying RNA class.
	RNA *Return_RNA();

	//Allocate a member of an array of RNA classes.  Requires a null terminated c string with the filename, .seq or .FASTA file.
	//Also requires a specificationa as to whether this is RNA or DNA folding.
	//Return an int that indicates whether an error occurred while reading the file.  0 = no error, otherwise there was an error.
	int Allocate(const char *filename, int type, bool IsRNA=true);

	// Allocate a copy of an RNA, using the template ONLY for thermodynamic data.
	int Allocate(const char* filenameOrSequence, RNAInputType type, RNA* copyThermo);

	//Destructor
	~RNA_Container();
	
private:

	//The underlying RNA class.
	RNA *rna;
	

};

RNA_Container* fromSaveFile(string filename, RNA* copyThermo=NULL);
RNA_Container* fromSeqFile(const char* filename, int n, RNA* copyThermo=NULL);
void toFile(string outfilename,RNA_Container* r,int numberofseqs);

#endif // !defined 
