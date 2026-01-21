


//#if !defined(STRUCTURE_CPP)
//#define STRUCTURE_CPP

//#include <stdlib.h>
#include "structure.h"
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "defines.h"
#include "pfunction_math.h"

using namespace std;

#define DBN_BRACKET_SYMBOLS "()<>{}[]AaBbCcDd" // List of characters that represent basepairs in dot-bracket notation. The symbols are listed in pairs: <open1><close1><open2><close2>....
#define DBN_UNPAIRED_SYMBOLS ".-," // List of characters accepted as unpaired bases in dot-bracket notation.

#if defined(WIN32) && !defined(__GNUC__)
	//This is compiling on Windows, where there is no definition for tgamma in cmath
	//Instead, use mathimf.h by including a separate file
	//Note: If compiling on Windows with GCC (e.g. MinGW) we do NOT need to include gamma.h

	#include "gamma.h"
	#define tgamma gammafunction
	//#define isnan _isnan
#endif

#ifdef COUNTING
	#define NO_COUNT(...) __VA_ARGS__._value
#else
	#define NO_COUNT(...) __VA_ARGS__
#endif //COUNTING

//Size the contained vectors in singlestructure for the length of the sequence:
singlestructure::singlestructure(int sequencelength) :
	basepr(sequencelength+1) //Add one to the sequence length so that the arrays can be one-indexed:
{
	//Initialize the energy of a new structure to 0.
	energy = 0;

}

//Establish an operator for comparing singlestructures bases on folding free energy change:
inline bool operator<(const singlestructure &a, const singlestructure &b) {

	return (a.energy < b.energy);
}


//Constructor
structure::structure(int structures)
{

	
	int i;

	allocatestructure(structures);

	allocated = false;

	
	
	intermolecular = false;
	
	templated = false;
	
	min_gu = 0;
	min_g_or_u = 0;
	nneighbors = 0;
	nregion = 0;
	nmicroarray=0;
	

	//initialize values for SHAPE slope and intercept as zero for both double and single stranded restraints.
	SHAPEslope_ss = 0;
	SHAPEintercept_ss = 0;

	SHAPEslope = 0;
	SHAPEintercept = 0;

	data = NULL;  //This is the pointer to the thermodynamic data, stored in a datatable.


	for (i=0;i<maxregions;i++) rnneighbors[i]=0;
	//for (i = 0; i < 100;i++) {
	//	for ( j = 0; j < 25; j++) {
	//		neighbors[i][j]=0;
	//	}
	//}

	stacking = false;
	limitdistance=false;//toogle to true to limit base pairing distance
	maxdistance=600;//default maximum distance between paired nucs
	shaped = false;//by default, a structure does not have SHAPE data associated with it
	SHAPE=SHAPEss=NULL;
	SHAPEss_region = NULL;//by default, the SHAPEss_region does not need to be allocated
	ssoffset = false;//by default, a structure does not have a single stranded offset read from disk
	experimentalPairBonusExists = false;//by default, no pairwise bonuses provided
	EX=NULL;

	constant = NULL;//by default, indicate that the array constant is not allocated


	numofbases = 0;//To start, set the length of the sequence to zero.  This indicates that no sequence has been read.
	distsread = false;

	sequencelabel="\n";//define a default sequence label

	lastErrorDetails = "";

	number_of_sequences = 1;  // number of sequences in sequences_from_alignment

}

// returns pointer to individual sequence i in vector<structure> 
structure* structure::get_individual_sequence(int i) 
{
	return &sequences_from_alignment[i];
}

//Get the label associated with a structure
string structure::GetCtLabel(int structurenumber) const {
	return arrayofstructures[structurenumber-1].ctlabel;
}

//Get an energy associated with a structure
int structure::GetEnergy(int structurenumber) const {

	return arrayofstructures[structurenumber-1].energy;


}

// Get the number of structures stored.
int structure::GetNumberofStructures() const {

	return (int) arrayofstructures.size();

}

//Get the pairing partner of i in structure structurenumber
int structure::GetPair(int i, int structurenumber) const {
	VERIFY_STRUCTURE_INDEX(structurenumber)
	VERIFY_NUC_INDEX(i)
	//Use structurenumber-1 so that the underlying zero indexing works with the 1 indexing used in the rest of the software.
	return arrayofstructures[structurenumber-1].basepr[i];

}

//Get the label associated with the sequence.
string structure::GetSequenceLabel() const {

	return sequencelabel;

}

const char* structure::GetSequence() const {
    if (GetSequenceLength()==0)
        return "";
    return nucs+1; // add one char position, so the first (unused) char is not included.
}

// Remove labels like "ENERGY = ...", "SCORE = ..." or "NED = ..." etc from CT labels.
void structure::RemoveEnergyLabels(const char* customLabel) {
	if (customLabel==NULL) customLabel = "ENERGY";
	for(int i = 1; i <= GetNumberofStructures(); i++) {
		string label(GetCtLabel(i));
		eraseEnergyLabel(label, customLabel);
		SetCtLabel(label, i);
	}
}

//Remove a pair from a specified structure
//Remove the pair for i and the pair for anything to which i was paired
void structure::RemovePair(int i, int structurenumber) {
	VERIFY_STRUCTURE_INDEX(structurenumber)
	VERIFY_NUC_INDEX(i)
	
	//Note that there is no index checking, for the sake of speed.
	//Refer to structurenumber-1 so that the underlying zero indexing works with the 1 indexing used in the rest of the software.
	if (arrayofstructures[structurenumber-1].basepr[i]!=0) {
		arrayofstructures[structurenumber-1].basepr[arrayofstructures[structurenumber-1].basepr[i]]=0;
		arrayofstructures[structurenumber-1].basepr[i] =0;
	}
}


// A stack useful for pseudoknot calculations. It stores intervals [i, j].
// An interval is pushed onto the stack with push(i', j')  (which does NOT affect the members i or j)
// An interval is popped off the stack with pop(), which sets the members i and j.
struct IntervalStack {
	vector<short> stack; unsigned short i, j; unsigned int pos;  /** i and j are the Current values */
	IntervalStack(int capacity = 8) :  stack(capacity), pos(0) {} 
    void push(int i, int j) {
		if (stack.size() < pos+2) stack.resize(pos+2); // ensure there is room for two more elements.
        stack[pos++] = (short)i;
        stack[pos++] = (short)j;
    }
    bool pop() {
        if (pos==0) return false;
        j = stack[--pos]; // pop in reverse order or push (j, then i)
        i = stack[--pos];
        return true;
    }
};

// Determine if the structure has one or more Pseudoknots (crossing bonds).
bool structure::HasPseudoknots(const int structurenumber) const {
	// call the global hasPseudoknots function that analyzes a vector of integers
	return ::hasPseudoknots(arrayofstructures[structurenumber-1].basepr);
}

// Splits the basepairs in the structure into two groups: normal pairs and pseudoknot pairs
void structure::FindPseudoknots(const int structurenumber, vector<int> *pseudoknotPairs /*default NULL*/, vector<int> *normalPairs /*default NULL*/) const {
	// call the global findPseudoknots function that analyzes a vector of integers
	return ::findPseudoknots(arrayofstructures[structurenumber-1].basepr, pseudoknotPairs, normalPairs);
}

// Fills the results output vector with the "pseudoknot rank" of each base pair in this structure.
void structure::GetPseudoknotRanks(vector<int> &results, const int structurenumber) const {
	const vector<int> &pairs = arrayofstructures[structurenumber-1].basepr;
	if (results.size() < pairs.size()) 
		results.resize(pairs.size());

	// `list` first contains ALL basepairs.  We then call findPseudoknots to fill it with just the pseudoknots.
	//  The pseudoknots found in each round become the "normal pairs" for each subsequent round of finding pseudoknots.
	vector<int> list(pairs.size());
	copy(pairs.begin(), pairs.end(), list.begin());

	for(unsigned int i=0;i<results.size();i++)
		results[i]=list[i]==0 ? 0 : 1;
	// results now has a 0 for all positions without basepairs and a 1 at each position with a pair (regardless of whether it is crossing or not)
	
	while(hasPseudoknots(list)) {
		findPseudoknots(list, &list);  // determine which are the cossing (pseudoknot) bonds. 
		// we've already accounted for the optimal set of non-crossing bonds.
		// so now just increment results[i] for each bond in the set of crossing bonds.
		for(unsigned int i=0;i<results.size();i++)
			if (list[i]!=0) results[i]++;
	}
}

//! Remove all basepairs that represent PseudoKnots (i.e. all bonds with a "crossing level" of 2 or higher).
//! This method is sequence and energy agnostic. 
void structure::BreakPseudoknots(int structurenumber, vector<int> *brokenPairs) {
	if (!HasPseudoknots(structurenumber)) return;
	// pass the singlestructure::basepr vector in as the OUTPUT vector for the 
	// NON-pseudoknot bonds. This will overwrite the basepr array with ONLY the normal bonds.
	// thereby removing the pseudoknots.
	FindPseudoknots(structurenumber, NULL, &arrayofstructures[structurenumber-1].basepr);
}

//Set the label for structure structurenumber
void structure::SetCtLabel(const string &label, const int structurenumber) {

	arrayofstructures[structurenumber-1].ctlabel = label;

}

//A version of lable setting that uses char pointer:
void structure::SetCtLabel(const char *label, const int structurenumber) {

	arrayofstructures[structurenumber-1].ctlabel = label;

}

//Assign the energy for a structure:
void structure::SetEnergy(int structurenumber, int energy) {
	
	arrayofstructures[structurenumber-1].energy = energy;

}


//Specify a basepair in a specific structure:
void structure::SetPair(int i, int j, int structurenumber) {

	arrayofstructures[structurenumber-1].basepr[i] = j;
	arrayofstructures[structurenumber-1].basepr[j] = i;

}


//Set a sequence label.
void structure::SetSequenceLabel(const string& label) {

	sequencelabel = label;

}


//Set base pairing cutoff
void structure::SetBasePairingCutoff(float i) {

	bp_cutoff = i;

}



void structure::FillPairingMatrix() {

	float cutoff = bp_cutoff;
	int n = number_of_sequences;
	int seq_len = get_individual_sequence(0)->GetSequenceLength();

	// i, j positions are 1-indexed
	pairing_matrix.resize(2*seq_len + 1, vector<bool>(2*seq_len + 1, false));

	for (int i = 1; i < seq_len + 1; i++) {
		for (int j = 1; j < seq_len + 1; j++) {
			bool flag = false;
			float count_pair = 0;			
			float count_gap_gap = 0;
			
			for (int p = 0; p < n; p++) {
				structure* ct_p = get_individual_sequence(p);
				if (ct_p->GetThermodynamicDataTable()->can_pair(i, j, ct_p->numseq)) {
					count_pair += 1;
				}
				// for RNAalifold style base pairing criteria
		//		if (ct_p->GetSequence()[i-1] == '-' && ct_p->GetSequence()[j-1] == '-') {
		//			count_gap_gap += 1;
		//		}
			}

			
			if ((count_pair / n) >= cutoff) {
				flag = true;
			}
			
			/*
			// base condition for RNAalifold
			if ((count_pair / (n - count_gap_gap)) >= bp_cutoff) {
				flag = true;
			}
			*/

			if (flag) {
				pairing_matrix[i][j] = true;
				pairing_matrix[i][j + seq_len] = true;
				pairing_matrix[i + seq_len][j] = true;
				pairing_matrix[i + seq_len][j+seq_len] = true;
			}
		}
	}
	/*
	for (int i = 1; i < 2*seq_len + 1; i++) {
		for (int j = 1; j < 2*seq_len + 1; j++) {
			cout << pairing_matrix[i][j] << "  ";
		}
		cout << endl;
	}
	*/
}


void structure::FillGapMatrix() 
{
	float cutoff = 0.5;
	
	int seq_len = get_individual_sequence(0)->GetSequenceLength();
	gap_matrix.resize(2*seq_len + 1, false);
	vector<int> gap_count(2*seq_len + 1, 0);

	for (int j = 0; j < number_of_sequences; j++) {
		const char* seq_j = get_individual_sequence(j)->GetSequence();		
		for (int i = 0; i < seq_len; i++) {
				if (seq_j[i] == '-') {
					gap_count[i+1] += 1;
				}
		}
	}

	for (int i = 1; i < seq_len + 1; i++) {
		if (gap_count[i] >= (number_of_sequences * cutoff)) {
			gap_matrix[i] = true;
			gap_matrix[seq_len + i] = true;
		}
	}
//	cout << gap_matrix << endl;
}



// Set the sequence for this structure.
// This includes the following operations:
//   - Allocate space for the bases.
//   - Verify each base exists in the alphabet.
//   - Setup the arrays nucs, numseq, and hnumber.
//   - Check to see if any nucleotide needs to be single-stranded and call AddSingle for them.
// The sequence can contain whitespace, which is ignored.
// If an error occurs (such as the data-table has not yet been loaded) the return value 
// will be an RNA error code (i.e. it corresponds to one of those defined in RNA::GetErrorMessage)
// Additionally lastErrorDetails may be set (which can be queried with GetErrorDetails())
int structure::SetSequence(const string& sequence) {
	if (!IsAlphabetLoaded()) return 30;
	int count = 0; string::const_iterator it, ite=sequence.end();
	//iterate through the string and count the number of non-whitespace characters
	for(it=sequence.begin();it!=ite;it++) if (!::isspace(*it)) count++; // skip whitespace

	allocate(count);
	count=0;
	nucs[0]='\0';
	hnumber[0]=0;
	// Setqup the CT sequence (i.e. arrays nucs, numseq, numofbases)
	for (unsigned int i=0;i<sequence.size();i++) {
		// record each nucleotide
		char base=sequence[i]; // sequence indices are 0-based while the CT arrays (and pos) are 1-based, so subtract 1 from pos.
		if (::isspace(base)) continue;
		nucs[++count]=base;
		int num = data->basetonum(base); 
		if (num==-1) {
			SetErrorDetails(sfmt("Invalid nucleobase %c at position %i.", base, i+1));
			return 28; // error reading sequence.
		}
		numseq[count]=num;
		hnumber[count]=count;
		//Check to see if any nucleotide needs to be single-stranded:
		//scan through the notpairing list for nucleotides that cannot pair
		for (unsigned int j=0;j<data->not_pairing.size();++j) {
			if (base==data->not_pairing[j]) {
				AddSingle(count); 
				break;
			}
		}
	}
	nucs[numofbases+1]='\0'; // add terminating null character.
	
	return 0; //success
}

int structure::SetSequenceWithoutGaps(const string& sequence) {
	string seq_without_gaps = sequence;
	seq_without_gaps.erase(remove(seq_without_gaps.begin(), seq_without_gaps.end(), '-'), seq_without_gaps.end());
	sequence_without_gaps->SetSequenceFast(seq_without_gaps + seq_without_gaps);
	return 0;
}

int structure::SetSequenceFast(const string& sequence) {
	int N = sequence.size();
	allocate(N);
	nucs[0] = '\0';
	hnumber[0] = 0;
	int count = 0;
	for (unsigned int i = 0; i < N; i++) {
		char base = sequence[i]; 
		nucs[++count] = base;
		int num = data->basetonum(base);
		numseq[count] = num;
		hnumber[count] = count;
	}
	nucs[numofbases + 1] = '\0'; 
	return 0; 
}


int structure::SetSequenceFast(const string& sequence, const short* numseq1) {
	int N = sequence.size();
	nucs[0] = '\0';
	hnumber[0] = 0;
	int count = 0;
	for (unsigned int i = 0; i < N; i++) {
		char base = sequence[i];
		nucs[++count] = base;
		numseq[count] = numseq1[count];
		hnumber[count] = count;
	}
	nucs[numofbases + 1] = '\0';
	return 0;
}

int structure::SetTmpSequence	(const string& sequence) {
	int N = sequence.size();
	nucs[0] = '\0';
	hnumber[0] = 0;
	int count = 0;
	for (unsigned int i = 0; i < N; i++) {
		char base = sequence[i];
		nucs[++count] = base;
		hnumber[count] = count;
	}
	nucs[numofbases + 1] = '\0';
	return 0;
}



// This function populates no_of_gaps_matrix
int structure::FillNoOfGapsMatrix() {
	
	string str = GetSequence();
	str = str + str;
	str = "x" + str;
	int seq_len = GetSequenceLength();

	no_of_gaps_matrix.resize(2*seq_len + 1, vector<int>(2*seq_len + 1, 0));

	for (int i = 1; i < 2*seq_len+1; i++) {
			for (int j = 1; j < 2*seq_len+1; j++) 
			{
				int count = 0;
				if (i == j) 
					no_of_gaps_matrix[i][j] = 0; 
				else if (i < j) { 
						for (int p = i+1; p < j; p++) {
							if (str[p] == '-') count += 1;
						}
					no_of_gaps_matrix[i][j] = count; 
				}
				else if (i > j) { 
					no_of_gaps_matrix[i][j] = no_of_gaps_matrix[j][i];
				}	
			}
	}
	return 0;
}

// This function populates gapped_to_ungapped_map
int structure::FillGappedtoUngappedMap() {
	
//	-AGAA--A-G
//	1123455566
//	gap will point to index of closest nucleotide to the right side
//	123456
//  AGAAAG

// -AGAA--A--
// 1123455566


	string seq = GetSequence();

	// making nucleotides 1-indexed
	//seq = "x" + seq;
	seq = seq + seq;
	seq = "x" + seq;
//	cout << seq << endl;

	int seq_length = GetSequenceLength();
	
	gapped_to_ungapped_map.resize(2*seq_length + 1, 0);
		
	for (int i = 1; i < 2*seq_length + 1; i++) {
		if		(i == 1)				     gapped_to_ungapped_map[i] = 1;
		else if (i > 1 && seq[i - 1] == '-') gapped_to_ungapped_map[i] = gapped_to_ungapped_map[i - 1];
		else if (i > 1 && seq[i - 1] != '-') gapped_to_ungapped_map[i] = gapped_to_ungapped_map[i - 1] + 1;
	}

	return 0;
}



//This allocates space in an array that is used for folding with phylogenetic data.
//	tem == template for allowed and disallowed pairs
void structure::allocatetem()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   tem = new bool *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	tem[i] = new bool [i+1];
   }
   templated = true;

   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i;j<=numofbases;j++) {
    		tem[j][i] = true;
		}
   }

}

// Add a nucleotide to the list of those that must pair.
void structure::AddDouble(int i) {

	doublestranded.push_back(i);
}

//Add a pair of nucleotides to the list of those not allowed to form.
void structure::AddForbiddenPair(int i, int j) {

	forbid5.push_back(i);
	forbid3.push_back(j);
}

//Add a nucleotide to the list of Us in GU pairs.
void structure::AddGUPair(int i) {

	GUpair.push_back(i);

}

//Add a nucleotide to the list of those accessible to traditional chemical modification.
void structure::AddModified(int i) {

	modified.push_back(i);

}

//Add a pair of nucleotides to the list of those that must form.
void structure::AddPair(int i, int j) {

	pair5.push_back(i);
	pair3.push_back(j);

}

//Add a nucleotide to the list of those not able to pair. 
void structure::AddSingle(int i) {

	singlestranded.push_back(i);

}

void structure::AddDomain(int i, int j){
	domains5.push_back(i);
	domains3.push_back(j);
}

//! Get a nucleotide that must be base paired.
int structure::GetDouble(int i) {

	return doublestranded[i];

}

//Get a nucleotide that must not be in a specific pair.
int structure::GetForbiddenPair5(int i) {

	return forbid5[i];

}


//Get a nucleotide that must not be in a specific pair.
int structure::GetForbiddenPair3(int i) {


	return forbid3[i];

}

//Get a nucleotide that must be a U in a GU pair.
int structure::GetGUpair(int i) {

	return GUpair[i];

}

//Get a nucleotide that is accessible to chemical modification.
int structure::GetModified(int i) {

	return modified[i];

}

//Get the number of nucleotides forced to be double-stranded.
int structure::GetNumberofDoubles() {

	return doublestranded.size();

}

//Get the number of pairs that are forbidden.
int structure::GetNumberofForbiddenPairs() {

	return forbid5.size();

}
		
//Get the number of Us forced to be in GU pairs.
int structure::GetNumberofGU() {

	return GUpair.size();

}
		
//Get the number of nucleotides that are accessible to chemical modification.
int structure::GetNumberofModified() {

	return modified.size();

}
		
//Get the number of nucleotides forced to be single-stranded.
int structure::GetNumberofSingles() {

	return singlestranded.size();

}
		
//Get the number of pairs that are constrained to occur.
int structure::GetNumberofPairs() {

	return pair5.size();

}

//Get the number of folding domains.
int structure::GetNumberofDomains() {

	return domains5.size();

}

//Get a nucleotide that must be in a specific pair.
int structure::GetPair5(int i) {

	return pair5[i];

}

//Get a nucleotide that must be in a specific pair.
int structure::GetPair3(int i) 
{
	return pair3[i];
}


//Get a nucleotide that defines a domain.
int structure::GetDomain5(int i) {

	return domains5[i];

}

//Get a nucleotide that defines a domain.
int structure::GetDomain3(int i) {

	return domains3[i];

}


//Get a nucleotide that must be single stranded.
int structure::GetSingle(int i) {

	return singlestranded[i];

}

//Get number of sequences to check if there is multiple sequence alignment 
int structure::GetNumberofSequences() {

	return number_of_sequences;

}

float structure::GetBasePairingCutoff() {

	return bp_cutoff;

}

//Check if i, j positions can base pair; different from can_pair in datatable struct
bool structure::can_pair(int i, int j, structure* ct) {
	if (ct->GetNumberofSequences() > 1) { return ct->pairing_matrix[i][j]; }
	else { return ct->data->can_pair(i, j, ct->numseq); }
}

//Check if base pairing between i, j is flanked by non-pairing nucleotides (not gaps)
bool structure::can_pair_isolated(int i, int j, structure* ct) {

	int ix = i, iy = i;
	int jx = j, jy = j;

	bool before = false, after = false;
	int seq_len = ct->GetSequenceLength();

	if (i > 1) {
		for (int m = i - 1; m > 0; m--) {
			if (gap_matrix[m] != true) {
				ix = m;
				break;
			}
		}
	}

	if (i < seq_len) {
		for (int m = i + 1; m < ct->GetSequenceLength() + 1; m++) {
			if (gap_matrix[m] != true) {
				iy = m;
				break;
			}
		}
	}

	if (j > 1) {
		for (int m = j - 1; m > 0; m--) {
			if (gap_matrix[m] != true) {
				jx = m;
				break;
			}
		}
	}

	if (j < seq_len) {
		for (int m = j + 1; m < ct->GetSequenceLength() + 1; m++) {
			if (gap_matrix[m] != true) {
				jy = m;
				break;
			}
		}
	}
	

	if (can_pair(ix, jy, ct)) before = true;
	if (ix == i || jx == j) before = false;

	if (can_pair(iy, jx, ct)) before = true;
	if (iy == i || jy == j) after = false;

	if (before || after) { return true; }
	else return false;
}


//Remove specified number of single stranded constraints from the end
//of the singlestranded vector
void structure::RemoveSingleStrandConstraints(int number) {
    
    for (int i=0; i<number; i++) singlestranded.pop_back();

}


//Reset, i.e. remove, all constraints
void structure::RemoveConstraints() {

	doublestranded.clear(); 
	singlestranded.clear();
	GUpair.clear();
	modified.clear();
	pair5.clear();
	pair3.clear();
	forbid5.clear();
	forbid3.clear();

}


//<<<<<<< structure.cpp
//Set a maximum distance between pairing partners
void structure::SetPairingDistance(int Maxdistance) {
	//store the distance
	maxdistance = Maxdistance;
//indicate that the distance is limited
	limitdistance = true;


}

// If the label contains "ENERGY = <NUMBER>" at the start of the string, it will be removed,
// leaving the remaining text intact. Leading whitespace is also removed.
// A custom label can be specified instead of "ENERGY", in which case, the full text to be removed 
// is "<customLabel> = <NUMBER>".
void eraseEnergyLabel(string &text, const char* const label /* default "ENERGY" */) {
	trimLeft(text); // defined in util.h
	unsigned int length=strlen(label);
	if (text.size() < length+3) return; // must contain at least "ENERGY = "
	 // verify that all the characters in ENERGY are present at the beginning of the text
	if (strncmp(text.c_str(), label, length) != 0) return;
	if (strncmp(text.c_str()+length, " = ", 3) != 0) return;
	trimLeft(text);
	// Now look for the end of the numeric energy value
	string::iterator it = text.begin()+length+3, end = text.end();
	while(it!=end&&!::isspace(*it)) it++; // read until we hit a space or the end of the string.
	text.erase(text.begin(), it);
	trimLeft(text);
}

//outputs a ct file (connection table)
//Provide a pointer to cstring with the filename
//	if append is true, the ct table is appended an existing file, otherwise a new file is created
//	append is false by default
//By default, the columns for indicies are only 5 characters wide, so this is a problem for sequences > 9,999 nucs.
//Now, when sequences are >9,999 nucs, columns are 6 characters wide.  (The code is written yet for sequences > 99,999 nucs.)
int structure::ctout(const char * const ctoutfile, const bool append, CTCommentProvider &commentProvider) const {
	int count,i;//length
	char line[2*ctheaderlength];//base[2]

	// if the filename (ctoutfile) is "-" then write to STDOUT.
	std::ostream output(std::cout.rdbuf());
	std::ofstream out_file; //only used if STDOUT is not used.
	if (!isStdIoFile(ctoutfile)) {
		//Open the file for writing:
		out_file.open(ctoutfile,append ? std::ios_base::app : std::ios_base::trunc);
		//Make sure the file opened.  If not, return an error indicator.
		if (!out_file.is_open()) {
			perror("Error opening ct output file");
			// not allowed because the method is marked const: SetErrorDetails(sfmt("Failed to open the file '%s'. Please verify the path and its permissions.", ctoutfile));
			return 2; // Error opening file
		}
		output.rdbuf(out_file.rdbuf());
	}

	for (count=1;count<=(GetNumberofStructures());count++) {
		line[0]='\0'; //initialize to ""
		if (numofbases>9999) sprintf(line,"%6i",numofbases);
		else sprintf(line,"%5i",numofbases);
   		strcat(line,"  ");
		// Get a comment like "ENERGY = ..." for this structure.
		string comment = commentProvider.getComment(this, count);
		if (!comment.empty()) { 
			strcat(line,comment.c_str());
			strcat(line,"  ");
		}
		string label = GetCtLabel(count);
		//cout << "Original CT Label: '" << label;
		trim(label); // remove leading and trailing whitepace (including trailing newline)
		strcat(line,label.c_str());
		//cout << "'  New Title-Line: '" << line << "'" <<endl;
		output << line << endl; // write out header line

		for (i=1;i<numofbases;i++) {
			//if (ct->stacking) {
			//	if (ct->numofbases>9999) sprintf(line,"%6i%2c%8i%6i%6i%6i%6i\n",
			//		i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->GetPair(basepr[count][i+ct->numofbases]);
			//	else sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
			//		i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
			//}
			//else {
			if (numofbases>9999)sprintf(line,"%6i%2c%8i%6i%6i%6i",
				i,nucs[i],(i-1),(i+1),GetPair(i,count),hnumber[i]);
			else sprintf(line,"%5i%2c%8i%5i%5i%5i",
				i,nucs[i],(i-1),(i+1),GetPair(i,count),hnumber[i]);
			//}
			output << line << endl;
		}

		//last nucleotide not connected--
		i = numofbases;
		//if (ct->stacking) {
		//	sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
		//		i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
		//}
		//else {
			if (numofbases>9999) sprintf(line,"%6i%2c%8i%6i%6i%6i",
				i,nucs[i],(i-1),0,GetPair(i,count),hnumber[i]); 
			
			else sprintf(line,"%5i%2c%8i%5i%5i%5i",
				i,nucs[i],(i-1),0,GetPair(i,count),hnumber[i]);
		//}
		output << line << endl;
	}
	if (output.fail())
		return 2; // Error opening file
	return 0;
}

//Open a Dot-Bracket (DBN) Structure File from disk and store all the information.
//Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable used to set the datatable pointer
//before this function is called.
// This can parse dot-bracket files in any of the following formats (see DotBracketFormat): 
//  DBN_FMT_SINGLE_TITLE,  DBN_FMT_SIDE_TITLES, or DBN_FMT_MULTI_TITLE
//  Note: DBN_FMT_MULTI_TITLE_AND_SEQ is NOT supported because a structure class can only contain a single sequence.
int structure::opendbn(const char *bracketFile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	// Create a variable that handles errors.
	long linenumber = 0;
	string line;
	string sequence;
	size_t start = string::npos;

	std::istream input(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(bracketFile)) { 
		if (!fileExists(bracketFile)) {
			SetErrorDetails(sfmt("The path '%s' is invalid or does not exist.", bracketFile));
			return 1; // file not found.
		}
		//Open the file for reading:
		file_in.open(bracketFile);
		//Make sure the file opened.  If not, return an error indicator.
		if (!file_in.is_open()) {
			SetErrorDetails(sfmt("Failed to open the file '%s'. Please verify the file and its permissions.", bracketFile));
			return 2; // Error opening file
		}
		input.rdbuf(file_in.rdbuf());
	}

	//scan for first relevant line. Ignore lines that are empty, contain only whitespace, or start with a semi-colon (;).
	#define WHITESPACE " \t\r"
	while (getline(input,line)) {
		++linenumber;
		if (line.length() == 0 || line[0] == '\r') continue; // line is empty
		if (line[0] == ';')
			continue;
		if (string::npos != (start = line.find_first_not_of(WHITESPACE))) // if the line has any non-whitespace character, exit the loop
			break;
	}
	if (input.bad()) {
		//Some kind of read error occurred.
		SetErrorDetails(sfmt("An error occured while reading the file '%s'.", bracketFile));
		return 2;
	}

	// Now past the comments. See if a valid line was found.
	if (start == string::npos || line[start] != '>') {
		//The end of the file was reached before finding the end of the comments, return with error
		SetErrorDetails("The dot-bracket file did not contain a sequence label starting with '>'.");
		return 29;
	}

	line.erase(0, start+1); // remove all characters up to and including start. (they were either whitespace or '>')
	SetSequenceLabel(line);

	// Read the next line, which is the sequence.
	if (getline(input,line).fail()) {
		SetErrorDetails("The dot-bracket file did not contain a sequence.");
		return 29;
	}
	linenumber++;

	//Process each character in the sequence
	for (unsigned int i=0;i<line.length();++i) {
		if (line[i] < 33) continue; // skip whitepace.
		// If the character is valid, add it. Otherwise return with an error.
		if (-1==data->basetonum(line[i])) {
			//This nucleotide wasn't recognized
			// report this error via the lastErrorDetails mechanism
			SetErrorDetails(sfmt("Invalid nucleobase '%c' at line %li column %i.", line[i], linenumber, i+1));
			return 28;
		} else {
			//Put the nucleotide in the sequence
			sequence+=line[i];
		}
	}//end of iteration over sequence characters

	if (sequence.size()==0) {
		SetErrorDetails("The file did not contain any nucleotides.");
		return 29;
	}

	SetSequence(sequence);
	
	// Read subsequent lines. Each is a new structure. e.g.: 
	//   >SequenceLabel
	//   AGUGACU...
	//   ((..(((...  (Str1 comment)
	//   ..(((((...  (Str2 comment)
	//   ((((..((..  (Str3 comment)
	int count=0; // number of structures
	string label("");
	while (!getline(input,line).fail()) {
		++linenumber;
		if (line.find_first_not_of(WHITESPACE)==string::npos) continue; // line is empty or whitespace

		if (line[0] == ';')
			continue; // a comment

		if (line[0] == '>') { // A structure label.
			label = line.substr(1);
			continue;
		}

		//Add a structure to which pairs can be placed.
		count++;
		AddStructure();
		SetCtLabel(label.empty()?GetSequenceLabel():label, count);
		label.clear();

		//Now parse the brackets.
		// Multiple bracket types can be used. Pseudoknots can be encoded by using different types of brackets: e.g.  ..<<<..(...>>>..)..
		struct Bracket { unsigned int id; unsigned int pos; };
		#define INVALID_BRACKET_ID (unsigned int)(-1)  // get the max unsigned int.
		int bracketCount = 0;
		Bracket *brackets= new Bracket[line.size()];
		
		//DEBUG: cout << "Starting bracket parsing." << endl;
		//DEBUG: cout << "Seq:" << sequence << endl;
		//DEBUG: cout << "Brk:" << line << endl;

		for(unsigned int i=0; i<line.size();i++) {
			char c = line[i];
			if (c=='\r'||findchr(DBN_UNPAIRED_SYMBOLS, c)!=string::npos) continue; // skip if it is a CR or a dot-symbol (. or -).
			if (c==' '||c=='\t') {
				//whitespace ends the structure. the remainder is a sequence label
				label=line.substr(i+1);
				cout << "Found side label: '" << label << "'" << endl;
				trim(label);
				if (!label.empty())
					SetCtLabel(label, count);
				break; // quit parsing the structure since the side-label has been found.
			}

			 size_t id = findchr(DBN_BRACKET_SYMBOLS, c);
			if (id==string::npos) {
				SetErrorDetails(sfmt("Invalid character in dot-bracket file: '%c' at line %i column %i.", c, linenumber, i+1));
				return 29;
			}
			if ((id&1)==0) {
				// it is a left-side (opening) bracket. Just add it to the list.
				//DEBUG: cout << "LeftBracket " << c << endl;
				brackets[bracketCount].id=id+1; // +1 because we store the closing bracket id, which is one greater than the opening bracket.
				brackets[bracketCount].pos=i;
				bracketCount++;
			} else {
				// it is a right-side (closing) bracket.
				// find its match and add the pair.
				//DEBUG: cout << "RightBracket " << c << endl;
				int found = -1;
				for(int j=bracketCount-1;j>-1;j--) { // look backwards through the pairs vector, searching for the matching bracket.
					if (brackets[j].id==id) {
						found = j; break;
					}
				}
				if (found==-1) {
					SetErrorDetails(sfmt("Unmatched bracket in dot-bracket file: '%c' at line %li column %i.", c, linenumber, i+1));
					return 29;
				}
				//DEBUG: cout << "SetPair " << (brackets[found].pos+1) << "," << (i+1) << endl;
				// Register the pair that was found.
				SetPair(brackets[found].pos+1,i+1, count);  // +1 because pos is the 1-based nucleotide position, while i is the 0-based string index.

				// If we are at the end of the array, decrement the count
				if (found==bracketCount-1) {
					bracketCount--;
					while(bracketCount>0 && brackets[bracketCount].id==INVALID_BRACKET_ID) // this loop removes all invalidated brackets that might have been introduced by a pseudoknot (i.e. when pairs of different-typed-brackets cross each other).
						bracketCount--;
				} else
					// this bracket is NOT at the end. It must be part of a pseudoknot. Just invalidate it for now.
					brackets[found].id = INVALID_BRACKET_ID;
			}
		} // end of parsing brackets
		delete[] brackets;
	} // end of getLine loop
	
	if (GetNumberofStructures()==0) {
		SetErrorDetails("The dot-bracket file did not contain any structures.");
		return 29;
	}
	return 0; // success
}

// Open a CT File from disk and store all the information.
// Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable 
// used to set the datatable pointer before this function is called.
// int structure::openct(const char *ctfile) {
// 	// if ctfile is "-", read from stdin.
// 	if (isStdIoFile(ctfile))
// 		return openct(std::cin, "STDIN");
// 	else {
// 		if (!fileExists(ctfile)) return 1;
// 		std::ifstream fin(ctfile);
// 		//Open the file for reading:
// 		fin.open(ctfile);
// 		if (!fin.is_open()) return 2;
// 		return openct(fin, ctfile);
// 		//Make sure the file opened.  If not, return an error indicator.
// 	}
// }

//Open a CT File from disk and store all the information.
//Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable used to set the datatable pointer
//before this function is called.
int structure::openct(const char *ctfile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	const int MAX_CT_BASES = 20000; // arbitrary number used mostly to prevent huge memory allocations if the CT file has text/garabage at the top that looks like a large number.
	int count, basepair;
	char base[2];
	string sequenceLabel;
	
	std::istream in(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(ctfile)) { 
		if (!fileExists(ctfile)) return 1;
		//Open the file for reading:
		file_in.open(ctfile);
		//Make sure the file opened.  If not, return an error indicator.
		if(!file_in.is_open()) return 2;
		in.rdbuf(file_in.rdbuf());
	}

	// detect whether this is a dot-bracket or a CT file.
	if (in.peek()=='>') {
		// if the first character is '>' this is a dot-bracket file
		file_in.close(); // close the input stream (but only if it was a file -- not stdin) before calling opendbn. Ideally we could just pass our ifstream to an overload of opendbn
		return opendbn(ctfile);
	}

	//First read the first item in the file.  If this is -100, then the file was flagged as a CCT-formatted file.  (A more compact format.)
	in >> count;
	if (!in) { // test the istream to make sure a valid number was read (and not just text/garbage etc)
		SetErrorDetails("Invalid character data at line 1 column 0. Expected a valid number of bases.");
		return 29;
	}

	if (count == -100) { 
		//this is a CCT formatted file:
		in >> numofbases;
		in >> count;
		getline(in, sequenceLabel);
		SetSequenceLabel(sequenceLabel);
		allocate(numofbases);
		for (int i=1;i<=numofbases;i++) {
   			in >> numseq[i];
			nucs[i]=GetThermodynamicDataTable()->numtobase(numseq[i]);
			// Todo: check for invalid base (-1)
			hnumber[i] = i;
		}
		for (int numofstructures=1;numofstructures<=count;numofstructures++) {
			AddStructure();
			SetCtLabel(sequenceLabel,numofstructures);
			int energy;
			in >> energy;
			SetEnergy(numofstructures,energy);
    		for (int i=1;i<=numofbases;i++) {
				in >> basepair;
				if (basepair>i) SetPair(i,basepair,numofstructures);
			}
		}
		nucs[numofbases+1]='\0'; // add terminating null character.
		return 0; //success (CCT format)
	} else if (count < 0 || count > MAX_CT_BASES) {
	    // Count would only be negative (other than -100) if the file is corrupt or it is not a true CT file. 
		// Report this error via the SetErrorDetails mechanism
		SetErrorDetails(count < 0 ? "Negative number of bases." : "Total number of bases exceeds maximum.");
		return 29; // invalid file format
	} else { // this is a ct file:
		
		//FOR NOW, DISABLE READING CT FILES WITH STACKING INFO: THIS SEEMS TO BE A PROBLEM WITH THE WEBSERVER BECAUSE THERE ARE MANY FORMATS IN CURRENT USE
		
		//first decide if it contains nucleotide stacking data at the far right
		//determine this based on the length of the line:
		//in.getline(temp,1000);
		//in.getline(temp,1000);

		//i = (int) strlen(temp);
		//if (i==35) {
		//	ct->stacking = true;
		//}
		//by default, ct->stacking is false


		//Allocate memory in the structure, the first thing read was the sequence length
		int position,linker=0; // linker is the number of intermolecular linkers found (i.e. linker+1 is the number of RNA molecules/strands)
		allocate(count);
		numofbases = count;
		long linenumber=1;//keep track of linenumber to report errors

		for (int numofstructures = 1;!in.eof();++(numofstructures))	{
			AddStructure();
			getline(in, sequenceLabel);
			SetCtLabel(sequenceLabel,numofstructures);
			//Also set the sequencelabel if this is structure #1
			if (numofstructures==1) SetSequenceLabel(sequenceLabel);
			
			for (count=1;count<=((numofbases));count++)	{
				++linenumber;//reading from a new line, so increment the line number
				in >> position;//base number

				//Check that the bases are incrementing correctly
				if (position!=count) {
					//There was a problem
					SetErrorDetails(sfmt("Invalid nucleobase index %i (expected %i) in structure %i at line %li.", position, count, numofstructures, linenumber));
					return 29;
				}

				in >> base;//read the base

				nucs[count]=base[0];

				numseq[count]=GetThermodynamicDataTable()->basetonum(base[0]);
				if (numseq[count]==-1) {
					//This means the bases was not recognized, and an error should be reported
					SetErrorDetails(sfmt("Invalid nucleobase '%c' in structure %i at line %li.", base[0], numofstructures, linenumber));
					return 29;
				}
				
				if (data->isLinker(numseq[count])&&numofstructures==1) {
      				intermolecular = true;
       				inter[linker++] = count;
				}
				in >> position;//read the previous connected nucleotide
				//Check the numbering
				if (position!=(count-1)&&position!=0) {
					//There was a problem
					SetErrorDetails(sfmt("Unexpected backbone connection between nucleotides %i and %i in structure %i at line %li.", count, position, numofstructures, linenumber));
					return 29;
				}
				in >> position;//read the next connected nucleotide
				//Check the numbering
				if (position!=(count+1)&&position!=0) {
					//There was a problem
					SetErrorDetails(sfmt("Unexpected backbone connection between nucleotides %i and %i in structure %i at line %li.", count, position, numofstructures, linenumber));
					return 29;
				}

				// basepair is the index that the current base is paired with, according to the CURRENT line in the file.
				// (basepair can be either before or after this one)
				in >> basepair;

				// Verify that the base is not paired with itself.
				if (basepair==count) {
					SetErrorDetails(sfmt("Base %i is paired with itself in structure %i at line %li.", count, numofstructures, linenumber));
					return 29;
				}

				// Get the index that THIS pair is ALREADY assigned to pair with.
				// This will be 0 unless the CURRENT nucleotide is the 3' base of a pair that was 
				// first mentioned on a PREVIOUS line in the file.
				int assigned = GetPair(count,numofstructures);

				// For base-pairing consistency one of the following must be true:
				//  1).  basepair == 0    AND assigned == 0         (This base is not paired to any others)
				//  2).  basepair > count AND assigned == 0         (This is the 5' base -- paired to one that comes later in the sequence)
				//  3).  basepair < count AND assigned == basepair  (This is the 3' base -- paired to one that came earlier in the sequence)
				//       Note that #1 is a specific case of #3, so if #3 is checked, #1 can be ignored.
				// As a corollary, here are some examples of inconsistencies:
				//  1)   basepair == 0      AND  assigned != 0         (A previous base said it was paired to this one, but this one says it is unpaired.)
				//  2)   basepair > count   AND  assigned != 0         (This base says it is paired to a later base, but the later base is already paired (which it shouldn't be until after THIS line is processed.)
				//  3)   basepair < count   AND  assigned != basepair  (A previous base said it was paired to this one, but this one says it is paired to a DIFFERENT base (not the previous one).

				if ((basepair<count&&assigned!=basepair) || (basepair>count&&assigned!=0)) {
					// For error reporting, we might have to swap some numbers. count is always THIS base. But either assigned or basepair could be 0. If so, swap them.
					if (assigned==0) assigned = basepair;
					SetErrorDetails(sfmt("Inconsistent base pairing information between nucleotides %i and %i in structure %i at line %li.", count, assigned, numofstructures, linenumber));
					return 29;
				}

				// One scenario that isn't covered above is that two bases are both paired to the same one (which comes AFTER the other two)
				// e.g.:  Pairs: 1->5  3->5  5->3   Steps: 1. SetPair(1,5);   2. Verify: GetPair(3)==0 and basepair>count; SetPair(3,5)  3. Verify: GetPair(5)==3==basepair.  (no errors detected so far)
				if (basepair>count) {
					assigned = GetPair(basepair,numofstructures); //here assigned is the base to which the TARGET base is assigned. 
					//assigned should be 0. If not, it means a previous base already listed it as being paired to it.
					if (assigned!=0&&assigned!=basepair) {
						SetErrorDetails(sfmt("Inconsistent base pairing information. Bases %i and %i are both paired to base %i in structure %i at line %li.", assigned, count, basepair, numofstructures, linenumber));
						return 29;
					}
				}

				if (basepair>count) SetPair(count,basepair,numofstructures);//set base pairing info
				in >> hnumber[count];//read historical numbering
				//if (stacking) in >> basepr[numofstructures][count+numofbases];
			}

			++linenumber;//read another line
			in >> count; //start on next structure and see whether the end of file is reached
		}
	}

	nucs[numofbases+1]='\0'; // add terminating null character.
	return 0; //success
}

// Write a sequence file
// seqFileType indicates the format of sequence file to write: 0=Plain Text (No label or delimiters etc), 1=SEQ, 2=FASTA (default)
// returns 0 on error or 1 without error.
int structure::writeseq(const char *seqfile, int seqFileType, bool append) {
	const int TEXT=0, SEQ=1, FASTA=2;
	ofstream file;

	if (seqFileType<TEXT||seqFileType>FASTA)
		SetErrorDetails("Unknown sequence file format.");
	else if (GetSequenceLength()==0)
		SetErrorDetails("The sequence has not yet been read or is empty.");
	else {
		file.open(seqfile, append ? std::ios_base::app : std::ios_base::trunc);
		if (!file.good()) {
			SetErrorDetails("The output file could not be opened.");
			file.close();
		}
	}
	if (GetErrorDetails().size()!=0) return 0;

	if (seqFileType == SEQ)
		file << ';' << endl << GetSequenceLabel() << endl;
	else if (seqFileType == FASTA)
		file << '>' << GetSequenceLabel() << endl;

	int pos = 1, end = GetSequenceLength() + 1;
	while(pos < end) {
		int sz = min(end - pos, 80);
		file.write(nucs+pos*sizeof(*nucs), sz);
		pos+=sz;
	}
	if (seqFileType == SEQ)
		file << '1' << endl;

	file.close();
	return 1;
}


//Open seq was originally designed for reading .seq files.
//It has now been extended to automatically identify and read FASTA files,
//	where the identity line starts with a ">".
//returns 0 on error or 1 without error.
// RMW 2016-12-20: Changed file type detection: Now allows comments (starting with ';') at the top of FASTA files. Also allows parsing plain TEXT sequences.
//   - All comment lines (starting with ';') and empty or whitespace-only lines are ignored.
//   - If the first non-empty, non-comment line starts with '>' it is FASTA.
//   - Otherwise it is SEQ as long as at least one comment line was found.
//   - Othwerwise it is a plain TEXT sequence.  (i.e. no comments lines were found and the first non-empty line does NOT start with '>') 
int structure::openseq (const char *seqfile) {
	int errnum = openseqx(seqfile);
	return errnum==0?1:0; // return 0 on error OR 1 on success (opposite of openseqx)
}
//# same as openseq except that the return value is 0 on success or an error code corresponding to RNA::GetErrorMessage(int)
int structure::openseqx (const char *seqfile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	string line;
	string sequence;
	size_t start = string::npos;
	long linenumber = 0;
	enum SeqFileType { TEXT, SEQ, FASTA } fileType = TEXT;  // TEXT means plain text (AUGC...) with no comments or other delimiters.

	std::istream input(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(seqfile)) { 
		if (!fileExists(seqfile)) return 1; // file not found.
		//Open the file for reading:
		file_in.open(seqfile);
		//Make sure the file opened.  If not, return an error indicator.
		if(!file_in.is_open()) return 2;
		input.rdbuf(file_in.rdbuf());
	}

	//Now identify the file type.
    //  Usually, a starting > means FASTA while a starting ; means .seq 
	//  However, the FASTA/Pearson format allows for comment lines starting with ; See: https://en.wikipedia.org/wiki/FASTA_format

	//scan for first relevant line. Ignore lines that are empty, contain only whitespace, or start with a semi-colon (;).
	while (!getline(input,line).fail()) {
		++linenumber;
		if (line.length() == 0 || line[0] == '\r') continue; // line is empty
		if (line[0] == ';') {
			fileType = SEQ; // assume it is a SEQ, but this will be change to FASTA if the first character of the sequence label is '>'
			continue;
		}
		if (string::npos != (start = line.find_first_not_of(" \t\r"))) // if the line has any non-whitespace character, exit the loop
			break;
	}

	// Now past the comments. See if a valid line was found.
	if (start == string::npos) {
		//The end of the file was reached before finding the end of the comments, return with error
		SetErrorDetails("The file did not contain any nucleotides.");
		return 28; //Error reading sequence
	}

	if (line[start] == '>') { fileType = FASTA; start++; } // FASTA files must start with a '>' (following any ;-comments)

	if (start != 0) line.erase(0, start); // remove all characters up to, but not including start. (they were either whitespace or '>')
	
	bool reuseLine = false; // used for plain-text sequences

	//Set the current line as the sequence label (if FASTA or SEQ)
	if (fileType == SEQ || fileType == FASTA)
		SetSequenceLabel(line);
	else {
		// The filetype is plain text. 
		// *** Uncomment the following to disable parsing of plain-text sequences: ***
		//	SetErrorDetails("The format of the sequence file was invalid. It did not appear to be a SEQ or FASTA file.");
		//	return 29;
		SetSequenceLabel(getFileName(seqfile, true)); // Use the file name as a sequence label.
		reuseLine = true; // I'd prefer to use seekg, but there seems to be a bug in it on Windows.
	}

	bool foundEndChar=false; // indicates that a character was found that indiates the termination of a sequence. (for SEQ this is a '1')

	//now read each line and process it:
	while ((reuseLine||getline(input,line))&&!foundEndChar) {
		if (reuseLine) 
			reuseLine = false;
		else
			++linenumber;

		//Process each character in the current line
		for (int i=0;i<line.length();++i) {
			if (line[i] < 33) continue; // skip whitepace.

			if ((fileType==SEQ && line[i]=='1') || (fileType==FASTA && line[i]=='>')) {
				//Done reading sequence
				foundEndChar=true;
				break;
			} 

			// If the character is valid, add it. Otherwise return with an error.
			if (-1==data->basetonum(line[i])) {
				//This nucleotide wasn't recognized
				// report this error via the lastErrorDetails mechanism
				SetErrorDetails(sfmt("Invalid nucleobase '%c' at line %li column %i.", line[i], linenumber, i+1));
				return 28; //Error reading sequence
			} else {
				//Put the nucleotide in the sequence
				sequence+=line[i];
			}
		}//end of iteration over characters in the line
	} //end of while getLine loop

	if (fileType == SEQ && !foundEndChar) { 
		//The end of the file was reached before finding the 1 that indicates the sequence end
		// report this error via the SetErrorDetails mechanism
		SetErrorDetails("The file was missing the required '1' (one) that indicates the end of a sequence in a SEQ file.");
		return 29; //Invalid file format
	}
	
	if (sequence.size()==0) {
		SetErrorDetails("The file did not contain any nucleotides.");
		return 28; //Error reading sequence.
	}

	//The sequence string is populated with the sequence, now allocate the needed space in structure, and record the sequence
	SetSequence(sequence);
	return 0; //success (no error code)
}

// opens alignment file in FASTA format and populates vector<structure> 
int structure::open_alignment(const char* alignment_file)
{
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.

	ifstream infile;

	//Error handling that mirrors open_seq
	if (!isStdIoFile(alignment_file)) 
	{
		if (!fileExists(alignment_file)) return 1; // file not found.
		//Open the file for reading:
		infile.open(alignment_file);
		//Make sure the file opened.  If not, return an error indicator.
		if (!infile.is_open()) return 2;
	}

	string line;
	vector<string> names;
	vector<string> sequences;

	int i = -1;
	while (getline(infile, line))
	{
		if (line[0] == '>') {
			names.push_back(&line[1]);
			i++;
			sequences.push_back("");
		}
		else if (line[0] != '>') {
			sequences[i] = sequences[i] + line;
		}
	}
	infile.close();

	number_of_sequences = names.size();

	sequences_from_alignment.resize(names.size());

	for (int i = 0; i < number_of_sequences; i++)
	{	
		datatable* dt = GetThermodynamicDataTable();
		sequences_from_alignment[i].SetThermodynamicDataTable(dt);
		sequences_from_alignment[i].SetSequenceLabel(names[i]);
		sequences_from_alignment[i].SetSequence(sequences[i]);
		sequences_from_alignment[i].FillNoOfGapsMatrix();
		sequences_from_alignment[i].FillGappedtoUngappedMap();
		sequences_from_alignment[i].sequence_without_gaps = new structure;
		sequences_from_alignment[i].sequence_without_gaps->SetThermodynamicDataTable(dt);
		sequences_from_alignment[i].SetSequenceWithoutGaps(sequences[i]);
	}	

	return 0;
}

//Write a dot-bracket file
int structure::writedotbracket(const char * const filename, const int structurenumber, 
							   const DotBracketFormat format /*default DBN_FMT_MULTI_TITLE*/,
							   CTCommentProvider &commentProvider, const bool append) const {

	// if the filename (ctoutfile) is "-" then write to STDOUT.
	std::ostream out(std::cout.rdbuf());
	std::ofstream out_file; //only used if STDOUT is not used.
	if (!isStdIoFile(filename)) {
		//Open the file for writing:
		out_file.open(filename,append ? std::ios_base::app : std::ios_base::trunc);
		//Make sure the file opened.  If not, return an error indicator.
		if (!out_file.is_open()) {
			// not allowed because the method is marked const: SetErrorDetails(sfmt("Failed to open the file '%s'. Please verify the path and its permissions.", ctoutfile));
			return 2; // Error opening file
		}
		out.rdbuf(out_file.rdbuf());
	}

	vector<int> ranks(GetSequenceLength()+1);
	const char*const brackets = DBN_BRACKET_SYMBOLS; // e.g. "()<>{}[]AaBbCcDd"
	const int MAX_LEVEL = strlen(brackets)/2; // highest pseudoknot level allowed.
	
	int count=0; // number of structures written.
	int start=1, end=GetNumberofStructures();
	if (structurenumber<-1||structurenumber>end) {
		// not allowed because the method is marked const: SetErrorDetails(sfmt("Invalid structure number: %i. Valid structure numbers are from 1 to %i. Enter -1 to write all structures.", structurenumber, end));
		return 3; //structure number out of range.
	}
	if (structurenumber>0) start=end=structurenumber; // accept -1 or 0 to indicate that ALL structures should be written.
	string label;
	for (int n=start;n<=end;n++) {
		
		// write out the TITLE LINE (including structure label) if this is the first structure or DBN_BIT_MULTI_TITLE is set
		if (count==0||format&DBN_BIT_MULTI_TITLE) {
			if (n==1&&structurenumber==-1&&format&DBN_BIT_SEQ_LABEL) { // We only want to use the sequence label if we are writing all structures and this is the first.
				label = GetSequenceLabel();
				eraseEnergyLabel(label);
			} else { // use the sequence label instead of the structure label if this format bit is set.
				string comment = commentProvider.getComment(this, n);
				label = GetCtLabel(n);
				if (!comment.empty()) {
					trim(label);
					label = comment + "  " + label;
				}
			}
			trim(label); // remove whitespace from around label, including the final newline
			out << ">" << label << endl;
		}

		// write out the SEQUENCE LINE
		if (count==0||format&DBN_BIT_MULTI_SEQ) {
			for (int i=1;i<=numofbases;i++) 
				out << nucs[i];
			out << endl;
		}

		// write out the STRUCTURE LINE
		GetPseudoknotRanks(ranks, n);
		//DEBUG: cout << "ranks: " << join(ranks) << endl;
		for (int i=1;i<=numofbases;i++) {
			// Normally pairs are encoded with "()" but psedoknots (and higher-order knots) are encoded with alternate brackets.
			int level = min(ranks[i], MAX_LEVEL) - 1;  // get the knot level. 0 = no pseudoknot, 1 = first-order pseudoknots, 2 = second-order knots (i.e. which crossed other pseudoknots)
			if (GetPair(i,n)>i) out << brackets[2*level]; // the opening bracket, e.g. '('
			else if (GetPair(i,n)==0) out << ".";
			else out << brackets[2*level+1]; // the closing bracket, e.g. ')'
		}

		// write out the SIDE COMMENTS (aka SIDE TITLE/LABEL) if a format-bit is set
		if (format&DBN_BIT_SIDE_TITLES) {
			string label(GetCtLabel(n)); 
			trim(label);
			out << "\t" << label;
		}
		out << endl;
		count++;
	}
	if (out.fail()) return 2; // Error opening file
	return 0;
}


//Add a new structure:
void structure::AddStructure() 
{
	arrayofstructures.push_back(singlestructure(numofbases));

	//If this is the first structure, go ahead and copy the sequence label to the structure label:
	if(arrayofstructures.size() == 1) 
	{
		arrayofstructures[0].ctlabel = sequencelabel;
	}
}

//Remove all pairs from a structure:
void structure::CleanStructure(int structurenumber) {
	int i;

	for (i=1;i<=numofbases;++i) arrayofstructures[structurenumber-1].basepr[i]=0;



}





//Remove the last structure:
void structure::RemoveLastStructure() {

	arrayofstructures.pop_back();

}

//Remove the last structure:
void structure::RemoveAllStructures() {

	arrayofstructures.clear();

}

//This function sizes the vectors ctlable and energy to the right size for the maximum expected number of structures.
//This function is called multiple times if the number of structures is going the exceed the allocatedstructures.
void structure::allocatestructure(int structures) {
	
	
	//Set the maximum number of structures, to aid the vector in finding memory once:
	arrayofstructures.reserve(structures+1);

	

}


//Remove a specific structure:
void structure::RemoveStructure(int structurenumber) {

	arrayofstructures.erase(arrayofstructures.begin()+(structurenumber-1));

}


//Get the pointer to the underlying thermodynamic data
datatable *structure::GetThermodynamicDataTable() {
	return data;
}

//Set the pointer to the underlying thermodynamic data 
void structure::SetThermodynamicDataTable(datatable *DataTablePointer) {

	data=DataTablePointer;
}

// Returns true if the data property has been set to a valid datatable and its alphabet has been successfully read.
bool structure::IsAlphabetLoaded() {
	return data!=NULL&&data->loadedAlphabet;
}
// Returns true if the data property has been set to a valid datatable and the thermodyanamic tables have been read.
bool structure::IsThermoDataLoaded() {
	return data!=NULL&&data->loadedTables;
}

//sort the structures from lowest to highest free energy
void structure::sort() {
	
	
	//Use sort, and refer to the energy sort function: 
		
	std::sort(arrayofstructures.begin(),arrayofstructures.end());

	

}

//! This function is for debugging.  It checks each pair in each structure to look for inconsistencies.  
bool structure::ProblemwithStructures() {


	//Loop over all structures
	for (int structures=1; structures<=GetNumberofStructures();++structures) {

		//Loop over nucleotides
		for (int nucleotides=1;nucleotides<=GetSequenceLength();++nucleotides) {\
		

			//If nucleotides is paired, make sure its pairing partner points back
			if (GetPair(nucleotides,structures)>0) {

				if (GetPair(GetPair(nucleotides,structures),structures)!=nucleotides) {


					int a = GetPair(nucleotides,structures);
					int b = GetPair(a,structures);
					int c = GetPair(b,structures);

					return true;

				}

			}


		}


	}

	return false;




}



//The destructor:
structure::~structure()
{
	int i;

	if (allocated) {
		delete[] numseq;

		delete[] hnumber;
		delete[] nucs;
	}
	if (templated) {
		for (i=0;i<=numofbases;i++) {
    		delete[] tem[i];
   		}

   		delete[] tem;
	}
	DeleteSHAPE();

	if ( experimentalPairBonusExists ) {
		delete[] EX;
	}
	if (constant!=NULL) {
		//delete the equilibrium constants
		for (i=0;i<=numofbases;i++) {
    		delete[] constant[i];
   		}

   		delete[] constant;

	}
	
}



void structure::allocate(int size)

{
	
	numofbases = size;

	//Try not to include the following, unless found necessary in the future:
	//if (allocated) {
	//		//Delete memory if already allocated
	//	delete[] numseq;
	//	delete[] hnumber;
	//	delete[] nucs;

	//}
	numseq = new short int [2*size+1];
	hnumber = new int [size+1];
	nucs = new char [size+2];
   
	allocated = true;

}


void structure::deallocate()
{
	numofbases = 0;

	delete[] numseq;
	delete[] hnumber;
	delete[] nucs;

	allocated = false;
}


//this allocates space in an array that is used for applying an equilibrium constant for missing specific pairs
void structure::allocateconstant()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   constant = new double *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	constant[i] = new double [i+1];
   }


   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i;j<=numofbases;j++) {
    		constant[j][i] = 1.0;
		}
   }

}

// Returns the textual name for the type of restraints (e.g. SHAPE)
// Used in warning messages etc.
const char* restraintTypeName(RestraintType modifier) {
	switch(modifier) {
		case RESTRAINT_SHAPE:return "SHAPE";
		case RESTRAINT_SHAPE_DIFF:return "diffSHAPE";
		case RESTRAINT_SHAPE_AC: return "SHAPE_AC";
		case RESTRAINT_SHAPE_GU: return "SHAPE_GU";
		case RESTRAINT_DMS: return "DMS";
		case RESTRAINT_CMCT:return "CMCT"; 
		// For unknown types, just return the simple description "restraint"
		default: return "restraint"; // e.g. "Warning: Error in >>restraint<< file ..."
	}
}

double structure::Gammadist(const double data, const double shape, const double loc, const double scale){
	return (1/scale)*pow((data - loc)*(1/scale), (shape - 1))*exp(-(1/scale)*(data - loc))/tgamma(shape);
}


/*
double structure::Potential(const double data, const std::vector< std::vector<double> > &params, const double kT){
	// params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of components 1 and 2 respectively.
	double pairedprob = params[0][6]*Gammadist(data, params[0][0], params[0][1], params[0][2]) + 
	                   params[0][7]*Gammadist(data, params[0][3], params[0][4], params[0][5]); 
	double unpairedprob = params[1][6]*Gammadist(data, params[1][0], params[1][1], params[1][2]) + 
	                     params[1][7]*Gammadist(data, params[1][3], params[1][4], params[1][5]);
	return -kT*log(pairedprob/unpairedprob);
}
*/

double structure::Potential(const double data, const std::vector< std::vector<double> > &params, const double kT, const int ntcode) {
    // params[i] is for paired, params[i+1] for unpaired...
    // params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 
    // j=6,7 for weights of components 1 and 2 respectively.
    
    int idx;
    switch(ntcode) {
        case 1: idx=0;  //A (or for default case without nuc-specific potentials)
                break;
        case 2: idx=2;  //C
                break;
        case 3: idx=4;  //G
                break;
        case 4: idx=6;  //U
                break;
        case 0: return 0.0; // nuc is X,N
        case 5: return 0.0; // nuc is I 
    }

	double pairedprob = params[idx][6]*Gammadist(data, params[idx][0], params[idx][1], params[idx][2]) + 
	                    params[idx][7]*Gammadist(data, params[idx][3], params[idx][4], params[idx][5]); 
	double unpairedprob = params[idx+1][6]*Gammadist(data, params[idx+1][0], params[idx+1][1], params[idx+1][2]) + 
	                      params[idx+1][7]*Gammadist(data, params[idx+1][3], params[idx+1][4], params[idx+1][5]);
    
    return -kT*log(pairedprob/unpairedprob);        
}



void structure::ReadProbabilisticPotentialParams() {
	string filedir(getDataPath()); // getDataPath will return DATAPATH if it exists or "." otherwise.
	filedir += "/dists/";

	string line;
	int start, end, i, j;
	int nparams = 8;
	// Initialize parameters
	std::vector<double> shapecol1;
	std::vector<double> shapecol2;
	for(i=0; i<nparams; i++){
		shapecol1.push_back(0.0);
		shapecol2.push_back(0.0);
	}
	SHAPE_params.push_back(shapecol1);
	SHAPE_params.push_back(shapecol2);
    
    std::vector<double> dmscol1;
	std::vector<double> dmscol2;
	for(i=0; i<nparams; i++){
		dmscol1.push_back(0.0);
		dmscol2.push_back(0.0);
	}
	DMS_params.push_back(dmscol1);
	DMS_params.push_back(dmscol2);

    for(i=0; i<8; i++) {
        std::vector<double> dmscol;
	    for(j=0; j<nparams; j++){
		    dmscol.push_back(0.0);
	    }
	    DMS_paramsnt.push_back(dmscol);
    }
    
	std::vector<double> cmctcol1;
	std::vector<double> cmctcol2;
	for(i=0; i<nparams; i++){
		cmctcol1.push_back(0.0);
		cmctcol2.push_back(0.0);
	}
	CMCT_params.push_back(cmctcol1);
	CMCT_params.push_back(cmctcol2);

	// Start reading distribution parameters....note that this code is ugly, can be collapsed into a single parameter array
	// SHAPE
	string tmp = filedir + "SHAPEdist.txt";
	char *filename = (char*)tmp.c_str();
	ifstream SHAPEfile;
	SHAPEfile.open(filename);
	if (SHAPEfile.good()){
		getline(SHAPEfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(SHAPEfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				SHAPE_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		SHAPEfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}

	// DMS
	tmp = filedir + "DMSdist.txt";
	filename = (char*)tmp.c_str();
	ifstream DMSfile;
	DMSfile.open(filename);
	if (DMSfile.good()){
		getline(DMSfile, line); // pop off the header
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(DMSfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				DMS_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		DMSfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
    
    
    // DMS NT
	tmp = filedir + "DMSdist_nt.txt";
	filename = (char*)tmp.c_str();
	ifstream DMSNTfile;
	DMSNTfile.open(filename);
	if (DMSNTfile.good()){
		getline(DMSNTfile, line);
		for(i=0;i<8;i++) {
			start = 0;
			end = 0;
			getline(DMSNTfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				DMS_paramsnt[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		DMSNTfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
    
	// CMCT
	tmp = filedir + "CMCTdist.txt";
	filename = (char*)tmp.c_str();
	ifstream CMCTfile;
	CMCTfile.open(filename);
	if (CMCTfile.good()){
		getline(CMCTfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(CMCTfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				CMCT_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		CMCTfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
}

// This function calculates the pseudoenergy for a given reactivity data. It changes the calculation
// depending on the modifier specified, giving either the log-likelihood-ratio of the unpaired/paired
// probabilities given a reactivity distribution per modifier, or the "classic" Deigan et al. bonus
// term when no modifier or an unrecognized modifier is provided. 
double structure::CalculatePseudoEnergy(const double data, const RestraintType modifier, const double slope, const double intercept, const int ntcode, const bool use_params_if_available ) {
	

    std::vector< std::vector<double> > *params; // use a pointer to avoid making a copy
    int ntcode2;
    double data2 = data;
    double kT = 5.904976983149999;  // kT at 24 C in 10ths of kcal/mol. Value used by Cordero et al.
    

    // slope==0 && intercept==0 & !use_params_if_available is triggered by SHAPEss calls to CalculatePseudoEnergy
    if( data <= -500 || (slope == 0 && intercept == 0 && !use_params_if_available )) {
		return 0;
    }
    
    switch(modifier) {
		case RESTRAINT_SHAPE_AC:
		case RESTRAINT_SHAPE_GU:
			// This is only applied if SHAPE_AC or SHAPE_GU is specified
			// By default, use the Deigan "default" method when the modifier is "SHAPE"
			params = &SHAPE_params;
            ntcode2 = 1;
			break;
		case RESTRAINT_DMS:
			params = &DMS_params;
            ntcode2 = 1;
			break;
        case RESTRAINT_DMSNT:
            params = &DMS_paramsnt;
            ntcode2 = ntcode;
            kT = 3.0816567; // 0.5*kT at 37 C in 10ths of kcal/mol; 0.5 adjust is for double counting of internal pairs
            break;
		case RESTRAINT_CMCT:
			params = &CMCT_params; 
            ntcode2 = 1;
			break;
		default: 
			// Calculate energies for SHAPE and diffSHAPE.
			// For negative values of data, just return the intercept (baseline).
			// Otherwise, return the calculated pseudo-energy.
			return data>0 ? log(data+1.0)*slope+intercept : intercept;
	}
	// If we get here, we are using DMS, DMSNT, CMCT, SHAPE_AC or SHAPE_GU.
    
    if( params->empty() ) 
        return 0;
    
    // SHAPE default returns `intercept` when data<0
    // For Gamma potentials, data<=0 is undefined, and is poorly behaved for data->0
    // For DMSNT params, 0.01 is well-behaved valuea, so evalute here as `intercept`
    // Behavior for earlier (Cordero et al 2012) parameters has not been benchmarked, 
    // so preserve `return 0` for these restraints. However, I think it would make sense
    // to consider also evaluting these potentials at data~0.01 to return the equivalent 
    // of an intercept.
    if (modifier == RESTRAINT_DMSNT && data < 0.005) {
        data2 = 0.005;
    } else if (data < 0) {
        return 0;
    }

	const double val  = Potential(data2, *params, kT, ntcode2);
    


 	// RMW: why not return `intercept` instead of 0?  If this function returns `intercept` for SHAPE data <= 0.
	return std::isnan(val) ? 0 : val;

}

// Allocate the SHAPE and SHAPEss arrays if they haven't already been created.
// Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
void structure::AllocateSHAPE() {
	if (shaped) return; // already created.
	//initializes an array to hold SHAPE data for double stranded segments
	SHAPE = new double [2*numofbases+1];
	//initializes an array to hold SHAPE data for single stranded segments
	SHAPEss = new double [2*numofbases+1];
	shaped = true;
	// Initialize all points to zero.
	for (int position=0; position<=2*numofbases; position++) {
		SHAPE[position] = 0;
		SHAPEss[position] = 0;
	}
	// initialize triangular 2-d array that stores ss SHAPE energies for loops. 1st index is ending location, 2nd index is starting location
	SHAPEss_region = new short int*[numofbases+1];
	for (int i=1; i<=numofbases; i++) SHAPEss_region[i] = new short int[i];
}

// Deletes the SHAPE, SHAPEss, and SHAPEss_region arrays if the have been allocated.
// Also sets the `shaped` field to false and sets the SHAPE, SHAPEss, and SHAPEss_region
// pointers to NULL.
void structure::DeleteSHAPE() {
	if (!shaped) return; // already created.
	if (SHAPE!=NULL)    delete[] SHAPE;
	if (SHAPEss!=NULL)  delete[] SHAPEss;
	if (SHAPEss_region!=NULL) {
		for (int i=1; i<=numofbases; i++)
			delete[] SHAPEss_region[i];
		delete[] SHAPEss_region;
	}
	shaped=false;
	SHAPE=SHAPEss=NULL;
	SHAPEss_region=NULL;
}

// This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
// calculatePseudoEnergies (default true) indicate whether these data are being read for folding.  
//  (false means the raw values need to be stored.)
int structure::ReadSHAPE(const char *filename, RestraintType modifier, bool calculatePseudoEnergies) {
	int position;
	double data;

	//Only read the dists if they need to be read.  These are not used for SHAPE or diffSHAPE.
	if( !distsread && !(modifier==RESTRAINT_SHAPE||modifier==RESTRAINT_SHAPE_DIFF) ){
		ReadProbabilisticPotentialParams();
		distsread=true;
	}

	// Allocate the SHAPE and SHAPEss arrays if they haven't already been created.
	// Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
	AllocateSHAPE();

	// temporary array to hold reactivity data for double-stranded segments
	double *SHAPEnew = new double [2*numofbases+1];
	// temporary array to hold reactivity data for single stranded segments
	double *SHAPEssnew = new double [2*numofbases+1];
	// useful for bootstrapping -- keep count of how many times the data for a given position is specified.
	int *num_data_points =  new int [numofbases+1];

	// the array pointers above will be automatically deleted if return is called before the end of the function.
	auto_delete<double,true> d1(SHAPEnew), d2(SHAPEssnew); auto_delete<int,true> d3(num_data_points);

	for (position=0; position <= 2*numofbases; position++) {
		SHAPEnew[position] = 0.0;
		SHAPEssnew[position] = 0.0;
	}
	for (position=0; position <= numofbases; position++) {
		num_data_points[position] = 0;
	}

	if (!fileExists(filename)) return 201; //file not found
	ifstream in(filename);
	if (!in.good()) return 202; // error opening file


	vector<int> invalidPositions;
	bool hasRepeats = false; //track whether any values have been repeated.
	while (!(in >> position >> data).fail()) { // this will ensure the last data point is read, even if there is no newline at the end of the file. (EOF can be set after a valid datapoint is read if there is no whitespace following it.)
		//read and parse all data
			//required format is rows with sequence position followed by reactivity
		if (position<1||position>numofbases) {
			invalidPositions.push_back(position);
			continue;
		}
			if (calculatePseudoEnergies) { 
				//store calculated pseudo-energy values.
				SHAPEnew[position]   += CalculatePseudoEnergy(data, modifier, SHAPEslope, SHAPEintercept, numseq[position], true /*use_params_if_available*/ );
				SHAPEssnew[position] += CalculatePseudoEnergy(data, modifier, SHAPEslope_ss, SHAPEintercept_ss, numseq[position],  false /*use_params_if_available*/ );
                
                //cout << position << " " << nucs[position] << " " << data << " " << SHAPEnew[position]/10 << "\n";

			} else { 
				//store the raw SHAPE values
				SHAPE[position] = data;
				SHAPEss[position] = data;
			}
			if (num_data_points[ position ]) hasRepeats = true;
			num_data_points[ position ]++;
	}
	in.close();
    

	if (!invalidPositions.empty())
		cwarn() << "Warning: Invalid nucleobase positions in " << restraintTypeName(modifier) << " file " << filename << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;

	if (calculatePseudoEnergies) {
		for (position=1;position<=numofbases;position++) {
			if(num_data_points[position] > 0){
				// SumShapeRepeats is true by default, but can be turned off by setting the environment variable AVG_SHAPE_REPEATS to 1.
				if (structure::SumShapeRepeats) { 
					// SHAPE/DMS/CMCT files can be ‘resampled with replacement’ (bootstrapped) to
					//  mock simulate having extra or no data at some positions — we need
					//  to *accumulate* the pseudo energies there. 
					SHAPE[position]+=SHAPEnew[position];
					SHAPEss[position]+=SHAPEssnew[position];
				} else {
					SHAPE[position]+=SHAPEnew[position]/num_data_points[position];
					SHAPEss[position]+=SHAPEssnew[position]/num_data_points[position];
				}
			}
		}

		// SHAPE values in ‘second set’ must exactly match SHAPE values in first set.
		for (position=1;position<=numofbases;position++) {
			SHAPE[position+numofbases]=SHAPE[position];
			SHAPEss[position+numofbases]=SHAPEss[position];
		}
	}

	if (hasRepeats && ShowWarnings!=0 && SumShapeRepeats) {
		// Show a warning about duplicate SHAPE values
		// unless warnings are suppressed or the user has expressly 
		// permitted multiple values by setting the AVG_SHAPE_REPEATS environment variable.
		ostream& warn = cwarn(); // (cwarn() gives the user some control over how warnings are shown.)
		warn << "Warning: The following nucleobase positions were repeated in " << restraintTypeName(modifier) << " file " << filename << ":";
		for (position=1;position<=numofbases;position++)
			if (num_data_points[position])
				warn << " " << position;
		warn << endl << "(This may be OK if you are bootstrapping -- resampling with replacement.)" << endl;
	}

	//fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	FillSHAPEssRegions();

	// (deletion of temporary array is handled by the auto_delete objects defined above)
	return 0;
}

// Return an exact copy of SHAPE data.
// The return value is a pointer to a dynamically allocated array of doubles.
// If includeSHAPEss is true, the array will contain both SHAPE and SHAPEss values and will
// have a size of 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and SHAPEss from arr[2N+1] to arr[4N+1])
// If includeSHAPEss is false, the array will only contain SHAPE data and will have a size of 
// 2*numofbases+1.
double* structure::CopySHAPE(const bool includeSHAPEss) {
	if (!shaped) return NULL;
	double* arr = new double[(includeSHAPEss?2:1)*(2*numofbases+1)]; // 2N+1 or 4N+2
	// Fill arr[0] to arr[2N] with SHAPE   (note that arr[0] is unused since nucleotide indexing is 1-based)
	for (int i=0; i<=2*numofbases; i++) 
		arr[i]=SHAPE[i];
	if (includeSHAPEss)
		// Fill arr[2N+1] to arr[4N+1] with SHAPEss  (note that arr[2N+1] is unused since nucleotide indexing is 1-based)
		for (int i=0; i<=2*numofbases; i++) 
			arr[2*numofbases+1+i]=SHAPEss[i];
	return arr;
}

// Load data from an array of doubles into the SHAPE array.
// This function copies the data (it does NOT take ownership of the passed-in array).
// If includeSHAPEss is true, the array must contain both SHAPE and SHAPEss values and must
// have a size of at least 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and SHAPEss from arr[2N+1] to arr[4N+1])
// If includeSHAPEss is false, the array need only contain SHAPE data and must have a size of 
// at least 2*numofbases+1.
// If the passed-in array is NULL, this structure's SHAPE data will be deleted.
void structure::LoadSHAPE(const double* shapeArray, const bool includeSHAPEss) {
	if (shapeArray==NULL)
		DeleteSHAPE();
	else {
		AllocateSHAPE(); // create arrays if not already created.
		// Fill arr[0] to arr[2N] with SHAPE   (note that arr[0] is unused since nucleotide indexing is 1-based)
		for (int i=0; i<=2*numofbases; i++) 
			SHAPE[i]=shapeArray[i];
		if (includeSHAPEss)
			// Fill arr[2N+1] to arr[4N+1] with SHAPEss  (note that arr[2N+1] is unused since nucleotide indexing is 1-based)
			for (int i=0; i<=2*numofbases; i++) 
				SHAPEss[i]=shapeArray[2*numofbases+1+i];
	}
}

//Read offset files
//The SSOffset provides free energies to add to nucleotides if they are single-stranded
//The DSOffset provides free energies to add to nucleotides if they are double-stranded
//This function must be called ofter reading SHAPE files, if read, because the same infrastructure is used.
//Either filename can be NULL, in which case that offset is not recorded.
int structure::ReadOffset(const char *SSOffset, const char *DSOffset) {
	int i,position;
	double data;

 	// Allocate the SHAPE and SHAPEss arrays if they haven't already been created.
	// Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
	AllocateSHAPE();

	vector<int> invalidPositions;
	if (SSOffset!=NULL) {
		ssoffset = true;//we are reading an offset
		if (!fileExists(SSOffset)) return 201; //file not found
		ifstream in(SSOffset);
		if (!in.good()) return 202; // error opening file

		while (!(in >> position >> data).fail()) {
			if (position<1||position>numofbases) {
				invalidPositions.push_back(position);
				continue;
			}
			//read and parse all data
				//required format is rows with sequence position followed by reactivity
				//multiply the specified free energy change by converstionfactor, the factor by which energies in kcal/mol are multipled
			SHAPEss[position]+=(data*conversionfactor);
			SHAPEss[position+numofbases]+=(data*conversionfactor);
		}
		in.close();
		if (!invalidPositions.empty())
			cwarn() << "Warning: Invalid nucleobase positions in SS Offset file " << SSOffset << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;
	}

	invalidPositions.clear();
	if (DSOffset!=NULL) {
		if (!fileExists(DSOffset)) return 201; //file not found
		ifstream in(DSOffset);
		if (!in.good()) return 202; // error opening file

		while (!(in >> position >> data).fail()) {
			//read and parse all data
				//required format is rows with sequence position followed by reactivity
				//multiply the specified free energy change by converstionfactor, the factor by which energies in kcal/mol are multipled
			if (position<1||position>numofbases) {
				invalidPositions.push_back(position);
				continue;
			}
			SHAPE[position]+=(data*conversionfactor);
			SHAPE[position+numofbases]+=(data*conversionfactor);
		}
		in.close();
		if (!invalidPositions.empty())
			cwarn() << "Warning: Invalid nucleobase positions in DS Offset file " << DSOffset << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;
	}

	//fills SHAPEss_region, 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	FillSHAPEssRegions();
	return 0;
}

// Fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using
//   single-stranded (ss) SHAPE parameters or offsets.
// The SHAPEss and SHAPEss_region arrays must already be allocated (by calling AllocateSHAPE)
//   and SHAPEss must already contain ss SHAPE pseudo-energies.
// This is called from e.g. ReadSHAPE and ReadOffset.
void structure::FillSHAPEssRegions() {
	//fills 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	//If SHAPE was previously read, this is a redo of the action
	for (int j = 2; j <= numofbases; j++) {
		SHAPEss_region[j][j - 1] = (short int)(SHAPEss[j] + SHAPEss[j-1]); //sets energy for "zero sized loop".  Acts as starting value to add onto below
		for (int i = j - 2; i >= 1; i--) {
			SHAPEss_region[j][i] = SHAPEss_region[j][i + 1] + (short int)(SHAPEss[i]); //adds energy for additional loop element
		}
	}
}

//this function returns a psuedo energy term for a single stranded nucleotide based on SHAPE data
short int structure::SHAPEss_give_value(int index) 
{
	if (shaped) 
	{
		if (index > numofbases) return (short int)SHAPEss[index - numofbases];
		else return (short int)SHAPEss[index];
	}
	else return 0;
}

//this is the function that will return a pseudo energy term for hairpin loops based off of SHAPE data
int structure::SHAPEss_calc(int index_i, int index_j) 
{
	if (shaped) 
	{
		//accounts for the SHAPEss_region array being only NxN, not 2Nx2N
		if (index_i > numofbases) index_i -= numofbases;
		if (index_j > numofbases) index_j -= numofbases;
		if (index_i > index_j) {
			int temp_index = index_i;
			index_i = index_j;
			index_j = temp_index;
		}
		return SHAPEss_region[index_j][index_i];
	} else return 0;  //if no shaped data is being used, return zero
}


// This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
// This function is largely deprecated by the pseudo-free energy approach.  It is still available for experimentation.
// Returns 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
int structure::ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold) {
	int position;
	float data;

	if (!fileExists(filename)) return 201;  // file not found
	ifstream in(filename);
	if (!in.good()) return 202; // error opening file.

	vector<int> invalidPositions;
	while (!(in >> position >> data).fail()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity
		if (position<1||position>numofbases) {
			invalidPositions.push_back(position);
			continue;
		}
		if (data>=SingleStrandThreshold) {
			AddSingle(position);
		}
		else if (data>=ModificationThreshold) {
			AddModified(position);
		}
	}
	in.close();
	if (!invalidPositions.empty())
		cwarn() << "Warning: Invalid nucleobase positions in SHAPE file " << filename << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;
	return 0;
}

// This function reads an experimental pair bonus file, similar to SHAPE, but just straightforward
// application as kcal bonuses.  As with SHAPE, bonus is applied at 0x, 1x, and 2x for
//  single stranded, edge base pairs, and internal base pairs.
// Returns 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
int structure::ReadExperimentalPairBonus(const char *filename, double const experimentalOffset, double const experimentalScaling ) {
	//int position;
	//float data;

	//initializing 2-d array that stores bonuses for base pairs.
	// actually this is already done in RNA.cpp. Thanks Pablo!
	//delete[] EX;
	EX = new double *[ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++) EX[i] = new double [ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++)
	  for (int j = 0; j < 2*numofbases+1; j++)
	    EX[i][j] = 0.0;

	for (int i = 1; i <= 2*numofbases; i++)
	  for (int j = 1; j <= 2*numofbases; j++)
	    EX[i][j] = experimentalOffset * conversionfactor;

	int i( 1 ), j( 1 ), count( 0 );
	double val;
    char header;
    
    int numofbases2 = numofbases*numofbases;


    // check that file exists 
    if (is_blank(filename) || !fileExists(filename)) return 201; // file doesn't exist
    
    ifstream in(filename);
    if (!in.good()) return 202; //error opening file
    
    // determine file format
    in.get(header);
    if (header == ';') {
        // format is columnar; expecting three columns: i(int) j(int) val(float)

        in.ignore(1000,'\n'); // pop off the complete header line
        
        while (count < numofbases2) {
            
            in >> i;
            in >> j;
            in >> val;
            
            if (!in.good()) break; // either end of file or bad value couldn't be coerced into i/j/val
            
            EX[i           ][j           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[i+numofbases][j           ] = EX[i][j];
	        EX[i           ][j+numofbases] = EX[i][j];
	        EX[i+numofbases][j+numofbases] = EX[i][j];
            
            EX[j           ][i           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[j+numofbases][i           ] = EX[j][i];
	        EX[j           ][i+numofbases] = EX[j][i];
	        EX[j+numofbases][i+numofbases] = EX[j][i];
       
            count++;
        }
        

        if (in.eof()) {
            cout << count << " columnar pairing restraints read...";
        } else {
            SetErrorDetails(sfmt("Experimental bonus file '%s' intrepreted as columnar format contains improper value or is incorrectly formatted", filename));
            return 203; // see RNA::GetErrormessage
        }


    // no header, so matrix format
    } else {
        in.unget(); //rewind file
    
	    while (!in.eof() && j <= numofbases) {
	        //read and parse all data
	        //required format is bonuses in square matrix
	        in >> val;

	        EX[i           ][j           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[i+numofbases][j           ] = EX[i][j];
	        EX[i           ][j+numofbases] = EX[i][j];
	        EX[i+numofbases][j+numofbases] = EX[i][j];
	        count++;

	        i++;
	        if ( i > numofbases ){
	            i = 1;
	            j++;
	        }
	    }
	
        if ( count != numofbases2 ){
		    SetErrorDetails(sfmt("Experimental bonus file '%s' intrepreted as matrix format but did not have expected number of values\nFound %i but expected %i.\nIf columnar format, first line needs to start with ';'", filename, count, numofbases*numofbases));
	        return 203; // see RNA::GetErrormessage

	    }
	}
 

    in.close();
	experimentalPairBonusExists = true;
	return 0;
}



// Write SHAPE and SHAPEss parameters out to file, exactly as they are currently stored.
// Currently this is just for debugging purposes.
// Returns 0 on success or 2 if the file cannot be opened for writing.
int structure::WriteSHAPE(const string &outfile, bool printHeaders) {
	ofstream out(outfile.c_str());
	if (!out.good()) return 2; // failed to open file.
	if (printHeaders) out << "# " << GetSequenceLabel() << endl;
	if (printHeaders) out << "# SHAPE" << endl;
	for(int i=1;i<=2*GetSequenceLength();i++)
		out << i << "\t" << SHAPE[i] << endl;
	if (printHeaders) out << endl << "# SHAPEss" << endl;
	for(int i=1;i<=2*GetSequenceLength();i++)
		out << i << "\t" << SHAPEss[i] << endl;
	out.close();
	return 0;
}

const string& structure::GetErrorDetails() {
	return lastErrorDetails;
}
void structure::SetErrorDetails(const string& details) {
	lastErrorDetails = details;
}

// returns 0, 1 or 2 based on the value of a c-string (ON=1, OFF=0, ERR=2)
int parse_OnOffErrFlag(const char*const cstr) {
	string warn=as_str(cstr); toUpper(warn);
	return (warn=="OFF"||warn=="0")?0:(warn=="ERR"||warn=="2")?2:1;
}
// Set default behavior regarding warnings. 1=ON (stdout) 1=ERR (stderr) 0=OFF
int structure::ShowWarnings(parse_OnOffErrFlag(getenv("RNA_WARNINGS")));

// affects how multiple datapoints for the same nucleboase are handled in SHAPE data.
// the default is to add them for the "resample with replacement" technique 
// but setting the environment variable AVG_SHAPE_REPEATS to 1 causes
// multiple values to be averaged instead.
bool structure::SumShapeRepeats(is_blank(getenv("AVG_SHAPE_REPEATS")));


// Return an appropriate output stream for warnings --
// one of cout, cerr or nullstream, depending on the value of ShowWarnings.
// TODO: Later this can be changed to allow a third option -- a stringstream to capture warnings for library clients (e.g. Java)
ostream& structure::cwarn() {
	return ShowWarnings==0?NullStream::Default:ShowWarnings==2?cerr:cout;
}



// This is the default implementation of a function to return a comment/label like 
// "ENERGY = ..." which is inserted before the existing structure label when writing 
// CT or dot-bracket files.
// If GetEnergy() is zero, this function returns "" to indicate that no comment
// should be inserted.
string CTComments::EnergyCommentProvider::getComment(const structure* ct, const int structurenumber) {
	const int energy = ct->GetEnergy(structurenumber);
	if (energy==0) return "";
	stringstream comment("ENERGY = "); comment.seekp(0, std::ios::end); // move to the end of the string so we append the rest.
	comment << std::fixed << std::setprecision(conversionprecision) << (energy/(float)conversionfactor);
	return comment.str();
}
// Create memory storage location for static members of CTComments.
CTComments::EnergyCommentProvider CTComments::Energy;
CTComments::NoCommentProvider     CTComments::None;

//comparison function used by qsort in structure::sort
int ecompare(const void *i, const void *j) {

	return (**((short int**)i)-**((short int**)j));

}


//read a ct file with sequence and structural information
#define linelength 20






//takes a nucleotide input and stores a numeric representation
//This is depracated because the alphabet needs to come from structure, which reads it from disk
/*void tonum(char *base,structure *ct,int count)	{
if (!strcmp(base,"A")) (ct->numseq[count] = 1);
else if(!strcmp(base,"B")) {
	(ct->numseq[count] = 1);

}
else if(!strcmp(base,"a")) {
	ct->numseq[count]=1;
	ct->AddSingle(count);

}
else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
else if(!strcmp(base,"Z")) {
	(ct->numseq[count] = 2);

}
else if(!strcmp(base,"c")) {
	ct->numseq[count] = 2;
	ct->AddSingle(count);
}
else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
else if(!strcmp(base,"H")) {
	(ct->numseq[count] = 3);

}
else if(!strcmp(base,"g")) {
	ct->numseq[count] = 3;
	ct->AddSingle(count);
}

else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
else if(!strcmp(base,"V")||!strcmp(base,"W")) {
	(ct->numseq[count] = 4);

}
else if(!strcmp(base,"u")||!strcmp(base,"t")) {
	ct->numseq[count] = 4;
	ct->AddSingle(count);
}

else if(!strcmp(base,"I")) {
	ct->numseq[count] = 5;
	ct->intermolecular= true;
}

else (ct->numseq[count]=0);  //this is for others, like X
return;
}*/






//Depracated because of new extended alphabet
/*char *tobase (int i)

{  //function is a look up table to convert the base
	// 	integer represention to a familiar character

	if (i==1) return "A";
	 else if (i==2) return "C";
	 else if (i==3) return "G";
	 else if (i==4) return "U";
	 else if (i==0) return "X";
    else if (i==5) return "I";
	 else return "?";   //to track down mistakes

}*/








//#endif


//add the factor from SHAPE calculation
//This pseudo-energy was calculated when the file was loaded (see structure.cpp).
//The pseudo energy is applied twice for each nuc in interior pair and once for each nuc in terminal pair.
inline integersize SHAPEend(int i, structure *ct) {
	if (ct->shaped) return (integersize) ct->SHAPE[i];
	return 0;
}



//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

integersize erg1(int i, int j, int ip, int jp, structure *ct, datatable *data)	{

//	if (ct->GetNumberofSequences() > 1) return multi_erg1(i, j, ip, jp, ct, data);
	
	integersize energy;
	if ((i==(ct->GetSequenceLength()))||(j==((ct->GetSequenceLength())+1))) {
      	//this is not allowed because n and n+1 are not covalently attached
		energy = INFINITE_ENERGY;
	}		
	else 
	{
		energy = data->stack[(ct->numseq[i])][(ct->numseq[j])][(ct->numseq[ip])][(ct->numseq[jp])] 
				+ data->eparam[1];
			
		energy += SHAPEend(i, ct);
		energy += SHAPEend(j, ct);
		energy += SHAPEend(ip, ct);
		energy += SHAPEend(jp, ct);

		if (ct->experimentalPairBonusExists) 
		{
			energy = energy 
					+ 0.5 * ( ct->EX[i][j] + ct->EX[j][i] ) 
					+ 0.5 * (ct->EX[ip][jp] + ct->EX[jp][ip]);						
		}
	}
	return energy;
}

//calculate the sum of energies of stacked base pairs across a multiple sequence alignment by calling erg1 for each sequence.
integersize multi_erg1(int i, int j, int ip, int jp, structure* ct, datatable* data) 
{
	integersize energysum = 0;
	int N = ct->GetNumberofSequences();
	for (int seqnumber = 0; seqnumber < N; seqnumber++)	{
		energysum += erg1(i, j, ip, jp, ct->get_individual_sequence(seqnumber), data);
	}
	return energysum;
}


// Hamming distance between (ai, aj) and (bi, bj)
int d(int ai, int aj, int bi, int bj)
{
	return (2 - ((ai == bi) ? 1 : 0) - ((aj == bj) ? 1 : 0));
}

// check if base pairing
// function already there; use that! can_pair it will be able to handle other alphabet
int pij(int ai, int aj, datatable * data)
{
	return (data->pairing[ai][aj]) ? 1 : 0;
}

// delta returns 1 if both caharacters are same else returns 0
int delta(int a, int b)
{
	return (a == b) ? 1 : 0;
}

// factorial
int factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// qij function (Equtation 5 RNAalifold 2002)
float qij(const vector<int>& i, const vector<int>& j, datatable * data)
{
	float sum = 0.0;
	int N = i.size();

	for (int alpha = 0; alpha < N; alpha++)
	{
		sum = sum + pij(i[alpha], j[alpha], data) + delta(i[alpha], data->basetonum('-')) * delta(j[alpha], data->basetonum('-'));
	}
    
	return (1.0 - ( sum / N) );
}

// Cij function (Equtation 4 RNAalifold 2002)
float Cij(const vector<int>& i, const vector<int>& j, datatable * data) 
{
	int N   = i.size();
	int NC2 = (factorial(N) / (factorial(N - 2) * factorial(2)));
	float sum = 0.0;

	for (int alpha = 0; alpha < N; alpha++)
	{
		for (int beta = 0; beta < N; beta++)
		{
			if (alpha < beta)
			{
				sum = sum + d(i[alpha], j[alpha], i[beta], j[beta]) 
					* pij(i[alpha], j[alpha], data) * pij(i[beta], j[beta], data);
			}
		}
	}

	return (sum / NC2);
}

// Calculate covariation energy term for column i and j in multiple sequence alignment
// Bij = Cij - phi_1 * qij where phi_1 = 1
integersize multi_erg_covar(int i, int j, structure* ct) 
{
	int N = ct->number_of_sequences;
	vector<int> column_i, column_j;
	
	for (int seqnumber = 0; seqnumber < ct->GetNumberofSequences(); seqnumber++) 
	{
		column_i.push_back(ct->get_individual_sequence(seqnumber)->numseq[i]);
		column_j.push_back(ct->get_individual_sequence(seqnumber)->numseq[j]);
	}
	
	datatable* data = ct->GetThermodynamicDataTable();
	float B = Cij(column_i, column_j, data) - qij(column_i, column_j, data);
	
	// making vector not cheap; integrated single function; we come back to this later.
	return (int) B*conversionfactor;
}



//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp

// integersize erg2_fast(int i, int j, int ip, int jp, structure* ct, datatable* data, char a, char b) {	
	//cout << ct->GetSequenceLength() << endl;
	/* size,size1,size2 = size of a loop
	energy = energy calculated	
	loginc = the value of a log used in large hairpin loops
	*/
/*
	if (((i <= (ct->GetSequenceLength())) && (ip > (ct->GetSequenceLength()))) || ((
		jp <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength())))) {
		//A loop cannot contain the ends of the sequence
		return INFINITE_ENERGY;
	}

	integersize energy;
	int size, size1, size2, loginc, lopsid, energy2, count, k, energy_option;

	// a typical internal or bulge loop:
	size1 = ip - i - 1;
	size2 = j - jp - 1;

	if (size1 == 0 || size2 == 0) {//bulge loop

		size = size1 + size2;
		if (size == 1) {
			count = 1;
			energy = data->stack[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[ip]][ct->numseq[jp]]
				+ data->bulge[size] + data->eparam[2];


			if (size1 == 1) {

				//count the number of alternative bulges that exist:
				
				k = i;
				while (ct->numseq[k] == ct->numseq[i + 1]) {
					count++;
					k--;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() || k == 0) break;
				}

				k = ip;
				while (ct->numseq[k] == ct->numseq[i + 1]) {
					count++;
					k++;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() + 1 || k > 2 * ct->GetSequenceLength()) break;
				}
				
				//give bonus to C adjacent to single C bulge
				if ((ct->IsNuc(i + 1, 'C') || ct->IsNuc(i + 1, 'c')) && count > 1) energy += data->singlecbulge;

			}

			else {
				//size2 == 1

				//count the number of alternative bulges that exist:
				
				k = jp;
				while (ct->numseq[k] == ct->numseq[jp + 1]) {
					count++;
					k--;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() || k == 0) break;
				}
				k = j;
				while (ct->numseq[k] == ct->numseq[jp + 1]) {
					count++;
					k++;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() + 1 || k > 2 * ct->GetSequenceLength()) break;
				}
				
				//give bonus to C adjacent to single C bulge
				if ((ct->IsNuc(j - 1, 'C') || ct->IsNuc(j - 1, 'c')) && count > 1) energy += data->singlecbulge;

			}
			//apply a correction for the number of equivalent states because
				//the bulge can move to adjacent sites
			 energy -= (int)(data->RT * conversionfactor * log((double)count));
		}
		else if (size > 30) {

			loginc = int((data->prelog) * log(double((size) / 30.0)));
			energy = data->bulge[30] + loginc + data->eparam[2];
			energy = energy + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);

		}
		else {
			energy = data->bulge[size] + data->eparam[2];
			energy = energy + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
		}
	}
	else {//internal loop
		size = size1 + size2;
		lopsid = abs(size1 - size2);

		if (size > 30) {

			loginc = int((data->prelog) * log((double((size)) / 30.0)));
			if (size1 == 1 || size2 == 1) {
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
					data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING			
				if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
				else {
					energy += lopsid * data->poppen[min(2, min(size1, size2))];
					data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
				}
#else
				energy += min(data->maxpen, lopsid * data->poppen[min(2, min(size1, size2))]);
#endif					
			}

			else {
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstki[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
					data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING						
				if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
				else {
					energy += lopsid * data->poppen[min(2, min(size1, size2))];
					data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
				}
#else
				energy += min(data->maxpen, lopsid * data->poppen[min(2, min(size1, size2))]);
#endif //COUNTING						

			}
		}
		else if ((size1 == 2) && (size2 == 2)) {//2x2 internal loop
			energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
				[ct->numseq[j]][ct->numseq[jp]]
				[ct->numseq[i + 1]][ct->numseq[i + 2]]
				[ct->numseq[j - 1]][ct->numseq[j - 2]];
		}
		else if ((size1 == 1) && (size2 == 2)) {//2x1 internal loop
			energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]]
				[ct->numseq[j - 1]][ct->numseq[jp + 1]][ct->numseq[ip]][ct->numseq[jp]];


		}
		else if ((size1 == 2) && (size2 == 1)) {//1x2 internal loop
			energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]]
				[ct->numseq[ip - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]];
		}
		else if (size == 2) {//a single mismatch
			energy = data->iloop11[ct->numseq[i]][ct->numseq[i + 1]][ct->numseq[ip]]
				[ct->numseq[j]][ct->numseq[j - 1]][ct->numseq[jp]];
		}
		else if (size1 == 1 || size2 == 1) { //this loop is lopsided
		//this is a 1xn loop:
			energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];
#ifdef COUNTING					
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
			}
#else	
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
#endif	
		}
		else if ((size1 == 2 && size2 == 3) || (size1 == 3 && size2 == 2)) {
			//this is a 2x3 loop
			energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];

#ifdef COUNTING					
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
			}
#else	
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
#endif					


		}
		else {
			energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];

			#ifdef COUNTING						
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;

			}
			#else
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
			#endif //COUNTING				
		}
	}

	return energy;
}
*/
/*
//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
integersize erg2(int i,int j,int ip,int jp,structure *ct, datatable *data,
	char a, char b)
{	
//	if (ct->GetNumberofSequences() > 1) return multi_erg2(i, j, ip, jp, ct, data, a, b); 

	if (((i <= (ct->GetSequenceLength())) && (ip > (ct->GetSequenceLength()))) || ((
		jp <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength())))) {
		//A loop cannot contain the ends of the sequence
		return INFINITE_ENERGY;
	}

	integersize energy;
	int size, size1, size2, loginc, lopsid, energy2, count, k, energy_option;
	
	size1 = ip - i - 1;
	size2 = j - jp - 1;

	if ((a > 0) || (b > 0)) {
		// the loop contains a nuc that
		// should be double stranded
		if ((a & DUBLE) || (b & DUBLE)) { 
			return INFINITE_ENERGY; 
		}   
		else if ((a & INTER)) {
			//the loop is actually between two strands (ie: intermolecular)
			if (size2 > 1) {//free energy is that of two terminal mismatches
			//and the intermolecular initiation
				energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]];
			}
			else if (size2 == 1) {//find the best terminal mismatch and terminal
			 //stack free energies combination

				energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]) +
					erg4_nc(jp, ip, ip - 1, 2, ct, data, false) + penalty_nc(jp, ip, ct, data);

				energy_option = 1;

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]) +
					erg4_nc(i, j, i + 1, 1, ct, data, false) + penalty_nc(i, j, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 2;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
				   //now consider if coaxial stacking is better:
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[jp + 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 3;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j - 1]][ct->numseq[ip - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[j - 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 4;
					energy = energy2;
				}

				switch (energy_option) {
				case 1:
					energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]
						+ erg4(jp, ip, ip - 1, 2, ct, data, false) + penalty(jp, ip, ct, data);
					break;

				case 2:
					energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]
						+ erg4(i, j, i + 1, 1, ct, data, false) + penalty(i, j, ct, data);
					break;

				case 3:
					energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]
						+ data->coaxstack[ct->numseq[jp + 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;

				case 4:
					energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j - 1]][ct->numseq[ip - 1]]
						+ data->coaxstack[ct->numseq[j - 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				}
			else if (size2 == 0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc(jp, ip, ip - 1, 2, ct, data, false) +
					erg4_nc(i, j, i + 1, 1, ct, data, false) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING                  
				if (energy < energy2) {
					energy = data->init + erg4(jp, ip, ip - 1, 2, ct, data, false) +
						erg4(i, j, i + 1, 1, ct, data, false) + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
				else {
					energy = data->init + data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
#else
				energy = min(energy, energy2);
#endif //COUNTING				  
				}
			return energy;
		}
		else if (b & INTER) {
			//the loop is actually between two strands (ie: intermolecular)
			if (size1 > 1) {//free energy is that of two terminal mismatches
			//and the intermolecular initiation
				energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]];
			}
			else if (size1 == 1) {//find the best terminal mismatch and terminal
			 //stack free energies combination
				energy_option = 1;

				energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]]) +
					erg4_nc(ip, jp, jp + 1, 1, ct, data, false) + penalty_nc(ip, jp, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]]) +
					erg4_nc(i, j, j - 1, 2, ct, data, false) + penalty_nc(i, j, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 2;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING
				//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
				//now consider if coaxial stacking is better:
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[i + 1]][ct->numseq[j - 1]]
						[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 3;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[ip - 1]][ct->numseq[j - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[ip - 1]][ct->numseq[j - 1]]
						[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 4;
					energy = energy2;
				}

				switch (energy_option) {
				case 1:
					energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
						erg4(ip, jp, jp + 1, 1, ct, data, false) + penalty(ip, jp, ct, data);
					break;

				case 2:
					energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
						erg4(i, j, j - 1, 2, ct, data, false) + penalty(i, j, ct, data);
					break;

				case 3:
					energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]
						+ data->coaxstack[ct->numseq[i + 1]][ct->numseq[j - 1]][ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;

				case 4:
					energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip - 1]][ct->numseq[j - 1]]
						+ data->coaxstack[ct->numseq[ip - 1]][ct->numseq[j - 1]][ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING


				}
			else if (size1 == 0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc(jp, ip, jp + 1, 1, ct, data, false) +
					erg4_nc(i, j, j - 1, 2, ct, data, false) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[j]][ct->numseq[i]]
					[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				
				if (energy < energy2) {
					energy = data->init + erg4(jp, ip, jp + 1, 1, ct, data, false) +
						erg4(i, j, j - 1, 2, ct, data, false) + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
				else {
					energy = data->init + data->coax[ct->numseq[j]][ct->numseq[i]]
						[ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
#else
				energy = min(energy, energy2);
#endif //COUNTING				  
				}
			return energy;
		}
	}

	//a typical internal or bulge loop:
	//size1 = ip-i-1;
	//size2 = j - jp - 1;
	//Introduces single stranded pseudoenergies from SHAPE data for interior/bulge loops

	int SHAPEss_energy = 0;

	if (size1 == 1) SHAPEss_energy += ct->SHAPEss_give_value(i + 1);
	else if (size1 != 0) SHAPEss_energy += ct->SHAPEss_calc(i + 1, ip - 1);

	if (size2 == 1) SHAPEss_energy += ct->SHAPEss_give_value(j - 1);
	else if (size2 != 0) SHAPEss_energy += ct->SHAPEss_calc(jp + 1, j - 1);

	energy = erg2_fast(i, j, ip, jp, ct, data, a, b);
	energy += SHAPEss_energy;
	return energy;
}
*/

//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
integersize erg2(int i, int j, int ip, int jp, structure* ct, datatable* data,
	char a, char b)
{

	integersize energy;
	int size, size1, size2, loginc, lopsid, energy2, count, k, energy_option;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
	if (((i <= (ct->GetSequenceLength())) && (ip > (ct->GetSequenceLength()))) || ((
		jp <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength())))) {
		//A loop cannot contain the ends of the sequence

		return INFINITE_ENERGY;
	}



	size1 = ip - i - 1;
	size2 = j - jp - 1;



	if ((a > 0) || (b > 0)) {
		if ((a & DUBLE) || (b & DUBLE)) return INFINITE_ENERGY;//the loop contains a nuc that
			//should be double stranded
		else if ((a & INTER)) {
			//the loop is actually between two strands (ie: intermolecular)

			if (size2 > 1) {//free energy is that of two terminal mismatches
			//and the intermolecular initiation
				energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]];
			}
			else if (size2 == 1) {//find the best terminal mismatch and terminal
			 //stack free energies combination

				energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]) +
					erg4_nc(jp, ip, ip - 1, 2, ct, data, false) + penalty_nc(jp, ip, ct, data);

				energy_option = 1;

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]) +
					erg4_nc(i, j, i + 1, 1, ct, data, false) + penalty_nc(i, j, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 2;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
				   //now consider if coaxial stacking is better:
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[jp + 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 3;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j - 1]][ct->numseq[ip - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[j - 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 4;
					energy = energy2;
				}

				switch (energy_option) {
				case 1:
					energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]
						+ erg4(jp, ip, ip - 1, 2, ct, data, false) + penalty(jp, ip, ct, data);
					break;

				case 2:
					energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]
						+ erg4(i, j, i + 1, 1, ct, data, false) + penalty(i, j, ct, data);
					break;

				case 3:
					energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]]
						+ data->coaxstack[ct->numseq[jp + 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;

				case 4:
					energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j - 1]][ct->numseq[ip - 1]]
						+ data->coaxstack[ct->numseq[j - 1]][ct->numseq[ip - 1]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

			}
			else if (size2 == 0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc(jp, ip, ip - 1, 2, ct, data, false) +
					erg4_nc(i, j, i + 1, 1, ct, data, false) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING                  
				if (energy < energy2) {
					energy = data->init + erg4(jp, ip, ip - 1, 2, ct, data, false) +
						erg4(i, j, i + 1, 1, ct, data, false) + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
				else {
					energy = data->init + data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
#else
				energy = min(energy, energy2);
#endif //COUNTING				  
			}


			return energy;
		}
		else if (b & INTER) {
			//the loop is actually between two strands (ie: intermolecular)

			if (size1 > 1) {//free energy is that of two terminal mismatches
			//and the intermolecular initiation
				energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]];
			}
			else if (size1 == 1) {//find the best terminal mismatch and terminal
			 //stack free energies combination
				energy_option = 1;

				energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]]) +
					erg4_nc(ip, jp, jp + 1, 1, ct, data, false) + penalty_nc(ip, jp, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]]) +
					erg4_nc(i, j, j - 1, 2, ct, data, false) + penalty_nc(i, j, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 2;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING
				//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
				//now consider if coaxial stacking is better:
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[i + 1]][ct->numseq[j - 1]]
						[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 3;
					energy = energy2;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[ip - 1]][ct->numseq[j - 1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[ip - 1]][ct->numseq[j - 1]]
						[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				 
				if (energy2 < energy) {
					energy_option = 4;
					energy = energy2;
				}

				switch (energy_option) {
				case 1:
					energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
						erg4(ip, jp, jp + 1, 1, ct, data, false) + penalty(ip, jp, ct, data);
					break;

				case 2:
					energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
						erg4(i, j, j - 1, 2, ct, data, false) + penalty(i, j, ct, data);
					break;

				case 3:
					energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]
						+ data->coaxstack[ct->numseq[i + 1]][ct->numseq[j - 1]][ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;

				case 4:
					energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip - 1]][ct->numseq[j - 1]]
						+ data->coaxstack[ct->numseq[ip - 1]][ct->numseq[j - 1]][ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
					break;
				}
#else				  
				energy = min(energy, energy2);
#endif //COUNTING


			}
			else if (size1 == 0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc(jp, ip, jp + 1, 1, ct, data, false) +
					erg4_nc(i, j, j - 1, 2, ct, data, false) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[j]][ct->numseq[i]]
					[ct->numseq[ip]][ct->numseq[jp]]) + penalty_nc(i, j, ct, data) + penalty_nc(jp, ip, ct, data);

#ifdef COUNTING				
				if (energy < energy2) {
					energy = data->init + erg4(jp, ip, jp + 1, 1, ct, data, false) +
						erg4(i, j, j - 1, 2, ct, data, false) + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
				else {
					energy = data->init + data->coax[ct->numseq[j]][ct->numseq[i]]
						[ct->numseq[ip]][ct->numseq[jp]] + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
				}
#else
				energy = min(energy, energy2);
#endif //COUNTING				  
			}


			return energy;

		}
	}


	//a typical internal or bulge loop:
	  //size1 = ip-i-1;
	  //size2 = j - jp - 1;

	  //Introduces single stranded pseudoenergies from SHAPE data for interior/bulge loops
	int SHAPEss_energy = 0;

	if (size1 == 1) SHAPEss_energy += ct->SHAPEss_give_value(i + 1);
	else if (size1 != 0) SHAPEss_energy += ct->SHAPEss_calc(i + 1, ip - 1);

	if (size2 == 1) SHAPEss_energy += ct->SHAPEss_give_value(j - 1);
	else if (size2 != 0) SHAPEss_energy += ct->SHAPEss_calc(jp + 1, j - 1);


	if (size1 == 0 || size2 == 0) {//bulge loop


		size = size1 + size2;
		if (size == 1) {
			count = 1;
			energy = data->stack[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[ip]][ct->numseq[jp]]
				+ data->bulge[size] + data->eparam[2];


			if (size1 == 1) {

				//count the number of alternative bulges that exist:

				k = i;
				while (ct->numseq[k] == ct->numseq[i + 1]) {
					count++;
					k--;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() || k == 0) break;
				}

				k = ip;
				while (ct->numseq[k] == ct->numseq[i + 1]) {
					count++;
					k++;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() + 1 || k > 2 * ct->GetSequenceLength()) break;
				}
				//give bonus to C adjacent to single C bulge
				if ((ct->IsNuc(i + 1, 'C') || ct->IsNuc(i + 1, 'c')) && count > 1) energy += data->singlecbulge;

			}

			else {
				//size2 == 1

				//count the number of alternative bulges that exist:

				k = jp;
				while (ct->numseq[k] == ct->numseq[jp + 1]) {
					count++;
					k--;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() || k == 0) break;
				}
				k = j;
				while (ct->numseq[k] == ct->numseq[jp + 1]) {
					count++;
					k++;

					//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
					if (k == ct->GetSequenceLength() + 1 || k > 2 * ct->GetSequenceLength()) break;
				}
				//give bonus to C adjacent to single C bulge
				if ((ct->IsNuc(j - 1, 'C') || ct->IsNuc(j - 1, 'c')) && count > 1) energy += data->singlecbulge;

			}
			//apply a correction for the number of equivalent states because
				//the bulge can move to adjacent sites
			energy -= static_cast<int> (round(data->RT * conversionfactor * log((double)count)));
		}
		else if (size > 30) {
			loginc = static_cast<int>(round((data->prelog) * log(double((size) / 30.0))));
			energy = data->bulge[30] + loginc + data->eparam[2];
			energy = energy + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);

		}
		else {
			energy = data->bulge[size] + data->eparam[2];
			energy = energy + penalty(i, j, ct, data) + penalty(jp, ip, ct, data);
		}
	}
	else {//internal loop
		size = size1 + size2;
		lopsid = abs(size1 - size2);

		if (size > 30) {

				loginc = static_cast<int>(round((data->prelog)*log((double ((size))/30.0))));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING			
				if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
				else {
					energy += lopsid * data->poppen[min(2, min(size1, size2))];
					data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
				}
#else
				energy += min(data->maxpen, lopsid * data->poppen[min(2, min(size1, size2))]);
#endif					
			}

			else {
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i + 1]][ct->numseq[j - 1]] +
					data->tstki[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
					data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING						
				if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
				else {
					energy += lopsid * data->poppen[min(2, min(size1, size2))];
					data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
				}
#else
				energy += min(data->maxpen, lopsid * data->poppen[min(2, min(size1, size2))]);
#endif //COUNTING						

			}
		}
		else if ((size1 == 2) && (size2 == 2))//2x2 internal loop
			energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
			[ct->numseq[j]][ct->numseq[jp]]
			[ct->numseq[i + 1]][ct->numseq[i + 2]]
			[ct->numseq[j - 1]][ct->numseq[j - 2]];


		else if ((size1 == 1) && (size2 == 2)) {//2x1 internal loop
			energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]]
				[ct->numseq[j - 1]][ct->numseq[jp + 1]][ct->numseq[ip]][ct->numseq[jp]];


		}
		else if ((size1 == 2) && (size2 == 1)) {//1x2 internal loop
			energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]]
				[ct->numseq[ip - 1]][ct->numseq[i + 1]][ct->numseq[j]][ct->numseq[i]];

		}

		else if (size == 2) //a single mismatch

			energy = data->iloop11[ct->numseq[i]][ct->numseq[i + 1]][ct->numseq[ip]]
			[ct->numseq[j]][ct->numseq[j - 1]][ct->numseq[jp]];
		else if (size1 == 1 || size2 == 1) { //this loop is lopsided
		//this is a 1xn loop:
			energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];
#ifdef COUNTING					
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
			}
#else	
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
#endif	
		}


		else if ((size1 == 2 && size2 == 3) || (size1 == 3 && size2 == 2)) {
			//this is a 2x3 loop
			energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];

#ifdef COUNTING					
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;
			}
#else	
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
#endif					


		}
		else {



			energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]] +
				data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp + 1]][ct->numseq[ip - 1]] +
				data->inter[size] + data->eparam[3];

#ifdef COUNTING						
			if (NO_COUNT(data->maxpen) < lopsid * NO_COUNT(data->poppen[min(2, min(size1, size2))])) energy += data->maxpen;
			else {
				energy += lopsid * data->poppen[min(2, min(size1, size2))];

				data->poppen[min(2, min(size1, size2))].get += lopsid - 1;

			}
#else
			energy += min(data->maxpen, (lopsid * data->poppen[min(2, min(size1, size2))]));
#endif //COUNTING				
		}
		//energy+=(SHAPEend(i,ct)+SHAPEend(j,ct)+SHAPEend(ip,ct)+SHAPEend(jp,ct));

	}

	//adds SHAPE pseudoenergy to energy value

	energy += SHAPEss_energy;

	return energy;
}





integersize multi_erg2(int i, int j, int ip, int jp, structure* ct, datatable* data,
	char a, char b) 
{
	integersize energysum = 0;
	int N = ct->GetNumberofSequences();
	for (int seqnumber = 0; seqnumber < N; seqnumber++) 
	{	
	//	a b c d e f g  h i j  k  l  m  n  o
	//	1 2 3 4 5 6 7  8 9 10 11 12 13 14 15
	//		i     ip       jp       j

		structure* ct_i = ct->get_individual_sequence(seqnumber);		
		int seq_length = ct_i->GetSequenceLength();
		
		int no_of_nuc_i = ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip];
		int no_of_nuc_j = j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j];

		// if no nucleotides in loop then treat as stack
		if (no_of_nuc_i == 0 && no_of_nuc_j == 0) {
			energysum += erg1(i, j, ip, jp, ct_i, data);
		}
		// if no gaps in internal loop
		else if (no_of_nuc_i == ip - i - 1 && no_of_nuc_j == j - jp - 1) {
			energysum += erg2(i, j, ip, jp, ct_i, data, a, b);
		}
		// if gaps in internal loop
		else if (no_of_nuc_i > 0 || no_of_nuc_j > 0) {	
			int i_  = 1;
			int ip_ = 1 + no_of_nuc_i + 1;
			int jp_ = ip_ + 1;
			int j_  = jp_ + no_of_nuc_j + 1;
			
			#ifdef SMP
			#pragma omp critical
			#endif
			{
				structure* m = remove_gaps(i, j, ip, jp, ct_i, ct->tmp_ct_i);
				//	energysum += erg2_fast(i_, j_, ip_, jp_, m, ct_i->GetThermodynamicDataTable(), a, b);
				//  added 1 to accomodate for the gap introduced by remove_gaps between ip and jp.
				// This will allow sequence to work with alternative bulge loop code in erg2
				energysum += erg2(i_, j_ + 1, ip_, jp_ + 1, m, ct_i->GetThermodynamicDataTable(), a, b);
				ct->tmp_ct_i->deallocate();
			}
			/*
			structure m;
			m.SetThermodynamicDataTable(data);
			m.SetSequence(ct_i->GetSequence());
			energysum += erg2(i_, j_, ip_, jp_, remove_gaps(i, j, ip, jp, &m) , ct_i->GetThermodynamicDataTable(), a, b);
			*/
		}
	}
	return energysum;
}

//calculate the energy of the exterior part of a internal loop
//which includes AU/GU penalty of the interior base pair + first mismatch bonus
integersize erg2ex(int i,int j,int size,structure *ct, datatable *data)
{
	integersize energy;
	int loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   			if (size>30) {

				loginc = static_cast<int>(round((data->prelog)*log((double ((size))/30.0))));
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->inter[30] + loginc;

						}

			else {



         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
					data->inter[size] ;
			}

		//energy+=(SHAPEend(i,ct)+SHAPEend(j,ct));
	return energy;
}

//calculate the energy of the interior part of a large internal loop(size>=6 and 1x(n-1) excluded)
//which includes AU/GU penalty of the interior base pair + first mismatch bonus + assym
//where i is paired to j; ip is paired to jp; ip > i; j > jp
integersize erg2in(int i,int j,int ip,int jp,structure *ct, datatable *data, char a, char b)
{

	integersize energy;
	int size1,size2,lopsid,energy2,energy_option;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
	size1 = ip-i-1;
	size2 = j - jp - 1;
	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return INFINITE_ENERGY;//the loop contains a nuc that
      		//should be double stranded
      	else if ((a&INTER)) {
         	//the loop is actually between two strands (ie: intermolecular)

             	if (size2>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size2==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination
				  energy_option=1;

                  energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]]
                  		[ct->numseq[i+1]][ct->numseq[j-1]]) +
                     	erg4_nc (jp,ip,ip-1,2,ct,data,false)+penalty_nc(jp,ip,ct,data);
                  energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  		[ct->numseq[jp+1]][ct->numseq[ip-1]]) +
                     	erg4_nc (i,j,i+1,1,ct,data,false)+penalty_nc(i,j,ct,data);

                  if (energy2 < energy){
                  	energy_option=2;
                  	energy = energy2;
                  }

                  //if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
                     //now consider if coaxial stacking is better:
                     energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]]
                     	[ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[jp+1]][ct->numseq[ip-1]]
                        [ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

                  if (energy2 < energy){
                  	energy_option=3;
                  	energy = energy2;
                  }

                     energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]]
                     	[ct->numseq[ip]][ct->numseq[j-1]][ct->numseq[ip-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[j-1]][ct->numseq[ip-1]]
                        [ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);


                  if (energy2 < energy){
                  	energy_option=4;
                  	energy = energy2;
                  }

                  switch (energy_option){
                  	case 1:
                  		energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  			[ct->numseq[i+1]][ct->numseq[j-1]] +
                     			erg4 (jp,ip,ip-1,2,ct,data,false)+penalty(jp,ip,ct,data);
                  		break;
                  	case 2:
                  		energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  			[ct->numseq[jp+1]][ct->numseq[ip-1]] +
                     			erg4 (i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data);
                     		break;
                  	case 3:
                  		energy = data->init + data->tstackcoax[ct->numseq[jp]]
                     			[ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]
                        		+ data->coaxstack[ct->numseq[jp+1]][ct->numseq[ip-1]]
                        		[ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                        	break;
                  	case 4:
                  		energy = data->init + data->tstackcoax[ct->numseq[jp]]
                     			[ct->numseq[ip]][ct->numseq[j-1]][ct->numseq[ip-1]]
                        		+ data->coaxstack[ct->numseq[j-1]][ct->numseq[ip-1]]
                        		[ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
							break;
		  		  }
               }
               else if (size2==0) {//just have dangling ends or flush stacking

				energy_option=1;

               	energy = NO_COUNT(data->init) + erg4_nc (jp,ip,ip-1,2,ct,data,false) +
                    	erg4_nc (i,j,i+1,1,ct,data,false)+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
                  energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[ip]][ct->numseq[jp]]
                  	[ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
#ifdef COUNTING                
                 if (energy2 < energy){
                 	energy = data->init + erg4 (jp,ip,ip-1,2,ct,data,false) +
                    	erg4 (i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                 }
                 else {
                 	energy = data->init + data->coax[ct->numseq[ip]][ct->numseq[jp]]
                  	[ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                 }
#else
               energy = min(energy,energy2);
#endif //COUNTING
               }


         		return energy;
      	}
         else if (b&INTER) {
                  	//the loop is actually between two strands (ie: intermolecular)

             	if (size1>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size1==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination

				  energy_option=1;
                  energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]]
                  		[ct->numseq[i+1]][ct->numseq[j-1]]) +
                        erg4_nc (ip,jp,jp+1,1,ct,data,false)+penalty_nc(ip,jp,ct,data);
                  energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  		[ct->numseq[jp+1]][ct->numseq[ip-1]]) +
                        erg4_nc (i,j,j-1,2,ct,data,false)+penalty_nc(i,j,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=2;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
                     //now consider if coaxial stacking is better:
                     
					 energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
                     	[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]]
                        [ct->numseq[ip]][ct->numseq[jp]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
					
#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=3;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

                     energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
                     	[ct->numseq[j]][ct->numseq[ip-1]][ct->numseq[j-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[ip-1]][ct->numseq[j-1]]
                        [ct->numseq[ip]][ct->numseq[jp]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=2;
					energy = energy2;
				  }		 
                  switch (energy_option) {
                  	case 1:
                  		energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  			[ct->numseq[i+1]][ct->numseq[j-1]] +
                        		erg4 (ip,jp,jp+1,1,ct,data,false)+penalty(ip,jp,ct,data);
                  		break;
                  	case 2:
                  		energy = data->init + data->tstackcoax[ct->numseq[i]]
                     			[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
                        		+ data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]]
                        		[ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  		break;
                  	case 3:
                  		energy = data->init + data->tstackcoax[ct->numseq[i]]
                     			[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
                        		+ data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]]
                        		[ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  		break;
                  	case 4:
                  		energy = data->init + data->tstackcoax[ct->numseq[i]]
                     			[ct->numseq[j]][ct->numseq[ip-1]][ct->numseq[j-1]]
                        		+ data->coaxstack[ct->numseq[ip-1]][ct->numseq[j-1]]
                        		[ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  		break;
                  }
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING				  
				  
          
                  //}
               }
               else if (size1==0) {//just have dangling ends or flush stacking
			   energy = NO_COUNT(data->init) + erg4_nc (jp,ip,jp+1,1,ct,data,false) +
                    	erg4_nc (i,j,j-1,2,ct,data,false)+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
                  energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[j]][ct->numseq[i]]
                  	[ct->numseq[ip]][ct->numseq[j]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING
                  if (energy < energy2) {
                  	energy = data->init + erg4 (jp,ip,jp+1,1,ct,data,false) +
                    		erg4 (i,j,j-1,2,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  }
                  else{
                  	energy = data->init + data->coax[ct->numseq[j]][ct->numseq[i]]
                  		[ct->numseq[ip]][ct->numseq[j]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  }
#else               
                  energy = min(energy,energy2);
#endif //COUNTING				  
               }


         		return energy;

         }
      }

   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return INFINITE_ENERGY;
      }





      //a typical internal , no bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;
		if (size1!=0&&size2!=0)

		{//internal loop
			//size = size1 + size2;
			lopsid = abs(size1-size2);
			energy = data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] + data->eparam[3];

		#ifdef COUNTING						
			if (NO_COUNT(data->maxpen) < (lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))]))){
				energy += data->maxpen;
			}
			else{
				energy =+ (lopsid*data->poppen[min(2,min(size1,size2))]);
				
        		data->poppen[min(2,min(size1,size2))].get+=lopsid-1;
			}
		#else
			energy += min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
		#endif				


 }
		//energy+=(SHAPEend(ip,ct)+SHAPEend(jp,ct));
		return energy;
}

/*
//calculate the energy of a hairpin loop:
integersize erg3_fast(int i, int j, structure* ct, datatable* data, char dbl) {
	integersize energy;
	int size, loginc, count, key, k;
	const int alphabetSize = data->alphabet.size();

	size = j - i - 1;
	if (size > 30)
	{
		loginc = int((data->prelog) * log((double((size)) / 30.0)));
		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[30] + loginc + data->eparam[4];
	}
	else if (size < 3)
	{
		energy = data->hairpin[size] + data->eparam[4];
		energy += penalty(i, j, ct, data);
	}
	else if (size == 4)
	{
		//key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*alphabetSize+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 6; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numoftloops; count++) {

			if (key == data->tloop[count][0]) {
				return data->tloop[count][1];
			}
		}

		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}
	else if (size == 3) {

		//key = (ct->numseq[j])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 5; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numoftriloops; count++) {
			if (key == data->triloop[count][0]) return data->triloop[count][1];
		}

		energy = data->hairpin[size] + data->eparam[4] + penalty(i, j, ct, data);
	}
	else if (size == 6) {
		//key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125
		//	+ (ct->numseq[i+4])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 8; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numofhexaloops; count++) {
			if (key == data->hexaloop[count][0]) {
				return data->hexaloop[count][1];
			}
		}

		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}
	else
	{
		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}

	//check for GU closeure preceded by GG
	if (ct->IsNuc(i, 'G') || ct->IsNuc(i, 'g')) {
		if (ct->IsNuc(j, 'U') || ct->IsNuc(j, 'u')) {
			if ((i > 2 && i < ct->GetSequenceLength()) || (i > ct->GetSequenceLength() + 2)) {
				if (ct->IsNuc(i - 1, 'G') || ct->IsNuc(i - 1, 'g')) {
					if (ct->IsNuc(i - 2, 'G') || ct->IsNuc(i - 2, 'g')) {
						energy = energy + data->gubonus;
					}
				}//end if (ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g'))
			}//end if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2))
		}//end if (ct->IsNuc(j,'U')||ct->IsNuc(j,'u')) 
	}//end if (ct->IsNuc(i,'G')||ct->IsNuc(i,'g'))

	//check for an oligo-c loop
	for (k = 1; (k <= size); k++) {
		//this is not an oligo-C loop	
		if (ct->numseq[i + k] != 2) return energy;
	}

	//this is a poly c loop so penalize
	if (size == 3) return (energy + data->c3);
	else return (energy + data->cint + size * data->cslope);
}
*/

/*
//calculate the energy of a hairpin loop:
integersize erg3(int i, int j, structure *ct, datatable *data, char dbl)
{
//	if (ct->GetNumberofSequences() > 1) return multi_erg3(i, j, ct, data, dbl);

	integersize energy;
	int n = ct->GetSequenceLength();
	//A hairpin cannot contain the ends of the sequence
	if ( (i<=n) && (j>n) ) {
         return INFINITE_ENERGY;
    }
	
	//the loop contains a base that should be double stranded
	if (dbl & DUBLE) return INFINITE_ENERGY;
	else if (dbl & INTER)
	{	//  intermolecular interaction
		//  intermolecular "hairpin" free energy is that of intermolecular
		//	initiation plus the stacked mismatch
		energy = data->init;
		#ifdef COUNTING
		if (NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]) < erg4_nc(i, j, i + 1, 1, ct, data, false)) {
			energy += data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]];
		}
		else {
			energy += erg4(i, j, i + 1, 1, ct, data, false);
		}
		energy += penalty(i, j, ct, data);
		#else
		energy += min(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]], erg4(i, j, i + 1, 1, ct, data, false)) + penalty(i, j, ct, data);
		#endif
		
		return energy;
	}

	energy =  erg3_fast(i, j, ct, data, dbl);

	//adds pseudo energy term for a hairpin loop based off of SHAPE data
	if (ct->shaped) energy += ct->SHAPEss_calc(i + 1, j - 1);

	return energy;
}
*/


//calculate the energy of a hairpin loop:
integersize erg3(int i, int j, structure* ct, datatable* data, char dbl)
{
	integersize energy;
	int size, loginc, count, key, k;
	const int alphabetSize = data->alphabet.size();


	if ((i <= (ct->GetSequenceLength())) && (j > (ct->GetSequenceLength()))) {
		//A hairpin cannot contain the ends of the sequence
		return INFINITE_ENERGY;
	}

	if (dbl & DUBLE) return INFINITE_ENERGY;//the loop contains a base that should be
											//double stranded

	else if (dbl & INTER) {//intermolecular interaction
		//intermolecular "hairpin" free energy is that of intermolecular
		 //	initiation plus the stacked mismatch
		energy = data->init;

#ifdef COUNTING
		if (NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]]) < erg4_nc(i, j, i + 1, 1, ct, data, false)) {
			energy += data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]];
		}
		else {
			energy += erg4(i, j, i + 1, 1, ct, data, false);
		}
		energy += penalty(i, j, ct, data);
#else
		energy += min(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i + 1]][ct->numseq[j - 1]], erg4(i, j, i + 1, 1, ct, data, false)) + penalty(i, j, ct, data);
#endif
		return energy;
	}




	size = j - i - 1;



	if (size > 30) {

			loginc = static_cast<int>(round((data->prelog)*log((double ((size))/30.0))));



		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[30] + loginc + data->eparam[4];
	}
	else if (size < 3) {
		energy = data->hairpin[size] + data->eparam[4];
		energy += penalty(i, j, ct, data);
	}
	else if (size == 4) {

		//key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*alphabetSize+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 6; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numoftloops; count++) {

			if (key == data->tloop[count][0]) {
				return data->tloop[count][1];
			}
		}

		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}
	else if (size == 3) {

		//key = (ct->numseq[j])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 5; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numoftriloops; count++) {
			if (key == data->triloop[count][0]) return data->triloop[count][1];
		}

		energy = data->hairpin[size] + data->eparam[4]
			+ penalty(i, j, ct, data);
	}
	else if (size == 6) {
		//key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125
		//	+ (ct->numseq[i+4])*625 +
		//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
		int digit_factor = 1;
		key = 0;
		for (count = 0; count < 8; ++count) {
			key += ct->numseq[i + count] * digit_factor;
			digit_factor = digit_factor * alphabetSize;
		}
		for (count = 0; count < data->numofhexaloops; count++) {
			if (key == data->hexaloop[count][0]) {
				return data->hexaloop[count][1];
			}
		}

		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}

	else {
		energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
			[ct->numseq[i + 1]][ct->numseq[j - 1]]
			+ data->hairpin[size] + data->eparam[4];
	}




	//check for GU closeure preceded by GG

	if (ct->IsNuc(i, 'G') || ct->IsNuc(i, 'g')) {
		if (ct->IsNuc(j, 'U') || ct->IsNuc(j, 'u')) {


			if ((i > 2 && i < ct->GetSequenceLength()) || (i > ct->GetSequenceLength() + 2)) {

				if (ct->IsNuc(i - 1, 'G') || ct->IsNuc(i - 1, 'g')) {
					if (ct->IsNuc(i - 2, 'G') || ct->IsNuc(i - 2, 'g')) {


						energy = energy + data->gubonus;


					}
				}//end if (ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g'))
			}//end if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2))
		}//end if (ct->IsNuc(j,'U')||ct->IsNuc(j,'u')) 
	}//end if (ct->IsNuc(i,'G')||ct->IsNuc(i,'g'))

  //adds pseudo energy term for a hairpin loop based off of SHAPE data
	if (ct->shaped) energy += ct->SHAPEss_calc(i + 1, j - 1);

	//check for an oligo-c loop

	for (k = 1; (k <= size); k++) {
		if (ct->numseq[i + k] != 2) return energy;//this is not an oligo-C loop
	}
	//this is a poly c loop so penalize
	if (size == 3) return (energy + data->c3);
	else return (energy + data->cint + size * data->cslope);


}


// Takes structure ct with multiple sequences containing gaps and 1-indexed loop indices [i,j] inclusive 
// returns sum of hairpin loop energy across all the sequences
// assuming i and j are base paired
integersize multi_erg3(int i, int j, structure* ct, datatable* data, char dbl)
{
	//A hairpin cannot contain the ends of the sequence
	if ((i <= (ct->get_individual_sequence(0)->GetSequenceLength())) && 
		(j > (ct->get_individual_sequence(0)->GetSequenceLength()))) {
		return INFINITE_ENERGY;
	}

	integersize energysum = 0;

	int n = ct->GetNumberofSequences();
	for (int seqnumber = 0; seqnumber < n; ++seqnumber) {

		structure* ct_i = ct->get_individual_sequence(seqnumber);
		int no_of_nuc = j - i - 1 - ct_i->no_of_gaps_matrix[i][j];
		// if number of nucleotides in hairpin loop < 3
		if (no_of_nuc < 3) { 
			// global variable hairpinlt3 = 8 * conversion factor in defines.h
			energysum += hairpinlt3;
		}
		// if no gaps
		else if (no_of_nuc == j - i - 1) { 
			energysum += erg3(i, j, ct_i, ct_i->GetThermodynamicDataTable(), dbl); 
		}
		// if gaps in hairpin loop
		else { 
			#ifdef SMP
			#pragma omp critical
			#endif
			{
				int new_i = 1;
				int new_j = 1 + no_of_nuc + 1;
				energysum += erg3(new_i, new_j, remove_gaps(i, j, ct_i, ct->tmp_ct_hp), data, dbl);
				ct->tmp_ct_hp->deallocate();
			}
			
			/*
			structure temp_ct;
			temp_ct.SetThermodynamicDataTable(ct->GetThermodynamicDataTable());
			temp_ct.SetSequence(ct_i->GetSequence());
			energysum += erg3_fast(new_i, new_j, remove_gaps(i, j, &temp_ct), data, dbl);
			*/
		}
	}
	return energysum;
}

/*
// Returns structure with sequence not containing gaps in substring [i,j] of sequence in structure ct
// i and j are 1-indexed
// i and j are base paired
structure* remove_gaps(int i, int j, structure* ct) {
	if (i > ct->GetSequenceLength()) { i = i - ct->GetSequenceLength(); j = j - ct->GetSequenceLength(); }
	string str = ct->GetSequence();
	string loop_str = str.substr(i, j - i - 1);
	loop_str.erase(remove(loop_str.begin(),  loop_str.end(), '-'), loop_str.end());
	string new_str = str[i - 1] + loop_str + str[j - 1];
	ct->deallocate();
	ct->SetSequenceFast(new_str);
	return ct;
}
*/

// Returns structure with sequence not containing gaps in substring [i,j] of sequence in structure ct
// i and j are 1-indexed
// i and j are base paired
structure* remove_gaps(int i, int j, structure* ct_i, structure* tmp_ct) {	
	
	if (i > ct_i->GetSequenceLength()) { i = i - ct_i->GetSequenceLength(); j = j - ct_i->GetSequenceLength(); }
	
	ct_i->loop_str = ct_i->sequence_without_gaps->GetSequence();
	ct_i->loop_str.insert(0, 1, 'x');
	string loop_str_i		= ct_i->loop_str.substr(ct_i->gapped_to_ungapped_map[i+1], j - i - 1 - ct_i->no_of_gaps_matrix[i][j]);
	
	int new_numseq_len		= j - i - 1 - ct_i->no_of_gaps_matrix[i][j] + 2 + 1;
	tmp_ct->allocate(new_numseq_len - 1);
    tmp_ct->numseq[1] = ct_i->numseq[i];
	tmp_ct->numseq[new_numseq_len - 1] = ct_i->numseq[j];
	std::copy_n(ct_i->sequence_without_gaps->numseq + ct_i->gapped_to_ungapped_map[i + 1], new_numseq_len - 2 - 1, tmp_ct->numseq + 2);
	
	const char* str = ct_i->GetSequence();
	tmp_ct->SetTmpSequence(str[i - 1] + loop_str_i + str[j - 1]);
	
	ct_i->loop_str.clear();
	return tmp_ct;
}

/*
// Returns structure with sequence not containing gaps in substring (i,ip) and (j, jp) of 
// sequence in structure ct
// i, j, ip and jp are 1-indexed
// removes gaps from internal loop closed by base pairs i-j and ip-jp
// i < ip < jp < j
// i, ip, j, jp not part of the loop region
structure* remove_gaps(int i, int j, int ip, int jp, structure* ct) {
	//	a b c d e f g  h i j  k  l  m  n  o
	//	1 2 3 4 5 6 7  8 9 10 11 12 13 14 15
	//		i     ip       jp       j

	string str = ct->GetSequence();
	str = str + str;
	
	string loop_str_i = "";
	if (ip - i > 1) {
		loop_str_i = str.substr(i, ip - i - 1);
		loop_str_i.erase(remove(loop_str_i.begin(), loop_str_i.end(), '-'), loop_str_i.end());
	}
	
	string loop_str_j = "";
	if (j - jp > 1) {
		loop_str_j = str.substr(jp, j - jp - 1);
		loop_str_j.erase(remove(loop_str_j.begin(), loop_str_j.end(), '-'), loop_str_j.end());
	}

	string new_str = str[i-1] + loop_str_i + str[ip - 1] + str[jp - 1] + loop_str_j + str[j - 1];

	ct->deallocate();
	ct->SetSequenceFast(new_str);

	return ct;
}
*/

/*
structure* remove_gaps(int i, int j, int ip, int jp, structure* ct_i, structure* tmp_ct_i) {
	//	a b c d e f g  h i j  k  l  m  n  o
	//	1 2 3 4 5 6 7  8 9 10 11 12 13 14 15
	//		i     ip       jp       j

	
	int new_numseq_len =  (ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 2) 
						+ (j - jp - 1 - ct_i->no_of_gaps_matrix[j][jp] + 2) 
						+ 1;
	
	tmp_ct_i->allocate(new_numseq_len-1);

	tmp_ct_i->numseq[1] = ct_i->numseq[i];
	tmp_ct_i->numseq[1 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1] = ct_i->numseq[ip];
	tmp_ct_i->numseq[1 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1] = ct_i->numseq[jp];
	tmp_ct_i->numseq[new_numseq_len - 1] = ct_i->numseq[j];
	
	int a, b, c, d;
	int seq_len = ct_i->GetSequenceLength();
	const char* str = ct_i->GetSequence();
	(i  > seq_len)  ? a = i  - seq_len : a = i;
	(ip > seq_len)  ? b = ip - seq_len : b = ip;
	(jp > seq_len)  ? c = jp - seq_len : c = jp;
	(j  > seq_len)  ? d = j  - seq_len : d = j;

	ct_i->loop_str = ct_i->sequence_without_gaps->GetSequence();
	ct_i->loop_str.insert(0, 1, 'x');

	string loop_str_j = "";
	if (j - jp > 1) {
		loop_str_j = ct_i->loop_str.substr(ct_i->gapped_to_ungapped_map[jp + 1], j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j]);
		std::copy_n(
			ct_i->sequence_without_gaps->numseq + ct_i->gapped_to_ungapped_map[jp + 1],
			j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j],
			tmp_ct_i->numseq + 2 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1);
	}

	string loop_str_i = "";
	if (ip - i > 1) {
		loop_str_i = ct_i->loop_str.substr(ct_i->gapped_to_ungapped_map[i + 1], ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip]);
		std::copy_n(
			ct_i->sequence_without_gaps->numseq + ct_i->gapped_to_ungapped_map[i + 1],
			ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip],
			tmp_ct_i->numseq + 2);
	}
	
	tmp_ct_i->SetTmpSequence(str[a - 1] + loop_str_i + str[b - 1] + str[c - 1] + loop_str_j + str[d - 1]);
	ct_i->loop_str.clear();

	return tmp_ct_i;
}
*/


structure* remove_gaps(int i, int j, int ip, int jp, structure* ct_i, structure* tmp_ct_i) {
	//	a b c d e f g  h i j  k  l  m  n  o
	//	1 2 3 4 5 6 7  8 9 10 11 12 13 14 15
	//		i     ip       jp       j

	int new_numseq_len = (ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 2)
		+ 1 // for linker I, to work with alternative bulges code in erg2
		+ (j - jp - 1 - ct_i->no_of_gaps_matrix[j][jp] + 2)
		+ 1;

	tmp_ct_i->allocate(new_numseq_len - 1);

	tmp_ct_i->numseq[1] = ct_i->numseq[i];
	tmp_ct_i->numseq[1 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1] = ct_i->numseq[ip];
	tmp_ct_i->numseq[1 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1] = 5; // linker I numseq
	tmp_ct_i->numseq[1 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1 + 1] = ct_i->numseq[jp];
	tmp_ct_i->numseq[new_numseq_len - 1] = ct_i->numseq[j];

	int a, b, c, d;
	int seq_len = ct_i->GetSequenceLength();
	const char* str = ct_i->GetSequence();
	(i > seq_len) ? a = i - seq_len : a = i;
	(ip > seq_len) ? b = ip - seq_len : b = ip;
	(jp > seq_len) ? c = jp - seq_len : c = jp;
	(j > seq_len) ? d = j - seq_len : d = j;

	ct_i->loop_str = ct_i->sequence_without_gaps->GetSequence();
	ct_i->loop_str.insert(0, 1, 'x');

	string loop_str_j = "";
	if (j - jp > 1) {
		loop_str_j = ct_i->loop_str.substr(ct_i->gapped_to_ungapped_map[jp + 1], j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j]);
		std::copy_n(
			ct_i->sequence_without_gaps->numseq + ct_i->gapped_to_ungapped_map[jp + 1],
			j - jp - 1 - ct_i->no_of_gaps_matrix[jp][j],
		//	tmp_ct_i->numseq + 2 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1);
			tmp_ct_i->numseq + 2 + ip - i - 1 - ct_i->no_of_gaps_matrix[ip][i] + 1 + 1 + 1);
	}

	string loop_str_i = "";
	if (ip - i > 1) {
		loop_str_i = ct_i->loop_str.substr(ct_i->gapped_to_ungapped_map[i + 1], ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip]);
		std::copy_n(
			ct_i->sequence_without_gaps->numseq + ct_i->gapped_to_ungapped_map[i + 1],
			ip - i - 1 - ct_i->no_of_gaps_matrix[i][ip],
			tmp_ct_i->numseq + 2);
	}

//	tmp_ct_i->SetTmpSequence(str[a - 1] + loop_str_i + str[b - 1] + str[c - 1] + loop_str_j + str[d - 1]);
	tmp_ct_i->SetTmpSequence(str[a - 1] + loop_str_i + str[b - 1] + "I" + str[c - 1] + loop_str_j + str[d - 1]);

	ct_i->loop_str.clear();

	return tmp_ct_i;
}


// returns number of gaps in range [i,j] in sequence of structure ct
// i and j are 1-indexed
// i and j are base paired
int no_of_gaps(int i, int j, structure* ct) {
	if (i > j) { return 0; }
	string str = ct->GetSequence();
	int count = 0;
	for (int p = i; p < j - 1; p++) {
		if (str[p] == '-') {
			count += 1;
		}
	}
	return count;
}

//calculate the energy of a dangling end:
integersize erg4(int i, int j, int ip, int jp, structure *ct, datatable *data, bool lfce)
{	
//	if (ct->GetNumberofSequences() > 1) return multi_erg4(i, j, ip, jp, ct, data, lfce); 

	integersize energy;

	//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle

      if (lfce) return INFINITE_ENERGY;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker

		energy = data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];

		//also add SAHPE single stranded contribution:
		energy+=ct->SHAPEss_give_value(ip);

		return energy;

}

integersize multi_erg4(int i, int j, int ip, int jp, structure* ct, datatable* data, bool lfce)
{
	int n = ct->GetNumberofSequences();	
	integersize energysum = 0;
	for (int seqnumber = 0; seqnumber < n; ++seqnumber) {
		energysum += erg4(i, j, ip, jp, ct->get_individual_sequence(seqnumber), ct->get_individual_sequence(seqnumber)->GetThermodynamicDataTable(), lfce);
	}
	return energysum;
}

integersize erg4_nc(int i,int j,int ip,int jp,structure *ct, datatable *data, bool lfce)
{
	integersize energy;

	#ifdef	COUNTING
	int tmp_count = data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp].get;
	#endif

	energy = erg4(i,j,ip,jp,ct,data,lfce);

	#ifdef	COUNTING
	data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp].get = tmp_count;
	#endif

	return energy;
}



//  When considering mismatch at the end of a helix, consult this 
//	function to check whether the nucs are required to pair
int checknp(bool lfce1, bool lfce2) 
{
	if (lfce1 || lfce2) return INFINITE_ENERGY;
	else return 0;
}

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
integersize ergcoax(int i, int j, int ip, int jp, int k, structure *ct, datatable *data) 
{


	if (ip == j+1) {
		//flush stacking

		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];

	}

	else if (k>0) {
		//coaxial stacking with an intervening mismatch
		if (k==i-1) {
			return data->tstackcoax[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[i-1]] +
				data->coaxstack[ct->numseq[j+1]][ct->numseq[k]][ct->numseq[ip]][ct->numseq[jp]];

		}
		else { //if (k==jp+1) {
			return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
				data->coaxstack[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[k]];
		}
		//else {
			//some error -- give message
			//errmsg(100,1);

		//}

	}
	else return INFINITE_ENERGY;



}
//these functions calculate the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k (determined by the function called) indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
integersize ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data) 
{
//	if (ct->GetNumberofSequences() > 1) return multi_ergcoaxflushbases(i, j, ip, jp, ct, data);
	// flush stacking, k==0
	// Remapped 10/11/2013 to match the order of the helical stack table
	return data->coax[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][ct->numseq[jp]];
}
integersize multi_ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data)
{
	integersize energy = 0;
	int N = ct->number_of_sequences;
	for (int n = 0; n < N; n++) 
	{	
		energy = energy + ergcoaxflushbases(i, j, ip, jp, ct->get_individual_sequence(n), ct->get_individual_sequence(n)->GetThermodynamicDataTable());
	}
	return energy;
}

integersize ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data) 
{
//	if (ct->GetNumberofSequences() > 1) return multi_ergcoaxinterbases1(i, j, ip, jp, ct, data);
	return data->tstackcoax[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[i-1]] 
				+ data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]][ct->numseq[ip]][ct->numseq[jp]] 
				+ ct->SHAPEss_give_value(j+1) 
				+ ct->SHAPEss_give_value(i-1);
}
integersize multi_ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data)
{
	integersize energy = 0;
	int N = ct->number_of_sequences;
	for (int n = 0; n < N; n++)
	{
		energy = energy + ergcoaxinterbases1(i, j, ip, jp, ct->get_individual_sequence(n), ct->get_individual_sequence(n)->GetThermodynamicDataTable());
	}
	return energy; 
}

integersize ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//k==jp+1
	//if (ct->GetNumberofSequences() > 1) return multi_ergcoaxinterbases2(i, j, ip, jp, ct, data);
	return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
				data->coaxstack[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[jp+1]] +
				ct->SHAPEss_give_value(jp+1) + ct->SHAPEss_give_value(ip-1);

}
integersize multi_ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data)
{
	integersize energy = 0;
	int N = ct->number_of_sequences;
	for (int n = 0; n < N; n++)
	{
		energy = energy + ergcoaxinterbases2(i, j, ip, jp, ct->get_individual_sequence(n), ct->get_individual_sequence(n)->GetThermodynamicDataTable());
	}
	return energy;
}

//	this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair

// Note that this form of the function takes the sequences, not the number of the nuc
// in i, j, ip, and jp.
integersize ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data) 
{
		// Look up the energy in the coax array.  
		// This was remapped 10/10/2013 to be consistent in order with the stack array.
		return data->coax[j][i][ip][jp];
}

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i-k jp
//i.e. a discontinuity between k and jp

//Note that this form of the function takes the sequences, not the number of the nuc
//	in i, j, ip, jp,k and l.
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data) 
{
	return data->tstackcoax[j][i][l][k] + data->coaxstack[l][k][ip][jp];
}

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i k-jp
//i.e. a discontinuity between i and jp


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,jp,k and l.
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data) {

		//coaxial stacking with an intervening mismatch

			return data->tstackcoax[jp][ip][k][l] +
				data->coaxstack[j][i][l][k];



}

//These functions are used by ergmulti to deconvolute the base pair code back into nucs
int decon1(int x,int alphabetsize) {

	return (int) (floor(((float)(x))/((float)(alphabetsize+1)))-1);

}
int decon2(int x,int alphabetsize) {

	return x - (decon1(x,alphabetsize)+1)*(alphabetsize+1)-1;
}

//This function will calculate the free energy of a multiloop starting at nuc i for
//	structure #st
//This uses a recursive algorithm

#define js true //this switch can turn on and off the logarithmic dep. on
				//unpaired nucleotides.  This is helpful for debugging

//simplemb, when true, indicates that the logarithmic dependence should be off.
//This sets the energy function equal to that used by the dynamic programming
//algorithms.

integersize ergmulti(int st, int ip, structure *ct, datatable *data, bool simplemb) {
	//short int *element,**energy;
	short int* element;
	//short int i,count,b,c,size,j,current,first,recent,biggest,minimum;
	short int i, count, b, c, j, current, first, recent, biggest, minimum;
	bool intermolecular;
	float average;
	integersize auendpenalties;//track the sum of AU heix end penalties

	const int alphabetSize = data->alphabet.size();

//	//array for multibranch loop calculation
	int **A;
	
	//trace out the loop to learn the size:
	count = 0;
	i = ip;

	while(i!=ip||count==0) {
		i++;
		if (ct->GetPair(i,st)!=0) {
			i=ct->GetPair(i,st);
		}

		count++;
	}

//	//Initialize array
	A = new int *[4];
	for (j=0;j<4;j++)
		A[j]=new int [count+1];

//	//Set array to zero
	for (i=0;i<count;i++) {
		for (j=0;j<4;j++)
			A[j][i]=0;
	}
	
	//allocate element to store this info:
	element = new short int [count+4];

	biggest = 0;
	average = 0;
	b=0;//keep track of the number of unpaired nucs
	c=0;//keep track of the number of helixes
	auendpenalties=0;//keep track of the number of terminal AU/GU pairs
	intermolecular = false; //catch whether this multiloop contains the
							//intermolecular linker
							//if so, there are no multiloop penalties

	//record the info:

	element[0] = (alphabetSize+1)*(ct->numseq[ct->GetPair(ip,st)]+1);
	element[0] = element[0] + ct->numseq[ip]+1;

	//if (ct->numseq[ct->basepr[st][ip]]==4) au++;
	current = 0;
	i = ip;
	count = 1;
	while (i!=ip||count==1) {
		i++;
		if (ct->GetPair(i,st)>0) {

			if(c>0) {
				if(abs(current-recent)>biggest) biggest=abs(current-recent);
				average = average + (float) (abs(current-recent));
			}

			else first = current;
			recent = current;
			current=0;



			c++;
			//store pairs in a code:
			element[count] = (alphabetSize+1)*(ct->numseq[i]+1);
			auendpenalties+=penalty(i,ct->GetPair(i,st),ct,data);
			i=ct->GetPair(i,st);
			element[count]=element[count]+ct->numseq[i]+1;
			//if (ct->numseq[i]==4) au++;





		}
		else {
			current++;
			b++;
			element[count] = ct->numseq[i];
			if (data->isLinker(element[count])) intermolecular=true;


		}

		count++;



	}

	if(abs(current-first)>biggest) biggest=abs(current-first);
	average = average + (float) (abs(first-recent));

	average = average/((float) (c));

	element[count] = element[1];
	element[count+1] = element[2];
	element[count+2] = element[3];
	
	count--;  //correct the count elements
				//bcs the first helix was counted twice
	
	//********************************************************************
	// Start Minimization
	//********************************************************************

	int k;

//	//Initialize first values of the array
	for (i=0;i<4;i++)
		A[i][0]=0;

//	//Fill array for each of the four phases
	for (i=0;i<4;i++)	{
		for (j=1;j<=count;j++){
			// By default set array member to be equal to previous member
			A[i][j]=A[i][j-1];

			// Shifted element position
			k=i+j-1;

			//if the element is an unpaired base
			if (element[k]<(alphabetSize+1)) {
				//Add 3' Dangle
				if (j>=2){
					if (element[k-1]>(alphabetSize+1)) 
						A[i][j]=min(A[i][j],A[i][j-2]+NO_COUNT(data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1]));
				}
				//Add terminal mismatch
				if (j>=3) {
					if (element[k-2]<(alphabetSize+1) && element[k-1]>(alphabetSize+1))
						A[i][j]=min(A[i][j],A[i][j-3]+NO_COUNT(data->tstkm[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]));
				}

				//Add mismatch-mediated coaxial stacking
				if (j>=4) {
					if (element[k-3]>(alphabetSize+1) && element[k-1]>(alphabetSize+1) && element[k-2]<(alphabetSize+1))
						A[i][j]=min(A[i][j],A[i][j-4]+NO_COUNT(data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]])+NO_COUNT(data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]));
				}
			}
			else { //if the element is in a base pair
				//Add 5' Dangle
				if (j>=2) {
					if (element[k-1]<(alphabetSize+1))
						A[i][j]=min(A[i][j],A[i][j-2]+NO_COUNT(data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2]));
				}
				//Add flush coaxial
				if (j>=2) {
					if (element[k-1]>(alphabetSize+1))
						A[i][j]=min(A[i][j],A[i][j-2]+NO_COUNT(data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)]));
				}
				//Add mismatch-mediated coaxial stacking
				if (j>=4) {
					if (element[k-2]>(alphabetSize+1) && element[k-3]<(alphabetSize+1) && element[k-1]<(alphabetSize+1))
						A[i][j]=min(A[i][j],A[i][j-4]+data->NO_COUNT(coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]])+NO_COUNT(data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]]));
				}
			}

		}
	}

	// Calculate minimum 
	minimum=min(A[0][count],A[1][count]);
	minimum=min(minimum,A[2][count]);
	minimum=min(minimum,A[3][count]);

#ifdef COUNTING
	//********************************************************************
	// Start Traceback. Traceback is not necessary, can be commented out
	//********************************************************************

	// Find the phase
	if (minimum==A[0][count])
		i=0;
	else if (minimum==A[1][count])
		i=1;
	else if (minimum==A[2][count])
		i=2;
	else
		i=3;

	// Start at the end of the array
	j=count;

	// Keep track of current energy
	short int next_trace = minimum;

	short int foo;

	// Used to limit one energy contribution per loop iteration, could also use break statements
	bool ntraced = true;

	while (j >0) {
		// Calculate element position
		k=i+j-1;

		//if the element is an unpaired base
		if (element[k]<10) {
			// Consider 3' Dangle
			if (j>=2){
				if (element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace == A[i][j-2]+data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1]._value) {
						next_trace=A[i][j-2];
						ntraced=false;
						j=j-2;
						foo=data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1];
					}
				}
			}
			// Consider terminal mismatch
			if (j>=3) {
				if (element[k-2]<(alphabetSize+1) && element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace==A[i][j-3]+data->tstkm[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]._value) {
						next_trace=A[i][j-3];
						ntraced=false;
						j=j-3;
						foo=data->tstkm[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]];
					}
				}
			}

			// Consider mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-3]>(alphabetSize+1) && element[k-1]>(alphabetSize+1) && element[k-2]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[i][j-4]+data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]]._value+data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]._value) {
						next_trace=A[i][j-4];
						ntraced=false;
						j=j-4;
						foo=data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]]+data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]];
					}
				}
			}
		}
		else { //if the element is in a base pair
			// Consider 5' Dangle
			if (j>=2) {
				if (element[k-1]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[i][j-2]+data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2]._value) {
						next_trace=A[i][j-2];
						ntraced=false;
						j=j-2;
						foo=data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2];
					}
				}
			}
			// Consider flush coaxial
			if (j>=2) {
				if (element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace==A[i][j-2]+data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)]._value) {
						next_trace=A[i][j-2];
						ntraced=false;
						j=j-2;
						foo=data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)];
					}
				}
			}
			// Consider mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-2]>(alphabetSize+1) && element[k-3]<(alphabetSize+1) && element[k-1]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[i][j-4]+data->coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]]._value+data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]]._value) {
						next_trace=A[i][j-4];
						ntraced=false;
						j=j-4;
						foo=data->coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]]+data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]];
					}
				}
			}
		}
		
		// If not other conditions are met, do default action
		if (ntraced) {
			next_trace=A[i][j-1];
			j=j-1;
		}

		// Reset ntraced
		ntraced=true;
	}
	//--------------------------------------------------------------------
	//  End Traceback
	//--------------------------------------------------------------------
#endif //COUNTING

	//deallocate memory use:
	for (i=0;i<4;i++) delete[] A[i];
	delete[] A;
	delete[] element;


	//return the energy:
	if (intermolecular)
		return minimum+auendpenalties+data->init;

	//do not add strin term for simplemb==true
	if (((c%2)!=0)&&((b==0)||(b==1))&&(!simplemb)) minimum = minimum+data->strain;

	//do not use the assymetry for simplemb==true
	if (simplemb) average = 0;
	else if (average>2) average =2;

	minimum = minimum+(short) (((float) data->mlasym)*average+0.5);
	//Correct the counts
#ifdef COUNTING
	data->mlasym.get+=average-1;
#endif	

	//note that this function, used by efn2, has a logarithmic dependence in
	//unpaired nucleotides after 8.  This is hardwired below.
	if (b>8&&js&&!intermolecular&&!simplemb){
//		cout << endl << minimum + data->efn2a + c*data->efn2c + 8*(data->efn2b)	+ int(11.*log(double(((double (b))/8.))) + 0.5) + au*data->auend << endl;

		minimum= minimum + data->efn2a + c*data->efn2c +
		8*(data->efn2b)
			+ int(11.*log(double(((double (b))/8.))) + 0.5) + auendpenalties;

		//Correct parameter usage
#ifdef COUNTING
		data->efn2c.get+=c-1;
		data->efn2b.get+=7;			
#endif		
	}
      else{
		minimum = minimum + data->efn2a + b*data->efn2b + c*data->efn2c
			+auendpenalties;
#ifdef COUNTING			
		data->efn2c.get+=c-1;
		data->efn2b.get+=b-1;			
#endif
	  }
	return minimum;
}

// Energyout: writes to file a list of energys calculated by efn2
void energyout(structure *ct,char *energyfile) {
int i;
ofstream out(energyfile);

for (i=1;i<=ct->GetNumberofStructures();i++)
	out << "Structure: "<<i<<"   Energy = "<<(float (ct->GetEnergy(i))/10)<<"   \n";

}


//This function will calculate the free energy of an exterior loop
//  for structure #st
//This uses a recursive algorithm

//Zuber ergex
integersize ergexterior(int st, structure *ct, datatable *data, int start/*=1*/, int stop/*=0*/) {
	//short int* element, ** energy;
	//short int i, count, size, j, minimum, helices;
	short int* element;
	short int i, count, j, minimum, helices;
	bool intermolecular;
	integersize auendpenalties;//track the sum of terminal AU penalties
	short int *A;
	int k;
	//trace out the loop to learn the size:
	count = 0;///DALIA: the count of nt in the external loop
	//start counting at 1
	i = start-1;
	if (stop==0)
		stop = ct->GetSequenceLength();

	const int alphabetSize = data->alphabet.size();

	intermolecular = false; //keep track as to whether the bimolecular linker is found
							//	if so, add the intermolecular initiation free energy

	helices  = 0;//DALIA: count of the numbr of helices encountered in external loop
	//DALIA: Walk through loop and count nt and helices
	while(i!=ct->GetSequenceLength() && i <stop) {
		i++;
		if (ct->GetPair(i,st)>0) {
			i=ct->GetPair(i,st);
			helices++;
		}
		count++;
	}

	//check for empty structure and return 0 if empty
	if (helices==0) return 0;//DALIA: dont waste time dong anything else, no structure, dG=0

	//allocate element to store this info:
	element = new short int [count];//DALIA: this is array with an element for each nt in the external loop
	
	auendpenalties=0;//keep track of the sum of terminal AU/GU pair penalties

	//record the info:
	i = start-1;
	count = 0;
	// walk through again the external loop, this time filling in the element array
	while (i!=ct->GetSequenceLength() && i < stop) {
		i++;

		if (ct->GetPair(i,st)>0) {
			//store pairs in a code:
			element[count] = (alphabetSize+1)*(ct->numseq[i]+1);
			auendpenalties+=penalty(i,ct->GetPair(i,st),ct,data);
			i=ct->GetPair(i,st);
			element[count]=element[count]+ct->numseq[i]+1;
		}
		else {
			element[count] = ct->numseq[i];
			if (data->isLinker(element[count])) intermolecular=true;
		}

		count++;
	}


	// Allocate the minimization array
	A = new short int [count+1];
	
	// Initialize first value of the array
	A[0]=0;

	//********************************************************************
	// Fill array
	//********************************************************************
	for (j=1;j<=count;j++){
		// Fill default value
		A[j]=A[j-1];

		// Shifted index for element array, not strictly necessary
		k=j-1;

		//if the element is an unpaired base
		if (element[k]<(alphabetSize+1)) {//DALIA: remember the paired bases were multiplied by alphabetSize+1
			//Add 3' Dangle
			if (j>=2){
				if (element[k-1]>(alphabetSize+1)) //Check for dangle requirements
					A[j]=min(A[j],A[j-2]+NO_COUNT(data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1]));
			}
			//Add terminal mismatch
			if (j>=3) {
				if (element[k-2]<(alphabetSize+1) && element[k-1]>(alphabetSize+1)) // Check for mismatch requirements
					A[j]=min(A[j],A[j-3]+NO_COUNT(data->tstack[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]));
			}

			//Add mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-3]>(alphabetSize+1) && element[k-1]>(alphabetSize+1) && element[k-2]<(alphabetSize+1)) // Check for mismatch-mediated coaxial stacking requirements
					A[j]=min(A[j],A[j-4]+NO_COUNT(data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]])+NO_COUNT(data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]));
			}
		}
		else { //if the element is in a base pair
			//Add 5' Dangle
			if (j>=2) { 
				if (element[k-1]<(alphabetSize+1)) // Check for dangle requirements
					A[j]=min(A[j],A[j-2]+NO_COUNT(data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2]));
			}
			//Add flush coaxial
			if (j>=2) { 
				if (element[k-1]>(alphabetSize+1)) //Check for flush coaxial stacking requirements
					A[j]=min(A[j],A[j-2]+NO_COUNT(data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)]));
					//cout << "flush coax = " << min(A[j],A[j-2]+data->coax[decon1(element[k-1],alphabetSize)][decon2(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)]) << "  " << flush;
			}
			//Add mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-2]>(alphabetSize+1) && element[k-3]<(alphabetSize+1) && element[k-1]<(alphabetSize+1)) // Check for mismatch-mediated coaxial stacking requirements
					A[j]=min(A[j],A[j-4]+NO_COUNT(data->coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]])+NO_COUNT(data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]]));
			}
		}

	}

	// Get the minimum energy
	minimum = A[count];

#ifdef COUNTING	
	//********************************************************************
	//Start Traceback. Traceback is not necessary, can be commented out
	//********************************************************************

	// Start at the end of the minimization array
	j=count;

	// Keep track of current energy
	short int next_trace = minimum;
	short int foo;
	// Used to limit one energy contribution per loop iteration, could also use break statements
	bool ntraced = true;

	while (j >0) {
		k=j-1;

		// Try and match the energy in the array element with possible contributions

		//if the element is an unpaired base
		if (element[k]<(alphabetSize+1)) {
			// Consider 3' Dangle
			if (j>=2){
				if (element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace == A[j-2]+data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1]._value) {
						next_trace=A[j-2];
						ntraced=false;
						j=j-2;
						foo=data->dangle[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][1];
					}
				}
			}
			// Consider terminal mismatch
			if (j>=3) {
				if (element[k-2]<(alphabetSize+1) && element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace==A[j-3]+data->tstack[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]._value) {
						next_trace=A[j-3];
						ntraced=false;
						j=j-3;
						foo=data->tstack[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]];
					}
				}
			}

			// Consider mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-3]>(alphabetSize+1) && element[k-1]>(alphabetSize+1) && element[k-2]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[j-4]+data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]]._value+data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]]._value) {
						next_trace=A[j-4];
						ntraced=false;
						j=j-4;
						foo=data->coaxstack[decon2(element[k-3],alphabetSize)][decon1(element[k-3],alphabetSize)][element[k-2]][element[k]]+data->tstackcoax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][element[k]][element[k-2]];
					}
				}
			}
		}
		else { //if the element is in a base pair
			// Consider 5' Dangle
			if (j>=2) {
				if (element[k-1]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[j-2]+data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2]._value) {
						next_trace=A[j-2];
						ntraced=false;
						j=j-2;
						foo=data->dangle[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-1]][2];
					}
				}
			}
			// Consider flush coaxial stacking
			if (j>=2) {
				if (element[k-1]>(alphabetSize+1) && ntraced) {
					if (next_trace==A[j-2]+data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)]._value) {
						next_trace=A[j-2];
						ntraced=false;
						j=j-2;
						foo=data->coax[decon2(element[k-1],alphabetSize)][decon1(element[k-1],alphabetSize)][decon1(element[k],alphabetSize)][decon2(element[k],alphabetSize)];
					}
				}
			}
			// Consider mismatch-mediated coaxial stacking
			if (j>=4) {
				if (element[k-2]>(alphabetSize+1) && element[k-3]<(alphabetSize+1) && element[k-1]<(alphabetSize+1) && ntraced) {
					if (next_trace==A[j-4]+data->coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]]._value+data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]]._value) {
						next_trace=A[j-4];
						ntraced=false;
						j=j-4;
						foo=data->coaxstack[decon2(element[k],alphabetSize)][decon1(element[k],alphabetSize)][element[k-3]][element[k-1]]+data->tstackcoax[decon2(element[k-2],alphabetSize)][decon1(element[k-2],alphabetSize)][element[k-1]][element[k-3]];
					}
				}
			}
		}
		
		// By default move to the next element
		if (ntraced) {
			next_trace=A[j-1];
			j=j-1;
		}
		ntraced=true;
	}
#endif //COUNTING	

//DALIA:after traceback which is not included

	if (intermolecular) minimum = minimum + data->init;

	//deallocate memory use:
	delete[] A;
	delete[] element;

	//return the energy:
	//cout << "ergex returns" << minimum+au*data->auend << " ";
	return minimum+auendpenalties;

}

//Write the structure class-specific items in a save file
void writestructuresave(ofstream *sav, structure *ct1) {

	int i,j;
	int local = ct1->GetNumberofPairs();
    write(sav,&(local));
    for (i=0;i<ct1->GetNumberofPairs();++i) {
		local = ct1->GetPair5(i);
      write(sav,&(local));
		local = ct1->GetPair3(i);
      write(sav,&(local));
    }
    for (i=0;i<=ct1->GetSequenceLength();++i) {
      write(sav,&(ct1->hnumber[i]));
      sav->write(&(ct1->nucs[i]),1);
    }

    for (i=0;i<=2*ct1->GetSequenceLength();i++) {
      write(sav,&(ct1->numseq[i]));
    }

	local = ct1->GetNumberofDoubles();
    write(sav,&(local));
    for (i=0;i<ct1->GetNumberofDoubles();i++) {
		local = ct1->GetDouble(i);
      write(sav,&(local));
    }

    write(sav,&(ct1->intermolecular));
    if (ct1->intermolecular) {
      for (i=0;i<3;i++) {
        write(sav,&(ct1->inter[i]));
      }
    }

	local = ct1->GetNumberofSingles();
    write(sav,&local);
    for (i=0;i<ct1->GetNumberofSingles();i++) {
		local = ct1->GetSingle(i);
      write(sav,&(local));
    }

	local = ct1->GetNumberofModified();
    write(sav,&(local));
    for (i=0;i<ct1->GetNumberofModified();i++) {
		local = ct1->GetModified(i);
      write(sav,&(local));
    }

	local = ct1->GetNumberofGU();
    write(sav,&(local));
    for (i=0;i<ct1->GetNumberofGU();i++) {
		local = ct1->GetGUpair(i);
      write(sav,&(local));
    }

	string label=ct1->GetSequenceLabel();
    write(sav,&label);

    write(sav,&(ct1->templated));
    if (ct1->templated) {
      for (i=0;i<=ct1->GetSequenceLength();i++) {
        for (j=0;j<=i;j++) {
          write(sav,&(ct1->tem[i][j]));
        }
      }
    }



}

//Read the structure class-specific items in a save file
void openstructuresave(ifstream *sav, structure *ct1) {

	int i, j, local,local1,local2;

  read(sav,&(local));
  
  //read pairing constraints
  for (i=0;i<local;i++) {
    read(sav,&(local1));
    read(sav,&(local2));
	ct1->AddPair(local1,local2);
  }

  //read the historical number, and the nucleotides
  for (i=0;i<=ct1->GetSequenceLength();i++) {

    read(sav,&(ct1->hnumber[i]));
    sav->read(&(ct1->nucs[i]),1);

  }

  //read the numseq info
  for (i=0;i<=2*ct1->GetSequenceLength();i++) read(sav,&(ct1->numseq[i]));

  //Read the constrainst for double-stranded
  read(sav,&(local));
  for (i=0;i<local;++i) {
	  read(sav,&(local1));
	  ct1->AddDouble(local1);
  }

  //read the intermolecular array
  read(sav,&(ct1->intermolecular));
  if (ct1->intermolecular) {
    for (i=0;i<3;i++) read(sav,&(ct1->inter[i]));

  }

  //Read the unpaired constraints
  read(sav,&(local));
  for (i=0;i<local;++i) {
	  read(sav,&(local1));
	  ct1->AddSingle(local1);
  }

  //read the modified nucleotides
  read(sav,&(local));
  for (i=0;i<local;++i) {
	  read(sav,&(local1));
	  ct1->AddModified(local1);
  }

  //Read Us in GU pairs
  read(sav,&(local));
  for (i=0;i<local;++i) {
	  read(sav,&(local1));
	  ct1->AddGUPair(local1);
  }


  //Read the sequence label
  string localstring;
  read(sav,&localstring);
  ct1->SetSequenceLabel(localstring);


  //Read the template information
  read(sav,&(ct1->templated));
  if (ct1->templated) {
    ct1->allocatetem();
    for (i=0;i<=ct1->GetSequenceLength();i++) {
      for (j=0;j<=i;j++) {
        read(sav,&(ct1->tem[i][j]));
      }
    }
  }

}

void writehelixfile(char *filename,structure *ct,int StructureNumber) {
	//write a helix file that can be read by XRNA
	int i,count;
	ofstream out;

	out.open(filename);

	i=1;
	while (i<=ct->GetSequenceLength()) {
		if (ct->GetPair(i,StructureNumber)>i) {
			//found a base pair
			out << i << " " << ct->GetPair(i,StructureNumber) << " ";

			//determine length of helix
			count = 1;
			while (ct->GetPair(i+1,StructureNumber)==ct->GetPair(i,StructureNumber)-1) {
				i++;
				count++;

			}
			out << count << "\n";
			i++;


		}
		else i++;


	}



}

// Generates constraint matrix that can be passed to partition-cuda.
int *structure::generate_constraint_matrix(){
	// Declare base can pair matrix as a linearized triangular array.  This is the format partition cuda uses for its dynamic programming arrays
	// Also note that the arrays in partition-cuda a 0-indexed, not 1-indexed.
	int *base_can_pair = new int [numofbases*(numofbases-1)/2];
	int i, j, n;

	// Initialize matrix values based on sequence (enforce minimum hairpin loop size, allowed base pairs, no lonely pairs)
	for (i=0; i< numofbases-1; i++){
		for (j=i+1; j< numofbases; j++){
	    	if (j-i < minloop+1)
				// Potential hairpin loop is too small to be allowed 
        		base_can_pair[(j*(j-1))/2 + i]=0;
			else {
				if (data->can_pair(i+1, j+1, numseq)) {
					if (data->can_pair(i+2, j, numseq)) 
						// Possible base pair stack
						base_can_pair[(j*(j-1))/2 + i]=1;
					else {
						if  (i>0 && j<numofbases-1)
							if (data->can_pair(i, j+2, numseq))
								// Possible base pair stack
								base_can_pair[(j*(j-1))/2 + i]=1;
							else
								// Lonely pairs are not allowed
								base_can_pair[(j*(j-1))/2 + i]=0;
						else
							// Lonely pairs are not allowed
							base_can_pair[(j*(j-1))/2 + i]=0;
					}
				}
				else
					// The sequence does not form a canonical base pair
					base_can_pair[(j*(j-1))/2 + i]=0;
			}
		}
	}

//	vector<int> doublestranded; //nucleotides that must be double stranded
//	vector<int> GUpair; //Us in GU pairs
//	vector<int>	modified; //nucleotides accessible to tradictional chemical modification agents

	// Constrain specific forbidden base pairs
//	cout << "Forbidden Base Pairs:  " << forbid3.size() << endl;

	for(n=0; n<forbid5.size(); n++){
		i = forbid5[n];
		j = forbid3[n];
		cout << i << "\t" << j << endl;
		base_can_pair[(j*(j-1))/2 + i]=0;
	}

	//nucleotides that must be single stranded
//	cout << "Forced Single Stranded:  " << singlestranded.size() << endl;

	for(n=0; n<singlestranded.size(); n++){
		// Forbid all possible 5' pairing partners
		for (i=0; i<singlestranded[n]; i++){
			base_can_pair[((singlestranded[n]-1)*(singlestranded[n]-2))/2 + i]=0;
		}

		// Forbid all possible 3' pairing partners
		for (j=singlestranded[n]; j<numofbases; j++){
			base_can_pair[(j*(j-1))/2 + singlestranded[n]]=0;
		}
	}

	// Forced base pairs
//	cout << "Forced Pairs:  " << pair5.size() << endl;

	for(n=0; n<pair5.size(); n++){
		i = pair5[n];
		j = pair3[n];
		base_can_pair[(j*(j-1))/2 + i]=0;
		for(i=0; i<pair5[n]; i++){
			for (j=pair5[n]; j<pair3[n]; j++){
				base_can_pair[(j*(j-1))/2 + i]=0;
			}
		}
		for(i=pair5[n]; i<pair3[n]; i++){
			for (j=pair3[n]; j<numofbases; j++){
				base_can_pair[(j*(j-1))/2 + i]=0;
			}
		}
	}

	// Domain constraints
//	cout << "Forced Domains:  " << domains5.size() << endl;

	for(n=0; n<domains5.size(); n++){
		for(i=0; i<domains5[n]-1; i++){
			for (j=domains5[n]; j<domains3[n]-1; j++){
				base_can_pair[(j*(j-1))/2 + i]=0;
			}
		}
		for(i=domains5[n]; i<domains3[n]-1; i++){
			for (j=domains3[n]; j<numofbases-1; j++){
				base_can_pair[(j*(j-1))/2 + i]=0;
			}
		}
		i = domains5[n]-1;
		j = domains3[n]-1;
		base_can_pair[(j*(j-1))/2 + i]=1;
		
	}//*/
	
	return base_can_pair;
}

// Determine if the structure has one or more Pseudoknots (crossing bonds).
bool hasPseudoknots(const vector<int> &pairs) {
	const int length = pairs.size();
    IntervalStack stack(min(8,length/4));
	const int FIRST_NUC = 1; // in RNAstructure the convention is that index 1 is the first nucleotide (instead of 0)
    stack.push(FIRST_NUC, length-1);
    while (stack.pop()) {
		int k; // represents the 3' end of the basepair for which the 5' end is at stack.i
        // increase stack.i until it is either == stack.j or stack.i points to a basepair.
		while(stack.i <= stack.j && 0==(k = pairs[stack.i])) // find out if there is a bond at position i (and make sure that i is the 5' end).
            stack.i++;

		// if no basepair was found between i and j, pop that interval. we are done looking inside it.
		if (stack.i > stack.j) 
            continue; // pop next item from stack

		// whenever a bond is encountered, two new intervals are pushed -- the interval INSIDE the bond and the interval AFTER it.
		// The 5' end is always encountered first (even for crossing bonds) because the lower interval is processed before the upper interval.
		// So stack.i will ALWAYS be the 5' position and bases[stack.i] will ALWAYS be the 3' position.
		if (k < stack.i) cerr << "Programming logic error. 5' end encountered in ::hasPseudoknots" << endl;

        if (k > stack.j) return true; // if the 3' end is outside of the interval, this represents a crossing bond. So return true.
        if (k + 1 <= stack.j)  // if there are any nucleotides between the 3'-end and j, push that interval
            stack.push(k+1, stack.j);
        if (stack.i + 1 <= k-1) // if there are any nucleotides between the 5' and 3' end, push that interval
            stack.push(stack.i+1, k-1);
	}
    return false;
}

// Given a vector containing base-pairing information, this function finds the pseudoknots using a dynamic programming algorithm.
// See the description in the header for more information.
void findPseudoknots(const vector<int> &pairs, vector<int> *knotted, vector<int> *optimal) {
	const unsigned int length = pairs.size();
	const unsigned int FIRST_NUC = 1; // in RNAstructure the convention is that index 1 is the first nucleotide (instead of 0)
	if (length==0) return;
    
	if (optimal==NULL && knotted==NULL) return; // at least one result vector should be passed in or there is no point in doing the calculation.

	// Create an upper-triangular matrix for both outer and trace.
	// Since the first base is at position 1, the 0th column and row are never used, so we can create smaller 
	//   arrays and do pointer math to make outer[1][1] actually point to the 0th element (i.e. first memory location).
	// Additionally note that for trace, we can make the array even smaller because i>j for all trace[i][j] so 
	//   we can delete the center diagonal. So trace[1][2] points to the 0th element. (i.e. rows are shifted by 1, columns by 2)
	// In `outer`: rows go from i=1 to N; cols from i to N  (where N=length-1 is the index of the last base)
	short** outer=(new short*[length-FIRST_NUC])-FIRST_NUC;  for(int i=FIRST_NUC; i<length; ++i) outer[i]=(new short[length-i])-i;  
    // In `trace`: // rows go from i=1 to N-1; cols from i+1 to N
	// The size of trace is length-2, but it is but it is only shifted by 1. This is because it is one row shorter than outer.
	bool**  trace=(new bool*[length-FIRST_NUC-1])-FIRST_NUC; for(int i=FIRST_NUC; i<length-1; ++i) trace[i]=(new bool[length-i-1])-i-1; 

	// sanity test for triangular array and pointer math: (these will NOT necessarily result in access violations or segfaults, even if the pointer math is wrong)
	// DEBUG: for(int i=FIRST_NUC; i<length; i++){for(int j=i; j<length; j++){outer[i][j]=i*j; cout<<outer[i][j]<<" ";}cout<<endl;}
	// DEBUG: for(int i=FIRST_NUC; i<length-1; i++){for(int j=i+1; j<length; j++){outer[i][j]=i*j; cout<<outer[i][j]<<" ";}cout<<endl;}

	for (int i=FIRST_NUC; i<length; i++) outer[i][i]=0; // initialize diagonal of `outer` to 0 (this ensures proper start conditions)

	// The center diagonal of `outer` has been set to 0.
	// Each iteration of `n` fills the next diagonal to the (upper) right.
	// The first iteration represents all intervals [i, i+1]. The second is [i, i+2] etc. so the size of the interval
	// increases with each iteration. (Note: "interval" means a segment of the nucleobase sequence. The "size" is simply `j - i`)
	// Iteration 1:  [1,2], [2,3], [3,4]...[N-2,N-1] [N-1,N]  (size: 1)       (where N = length-1 is the last valid index in `bases`.)
	// Iteration 2:  [1,3], [2,4], [3,5]...[N-2, N]           (size: 2)
	// Ieration N-1: [1,N-1],[2,N]                            (size: N-1)
	// Iteration N:  [1,N]                                    (size: N)
	// The value at outer[i, j] represents the maximum number of non-crossing bonds **fully contained** within the sequence-interval [i, j].
	// (Fully contained means  bond.i >=i and bond.j <= j for every bond)
	// trace[i, j] is true if the bond at bases[i] should be included (i.e. it maximizes the number of bonds in the interval [i,j])
	
	for (int n = 1; n < length; n++) {
		for (int i = FIRST_NUC; i < length - n; i++) {
            int j = i + n;
            // First assume that either no bond starts here. So the total non-crossing bonds within the interval [i, j] is unchanged compared 
			// to the (smaller) interval [i+1, j] which would have been set in the previous iteration.
			outer[i][j] = outer[i + 1][j]; 
            trace[i][j] = false; // default value assume there is no bond here.
            int k = pairs[i];   // k == 0 if there is no basepair. 
			                    // if k > i, then i is the 5' end, and k is the 3' end of the pair. 
			                    // if k < i then k is the 5' end and i is the 3' end.
			if (k != 0 && k > i && k <= j) { // i.e. If this is the 3' end of a basepair and the 5' end (k) is inside the current interval [i, j]
                int tmp = 1; // add 1 for THIS bond
                
				if (i + 1 <= k - 1)
                    tmp += outer[i + 1][k - 1]; // add all the bonds **fully contained** by this one (i.e. which have i' > i and k' < k)

                if (k + 1 <= j)
                    tmp += outer[k + 1][j];  // add in all the bonds AFTER this one but still **fully contained** by the interval [k+1, j]

                if (tmp >= outer[i][j]) {     // if inclusion of this bond maximizes the total number of non-crossing bonds in the interval [i, j]
                    outer[i][j] = (short) tmp;
                    trace[i][j] = true;
                }
			}
        }
    }
	// DEBUG:
	//cout<<endl<<"pairs: ";  for(int i=FIRST_NUC; i<length; i++) cout<<pairs[i]<<"\t";  cout<<endl;
	//for(int i=FIRST_NUC; i<length; i++) {
	//	for(int j=FIRST_NUC; j<=i; j++) cout<<"\t";
	//	for(int j=i+1; j<length; j++)
	//		cout<<outer[i][j]<<(trace[i][j]?"**\t":"\t");
	//	cout<<endl;	}

	// We no longer need the "outer" data, so we reuse it here:
	// Note that the user may pass in the same vector for pairs AND either optimal or knotted, 
	// so we have to make a copy of the pairing data in case the pairs array is modified when setting values in `optimal` or `knotted`.
	short *results = outer[FIRST_NUC]; for(int i = FIRST_NUC; i < length; i++) results[i]=pairs[i]; // fill in results with the "pairs" information.
    // Backtrace
    IntervalStack stack(min(8,length / 4));
    stack.push(1, length - 1);
    while (stack.pop()) {
		//DEBUG: cout << "stack: " << stack.i << "," << stack.j << endl;
		while(stack.i < stack.j && !trace[stack.i][stack.j]) // trace[i][j] is true if the pair with its 5' end at i is in the optimal set.
            stack.i++;
        if (stack.i >= stack.j)
            continue; // pop next item from stack
		//DEBUG: cout << "found: " << stack.i << "," << stack.j << ": " << trace[stack.i][stack.j] << endl;
        
		// here k != -1 and i < j
		int k = pairs[stack.i];
        results[stack.i] = -k; // make the 5' end negative to indicate it's in the optimal set.
        results[k] = -pairs[k]; // make the 3' end negative
        if (stack.i + 1 < k-1)
            stack.push(stack.i+1, k-1);
        if (k + 1 < stack.j)
            stack.push(k+1, stack.j);
    }
	//DEBUG: cout<<endl<<"Results: ";for(int i=FIRST_NUC; i<length; i++) cout<<results[i]<<"   ";cout<<endl;

	if (optimal != NULL) {
		if (optimal->size() < length) optimal->resize(length);
		for (int i = FIRST_NUC; i < length; i++)
			(*optimal)[i]=results[i]<0 ? -results[i] : 0; // fill in pair information if this pair IS in the optimal results list.
	}
	if (knotted != NULL) {
		if (knotted->size() < length) knotted->resize(length);
		for (int i = FIRST_NUC; i < length; i++)
			(*knotted)[i]=results[i]>0 ? results[i] : 0; // fill in pair information if this pair is NOT in the optimal results list.
	}

	// Note: This should perform the inverse of the pointer math done at allocation.
	for(int i=FIRST_NUC;i<length;i++)   delete[] (outer[i]+i);   delete[] (outer+FIRST_NUC);	 
	for(int i=FIRST_NUC;i<length-1;i++) delete[] (trace[i]+i+1); delete[] (trace+FIRST_NUC);
	 // FYI - do not delete `results` because it is just a pointer to outer[0].
}

// Returns structure * from multiple_sequence_alignment[i] 
//structure* structure::sequence_from_alignment(int i)
//{
//	return &(*multiple_sequence_alignment)[i];
//}