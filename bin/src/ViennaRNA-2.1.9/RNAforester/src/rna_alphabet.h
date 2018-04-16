#ifndef _RNA_ALPHABET_
#define _RNA_ALPHABET_

#include <string>

using namespace std;

/*
not every c++ compiler allows that enum types are used as
template arguments

typedef enum
  {
    ALPHA_UNDEFINED=-1,
    ALPHA_BASE_A=0,
    ALPHA_BASE_C,
    ALPHA_BASE_G,
    ALPHA_BASE_U,
    ALPHA_BASEPAIR,
    ALPHA_GAP,
    ALPHA_BASE,
    RNA_ALPHABET_SIZE,
  } RNA_Alphabet;
*/

typedef char RNA_Alphabet;

/** pair alphabet of RNA_Alphabet */
typedef struct {RNA_Alphabet a; RNA_Alphabet b;} RNA_AlphaPair;


//const int ALPHA_UNDEFINED=-1;
/*
const int ALPHA_BASE_A=0;
const int ALPHA_BASE_C=1;
const int ALPHA_BASE_G=2;
const int ALPHA_BASE_U=3;
const int ALPHA_BASEPAIR=4;
const int ALPHA_GAP=5;
const int ALPHA_BASE=6;*/


const int ALPHA_BASE_A='a';
const int ALPHA_BASE_C='c';
const int ALPHA_BASE_G='g';
const int ALPHA_BASE_U='u';
const int ALPHA_BASEPAIR='P';
const int ALPHA_GAP='-';
const int ALPHA_BASE='B';

const int RNA_ALPHABET_SIZE=7;

/*const int ALPHABET_SIZE = 256;*/

/*
class RNA_AlphabetTransformation
{
private:
	const static int m_alphabetsize;
	const static int m_rna_alphabetsize;

	RNA_Alphabet m_alpha2RNA_Alpha[ALPHABET_SIZE];
	char m_RNA_Alpha2alpha[RNA_ALPHABET_SIZE];

	void setMapEntry(RNA_Alphabet a, char c);
public:
	RNA_AlphabetTransformation();

	inline RNA_Alphabet RNA_AlphabetTransformation::alpha2RNA_Alpha(char c)
	{
		return m_alpha2RNA_Alpha[(unsigned int)c];
	};

	inline char RNA_AlphabetTransformation::RNA_Alpha2alpha(RNA_Alphabet a)
	{
		return m_RNA_Alpha2alpha[(unsigned int)a];
	};

};*/

const int ALPHA_PRO_BASE_A=0;
const int ALPHA_PRO_BASE_C=1;
const int ALPHA_PRO_BASE_G=2;
const int ALPHA_PRO_BASE_U=3;
const int ALPHA_PRO_BASEPAIR=4;
const int ALPHA_PRO_GAP=5;
const int ALPHA_PRO_BASE=6;

struct RNA_Alphabet_Profile
{
	RNA_Alphabet_Profile() {};

	double p[RNA_ALPHABET_SIZE];
	string columnStr;
};

int alpha2RNA_Alpha(char c);

#endif
