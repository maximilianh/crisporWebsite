#ifndef _SEQUENCE_
#define _SEQUENCE_

#include <string>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "SafeVector.h"

using namespace std;

//! Sequence Class
/*!
    Class for storing sequence information.
*/
class Sequence {

    bool isValid;               // a boolean indicating whether the sequence data is valid or not
    string header;              // string containing the comment line of the FASTA file
    SafeVector<char> *data;     // pointer to character data
    int length;                 // length of the sequence
    int sequenceLabel;          // integer sequence label, typically to indicate the ordering of sequences
                                //   in a Multi-FASTA file
    int inputLabel;             // position of sequence in original input

    //! Default constructor.
    Sequence();

public:

    //! Constructor. Builds a sequence from existing data.
    Sequence(SafeVector<char> *data, string header, int length, int sequenceLabel, int inputLabel);

    //! Destructor.
    ~Sequence();

    //! Return the string comment associated with this sequence.
    string GetHeader() const ;

    //! Return the first word of the string comment associated with this sequence.
    string GetName() const ;

    //! Return the iterator to data associated with this sequence.
    SafeVector<char>::iterator GetDataPtr();
    //! Return the character at position i. 
    char GetPosition(int i) const ;

    //! Sets the sequence label to i.
    void SetLabel (int i);

    //! Sets the sequence sorting label to i.
    void SetSortLabel (int i);

    //! Retrieves the input label.
    int GetLabel() const ;

    //! Retrieves the sorting label.
    int GetSortLabel() const ;

    //! Checks to see if the sequence successfully loaded.
    bool Fail() const ;

    //! Returns the length of the sequence.
    int GetLength () const ;

    //! Writes the sequence to outfile in Fasta format. 
    void WriteMFA(ostream &outfile, int numColumns) const ;

    //! Returns a new deep copy of the seqeuence.
    Sequence *Clone () const ;

    //! Given an SafeVector<char> containing the skeleton for an alignment and the identity of the current character.
    //! Create a new sequence with all necesssary gaps added.
    //! For example: Given alignment = "XXXBBYYYBBYYXX", the new sequence is "ATGCC---GT--CA".
    //!                                                                      (XXXBBYYYBBYYXX)

    Sequence *AddGaps (SafeVector<char> *alignment, char id);

    //! Returns a SafeVector<int> containing the indices of every character in the sequence.
    //! For instance, if the data is "ATGCC---GT--CA", the method returns {1,2,3,4,5,9,10,13,14}.
    SafeVector<int> *GetMapping () const ;
};

#endif