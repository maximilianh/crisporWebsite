// #ifndef _SEQUENCE_
// #define _SEQUENCE_

#include "Sequence.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>

//! Sequence Class
/*!
    Class for storing sequence information.
*/

//! Default constructor.
Sequence::Sequence() : isValid (false), header (""), data (NULL), length (0), sequenceLabel (0), inputLabel (0) {}

//! Constructor. Builds a sequence from existing data.
Sequence::Sequence(SafeVector<char> *data, string header, int length, int sequenceLabel, int inputLabel) :
    isValid (data != NULL), header(header), data(data), length (length), sequenceLabel (sequenceLabel), inputLabel (inputLabel) {
    assert (data); 
    assert ((*data)[0] == '@');
}

//! Destructor.
Sequence::~Sequence(){
    if (data){
        assert (isValid);
        delete data;
        data = NULL;
        isValid = false;
    }
}

//! Return the string comment associated with this sequence.
string Sequence::GetHeader() const {
    return header;
}

//! Return the first word of the string comment associated with this sequence.
string Sequence::GetName() const {
    char name[1024];
    sscanf (header.c_str(), "%s", name);
    return string(name);
}

//! Return the iterator to data associated with this sequence.
SafeVector<char>::iterator Sequence::GetDataPtr(){
    assert (isValid);
    assert (data);
    return data->begin();
}

//! Return the character at position i. 
char Sequence::GetPosition(int i) const {
    assert (isValid);
    assert (data);
    assert (i >= 1 && i <= length);
    return (*data)[i];
}

//! Sets the sequence label to i.
void Sequence::SetLabel (int i){
    assert (isValid);
    sequenceLabel = i;
    inputLabel = i;
}

//! Sets the sequence sorting label to i.
void Sequence::SetSortLabel (int i){
    assert (isValid);
    sequenceLabel = i;
}

//! Retrieves the input label.
int Sequence::GetLabel() const {
    assert (isValid);
    return inputLabel;
}

//! Retrieves the sorting label.
int Sequence::GetSortLabel() const {
    assert (isValid);
    return sequenceLabel;
}

//! Checks to see if the sequence successfully loaded.
bool Sequence::Fail() const {
    return !isValid;
}

//! Returns the length of the sequence.
int Sequence::GetLength () const {
    assert (isValid);
    assert (data);
    return length;
}

//! Writes the sequence to outfile in Fasta format. 
void Sequence::WriteMFA(ostream &outfile, int numColumns) const {
    assert (isValid);
    assert (data);
    assert (!outfile.fail());

    // print out heading
    outfile << ">" << header << endl;

    // print out character data
    int ct = 1;
    for (; ct <= length; ct++){
        outfile << (*data)[ct];
        if (ct % numColumns == 0) outfile << endl;
    }
    if ((ct-1) % numColumns != 0) outfile << endl;
}

//! Returns a new deep copy of the seqeuence.
Sequence * Sequence::Clone () const {
    Sequence *ret = new Sequence();
    assert (ret);

    ret->isValid = isValid;
    ret->header = header;
    ret->data = new SafeVector<char>; 
    assert (ret->data);
    *(ret->data) = *data;
    ret->length = length;
    ret->sequenceLabel = sequenceLabel;
    ret->inputLabel = inputLabel;

    return ret;
}

//! Given an SafeVector<char> containing the skeleton for an alignment and the identity of the current character.
//! Create a new sequence with all necesssary gaps added.
//! For example: Given alignment = "XXXBBYYYBBYYXX", the new sequence is "ATGCC---GT--CA".
//!                                                                      (XXXBBYYYBBYYXX)

Sequence * Sequence::AddGaps (SafeVector<char> *alignment, char id){
    Sequence *ret = new Sequence();
    assert (ret);

    ret->isValid = isValid;
    ret->header = header;
    ret->data = new SafeVector<char>; assert (ret->data);
    ret->length = (int) alignment->size();
    ret->sequenceLabel = sequenceLabel;
    ret->inputLabel = inputLabel;
    ret->data->push_back ('@');

    SafeVector<char>::iterator dataIter = data->begin() + 1;
    for (SafeVector<char>::iterator iter = alignment->begin(); iter != alignment->end(); ++iter){
        if (*iter == 'B' || *iter == id){
            ret->data->push_back (*dataIter);
            ++dataIter;
        }
        else
            ret->data->push_back ('-');
    }

    return ret;
}

//! Returns a SafeVector<int> containing the indices of every character in the sequence.
//! For instance, if the data is "ATGCC---GT--CA", the method returns {1,2,3,4,5,9,10,13,14}.
SafeVector<int> * Sequence::GetMapping () const {
    SafeVector<int> *ret = new SafeVector<int>(1, 0);
    for (int i = 1; i <= length; i++){
        if ((*data)[i] != '-') ret->push_back (i);
    }
    return ret;
}

// #endif