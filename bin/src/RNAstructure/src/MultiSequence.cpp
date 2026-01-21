#include "MultiSequence.h"

#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include <cstdlib>


//! MultiSequence Class
/*!
    The MultiSequence Class provides utilities for reading and writing multiple sequence data. 
*/


//! Default Constructor.
MultiSequence::MultiSequence() : sequences (NULL) {}

//! Destructor.
MultiSequence::~MultiSequence()
{
    // if sequences allocated
    if (sequences){
        // free all sequences
        for (SafeVector<Sequence *>::iterator iter = sequences->begin(); iter != sequences->end(); ++iter){
            assert (*iter);
            delete *iter;
            *iter = NULL;
        }

        // free sequence vector
        delete sequences;
        sequences = NULL;
    }
}

//! Add another sequence to an existing sequence list.
void MultiSequence::AddSequence(Sequence *sequence){
    assert(sequence);
    assert(!sequence->Fail());
    // Add sequence.
    if(!sequences)
        sequences = new SafeVector<Sequence *>;
    sequences->push_back(sequence);
}

//! Write MFA to output file. 
//! Allowing user to specify the number of columns.
void MultiSequence::WriteMFA(ostream &outfile, int numColumns){
    if(!sequences)
        return;
    // Loop through all sequences and write.
    for(SafeVector<Sequence *>::iterator iter = sequences->begin(); iter != sequences->end(); ++iter){
        (*iter)->WriteMFA(outfile, numColumns);
    }
}

//! Write ALN to output file.
//! Allowing user to specify the number of columns.
void MultiSequence::WriteALN(ostream &outfile, int numColumns , bool isClustal ){
    if(!sequences)
        return;

    int longestComment = 0;
    SafeVector<SafeVector<char>::iterator> ptrs(GetNumSequences());
    SafeVector<int> lengths(GetNumSequences());
    for(int i = 0; i < GetNumSequences(); i++){
        ptrs[i] = GetSequence(i)->GetDataPtr();
        lengths[i] = GetSequence(i)->GetLength();
        longestComment = max(longestComment, (int) GetSequence(i)->GetName().length());
    }
    longestComment += 4;

    if( isClustal ){
        int writtenChars = 0;
        bool allDone = false;
        while(!allDone){
            allDone = true;

            // Clustal format
            for(int i = 0; i < GetNumSequences(); i++){
                if(writtenChars < lengths[i]){
                    outfile << GetSequence(i)->GetName();
                    for(int j = 0; j < longestComment - (int) GetSequence(i)->GetName().length(); j++)
                        outfile << ' ';

                    for(int j = 0; j < numColumns; j++){
                        if(writtenChars + j < lengths[i])
                            outfile << ptrs[i][writtenChars + j + 1];
                        else
                            break; 
                    }
                    outfile << endl;

                    if(writtenChars + numColumns < lengths[i])
                        allDone = false;
                }
            }
            outfile << endl;
            writtenChars += numColumns;
        }
    }

    if( !isClustal ){
        for(int i = 0; i < GetNumSequences(); i++){
            int writtenChars = 0;
            bool allDone = false;
            outfile << '>' << GetSequence(i)->GetName() << endl;
            while(!allDone){
                allDone = true;
                if(writtenChars < lengths[i]){
                    for(int j = 0; j < numColumns; j++){
                        if(writtenChars + j < lengths[i])
                            outfile << ptrs[i][writtenChars + j + 1];
                        else
                            break;
                    }
                    outfile << endl;

                    if(writtenChars + numColumns < lengths[i])
                        allDone = false;
                }
                //outfile << endl;
                writtenChars += numColumns;
            }
            outfile << endl;
        }
    }
}

//! Retrieve a sequence from MultiSequence object.
Sequence* MultiSequence::GetSequence(int i){
    assert(sequences);
    assert(0 <= i && i < (int) sequences->size());
    return(*sequences)[i];
}

//! Retrieve a sequence from MultiSequence object.(const version)
Sequence* MultiSequence::GetSequence(int i) const {
    assert(sequences);
    assert(0 <= i && i < (int) sequences->size());

    return(*sequences)[i];
}

//! Returns the number of sequences in the MultiSequence.
int MultiSequence::GetNumSequences() const {
    if(!sequences) 
        return 0;
    return(int) sequences->size();
}

//! Organizes the sequences according to their sequence labels in ascending order.
void MultiSequence::SortByLabel () {
    assert (sequences);

    // a quick and easy O(n^2) sort
    for (int i = 0; i < (int) sequences->size()-1; i++){
        for (int j = i+1; j < (int) sequences->size(); j++){
            if ((*sequences)[i]->GetSortLabel() > (*sequences)[j]->GetSortLabel())
            swap ((*sequences)[i], (*sequences)[j]);
        }
    }
}

//! Relabels sequences so as to preserve the current ordering.
void MultiSequence::SaveOrdering () {
    assert (sequences);

    for (int i = 0; i < (int) sequences->size(); i++)
        (*sequences)[i]->SetSortLabel (i);
}

//! Extracts all sequences from MultiSequence object whose index is given by a set. 
//! Projects the multiple sequences to subset and returns as a new MultiSequence object.
MultiSequence * MultiSequence::Project(const set<int> &indices){
    SafeVector<SafeVector<char>::iterator> oldPtrs(indices.size());
    SafeVector<SafeVector<char> *> newPtrs(indices.size());
    assert (indices.size() != 0);

    int i = 0;
    for(set<int>::const_iterator iter = indices.begin(); iter != indices.end(); ++iter){
        oldPtrs[i++] = GetSequence(*iter)->GetDataPtr();
    }

    // Computes new length.
    int oldLength = GetSequence(*indices.begin())->GetLength();
    int newLength = 0;
    for(i = 1; i <= oldLength; i++){
        // Counts the sequence length by non-gap columns.
        bool found = false;
        for(int j = 0; !found && j < (int) indices.size(); j++)
            found = (oldPtrs[j][i] != '-');
        if(found) 
            newLength++;
    }

    // Builds new alignments.
    for(i = 0; i < (int) indices.size(); i++){
        newPtrs[i] = new SafeVector<char>(); 
        assert(newPtrs[i]);
        newPtrs[i]->push_back('@');
    }

    for(i = 1; i <= oldLength; i++){
        // Checks there is not gap in sequences.
        bool found = false;
        for(int j = 0; !found && j < (int) indices.size(); j++)
            found = (oldPtrs[j][i] != '-');
        if(found){
            for(int j = 0; j < (int) indices.size(); j++)
                newPtrs[j]->push_back(oldPtrs[j][i]);
        }
    }

    MultiSequence *ret = new MultiSequence();
    i = 0;
    for(set<int>::const_iterator iter = indices.begin(); iter != indices.end(); ++iter){
        ret->AddSequence(new Sequence(newPtrs[i++], GetSequence(*iter)->GetHeader(), newLength, 
                                      GetSequence(*iter)->GetSortLabel(), GetSequence(*iter)->GetLabel()));
    }

    return ret;
}
