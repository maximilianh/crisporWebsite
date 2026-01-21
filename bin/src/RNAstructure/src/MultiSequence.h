#ifndef _MULTISEQUENCE_
#define _MULTISEQUENCE_

#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include <cstdlib>
#include "SafeVector.h"
#include "Sequence.h"

#ifdef TURBOHOMOLOGY
	#define MultiSequence MultiSequenceturbohomology
#endif

using namespace std;

//! MultiSequence Class
/*!
    The MultiSequence Class provides utilities for reading and writing multiple sequence data. 
*/

class MultiSequence{

    SafeVector<Sequence *> *sequences;
    vector<string> sequence_names;

public:
    //! Default Constructor.
    MultiSequence();

    //! Destructor.
    ~MultiSequence();

    //! Add another sequence to an existing sequence list.
    void AddSequence(Sequence *sequence);

    //! Write MFA to output file. 
    //! Allowing user to specify the number of columns.
    void WriteMFA(ostream &outfile, int numColumns = 60);

    //! Write ALN to output file.
    //! Allowing user to specify the number of columns.
    void WriteALN(ostream &outfile, int numColumns = 60, bool isClustal = true);

    //! Retrieve a sequence from MultiSequence object.
    Sequence* GetSequence(int i);

    //! Retrieve a sequence from MultiSequence object.(const version)
    Sequence* GetSequence(int i) const ;

    //! Returns the number of sequences in the MultiSequence.
    int GetNumSequences() const ;

    //! Organizes the sequences according to their sequence labels in ascending order.
    void SortByLabel () ;

    //! Relabels sequences so as to preserve the current ordering.
    void SaveOrdering ();

    //! Extracts all sequences from MultiSequence object whose index is given by a set. 
    //! Projects the multiple sequences to subset and returns as a new MultiSequence object.
    MultiSequence *Project(const set<int> &indices);
};

#endif