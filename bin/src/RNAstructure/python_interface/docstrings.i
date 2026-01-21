
// File: index.xml

// File: structbp.xml
%feature("docstring") bp "C++ includes: ProbScan.h ";


// File: class_dynalign__object.xml
%feature("docstring") Dynalign_object "

Dynalign_object Class.

The Dynalign_object class provides an entry point for the Dynalign
algorithm. The class is inherited from the TwoRNA class, which itself
contains two instances to the class RNA.

C++ includes: Dynalign_object.h ";

%feature("docstring")  Dynalign_object::Dynalign_object "Dynalign_object::Dynalign_object()

Constructor This is a default constructor that calls the TwoRNA
default constructor. ";

%feature("docstring")  Dynalign_object::Dynalign_object "Dynalign_object::Dynalign_object(const char sequence1[], const char
sequence2[], const bool IsRNA=true)

This constuctor is available to reach the TwoRNA constructors, from
which this class is inherited. This constructor uses two cstrings to
provide sequences. IsRNA is true for RNA folding and flase for DNA.
Constructor This constuctor is available to reach the TwoRNA
constructors, from which this class is inherited. This constructor
uses two cstrings to provide sequences. IsRNA is true for RNA folding
and false for DNA. Input sequences should contain
A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. T=t=u=U.
If IsRNA is true, the backbone is RNA, so U is assumed. If IsRNA is
false, the backbone is DNA, so T is assumed. x=X= nucleotide that
neither stacks nor pairs. For now, any unknown nuc is considered 'X'.
Both sequences are passed to underlying RNA classes for each sequence.

Parameters:
-----------

sequence1:  is a NULL terminated c string for sequence 1.

sequence2:  is a NULL terminated c string for sequence 2.

IsRNA:  is a bool that indicates whether these sequences are RNA or
DNA. true= RNA. false=DNA. Default is true. Both sequences must have
the same backbone. ";

%feature("docstring")  Dynalign_object::Dynalign_object "Dynalign_object::Dynalign_object(const char filename1[], const int
type1, const char filename2[], const int type2, const bool IsRNA=true)

Constructor.

This constuctor is available to reach the TwoRNA constructors, from
which this class is inherited. This constructor uses two ctsirngs as
filenames, accompanied by integers to set the file type. IsRNA is true
for RNA folding and false for DNA. The existing files, specified by
filenames, can either be a ct file, a sequence, or an RNAstructure
save file. Therefore, the user provides a flag for the file type: type
= 1 => .ct file, type = 2 => .seq file, type = 3 => partition function
save (.pfs) file, type = 4 => folding save file (.sav). The file
opening is performed by the constructors for the RNA classes that
underlie each sequence. This constructor generates internal error
codes that can be accessed by GetErrorCode() after the constructor is
called. 0 = no error. The errorcode can be resolved to a c string
using GetErrorMessage. Note that the contructor needs to be explicitly
told, via IsRNA, what the backbone is because files do not store this
information. Note also that save files explicitly store the
thermodynamic parameters, therefore changing the backbone type as
compaared to the original calculation will not change structure
predictions.

Parameters:
-----------

filename1:  is a null terminated c string and refers to sequence 1.

filename2:  is a null terminated c string and refers to sequence 2.

type1:  is an integer that indicates the file type for sequence 1.

type2:  is an integer that indicates the file type for sequence 2.

IsRNA:  is a bool that indicates whether these sequences are RNA or
DNA. true= RNA. false=DNA. Default is true. Only one backbone is
allowed for both sequences. ";

%feature("docstring")  Dynalign_object::Dynalign_object "Dynalign_object::Dynalign_object(const char filename[])

Constructor.

This constructor allows the user to read a dynalign save file (.dsv)
to get base pairing information. This constructor generates internal
error codes that can be accessed by GetErrorCode() after the
constructor is called. 0 = no error. The errorcode can be resolved to
a c string using GetErrorMessage.

Parameters:
-----------

filename[]:  is a cstring that indicates the filename of a .dsv file.
";

%feature("docstring")  Dynalign_object::Dynalign_object "Dynalign_object::Dynalign_object(const char *filename, const short
maxtrace, const short bpwin, const short awin, const short percent)

Constructor This constructor is used to perform Dynaligh refolding.
This does not allow any changes in constraints, but does allow the
creation of different set of suboptimal structures. This constructor
generates internal error codes that can be accessed by GetErrorCode()
after the constructor is called. 0 = no error. The errorcode can be
resolved to a c string using GetErrorMessage.

Parameters:
-----------

filename:  is the name of a Dynalign save file name (.dsv).

maxtrace:  is the maximum number of common structures to be
determined. The recommended default is 20.

bpwin:  the the base pair window parameter, where 0 allows the
structures to have similar pairs and larger windows make the
structures more diverse. The recommended default is 5.

awin:  is the alignment window parameter, where 0 allows the
alignments to be similar and larger values make the alignments more
diverse. The recommended default is 1.

percent:  is the maximum percent difference in total folding free
energy change above the lowest for suboptimal common structures. The
recommended default is 20. ";

%feature("docstring")  Dynalign_object::Dynalign "int
Dynalign_object::Dynalign(const short int maxtrace=20, const short int
bpwin=5, const short int awin=1, const short int percent=20, const
short int imaxseparation=-99, const float gap=0.4, const bool
singleinsert=true, const char savefile[]=NULL, const bool
optimalonly=false, const short int singlefold_subopt_percent=30, const
bool local=false, const short int numProcessors=1, const int
maxpairs=-1)

Predict the lowest free energy structure common to two sequences and
suboptimal solutions with the Dynalign algorithm.

In case of error, the function returns a non-zero that can be parsed
by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

maxtrace:  is the maximum number of common structures to be
determined. The defaults is 20.

bpwin:  the the base pair window parameter, where 0 allows the
structures to have similar pairs and larger windows make the
structures more diverse. The default is 5.

awin:  is the alignment window parameter, where 0 allows the
alignments to be similar and larger values make the alignments more
diverse. The default is 1.

percent:  is the maximum percent difference in total folding free
energy change above the lowest for suboptimal common structures. The
defaults is 20.

imaxseparation:  is the maximum separation between aligned
nucleotides. Values >= 0 are the traditional parameter, those below
zero trigger the HMM alignment method, which is now prefered.

gap:  is the cost of adding gap nucleotides in the alignment in
kcal/mol.

singleinsert:  is whether single basepair inserts are allowed in one
sequence vs the other.

savefile:  is c-string with the name of a dynalign savefile (*.dsv) to
be created.

optimalonly:  can be used to turn on a calculation of only the energy
(when true) and not the structures.

singlefold_subopt_percent:  is the maximum % difference of folding
energy above the lowest free energy structure for pairs in single
sequence folding that will be allowed in the dynalign calculation.

local:  is whether Dynalign is being run in local (true) or global
mode (false).

numProcessors:  is the number of processors to use for the
calculation. This requires a compilation for SMP.

maxpairs:  is under development for multiple sequence folding. Use -1
(default) for now.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error
occurred). ";

%feature("docstring")  Dynalign_object::WriteAlignment "void
Dynalign_object::WriteAlignment(const char filename[])

Write the alignment to disk.

This function should be called after loading a dynalign save file or
after a dynalign calculation has been performed. This function
generates no error flag. Nothing can go wrong...

Parameters:
-----------

filename:  is the file to which the alignment should be written. ";

%feature("docstring")  Dynalign_object::ForceAlignment "int
Dynalign_object::ForceAlignment(const int i, const int k)

Force an alignment during a Dynalign calculation).

Nucleotide i from sequence 1 will be aligned to nucleotide k in
sequence 2 in subsequent Dynalign calculation. The function returns 0
with no error and a non-zero otherwise that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

i:  is the index of nucleotide from sequence 1.

k:  is the index of nucleotide from sequence 2.

Returns:
--------

An integer that indicates an error code (0 = no error, 100 =
nucleotide i out of range, 101 = nucleotide k out of range). ";

%feature("docstring")  Dynalign_object::GetForcedAlignment "int
Dynalign_object::GetForcedAlignment(const int i, const int seq)

Get an alignment constraint.

Parameters:
-----------

i:  is the nucleotide number.

seq:  is the sequence (1 or 2) from which i is derived.

Returns:
--------

An integer that indicates the nucleotide to which i is forced to be
aligned, where 0 indicates no alignment. ";

%feature("docstring")  Dynalign_object::ReadAlignmentConstraints "int
Dynalign_object::ReadAlignmentConstraints(const char filename[])

Read alignment constraints from disk.

The file format is: i1 k1 i2 k2 -1 -1 Where each line gives a aligned
pair (i from sequence 1 and k from sequence 2). The file terminates
with -1 -1 to indicate the file end. The function returns 0 with no
error and a non-zero otherwise that can be parsed by GetErrorMessage()
or GetErrorMessageString().

Parameters:
-----------

filename:  is a c string that is the file name to be read.

Returns:
--------

An integer that indicates an error code (0 = no error, 102 = file not
found, 103 = error reading constraint file). ";

%feature("docstring")  Dynalign_object::Templatefromct "int
Dynalign_object::Templatefromct(const char ctfilename[])

Read a ct file to determine what pairs will be allowed for sequence 1
in a subsequent dynalign calculation.

This results in all pairs but those in the ct being disallowed.

Parameters:
-----------

ctfilename:  is the name of the ct file to be read to provide the
template.

Returns:
--------

An integer that indicates an error code (0=no error, 104=file not
found, 105=template is already specified) ";

%feature("docstring")  Dynalign_object::Templatefromdsv "int
Dynalign_object::Templatefromdsv(const char dsvfilename[], const float
maxdsvchange)

This reads a dsv file and only allows pairs with folding free energy
change between the lowest and lowest + maxdsvchange in a subsequent
dynalign calculation.

Parameters:
-----------

dsvfilename:  is the name of the ct file to be read to provide the
template.

maxdsvchange:  in a float that gives a percent difference in free
energy above the lowest free energy change.

Returns:
--------

An integer that indicates an error code (0=no error, 106=file not
found, 105=template is already specified) ";

%feature("docstring")  Dynalign_object::GetBestPairEnergy "double
Dynalign_object::GetBestPairEnergy(const int sequence, const int i,
const int j)

Report the best energy for pair i-j from sequence number sequence (1
or 2).

This function reports the lowest ffolding free energy for any pairs
between i-j in sequence number sequence (1 or 2). This requires a
search over all possible pairs in the second sequence. NOTE: This
function ONLY works after reading a Dynalign save file (.dsv) using
the constructor. This is because the Dynalign energies are not
normally stored after calling Dynalign. This function generates
internal error codes that can be accessed by GetErrorCode() after the
constructor is called. 0 = no error, 107 = Data not available, 108 =
nucleotide out of range. The errorcode can be resolved to a c string
using GetErrorMessage.

Parameters:
-----------

sequence:  is an integer indicating the sequence # (must be 1 or 2).

i:  is the 5' nucleotide in a pair.

j:  is the 3' nucleotide in a pair.

Returns:
--------

A double that gives an energy in kcal/mol. ";

%feature("docstring")  Dynalign_object::GetLowestEnergy "double
Dynalign_object::GetLowestEnergy()

Report the lowest total free energy change from a Dynalign
calculation.

NOTE: This function ONLY works after reading a Dynalign save file
(.dsv) using the constructor. This is because the Dynalign energies
are not normally stored after calling Dynalign. This function
generates internal error codes that can be accessed by GetErrorCode()
after the constructor is called. 0 = no error, 107 = Data not
available, 108 = nucleotide out of range. The errorcode can be
resolved to a c string using GetErrorMessage.

Returns:
--------

a double that gives an energy in kcal/mol. ";

%feature("docstring")  Dynalign_object::GetErrorMessage "char *
Dynalign_object::GetErrorMessage(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

0 = no error 100-999 = Error associated with Dynalign, to be handled
here. >=1000 = Errors for underlying sequence, get message from TwoRNA
base class. Current errors handled here are: 100 \"Nucleotide from
sequence 1 is out of range.\\\\n\"; 101 \"Nucleotide from sequence 2
is out of range.\\\\n\"; 102 \"Alignment constraint file not
found.\\\\n\"; 103 \"Error reading alignment constraint file.\\\\n\";
104 \"CT file not found.\\\\n\"; 105 \"A template has already been
specified; only one is allowed.\\\\n\"; 106 \"DSV file not
found.\\\\n\"; 107 \"Data not available to calculate energy.\\\\n\"
108 \"Nucleotide out of range.\\\\n\"; 109 \"Value of maxpairs is too
large to be achievable.\\\\n\" 110 \"Error reading thermodynamic
parameters.\"

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode() or by a
call to a function that returns an error.

Returns:
--------

A pointer to a c string that provides an error message or from other
functions that return integer error codes. ";

%feature("docstring")  Dynalign_object::SetProgress "void
Dynalign_object::SetProgress(TProgressDialog &Progress)

Provide a TProgressDialog for following calculation progress. A
TProgressDialog class has a public function void update(int percent)
that indicates the progress of a long calculation.

Parameters:
-----------

Progress:  is a TProgressDialog class. ";

%feature("docstring")  Dynalign_object::StopProgress "void
Dynalign_object::StopProgress()

Provide a means to stop using a TProgressDialog. StopProgress tells
the RNA class to no longer follow progress. This should be called if
the TProgressDialog is deleted, so that this class does not make
reference to it. ";

%feature("docstring")  Dynalign_object::~Dynalign_object "Dynalign_object::~Dynalign_object() ";


// File: structhp.xml
%feature("docstring") hp "C++ includes: ProbScan.h ";


// File: class_hybrid_r_n_a.xml
%feature("docstring") HybridRNA "

HybridRNA Class.

The HybridRNA class provides an entry point for all the bimolecular
structure prediction routines of RNAstructure. The class is inherited
from the RNA class and contains an instance of TwoRNA, which itself
contains two instances to the class RNA.

C++ includes: HybridRNA.h ";

%feature("docstring")  HybridRNA::HybridRNA "HybridRNA::HybridRNA(const char sequence1[], const char sequence2[],
const bool IsRNA=true)

Constructors.

The two constuctors are available to reach the TwoRNA constructors,
from which this class is built. Please see the TwoRNA constructor
entries for documentation. ";

%feature("docstring")  HybridRNA::HybridRNA "HybridRNA::HybridRNA(const char filename1[], const int type1, const
char filename2[], const int type2, const bool IsRNA=true) ";

%feature("docstring")  HybridRNA::AccessFold "int
HybridRNA::AccessFold(const double gamma=0.4, const float percent=50,
const int maximumstructures=20, const int window=0, const int
maxinternalloopsize=30)

Predict the lowest free energy secondary structure for two interacting
strands and generate suboptimal structures. Thuis method does not
allow intramolecular pairs. It considers accessibility with a
heuristic that uses tge partition function. If the temperature has not
been specified using the RNA base class SetTemperature and no free
energies have been calculated, the thermodynamic parameters have not
been read and therefore they will be read by this function call. The
parameter files should be located in the directory specified by the
environment variable $DATAPATH of the pwd. In case of error, the
function returns a non- zero that can be parsed by GetErrorMessage()
or GetErrorMessageString().

Parameters:
-----------

gamma:  is a scaling factor that weights accessibility. The defaults
is 0.4

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 50.

maximumstructures:  is the maximum number of suboptimal structures to
generate. The fefault is 20.

window:  is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). The default is 0.

maxinternalloopsize:  is the maximum number of unpaired nucleotides in
bulge and internal loops. This is used to accelerate the prediction
speed. The default is 30.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files, 14 = traceback error). ";

%feature("docstring")  HybridRNA::FoldBimolecular "int
HybridRNA::FoldBimolecular(const float percent=10, const int
maximumstructures=20, const int window=0, const char savefile[]=\"\",
const int maxinternalloopsize=30)

Predict the lowest free energy secondary structure and generate
suboptimal structures using a heuristic.

This function predicts the lowest free energy structure and suboptimal
structures. If the temperature has not been specified using the RNA
base class SetTemperature and no free energies have been calculated,
the thermodynamic parameters have not been read and therefore they
will be read by this function call. The parameter files should be
located in the directory specified by the environment variable
$DATAPATH of the pwd. In case of error, the function returns a non-
zero that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 10.

maximumstructures:  is the maximum number of suboptimal structures to
generate. The defaults is 20.

window:  is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). The default is 0.

savefile:  is c string containing a file path and name for a savefile
(.sav)that can be used to generate energy dot plots and to refold the
secondary structure using different suboptimal structure parameters.
The default is \"\", which results in no save file written.

maxinternalloopsize:  is the maximum number of unpaired nucleotides in
bulge and internal loops. This is used to accelerate the prediction
speed. The default is 30.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files, 14 = traceback error). ";

%feature("docstring")  HybridRNA::FoldDuplex "int
HybridRNA::FoldDuplex(const float percent=40, const int
maximumstructures=10, const int window=0, const int
maxinternalloopsize=30)

Predict the lowest free energy secondary structure for two strands
that cannot form intramolecular pairs and generate suboptimal
structures using a heuristic.

This function predicts the lowest free energy bimolecular structure
and suboptimal structures. This function does not allow any folding
constraints. If the temperature has not been specified using
SetTemperature and no free energies have been calculated, the
thermodynamic parameters have not been read and therefore they will be
read by this function call. The parameter files should be located in
the directory specified by the environment variable $DATAPATH of the
pwd. In case of error, the function returns a non-zero that can be
parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 40.

maximumstructures:  is the maximum number of suboptimal structures to
generate. The default is 10.

window:  is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). The default is 0.

maxinternalloopsize:  is the maximum number of unpaired nucleotides in
bulge and internal loops. This is used to accelerate the prediction
speed. The default is 30.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files, 14 = traceback error). ";

%feature("docstring")  HybridRNA::PartitionFunctionBimolecular "int
HybridRNA::PartitionFunctionBimolecular(const char savefile[]=\"\")

Predict the bimolecular partition function for a sequence (with no
intramolecular pairs).

This function must be called to predict base pair probabilities,
perform stochastic traceback, or for maximizing expected accuracy.
This predicts the partition function without intramolecular pairs. If
the temperature has not been specified using SetTemperature and no
free energies have been calculated, the thermodynamic parameters have
not been read and therefore they will be read by this function call.
The parameter files should be located in the directory specified by
the environment variable $DATAPATH of the pwd. In case of error, the
function returns a non-zero that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

savefile:  is a c string that contains the path and filename for
creating a save file. This defaults to \"\", which indicates no file
is to be written.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files). ";

%feature("docstring")  HybridRNA::GetRNA1 "RNA * HybridRNA::GetRNA1()

Access the underlying RNA class from an instance of TwoRNA. This is
provided for use with two sequence methods. Generally, there is no
need for end users to use this function.

Returns:
--------

A pointer to the underlying RNA class for sequence 1. ";

%feature("docstring")  HybridRNA::GetRNA2 "RNA * HybridRNA::GetRNA2()

Access the underlying RNA class from an instance of TwoRNA. This is
provided for use with two sequence methods. Generally, there is no
need for end users to use this function.

Returns:
--------

A pointer to the underlying RNA class for sequence 2. ";

%feature("docstring")  HybridRNA::GetErrorCode "int
HybridRNA::GetErrorCode()

Return an error code, where a return of zero is no error.

This function returns and error flag that is generated during
construction by RNA(const char &filename, const int type, const bool
IsRNA=true) or from CalculateFreeEnergy(). An error of zero is always
no error. Other codes are errors and a c-string can be fetched for the
error with GetErrorMessage().

Returns:
--------

An integer that provides the error code. ";

%feature("docstring")  HybridRNA::GetErrorMessage "char *
HybridRNA::GetErrorMessage(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

0 = no error <100 = Error to be fetched from RNA base class. 100-999 =
Error associated with bimolecular folding, to be handled here. >=1000
= Errors for underlying sequence, get message from TwoRNA base class.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode().

Returns:
--------

A pointer to a c string that provides an error message or from other
functions that return integer error codes. ";

%feature("docstring")  HybridRNA::GetForbidIntramolecular "bool
HybridRNA::GetForbidIntramolecular()

Get whether intramolecular pairs are allowed.

Returns:
--------

A bool that indicates whether intramolecular pairs are forbidden (true
= forbidden, false = not). ";

%feature("docstring")  HybridRNA::SetForbidIntramolecular "void
HybridRNA::SetForbidIntramolecular(const bool forbid)

Set whether intramolecular pairs are allowed.

If true is passed to this function, intramolecular pairs will be
forbidden in FoldBimolecular.

Parameters:
-----------

forbid:  is a bool that indicates whether intramolecular pairs are
forbid. ";

%feature("docstring")  HybridRNA::SetProgress "void
HybridRNA::SetProgress(TProgressDialog &Progress)

Provide a TProgressDialog for following calculation progress. A
TProgressDialog class has a public function void update(int percent)
that indicates the progress of a long calculation.

Parameters:
-----------

Progress:  is a TProgressDialog class. ";

%feature("docstring")  HybridRNA::StopProgress "void
HybridRNA::StopProgress()

Provide a means to stop using a TProgressDialog. StopProgress tells
the RNA class to no longer follow progress. This should be called if
the TProgressDialog is deleted, so that this class does not make
reference to it. ";

%feature("docstring")  HybridRNA::~HybridRNA "HybridRNA::~HybridRNA()
";


// File: structil.xml
%feature("docstring") il "C++ includes: ProbScan.h ";


// File: structmb.xml
%feature("docstring") mb "C++ includes: ProbScan.h ";


// File: classmb__element.xml
%feature("docstring") mb_element "C++ includes: ProbScan.h ";

%feature("docstring")  mb_element::mb_element "mb_element::mb_element(std::pair< int, int > h) ";

%feature("docstring")  mb_element::mb_element "mb_element::mb_element(int nuc) ";


// File: class_multifind__object.xml
%feature("docstring") Multifind_object "

Multifind_object Class.

The Multifind_object class provides an entry point for the Multifind
algorithm.

C++ includes: Multifind_object.h ";

%feature("docstring")  Multifind_object::Multifind_object "Multifind_object::Multifind_object(const string &outputmultifind,
const vector< string > &ctfiles, const vector< string >
&inputalignment, const vector< string > &inputsequences, const int
&processors, TProgressDialog *progress=NULL)

Constructor:

Parameters:
-----------

outputmultifind:  is the name of the Multifind output file to which
the output is written to.

ctfiles:  is a vector of strings storing the names of the ct files to
which the output structures are written to.

inputalignment:  is a vector of strings storing the input sequences in
the alignment (with gaps).

inputsequences:  is a vector of strings storing the input sequences in
the alignment (without gaps).

processors:  is a interger indicating the number of processors
required by Multifind in smp calculations.(only applicable in smp
version)

progress:  is a TProgressDialog for reporting progress of the
calculation to the user. The default value of NULL means that no
communication is provided. ";

%feature("docstring")  Multifind_object::Multifind_Predict "int
Multifind_object::Multifind_Predict()

The core function doing Multilign calculation and SVM prediction. ";


// File: class_multilign__object.xml
%feature("docstring") Multilign_object "

Multilign_object Class.

The Multilign_object class provides an entry point for the Multilign
algorithm.

C++ includes: Multilign_object.h ";

%feature("docstring")  Multilign_object::Multilign_object "Multilign_object::Multilign_object()

Default constructor: ";

%feature("docstring")  Multilign_object::Multilign_object "Multilign_object::Multilign_object(const vector< vector< string > >
&inputlist, const bool isrna=true, TProgressDialog *progress=NULL)

Constructor:

Parameters:
-----------

inputlist:  is a vector of vectors of strings storing the name of the
filenames. Currently, inputList is an matrix of 4 columns: col 1 is
the input seq filename; col 2 is the output ct filename; col 3 is the
input constraint filename; col 4 is the input SHAPE filename. Empty
string is not allowed for Col 1 and Col 2; If no SHAPE or folding
constraints are given, Col 3 and Col 4 for the correponding sequence
are empty strings.

isrna:  is a bool indicating the sequences are RNA or DNA. The default
of true indicates RNA.

progress:  is a TProgressDialog for reporting progress of the
calculation to the user. The default value of NULL means that no
communication is provided. ";

%feature("docstring")  Multilign_object::Multilign_object "Multilign_object::Multilign_object(const bool Multifind, const string
&outputmultifind, const vector< string > &ctfiles, TProgressDialog
*progress=NULL, const bool isrna=true) ";

%feature("docstring")  Multilign_object::~Multilign_object "Multilign_object::~Multilign_object() ";

%feature("docstring")  Multilign_object::CountBP "int
Multilign_object::CountBP(const int i=0, const int j=0, const double
percent=0.8) const

count the number of basepairs with the lowest free energies below the
percent of the miminal free energy. Dsv file is used for the counting.
By default, the first dsv file in the progressive dynalign
calculations is used, i.e. i = 0 and j = 0.

Parameters:
-----------

i:  is an int value indicating which one in the iteration.

j:  is an int value indicating which iteration.

percent:  is threshold of double value in percentage.

Returns:
--------

an int of the number of basepairs counted. ";

%feature("docstring")  Multilign_object::ProgressiveMultilign "int
Multilign_object::ProgressiveMultilign(const short int
numProcessors=1, const bool Dsv=1, const bool Ali=1, const short int
maxtrace=750, const short int bpwin=2, const short int awin=1, const
short int percent=20, const short int imaxseparation=-99, const float
gap=0.4, const bool singleinsert=true, const short int
singlefold_subopt_percent=30, const bool local=false)

The core function doing dynalign calculation and templating In case of
error, the function returns a non-zero that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

Dsv:  is a boolean value indicating to output pairwise dsv files or
not. It has to be set to true currently.

Ali:  is a boolean value indicating to output pairwise ali files or
not. It has to be set to true currently.

maxtrace:  is the maximum number of common structures to be
determined.

bpwin:  the the base pair window parameter, where 0 allows the
structures to have similar pairs and larger windows make the
structures more diverse.

awin:  is the alignment window parameter, where 0 allows the
alignments to be similar and larger values make the alignments more
diverse.

percent:  is the maximum percent difference in total folding free
energy change above the lowest for suboptimal common structures.

imaxseparation:  is the maximum separation between aligned
nucleotides. Values >= 0 are the traditional parameter, those below
zero trigger the HMM alignment method, which is now prefered.

gap:  is the cost of adding gap nucleotides in the alignment in
kcal/mol.

singleinsert:  is whether single basepair inserts are allowed in one
sequence vs the other.

singlefold_subopt_percent:  is the maximum % difference of folding
energy above the lowest free energy structure for pairs in single
sequence folding that will be allowed in the dynalign calculation.

local:  is whether Dynalign is being run in local (true) or global
mode (false).

numProcessors:  is the number of processors to use for the
calculation. This requires a compilation for SMP.

Returns:
--------

an int that indicates an error code (0 = no error, non-zero = error
occurred). ";

%feature("docstring")  Multilign_object::MultiTempMultilign "int
Multilign_object::MultiTempMultilign() ";

%feature("docstring")  Multilign_object::WriteAlignment "int
Multilign_object::WriteAlignment(const string allali=\"all.ali\")
const

calculate and output multiple alignment

Parameters:
-----------

allali:  is the output filename of multiple alignment

Returns:
--------

an int of error code. ";

%feature("docstring")  Multilign_object::GetErrorCode "int
Multilign_object::GetErrorCode() const

Return an error code, where a return of zero is no error.

This function returns and error flag that is generated during
construction by RNA(const char &filename, const int type, const bool
IsRNA=true) or from CalculateFreeEnergy(). An error of zero is always
no error. Other codes are errors and a c-string can be fetched for the
error with GetErrorMessage().

Returns:
--------

An integer that provides the error code. ";

%feature("docstring")  Multilign_object::GetErrorMessage "string
Multilign_object::GetErrorMessage(const int error) const

Return error messages based on code from GetErrorCode and other error
codes. 0 = no error 1000 = Error associated with sequence 1 or with a
procedure, function will get message from sequence 1 (the inherited
RNA class). 2000 = Error associated with sequence 2, function will get
message from sequence 2 (the RNA2 class). 3000 = Errors with each
sequence, function will get messages from each.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode().

Returns:
--------

A string that provides an error message or from other functions that
return integer error codes. ";

%feature("docstring")  Multilign_object::SetMaxPairs "int
Multilign_object::SetMaxPairs(const int maxpairs=-1)

Parameters:
-----------

maxpairs:  is int value defining how the MaxPairs will be set. By
default it is set to be -1, meaning the average length of all the
sequences.

Returns:
--------

an errorcode. ";

%feature("docstring")  Multilign_object::GetMaxPairs "int
Multilign_object::GetMaxPairs() const

get the value of MaxPairs

Returns:
--------

the value of MaxPairs ";

%feature("docstring")  Multilign_object::AverageLength "int
Multilign_object::AverageLength() const

get the average length of the input sequences

Returns:
--------

the average length of the input sequences. ";

%feature("docstring")  Multilign_object::SetIterations "int
Multilign_object::SetIterations(const int it=2)

set the value of iterations

Parameters:
-----------

it:  is an value of int assigned to iterations. By default it is set
to 2.

Returns:
--------

an errorcode ";

%feature("docstring")  Multilign_object::GetIterations "int
Multilign_object::GetIterations() const

get the value of iterations

Returns:
--------

the value of iterations. ";

%feature("docstring")  Multilign_object::SetMaxDsv "int
Multilign_object::SetMaxDsv(const float maxdsvchange=1)

set the value of MaxDsv/maxdsvchange

Parameters:
-----------

maxdsvchange:  is a value of float assigned to MaxDsv. By default it
is set to 1.

Returns:
--------

an errorcode ";

%feature("docstring")  Multilign_object::GetMaxDsv "float
Multilign_object::GetMaxDsv() const

get the value of MaxDsv/maxdsvchange

Returns:
--------

the value of MaxDsv/maxdsvchange. ";

%feature("docstring")  Multilign_object::GetSequenceNumber "int
Multilign_object::GetSequenceNumber() const

get the sequence number

Returns:
--------

the number of input sequences ";

%feature("docstring")  Multilign_object::SetIndexSeq "int
Multilign_object::SetIndexSeq(size_t indexSeq=1)

set the Index Sequence for Multilign calculation.

Parameters:
-----------

indexSeq:  is a size_t value indicating which sequence is the index
sequence; by default it is the 1st one.

Returns:
--------

a int value of ErrorCode ";

%feature("docstring")  Multilign_object::SetIndexSeq "int
Multilign_object::SetIndexSeq(const string seqname)

an overloaded function accepting a string as its parameter.

Parameters:
-----------

seqname:  is the seq filename that will be set as the index sequence.

Returns:
--------

an int value of ErrorCode ";

%feature("docstring")  Multilign_object::GetIndexSeq "string
Multilign_object::GetIndexSeq() const

return the filename of the index seq.

Returns:
--------

a string of index seq filename. ";

%feature("docstring")  Multilign_object::Randomize "void
Multilign_object::Randomize()

randomize the order of inputList. ";

%feature("docstring")  Multilign_object::AddOneInput "int
Multilign_object::AddOneInput(const string seq, const string ct, const
string constraint=\"\", const string shape=\"\")

add one entry into inputList.

Parameters:
-----------

seq:  is a string value of sequence filename to be appended1

ct:  is a string value of corresponding ct filename

constraint:  is a string value of corresponding constraint filename.
By default, it is empty, meaning no folding constraint exists

shape:  is string value of corresponding SHAPE filename. By default,
it is empty, meaning no SHAPE exists.

Returns:
--------

a is int value of ErrorCode ";

%feature("docstring")  Multilign_object::RemoveOneInput "int
Multilign_object::RemoveOneInput(const string seq)

remove one entry from inputList.

Parameters:
-----------

seq:  is a string value of sequence filename of which the entry in
inputList will be removed

Returns:
--------

a int value of ErrorCode ";

%feature("docstring")  Multilign_object::SetSHAPESlope "void
Multilign_object::SetSHAPESlope(const double slope=1.8)

set the slope parameter for SHAPE

Parameters:
-----------

slope:  is a double value assigned to SHAPESlope. By default, it is
set to 1.8. ";

%feature("docstring")  Multilign_object::GetSHAPESlope "double
Multilign_object::GetSHAPESlope() const

get the SHAPESlope

Returns:
--------

a SHAPESlope of double value. ";

%feature("docstring")  Multilign_object::SetSHAPEIntercept "void
Multilign_object::SetSHAPEIntercept(const double intercept=-0.6)

set the intercept parameter for SHAPE.

Parameters:
-----------

intercept:  is a double value assigned to SHAPEIntercept. By default,
it is set to -0.6. ";

%feature("docstring")  Multilign_object::GetSHAPEIntercept "double
Multilign_object::GetSHAPEIntercept() const

get the SHAPEIntercept.

Returns:
--------

SHAPEIntercept of double value. ";

%feature("docstring")  Multilign_object::SetTemperature "void
Multilign_object::SetTemperature(const double temp=310.15)

set the temperature to fold the sequences.

Parameters:
-----------

temp:  is a double value of temperature; by default it is set to
310.15K ";

%feature("docstring")  Multilign_object::GetTemperature "double
Multilign_object::GetTemperature() const

get the temperature to fold the sequences

Returns:
--------

a double value of the set temperature. ";

%feature("docstring")  Multilign_object::SetNucType "void
Multilign_object::SetNucType(const bool isrna=true)

set the flag isRNA to be true or false. By default it is true. When it
is true, RNA nearest neighbor parameters are used.

Parameters:
-----------

isrna:  ";

%feature("docstring")  Multilign_object::GetNucType "bool
Multilign_object::GetNucType() const

get the type fo nucleic acid

Returns:
--------

return true when it is of RNA prediction; otherwise, false. ";

%feature("docstring")  Multilign_object::CleanupIntermediateFiles "int Multilign_object::CleanupIntermediateFiles() const

delete intermediate pairwise dsv and aout files

Returns:
--------

an int value of error code. ";

%feature("docstring")  Multilign_object::SetProgress "void
Multilign_object::SetProgress(TProgressDialog *Progress=NULL)

Provide a TProgressDialog for following calculation progress.

Parameters:
-----------

Progress:  is a pointer to TProgressDialog ";

%feature("docstring")  Multilign_object::StopProgress "void
Multilign_object::StopProgress()

Provide a means to stop using a TProgressDialog by assigning NULL to
progress pointer. ";

%feature("docstring")  Multilign_object::GetProgress "TProgressDialog
* Multilign_object::GetProgress() const

get the progress

Returns:
--------

the pointer to TProgressDialog ";

%feature("docstring")  Multilign_object::GetInputFilenames "void
Multilign_object::GetInputFilenames()

For diagnostic purpose only. Output the input sequence, ct,
constraints, and SHAPE filenames to stdout.

The following functions are used for Diagnostic purpose only
//////////// Generally not needed, but for debugging input. ";

%feature("docstring")  Multilign_object::GetPairs "void
Multilign_object::GetPairs()

For diagnostic purpose only. Output the paired sequence filenames to
stdout. ";

%feature("docstring")  Multilign_object::get_energies "vector<float>
Multilign_object::get_energies() ";

%feature("docstring")  Multilign_object::get_dGIndex "vector<float>
Multilign_object::get_dGIndex() ";

%feature("docstring")  Multilign_object::get_pair_alignments "vector<vector<string> > Multilign_object::get_pair_alignments() ";


// File: class_oligowalk__object.xml
%feature("docstring") Oligowalk_object "

Oligowalk_object Class.

The Oligowalk_class inhereits from RNA and provides the OligoWalk
functionality. Additionally, it provides OligoScreen.

C++ includes: Oligowalk_object.h ";

%feature("docstring")  Oligowalk_object::Oligowalk_object "Oligowalk_object::Oligowalk_object(const char sequence[])

Constructor - user provides a sequence as a c string.

This constructor simply calls the underlying RNA(const char
sequence[], const bool IsRNA=true) constructor and sets IsRNA=true.
Input sequence should contain A,C,G,T,U,a,c,g,t,u,x,X. Capitalization
makes no difference. T=t=u=U. U is assumed by OligoWalk. x=X=
nucleotide that neither stacks nor pairs. For now, any unknown nuc is
considered 'X'. Note that sequences will subsequently be indexed
starting at 1 (like a biologist), so that the 0th position in the
sequence array will be nucleotide 1.

Parameters:
-----------

sequence:  is a NULL terminated c string. ";

%feature("docstring")  Oligowalk_object::Oligowalk_object "Oligowalk_object::Oligowalk_object(const char filename[], const int
type)

Constructor - user provides a filename for existing file as a c
string.

This constructor simply calls the underlying RNA(const char
filename[], const int type, const bool IsRNA=true) constructor and
sets IsRNA=true. The existing file, specified by filename, can either
be a ct file, a sequence, or an RNAstructure save file. Therefore, the
user provides a flag for the file: type = 1 => .ct file, type = 2 =>
.seq file, type = 3 => partition function save (.pfs) file, type = 4
=> folding save file (.sav). For OligoWalk calculations, types 1 and 2
are the relevant types. This constructor generates internal error
codes that can be accessed by GetErrorCode() after the constructor is
called. 0 = no error. The errorcode can be resolved to a c string
using GetErrorMessage. Note also that save files explicitly store the
thermodynamic parameters, therefore changing the backbone type as
compaared to the original calculation will not change structure
predictions.

Parameters:
-----------

filename:  is null terminated c string.

type:  is an integer that indicates the file type. ";

%feature("docstring")  Oligowalk_object::Oligowalk_object "Oligowalk_object::Oligowalk_object(const bool IsRNA=true)

Default Constructor - user provides nothing.

This constructor calls the underlying RNA(IsRNA=true) constructor.
This basic constructor is provided for performing OligoScreen
calculations, which do not need an input sequence. This constructor
needs to be called with RNA=true for Oligos to be RNA and RNA=false
for Oligos to be DNA.

Parameters:
-----------

IsRNA:  is a bool where true= RNA and false=DNA. ";

%feature("docstring")  Oligowalk_object::Oligowalk "int
Oligowalk_object::Oligowalk(const int oligo_length, const bool isDNA,
const int option, const double oligo_concentration, const int usesub,
const int start, const int stop)

Perform an OligoWalk calculation. Note that this can only be performed
once.

Parameters:
-----------

oligo_length:  is an int that gives the length of the
oligonucleotides.

isDNA:  is a bool that indicates the oligonucleotide chemistry, where
true = DNA and false = RNA.

option:  is an int that gives the calculation type, where option 1 =
break local target structure to bind oligo, option 2 = refold target
RNA after oligo binding, and option 3 = no target structure
considered.

oligo_concentration:  is a double that indicates the oligonucleotide
concentration in M.

usesub:  is an int that indicates whether suboptimal structures are to
be used, where 0 = none and 3 = use heuristic method.

start:  is an int that indicates the starting location of the walk.

stop:  is an int that indicates the ending location of the walk.

Returns:
--------

An integer that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  Oligowalk_object::GetBreakTargetDG "double
Oligowalk_object::GetBreakTargetDG(const int index)

Get the breaking target DG for a given nucleotide. This can only be
called after performing a valid OligoWalk calculation. In case of
error, the function returns a free energy change of zero. Note!: That
a free energy change of zero is also a valid folding free energy
change. Errors will also generate an internal error code, accessible
with GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the free energy change in kcal/mol. ";

%feature("docstring")  Oligowalk_object::GetDuplexDG "double
Oligowalk_object::GetDuplexDG(const int index)

Get the duplex DG for a given nucleotide. This can only be called
after performing a valid OligoWalk calculation. In case of error, the
function returns a free energy change of zero. Note!: That a free
energy change of zero is also a valid folding free energy change.
Errors will also generate an internal error code, accessible with
GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the free energy change in kcal/mol. ";

%feature("docstring")  Oligowalk_object::GetOligoOligoDG "double
Oligowalk_object::GetOligoOligoDG(const int index)

Get the bimolecular oligo-oligo DG for a given nucleotide. This can
only be called after performing a valid OligoWalk calculation. In case
of error, the function returns a free energy change of zero. Note!:
That a free energy change of zero is also a valid folding free energy
change. Errors will also generate an internal error code, accessible
with GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the free energy change in kcal/mol. ";

%feature("docstring")  Oligowalk_object::GetOligoSelfDG "double
Oligowalk_object::GetOligoSelfDG(const int index)

Get the oligo-self DG for a given nucleotide. This can only be called
after performing a valid OligoWalk calculation. In case of error, the
function returns a free energy change of zero. Note!: That a free
energy change of zero is also a valid folding free energy change.
Errors will also generate an internal error code, accessible with
GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the free energy change in kcal/mol. ";

%feature("docstring")  Oligowalk_object::GetOverallDG "double
Oligowalk_object::GetOverallDG(const int index)

Get the overall DG for a given nucleotide. This can only be called
after performing a valid OligoWalk calculation. In case of error, the
function returns a free energy change of zero. Note!: That a free
energy change of zero is also a valid folding free energy change.
Errors will also generate an internal error code, accessible with
GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the free energy change in kcal/mol. ";

%feature("docstring")  Oligowalk_object::GetTm "double
Oligowalk_object::GetTm(const int index)

Get the Tm for a given nucleotide. This can only be called after
performing a valid OligoWalk calculation. In case of error, the
function returns a free energy change of zero. Note!: That a free
energy change of zero is also a valid folding free energy change.
Errors will also generate an internal error code, accessible with
GetErrorCode().

Parameters:
-----------

index:  is an int that specifies the 5' end of an oligonucleotide
binding site on the target.

Returns:
--------

A double that is the Tm in degrees C. ";

%feature("docstring")  Oligowalk_object::WriteReport "int
Oligowalk_object::WriteReport(const char outputfilename[], const int
oligo_length, const bool isDNA, const int option, const double
oligo_concentration, const int usesub, const int start, const int
stop)

Write a report for an OligoWalk calculation. This must be called after
performing a valid OligoWalk calculation.

Parameters:
-----------

outputfilename:  is a c-string that provides a filename to which the
report will be written.

oligo_length:  is an int that gives the length of the
oligonucleotides.

isDNA:  is a bool that indicates the oligonucleotide chemistry, where
true = DNA and false = RNA.

option:  is an int that gives the calculation type, where option 1 =
break local target structure to bind oligo, option 2 = refold target
RNA after oligo binding, and option 3 = no target structure
considered.

oligo_concentration:  is a double that indicates the oligonucleotide
concentration in M.

usesub:  is an int that indicates whether suboptimal structures are to
be used, where 0 = none and 3 = use heuristic method.

start:  is an int that indicates the starting location of the walk.

stop:  is an int that indicates the ending location of the walk.

Returns:
--------

An integer that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  Oligowalk_object::OligoScreen "int
Oligowalk_object::OligoScreen(const char infilename[], const char
outfilename[])

This function runs OligoScreen.

Read a list of oligonucleotides in infilename and output thermodynamic
characteristics in outfilename. The backbone type is set when the
constructor is called. Note that oligoscreen has no target sequence,
so the default constructor OligoWalk(bool IsRNA) should be used.

Parameters:
-----------

infilename:  is a c-string that indicates the filename to be read.

outfilename:  is a c-string that indicates the name of an output file
to be written with report.

Returns:
--------

An int that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  Oligowalk_object::GetErrorMessage "char *
Oligowalk_object::GetErrorMessage(const int error)

Return error messages based on code from GetErrorCode and other error
codes. 100 = no OligoWalk data present 101 = OligoWalk has already
been performed other codes are handled by the RNA base class function
GetErrorMessage.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode()/

Returns:
--------

A pointer to a c string that provides an error message or from other
functions that return integer error codes. ";

%feature("docstring")  Oligowalk_object::~Oligowalk_object "Oligowalk_object::~Oligowalk_object()

This is the destructor. ";


// File: class_prob_scan.xml
%feature("docstring") ProbScan "C++ includes: ProbScan.h ";

%feature("docstring")  ProbScan::ProbScan "ProbScan::ProbScan(const
char sequence[], bool isRNA=true)

Constructor - user provides a sequence as a c string.

The partition function will be calculated. If the sequence is long,
this may take some time. Input sequence should contain
A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. T=t=u=U.
If IsRNA is true, the backbone is RNA, so U is assumed. If IsRNA is
false, the backbone is DNA, so T is assumed. x=X= nucleotide that
neither stacks nor pairs. For now, any unknown nuc is considered 'X'.
Note that sequences will subsequently be indexed starting at 1 (like a
biologist), so that the 0th position in the sequence array will be
nucleotide 1.

Parameters:
-----------

sequence:  is a NULL terminated c string containing the nucleotide
sequence.

isRNA:  is a bool that indicates whether this sequence is RNA or DNA.
true= RNA. false=DNA. Default is true. ";

%feature("docstring")  ProbScan::ProbScan "ProbScan::ProbScan(const
char filename[], bool from_sequence_file, bool isRNA=true)

Constructor - user provides a filename for existing file as a c
string.

The existing file, specified by filename, can either be a ct file, a
sequence, or an RNAstructure save file. Therefore, the user provides a
flag for the file: type = 1 => .ct file, type = 2 => .seq file, type =
3 => partition function save (.pfs) file, type = 4 => folding save
file (.sav). If the input file is ont a partition function save file,
the partition function will be calculated. If the sequence is long,
this may take some time. This constructor generates internal error
codes that can be accessed by GetErrorCode() after the constructor is
called. 0 = no error. The errorcode can be resolved to a c string
using GetErrorMessage. Note that the contructor needs to be explicitly
told, via IsRNA, what the backbone is because files do not store this
information. Note also that save files explicitly store the
thermodynamic parameters, therefore changing the backbone type as
compaared to the original calculation will not change structure
predictions.

Parameters:
-----------

filename:  is null terminated c string containing the path to the
input file.

from_sequence_file:  is a bool which tells the constructor whether we
are initializing from a sequence file, in which case the partition
function must be calculated

isRNA:  is a bool that indicates whether this sequence is RNA or DNA.
true= RNA. false=DNA. Default is true. ";

%feature("docstring")  ProbScan::probability_of_hairpin "double
ProbScan::probability_of_hairpin(int i, int j)

Returns probability of a hairpin closed at a specific position.

Parameters:
-----------

i:  The 5' nucleotide closing the hairpin

j:  The 3' nucleotide closing the hairpin

Returns:
--------

A double containing the probability of the hairpin ";

%feature("docstring")  ProbScan::probability_of_all_hairpins "vector<
hairpin_t > ProbScan::probability_of_all_hairpins(int min, int max,
double threshold)

Calculates the probabilities of all possible hairpins in this
sequence.

Parameters:
-----------

min:  The minimum size of a hairpin

max:  The maximum size of a hairpin

threshold:  The minimum probability for candidate hairpins

Returns:
--------

A vector of hairpin objects, containing the positions of the hairpins
and their probabilities ";

%feature("docstring")  ProbScan::probability_of_internal_loop "double
ProbScan::probability_of_internal_loop(int i, int j, int k, int l)

Returns probability of an internal loop or bulge loop closed at a
specific position.

Parameters:
-----------

i:  The 5' nucleotide closing the loop on the exterior

j:  The 3' nucleotide closing the loop on the exterior

k:  The 5' nucleotide closing the loop on the interior

l:  The 3' nucleotide closing the loop on the interior

Returns:
--------

A double containing the probability of the internal loop ";

%feature("docstring")  ProbScan::probability_of_all_internal_loops "vector< internal_loop_t >
ProbScan::probability_of_all_internal_loops(double threshold,
std::string mode=std::string(\"both\"))

Calculates the probabilities of all possible internal loops and/or
bulge loops in this sequence.

Parameters:
-----------

threshold:  the minimum probability of candidate loops

mode:  a string which indicates what type of loops should be searched
for. Allowed values are \"internal\", \"bulge\", and \"both\"

Returns:
--------

A vector of internal loop objects, containing the positions of the
loops and their probabilities ";

%feature("docstring")  ProbScan::probability_of_stack "double
ProbScan::probability_of_stack(int i, int j)

Calculates probability of a base pair stack closed at a specific
position Note that this is a special case of probability_of_helix
where the size is set to 1

Parameters:
-----------

i:  The 5' nucleotide closing the stack

j:  The 3' nucleotide closing the stack

Returns:
--------

A double containing the probability of the stack ";

%feature("docstring")  ProbScan::probability_of_helix "double
ProbScan::probability_of_helix(const int i, const int j, const int
how_many_stacks)

Calculates probability of an helix at a specific position.

Parameters:
-----------

i:  The 5' nucleotide closing the helix on the exterior

j:  The 3' nucleotide closing the helix on the exterior

how_many_stacks:  The number of base pair STACKS in the helix (this is
the number of pairs minus 1)

Returns:
--------

A double containing the probability of the helix ";

%feature("docstring")  ProbScan::probability_of_all_helices "std::vector< basestack_t > ProbScan::probability_of_all_helices(double
threshold, int length)

Calculates the probabilities of all possible helices in this sequence
of a specific length.

Parameters:
-----------

threshold:  the minimum probability of candidate helices

length:  the number of base pair stacks to search for

Returns:
--------

A vector of helix objects, containing the positions of the helices and
their probabilities ";

%feature("docstring")  ProbScan::probability_of_multibranch_loop "double ProbScan::probability_of_multibranch_loop(const
multibranch_loop_t &mb)

Calculates probability of a multibranch loop at a specific position.

Parameters:
-----------

mb:  A multibranch loop object, containing a vector of pairs
describing the multibranch loop. These can be created with the
multibranch_loop function. See the text interface for the ProbScan
program for an example of usage.

Returns:
--------

A double containing the probability of the multibranch loop ";


// File: class_r_n_a.xml
%feature("docstring") RNA "

RNA Class.

The RNA class provides an entry point for all the single sequence
operations of RNAstructure.

C++ includes: RNA.h ";

%feature("docstring")  RNA::RNA "RNA::RNA(const char sequence[],
const bool IsRNA=true)

Constructor - user provides a sequence as a c string.

Input sequence should contain A,C,G,T,U,a,c,g,t,u,x,X. Capitalization
makes no difference. T=t=u=U. If IsRNA is true, the backbone is RNA,
so U is assumed. If IsRNA is false, the backbone is DNA, so T is
assumed. x=X= nucleotide that neither stacks nor pairs. For now, any
unknown nuc is considered 'X'. Note that sequences will subsequently
be indexed starting at 1 (like a biologist), so that the 0th position
in the sequence array will be nucleotide 1.

Parameters:
-----------

sequence:  is a NULL terminated c string.

IsRNA:  is a bool that indicates whether this sequence is RNA or DNA.
true= RNA. false=DNA. Default is true. ";

%feature("docstring")  RNA::RNA "RNA::RNA(const char filename[],
const int type, const bool IsRNA=true)

Constructor - user provides a filename for existing file as a c
string.

The existing file, specified by filename, can either be a ct file, a
sequence, or an RNAstructure save file. Therefore, the user provides a
flag for the file: type = 1 => .ct file, type = 2 => .seq file, type =
3 => partition function save (.pfs) file, type = 4 => folding save
file (.sav). This constructor generates internal error codes that can
be accessed by GetErrorCode() after the constructor is called. 0 = no
error. The errorcode can be resolved to a c string using
GetErrorMessage. Note that the contructor needs to be explicitly told,
via IsRNA, what the backbone is because files do not store this
information. Note also that save files explicitly store the
thermodynamic parameters, therefore changing the backbone type as
compaared to the original calculation will not change structure
predictions.

Parameters:
-----------

filename:  is null terminated c string.

type:  is an integer that indicates the file type.

IsRNA:  is a bool that indicates whether this sequence is RNA or DNA.
true= RNA. false=DNA. Default is true. ";

%feature("docstring")  RNA::RNA "RNA::RNA(const bool IsRNA=true)

Default Constructor - user provides nothing. This basic constructor is
provided for bimolecular folding and should not generally need to be
accessed by end users of the RNA class.

Parameters:
-----------

IsRNA:  is a bool that indicates whether this sequence is RNA or DNA.
true= RNA. false=DNA. Default is true. ";

%feature("docstring")  RNA::GetErrorCode "int RNA::GetErrorCode()

Return an error code, where a return of zero is no error.

This function returns and error flag that is generated during
construction by RNA(const char &filename, const int type, const bool
IsRNA=true) or from CalculateFreeEnergy(). An error of zero is always
no error. Other codes are errors and a c-string can be fetched for the
error with GetErrorMessage().

Returns:
--------

An integer that provides the error code. ";

%feature("docstring")  RNA::GetErrorMessage "char *
RNA::GetErrorMessage(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

0 = no error 1 = input file not found 2 = error opening file 3 =
structure number out of range 4 = nucleotide number out of range 5 =
error reading thermodynamic parameters 6 = pseudoknot formation 7 =
non-canonical pair 8 = too many restraints specified 9 = same
nucleotide in conflicting restraint 10 = no structures to write 11 =
nucleotide not a U (caused by ForceFMNCleavage() 12 = distance too
short 13 = error reading constraint file 14 = traceback error 15 = no
partition function data present 16 = incorrect save file version used
17 = cannot be performed without having read a save file (.sav) 18 =
threshold is too low to be valid 19 = drawing coordinates have not
been determined 20 = no sequence has been read 21 = over 1 probability
error on stochastic traceback 22 = programming error, unrecognized
input to constructor 23 = no structures present 24 = too few
iterations 25 = index (for drawing) is not a multiple of 10

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode() or from
other functions that return integer error codes.

Returns:
--------

A pointer to a c string that provides an error message. ";

%feature("docstring")  RNA::GetErrorMessageString "std::string
RNA::GetErrorMessageString(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

Although RNA generally uses c strings, this member function returns a
string that is suitable for interfacing with JAVA, etc. See the error
list in the GetErrorMessage() entry.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode() or from
other functions that return integer error codes.

Returns:
--------

A string that provides an error message. ";

%feature("docstring")  RNA::SpecifyPair "int RNA::SpecifyPair(const
int i, const int j, const int structurenumber=1)

Specify a base pair between nucleotides i and j.

The base pair is in structure number structurenumber, which is assumed
to be structure 1. Return 0 if there is no problem, otherwise return
an error code: error = 3 -> structurenumber out of range. error = 4 ->
nucleotide number out of range. A c string or string description of
the error are available using GetErrorMessage() or
GetErrorMessageString(). Note!: Sequences with the 5' end = nucleotide
1. Note!: Structures start at structure 1.

Parameters:
-----------

i:  is an integer for the position of the first nucleotide in the
pair.

j:  in an integer for the position of the second nucleotide in the
pair.

structurenumber:  is the structure that has the pair. This defaults to
1.

Returns:
--------

An integer that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  RNA::RemovePairs "int RNA::RemovePairs(const
int structurenumber=1)

Remove all the current base pairs in a specified structure.

Return 0 if there is no error. Return 5 if structurenumber never had
pairs specified.

Parameters:
-----------

structurenumber:  is an integer specifying the structure from which to
remove the pairs.

Returns:
--------

An integer that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  RNA::RemoveBasePair "int
RNA::RemoveBasePair(const int i, const int structurenumber=1)

Remove a specified pair in a specified structure.

Break the pair between i and i's pairing partner Return 0 if there is
no error. Return 3 if structurenumber out of range. Return 4 if
nucleotide number out of range.

Parameters:
-----------

i:  is the index of a nucleotide in a pair that will be broken.

structurenumber:  is an integer specifying the structure from which to
remove the pairs.

Returns:
--------

An integer that indicates an error code that can be parsed by
GetErrorMessage() or GetErrorMessageString(), 0 = no error. ";

%feature("docstring")  RNA::CalculateFreeEnergy "double
RNA::CalculateFreeEnergy(const int structurenumber=1, const bool
UseSimpleMBLoopRules=false)

Return the predicted Gibb's free energy change for structure #
structurenumber, defaulted to 1.

Free energies are in kcal/mol. The first time this is called, if no
other free energy calculation has been performed and the folding
temperature has not been specifed, thermodynamic parameter files
(.dat) files will be read from disk. The parameter files should be
located in the directory specified by environment variable $DATAPATH,
or the pwd. In case of error, the function returns a free energy
change of zero. Note!: That a free energy change of zero is also a
valid folding free energy change. Errors will also generate an
internal error code, accessible with GetErrorCode(). GetErrorCode()
will return 0 when there is no error and other codes can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

structurenumber:  is an integer that refers to the index of the
structure for which to calculate the folding free energy change. This
defaults to 1.

UseSimpleMBLoopRules:  is a bool that indicates what energy rules to
use. The default, false, uses the complete nearest neighbor model for
multibranch loops. When true is passed, the energy model is instead a
simplified model that is the one used by the dynamic programming
algorithms.

Returns:
--------

A double which is the folding free energy change in kcal/mol. ";

%feature("docstring")  RNA::WriteThermodynamicDetails "int
RNA::WriteThermodynamicDetails(const char filename[], const bool
UseSimpleMBLoopRules=false)

Calculate the folding free energy change for all structures and write
the details of the calculation to a file.

Free energies are in kcal/mol. The first time this is called, if no
other free energy calculation has been performed and the folding
temperature has not been specifed, thermodynamic parameter files
(.dat) files will be read from disk. The parameter files should be
located in the directory specified by environment variable $DATAPATH,
or the pwd. In case of error, the function returns a non-zero.

Parameters:
-----------

filename:  is a NULL terminated c string that provides the name of the
output file to be written.

UseSimpleMBLoopRules:  is a bool that indicates what energy rules to
use. The default, false, uses the complete nearest neighbor model for
multibranch loops. When true is passed, the energy model is instead a
simplified model that is the one used by the dynamic programming
algorithms.

Returns:
--------

An int that indicates whether an error occurred (0 = no error; 5 =
error reading parameter files). ";

%feature("docstring")  RNA::FoldSingleStrand "int
RNA::FoldSingleStrand(const float percent=20, const int
maximumstructures=20, const int window=5, const char savefile[]=\"\",
const int maxinternalloopsize=30, bool mfeonly=false)

Predict the lowest free energy secondary structure and generate
suboptimal structures using a heuristic.

This function predicts the lowest free energy structure and suboptimal
structures. If the temperature has not been specified using
SetTemperature and no free energies have been calculated, the
thermodynamic parameters have not been read and therefore they will be
read by this function call. The parameter files should be located in
the directory specified by the environment variable $DATAPATH of the
pwd. In case of error, the function returns a non-zero that can be
parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 20.

maximumstructures:  is the maximum number of suboptimal structures to
generate. The default is 20.

window:  is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). The defaults is 5,
but this should be customized based on sequence length.

savefile:  is c string containing a file path and name for a savefile
(.sav)that can be used to generate energy dot plots and to refold the
secondary structure using different suboptimal structure parameters.
The default is \"\", which results in no save file written.

maxinternalloopsize:  is the maximum number of unpaired nucleotides in
bulge and internal loops. This is used to accelerate the prediction
speed. The default is 30.

mfeonly:  is a bool that indicates whether only the minimum free
energy structure will be generated. This saves half the calculation
time, but no save file can be generated. Default is false.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files, 14 = traceback error). ";

%feature("docstring")  RNA::GenerateAllSuboptimalStructures "int
RNA::GenerateAllSuboptimalStructures(const float percent=5, const
double deltaG=0.6)

Predict the lowest free energy secondary structure and generate all
suboptimal structures.

This function predicts the lowest free energy structure and suboptimal
structures. If the temperature has not been specified using
SetTemperature and no free energies have been calculated, the
thermodynamic parameters have not been read and therefore they will be
read by this function call. The parameter files should be located in
the directory specified by the environment variable $DATAPATH of the
pwd. In case of error, the function returns a non-zero that can be
parsed by GetErrorMessage() or GetErrorMessageString(). Two controls
are available for limiting the number of structures, the maximum %
difference in energy (percent) and the maximum absolute change in
energy (deltaG). The smaller of the two will be used as the limit.

Parameters:
-----------

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 5.

deltaG:  is the maximum difference in free energy change above the
lowest free energy structure (in kcal/mol). The defaults is 0.6
kcal/mol.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::MaximizeExpectedAccuracy "int
RNA::MaximizeExpectedAccuracy(const double maxPercent=20, const int
maxStructures=20, const int window=1, const double gamma=1.0)

Predict the structure with maximum expected accuracy and suboptimal
structures.

This function predicts structures composed of probable base pairs and
single-srtranded nucleotide, weighted by gamma. The score for a
structure is = gamma * 2 * (sum of pairing probabilities for pairs) +
(sum of unpairing probabilities for single stranded nucleotides). This
function requires partition function data from either a previous
partition function calculations or from having read a partition
function save file during construction of the class. In case of error,
the function returns a non-zero that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

maxPercent:  is the maximum percent difference is score in generating
suboptimal structures. The default is 20.

maxStructures:  is the maximum number of suboptimal structures to
generate. The default is 20.

window:  is the window parameter, where a higher value generates
suboptimal structures that are more different from each other. The
default is 1.

gamma:  is the weight given to base pairs

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::PartitionFunction "int
RNA::PartitionFunction(const char savefile[]=\"\", double
temperature=-10.0)

Predict the partition function for a sequence.

This function must be called to predict base pair probabilities,
perform stochastic traceback, or for maximizing expected accuracy. If
the temperature has not been specified using SetTemperature and no
free energies have been calculated, the thermodynamic parameters have
not been read and therefore they will be read by this function call.
The parameter files should be located in the directory specified by
the environment variable $DATAPATH of the pwd. In case of error, the
function returns a non-zero that can be parsed by GetErrorMessage() or
GetErrorMessageString(). Note that the parameter temperature is used
when calculating equilibrium constants, but does not change the
temperature at which the free energies are determined. SetTemperature,
from the underlying base class Thermodynamics, should be used to
change the temperature for most calculations. This parameter should
generally not be used. The default is -10.0 and values below zero
cause this parameter to be ignored (the correct default behavior).
Note also that if SetTemperature is not used, the temperature defaults
to 310.15 K (37 deg. C), which is the desired behavior for most
purposes.

Parameters:
-----------

savefile:  is a c string that contains the path and filename for
creating a save file. This defaults to \"\", which indicates no file
is to be written.

temperature:  is a double that indicates a pseudo-temperature for
calculating equilibrium constants from free energies at fixed
temperature previously specified.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files). ";

%feature("docstring")  RNA::PredictProbablePairs "int
RNA::PredictProbablePairs(const float probability=0)

Predict structures containing highly probable pairs.

This function predicts structures composed of probable base pairs.
This function requires partition function data from either a previous
partition function calculations or from having read a partition
function save file during construction of the class. In case of error,
the function returns a non-zero that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

probability:  is the pairing probability threshold, where pairs will
be predicted if they have a higher probability. Note that a value of
less than 0.5 (50%), will cause an error. The default value of zero
will trigger the creation of 8 structures, with thresholds of >=0.99,
>=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::ProbKnot "int RNA::ProbKnot(int
iterations=1, int MinHelixLength=1)

Predict maximum expected accuracy structures that contain pseudoknots
from either a sequence or a partition function save file.

This function uses base pair probabilities to predict structures that
contains pseudoknots. This function requires partition function data
from either a previous partition function calculations or from having
read a partition function save file during construction of the class.
In case of error, the function returns a non-zero that can be parsed
by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

iterations:  is the number of iterations of pair selection that are
performed. The default and recommended value is 1.

MinHelixLength:  is the shortest helix that is allowed. If this is set
>1, a post- processing step is performed to remove short helices.
Default = 1, i.e. no post-processing.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::ProbKnotFromSample "int
RNA::ProbKnotFromSample(int iterations=1, int MinHelixLength=1)

Predict maximum expected accuracy structures that contain pseudoknots
from a file containing ensemble of structures.

This function uses base pair probabilities to predict structures that
contains pseudoknots. This function requires a file with ensemble of
structures. This function processes the file to calculate pair
probabilities. In case of error, the function returns a non-zero that
can be parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

iterations:  is the number of iterations of pair selection that are
performed. The default and recommended value is 1.

MinHelixLength:  is the shortest helix that is allowed. If this is set
>1, a post- processing step is performed to remove short helices.
Default = 1, i.e. no post-processing.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::ReFoldSingleStrand "int
RNA::ReFoldSingleStrand(const float percent=20, const int
maximumstructures=20, const int window=5)

Re-predict the lowest free energy secondary structure and generate
suboptimal structures using a heuristic.

This function predicts the lowest free energy structure and suboptimal
structure after a save file (.sav) was specified to the constructor.
The step of predicting structures from the save file is rapid, so this
is laregely a method to quickly generate a different set of suboptimal
structures. Refolding can only be performed if the RNA constructor was
called with a save file name. (That is, you cannot call this after
calling fold single strand, without loading the data from disk with a
new instance of RNA. This is for historical reasons.) In case of
error, the function returns a non-zero that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

percent:  is the maximum % difference in free energy in suboptimal
structures from the lowest free energy structure. The default is 20.

maximumstructures:  is the maximum number of suboptimal structures to
generate. The default is 20.

window:  is a parameter that specifies how different the suboptimal
structures should be from each other (0=no restriction and larger
integers require structures to be more different). The default is 5.

Returns:
--------

An int that indicates an error code (0 = no error, 5 = error reading
thermodynamic parameter files, 14 = traceback error). ";

%feature("docstring")  RNA::Stochastic "int RNA::Stochastic(const int
structures=1000, const int seed=1)

Sample structures from the Boltzman ensemable.

This function requires partition function data from either a previous
partition function calculations or from having read a partition
function save file during construction of the class. In case of error,
the function returns a non-zero that can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

structures:  is the number of structures to be sampled. The default is
1000.

seed:  is an integer that seeds the random number generator that is
required for sampling, which defaults to 1.

Returns:
--------

An int that indicates an error code (0 = no error, non-zero = error).
";

%feature("docstring")  RNA::ForceDoubleStranded "int
RNA::ForceDoubleStranded(const int i)

Force a nucleotide to be double stranded (base paired).

This function indicates a nucleotide that is double stranded (paired).
In subsequent structure prediction, this nucleotide will be double
stranded. The function returns 0 with no error and a non-zero
otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

i:  is the index of the paired nucleotide.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 8 = too many restraints specified, 9 = same nucleotide
in conflicting restraint). ";

%feature("docstring")  RNA::ForceFMNCleavage "int
RNA::ForceFMNCleavage(const int i)

Indicate a nucleotide that is accessible to FMN cleavage (a U in GU
pair).

In subsequent structure prediction, this nucleotide will be in a GU
pair. The function returns 0 with no error and a non-zero otherwise
that can be parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

i:  is the index of the FMN-cleaved nucleotide.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 8 = too many restraints specified, 9 = same nucleotide
in conflicting restraint, 11 = nucleotide not U). ";

%feature("docstring")  RNA::ForceMaximumPairingDistance "int
RNA::ForceMaximumPairingDistance(const int distance)

Force a maximum distance between apired nucleotides.

In a subsequent structure prediction, there will be no pairs allowed
between nucleotides more distant than distance, i.e. |j-i| < distance
for i to pair to j. The function returns and error code; 0==no error,
12== too long or too short distance.

Parameters:
-----------

distance:  is the maximum pairing distance.

Returns:
--------

An integer that indicates an error code (0 = no error, 12 = too
short). ";

%feature("docstring")  RNA::ForceModification "int
RNA::ForceModification(const int i)

Force modification for a nucleotide.

This function indicates a nucleotide that is accessible to chemical
modification. In subsequent structure prediction, this nucleotide will
be single stranded, at the end of a helix, or in or adjacent to a GU
pair. The function returns 0 with no error and a non-zero otherwise
that can be parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

i:  is the index of the nucleotide accessible to chemical
modification.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 8 = too many restraints specified). ";

%feature("docstring")  RNA::ForcePair "int RNA::ForcePair(const int
i, const int j)

Force a pair between two nucleotides.

This function forces a pair between two nucleotides in subsequent
structure predictions. When multiple pairs are specified, the pairs
must not force a pseudoknot. The function returns 0 with no error and
a non-zero otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

i:  is the index of one nucleotide in the pair.

j:  is the index of the second nucleotide in the pair.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 6 = pseudoknot formation, 7 = non-canonical pair, 8 =
too many restraints specified, 9 = same nucleotide in conflicting
restraint). ";

%feature("docstring")  RNA::ForceProhibitPair "int
RNA::ForceProhibitPair(const int i, const int j)

Prohibit a pair between two nucleotides.

This function prevents a pair between two nucleotides in subsequent
structure predictions. The function returns 0 with no error and a non-
zero otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

i:  is the index of one nucleotide in the pair.

j:  is the index of the second nucleotide in the pair.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 8 = too many restraints specified, 9 = nucleotide in
conflicting restraint). ";

%feature("docstring")  RNA::ForceSingleStranded "int
RNA::ForceSingleStranded(const int i)

Force a nucleotide to be single stranded.

This function indicates a nucleotide that is single stranded. In
subsequent structure prediction, this nucleotide will be single
stranded. The function returns 0 with no error and a non-zero
otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

i:  is the index of the nucleotide that is single stranded.

Returns:
--------

An integer that indicates an error code (0 = no error, 4 = nucleotide
out of range, 8 = too many restraints specified, 9 = same nucleotide
in conflicting restraint). ";

%feature("docstring")  RNA::GetForcedDoubleStranded "int
RNA::GetForcedDoubleStranded(const int constraintnumber)

Return a nucleotide that is forced double stranded.

This function returns a nucleotide that is constrainted to be paired.
Constraints are numbered from zero to
GetNumberofForcedDoubleStranded()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetForcedFMNCleavage "int
RNA::GetForcedFMNCleavage(const int constraintnumber)

Return a nucleotide that is accessible to FMN cleavage.

This function returns a nucleotide that is constrainted to be
accessible to FMN cleavage (a U in a GU pair). Constraints are
numbered from zero to GetNumberofForcedFMNCleavages()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetForcedModification "int
RNA::GetForcedModification(const int constraintnumber)

Return a nucleotide that is accessible to modification.

This function returns a nucleotide that is constrainted to be
accessible to chemical modification. Constraints are numbered from
zero to GetNumberofModifications()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetForcedPair "int
RNA::GetForcedPair(const int constraintnumber, const bool fiveprime)

Return a nucleotide in a forced pair.

This function returns either the five prime or three prime nucleotide
in a forced pair constraint, depending on the value of fiveprime.
Constraints are numbered from zero to GetNumberofForcedPairs()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

fiveprime:  determines if the nucleotide is the five prime or the
three prime nucleotide in the constraint. true = five prime
nucleotide.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetForcedProhibitedPair "int
RNA::GetForcedProhibitedPair(const int constraintnumber, const bool
fiveprime)

Return a nucleotide in a prohibited pair.

This function returns either the five prime or three prime nucleotide
in a prohibited pair constraint, depending on the value of fiveprime.
Constraints are numbered from zero to GetNumberofForcedProhibited()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

fiveprime:  determines if the nucleotide is the five prime or the
three prime nucleotide in the constraint. true = five prime
nucleotide.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetForcedSingleStranded "int
RNA::GetForcedSingleStranded(const int constraintnumber)

Return a nucleotide that is forced single stranded.

This function returns a nucleotide that is constrainted to be single
stranded. Constraints are numbered from zero to
GetNumberofForcedSingleStranded()-1.

Parameters:
-----------

constraintnumber:  is the index to the constraint number.

Returns:
--------

An integer that is the nucleotide index. If the constraintnumber is
for a constraint that does not exist, zero is returned. ";

%feature("docstring")  RNA::GetMaximumPairingDistance "int
RNA::GetMaximumPairingDistance()

Return the maximum pairing distance.

return An integer that indicates the maximum distance allowed between
paired nucleotides, where -1 indicates that the maximum distance is
not set. ";

%feature("docstring")  RNA::GetNumberOfForcedDoubleStranded "int
RNA::GetNumberOfForcedDoubleStranded()

Return the number of nucletides forced to be paired.

Returns:
--------

An integer that indicates the number of nucleotides that are forced
pair. ";

%feature("docstring")  RNA::GetNumberOfForcedFMNCleavages "int
RNA::GetNumberOfForcedFMNCleavages()

Return the number of nucleotides accessible to FMN cleavage.

Returns:
--------

An integer that indicates the number of FMN cleavage nucleotides (Us
in GU pairs). ";

%feature("docstring")  RNA::GetNumberOfForcedModifications "int
RNA::GetNumberOfForcedModifications()

Return the number of nucleotides accessible to chemical modification.

Returns:
--------

An integer that indicates the number of modified nucleotides. ";

%feature("docstring")  RNA::GetNumberOfForcedPairs "int
RNA::GetNumberOfForcedPairs()

Return the number of forced base pairs.

Returns:
--------

An integer that indicates the number of forced pairs. ";

%feature("docstring")  RNA::GetNumberOfForcedProhibitedPairs "int
RNA::GetNumberOfForcedProhibitedPairs()

Return the number of prohibited base pairs.

Returns:
--------

An integer that indicates the number of pairs that are prohibited. ";

%feature("docstring")  RNA::GetNumberOfForcedSingleStranded "int
RNA::GetNumberOfForcedSingleStranded()

Return the number of nucleotides that are not allowed to pair.

Returns:
--------

An integer that indicates the number of nucleotides not allowed to
pair. ";

%feature("docstring")  RNA::ReadConstraints "int
RNA::ReadConstraints(const char filename[])

Read a set of folding constraints to disk in a plain text file.

The file format for constraints is that generated by the
WriteConstraints() function. The function returns 0 with no error and
a non-zero otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString(). Note that calling ReadConstraints() will
erase previously defined constraints (except for SHAPE pseudoenergy
restraints).

Parameters:
-----------

filename:  is a c string that is the file name to be read.

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = file not
found, 13 = error reading constraint file). ";

%feature("docstring")  RNA::ReadSHAPE "int RNA::ReadSHAPE(const char
filename[], const double parameter1, const double parameter2,
std::string modifier=\"SHAPE\", const bool IsPseudoEnergy=true)

Read SHAPE data from disk.

The SHAPE data is used to constrain structure prediction on subsequent
structure predictions. The function returns 0 with no error and a non-
zero otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString(). Pseudo folding free energy change parameters
should be in units of kcal/mol.

Parameters:
-----------

filename:  is a c string that indicates a file that contains SHAPE
data.

IsPseudoEnergy:  indicates whether this is the pseudo folding free
energy constraint (the preferred method). This defaults to true.

parameter1:  is the slope when IsPseudoEnergy=true and is a threshold
above which nucleotides are forced single stranded otherwise.

parameter2:  is the intercept when IsPseudoEnergy=true and is a
threshold above which a nucleotide is considered chemically modified
otherwise.

modifier:  is the type of chemical modification probe that was used
(currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT).
Defaults to SHAPE.

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = input file
not found). ";

%feature("docstring")  RNA::ReadSHAPE "int RNA::ReadSHAPE(const char
filename[], const double parameter1, const double parameter2, const
double ssm, const double ssb, std::string modifier=\"SHAPE\")

Read SHAPE data from disk including single-stranded SHAPE pseudo free
energys.

The SHAPE data is used to constrain structure prediction on subsequent
structure predictions. This version of the overloaded function
includes a single-stranded pseudo free energy change. The function
returns 0 with no error and a non-zero otherwise that can be parsed by
GetErrorMessage() or GetErrorMessageString(). Pseudo folding free
energy change parameters should be in units of kcal/mol.

Parameters:
-----------

filename:  is a c string that indicates a file that contains SHAPE
data.

parameter1:  is the double-stranded slope.

parameter2:  is the double-stranded intercept.

modifier:  is the type of chemical modification probe that was used
(currently accepted values are SHAPE, DMS, and CMCT). Defaults to
SHAPE.

ssm:  is the single-stranded slope.

ssb:  is the single-stranded intercept.

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = input file
not found). ";

%feature("docstring")  RNA::ReadDSO "int RNA::ReadDSO(const char
filename[])

Read double strand offset data from disk.

The double strand offset is data that is used to constrain structure
prediction on subsequent structure predictions. This is a free energy
in kcal/mol that is added to a specific nucleotide that is double
stranded. The function returns 0 with no error and a non-zero
otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

filename:  is a c string that indicates a file that contains data, in
a raw format with nucleotide index and offset (one set per line).

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = input file
not found). ";

%feature("docstring")  RNA::ReadSSO "int RNA::ReadSSO(const char
filename[])

Read single strand offset data from disk.

The single strand offset is data that is used to constrain structure
prediction on subsequent structure predictions. This is a free energy
in kcal/mol that is added to a specific nucleotide that is single
stranded. The function returns 0 with no error and a non-zero
otherwise that can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

filename:  is a c string that indicates a file that contains data, in
a raw format with nucleotide index and offset (one set per line).

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = input file
not found). ";

%feature("docstring")  RNA::ReadExperimentalPairBonus "int
RNA::ReadExperimentalPairBonus(const char filename[], double const
experimentalOffset, double const experimentalScaling)

Read experimental pair bonuses from disk.

This is a quantity that results in a bonus added to a specific pair,
once per stack, so that pairs in the middle of a helix get the bonus
twice and those at the end of a helix get the bonus once. The bonus is
in the form of experimentalScaling*value + experimentalOffset. The
data is formatted using a simple square matrix of values and no
headers. The format requires that there be N^2 entries for a sequence
of N nucleotides.

Parameters:
-----------

filename:  is a c string that indicates a file that contains data.

experimentalOffset:  is a double that is added to each value.

experimentalScaling:  is a double by which each value is multiplied.

Returns:
--------

An integer that indicates an error code (0 = no error, 1 = input file
not found). ";

%feature("docstring")  RNA::RemoveConstraints "void
RNA::RemoveConstraints()

Remove all folding constraints.

This function strips all previously assigned folding constraints. Note
that this function does not delete SHAPE constraints or pseudo free
energies. ";

%feature("docstring")  RNA::SetExtrinsic "int RNA::SetExtrinsic(int
i, int j, double k)

Add extrinsic restraints for partition function calculations.

This function multiplies the equilibrium constant for structures
including the i-j basepair by k. This applies only to partition
functions and to stochastic traceback. If k>1, then the i-j pair is
favored and if k<1, the i-j pair is disfavored. k should always be >=
0. In case of error, the function returns a non-zero that can be
parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

i:  is the index of a nucleotide in the i-j pair.

j:  is the index of the other nucleotide in the i-j pair.

k:  is an equilibrium constant that is >= 0.

Returns:
--------

An integer that indicates an error code (0 = no error, >0 indicates an
error). ";

%feature("docstring")  RNA::WriteConstraints "int
RNA::WriteConstraints(const char filename[])

Write the current set of folding constraints to disk in a plain text
file.

This function does not write SHAPE pseudo energies.

Parameters:
-----------

filename:  is a c string that is the file name to be written.

Returns:
--------

An integer that indicates an error code (0 = no error). Currently,
this function does not generate errors, but the return is provided to
add error handling in the future. ";

%feature("docstring")  RNA::AddComment "int RNA::AddComment(const
char comment[], const int structurenumber=1)

Add a comment associated with a structure.

This comment will appear in a written .ct file. The comment is
appended to any existing comments, like titles read from .seq files.
This function is especially useful if the constructor is used in which
a character array is provided with the sequence. In that case, there
is no sequence title read. The function returns 0 in the case of no
errors, or 3 if the structurenumber is invalid. An error message can
be retrieved using GetErrorMessage() called with the errorcode.

Parameters:
-----------

comment:  is a character array that contains a null terminated
c-string with the comment to be registered.

structurenumber:  is an integer that specifies to which structure the
comment should be added.

Returns:
--------

An integer that contains an error code, where 0 is no error and non-
zero is an error. ";

%feature("docstring")  RNA::WriteCt "int RNA::WriteCt(const char
filename[], bool append=false)

Write a ct file of the structures.

Return 0 if no error and non-zero errors can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

filename:  is a NULL terminated c string that specifies the name of
the ct file to be written.

append:  is a bool that indiactes whether the ct data should be
appended to an existing file. If true, data will be appended if the
file exists, or a new file created if the file does not exist. If
false, any esiting file is overwritten. This is false by default.

Returns:
--------

An integer that provides an error code. 0 = no error, 10 = no
structure to write. ";

%feature("docstring")  RNA::WriteDotBracket "int
RNA::WriteDotBracket(const char filename[])

Write dot-bracket file of structures.

Return 0 if no error and non-zero errors can be parsed by
GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

filename:  is a NULL terminated c string that specified the name of
the file to be written.

Returns:
--------

An integer that provides an error code. 0 = no error. ";

%feature("docstring")  RNA::BreakPseudoknot "int
RNA::BreakPseudoknot(const bool minimum_energy=true, const int
structurenumber=0)

Break any pseudoknots that might be in a structure.

This function uses the method of Smit et al. to break pseudoknots by
running a dynamic programming algorithm that cannot predict
pseudoknots while only allowing the pairs that already exist in the
structure. When minimum_energy = true (the default), this function
predicts the lowest free energy structure that has no pseudoknots.
Note that when minimum_energy = true, this function might additionally
break pairs that are not pseudoknotted if the pairs increase the
folding free energy change or are forbidden (as in an isolated pair).
Also note that when minum_energy=true, this function uses the
GenerateAllSuboptimalStructures methodology behind the scenes, so
large internal loops in the input would lead to a loss of pairs. When
minumum_energy is set to false, this function maximizes the number of
base pairs in the pseudoknot free structure. Return 0 if no error and
non-zero errors can be parsed by GetErrorMessage() or
GetErrorMessageString().

Parameters:
-----------

minimum_energy:  is a bool thgat indicates where the structure should
be minimum in free energy (true) or maximize pairs (false).

structurenumber:  is an int that indicates a specific structure for
which to break pseudoknots (indexed from 1). The default value, 0,
indicates that all structures should have pseudoknots broken.

Returns:
--------

an int that provides an error code. 0 = no error. ";

%feature("docstring")  RNA::ContainsPseudoknot "bool
RNA::ContainsPseudoknot(const int structurenumber)

Report if there are any pseudoknots in a structure.

This method checks for any \"crossing pairs,\" i.e. i-j and i'-j' s.t.
i < i' < j < jp. If there is at least one crossing pair set, then
there is a pseudoknot and the function returns true. This function
generates internal error codes that can be accessed by GetErrorCode()
after the constructor is called: 0 = no error, nonzero = error. The
errorcode can be resolved to a c string using GetErrorMessage.

Parameters:
-----------

structurenumber:  is an int that indicates the structure number to
check. Note that indexing of structures starts with structure #1.

Returns:
--------

A bool that indicates whether there is a pseudoknot. ";

%feature("docstring")  RNA::GetEnsembleEnergy "double
RNA::GetEnsembleEnergy()

Get the ensemble folding free energy change.

Returns the ensemble folding free energy change as determined by the
partition function. This is a handy way of getting the size of the
partition function Q, which itself is too large to fit in a double for
all but the shortest sequences. The ensemble folding free energy
change = -RT ln (Q). Function requires that the partition function
data be present either because PartitionFunction() has been called or
the constructor that reads a partition function save was used. This
function generates internal error codes that can be accessed by
GetErrorCode() after the constructor is called: 0 = no error, nonzero
= error. The errorcode can be resolved to a c string using
GetErrorMessage.

Returns:
--------

A double that is the ensemble folding free energy change in kcal/mol.
";

%feature("docstring")  RNA::GetFreeEnergy "double
RNA::GetFreeEnergy(const int structurenumber)

Get the folding free energy change for a predicted structure.

Returns the folding free energy change of structure i as determined by
a previous folding calculation. Function requires that the structure
be predicted by a structure prediction method. This function generates
internal error codes that can be accessed by GetErrorCode() after the
constructor is called: 0 = no error, nonzero = error. The errorcode
can be resolved to a c string using GetErrorMessage.

Parameters:
-----------

structurenumber:  is an integer indicating the predicted structure
number.

Returns:
--------

A double that is the folding free energy change in kcal/mol. ";

%feature("docstring")  RNA::GetPair "int RNA::GetPair(const int i,
const int structurenumber=1)

Get the nucleotide to which the specified nucleotide is paired.

Returns the pairing partner of the ith nucleotide in structure number
structurenumber. Zero means the nucleotide is unpaired. This function
generates internal error codes that can be accessed by GetErrorCode()
after the constructor is called: 0 = no error, nonzero = error. The
errorcode can be resolved to a c string using GetErrorMessage.

Parameters:
-----------

i:  is an int that indicates the nucleotide to which the pairing
partner is being queried.

structurenumber:  is an int that indicates the structure number, where
the default is 1.

Returns:
--------

An int that indicates the other nucleotide in pair, where 0 is no
paired. ";

%feature("docstring")  RNA::GetPairEnergy "double
RNA::GetPairEnergy(const int i, const int j)

Get the lowest folding free energy possible for a structure containing
pair i-j.

Returns a folding free energy change in kcal/mol for use in energy dot
plots. This function requires that the RNA constructor be called with
a save file (.sav) name. (That is, for historical resaons, this cannot
be called after FoldSingleStrand without writing the data to disk.)
This function generates internal error codes that can be accessed by
GetErrorCode() after the constructor is called: 0 = no error, nonzero
= error. The errorcode can be resolved to a c string using
GetErrorMessage. param i and j are ints that provide indexes the 5'
and 3' nucleotides, respectively, in a pair. return A double that is
the folding free energy change in kcal/mol. ";

%feature("docstring")  RNA::GetPairProbability "double
RNA::GetPairProbability(const int i, const int j)

Get a base pair probability.

Returns the base pair probability for the pair between i and j.
Function requires that the partition function data be present either
because PartitionFunction() has been called or the constructor that
reads a partition function save was used. This function generates
internal error codes that can be accessed by GetErrorCode(): 0 = no
error, nonzero = error. The errorcode can be resolved to a c string
using GetErrorMessage.

Parameters:
-----------

i:  provides the 5' nucleotide in a pair.

j:  provides the 3' nucleotides in a pair.

Returns:
--------

A double that is the base pair probability. If i and j cannot pair,
0.0 is returned. If an error occurs, 0.0 is returned. ";

%feature("docstring")  RNA::GetStructureNumber "int
RNA::GetStructureNumber()

Get the total number of specified or predicted structures.

Returns:
--------

An integer specify the total number of structures. ";

%feature("docstring")  RNA::DetermineDrawingCoordinates "int
RNA::DetermineDrawingCoordinates(const int height, const int width,
const int structurenumber=1)

Determine the coordinates for drawing a secondary structure.

This function determines drawing coordinates for all nucleotides in
structure number structurenumber. User must specify the height and
width of a character (use the largest of nucleotides). The coordinates
are in an abstract palette; the user must determine the minimum and
maximum coordinate in both the x and y direction. The actual
coordinates are fetched using GetNucleotideXCoordinate(int i) and
GetNucleotideYCoordinate(int i). This function returns are error code,
where 0=no error and other messages can be resolved to a c string
using GetErrorMessage. The structure to be drawn must be free of
pseudoknots.

Parameters:
-----------

height:  is an integer that refers to the height of a nucleotide.

width:  is an integer that refers to the width of a nucleotide, where
the largest nucleotide should be provided or a non-proportional font
should be used.

structurenumber:  is an integer that refers to the structure to be
drawn.

Returns:
--------

An int that provides an error code, 0 = no error and other errors can
be resolved to a c string using GetErrorMessage. ";

%feature("docstring")  RNA::GetCommentString "std::string
RNA::GetCommentString(const int structurenumber=1)

Provide the comment from the ct file as a string.

This function provides the comment from the CT file for a structure as
a string. This function generates internal error codes that can be
accessed by GetErrorCode() after the function is called: 0 = no error,
nonzero = error. The errorcode can be resolved to a c string using
GetErrorMessage.

Parameters:
-----------

structurenumber:  is the structure for which the comment is to be
provided.

Returns:
--------

A string that provides the comment. ";

%feature("docstring")  RNA::GetNucleotideXCoordinate "int
RNA::GetNucleotideXCoordinate(const int i)

Get the X coordinate for nucleotide i for drawing a structure.

This function gets the X coordinate for placing the nucleotide
specified by i. The user needs to have determined the coordinates for
a complete structure using DetermineDrawingCoordinates prior to making
this call. This function generates internal error codes that can be
accessed by GetErrorCode(): 0 = no error, nonzero = error. The
errorcode can be resolved to a c string using GetErrorMessage. Zero is
returned in case of error, but note that zero is also a valid
coordinate.

Parameters:
-----------

i:  is an integer refering to the nucleotide to be drawn.

Returns:
--------

An int that gives the X coordinate. ";

%feature("docstring")  RNA::GetNucleotideYCoordinate "int
RNA::GetNucleotideYCoordinate(const int i)

Get the Y coordinate for nucleotide i for drawing a structure.

This function gets the Y coordinate for placing the nucleotide
specified by i. The user needs to have determined the coordinates for
a complete structure using DetermineDrawingCoordinates prior to making
this call. This function generates internal error codes that can be
accessed by GetErrorCode(): 0 = no error, nonzero = error. The
errorcode can be resolved to a c string using GetErrorMessage. Zero is
returned in case of error, but note that zero is also a valid
coordinate.

Parameters:
-----------

i:  is an integer refering to the nucleotide to be drawn.

Returns:
--------

An int that gives the Y coordinate. ";

%feature("docstring")  RNA::GetLabelXCoordinate "int
RNA::GetLabelXCoordinate(const int i)

Get the X coordinate for placing the nucleotide index label specified
by i.

This function gets the X coordinate for placing the nucleotide index
label specified by i. The user needs to have determined the
coordinates for a complete structure using DetermineDrawingCoordinates
prior to making this call. This function generates internal error
codes that can be accessed by GetErrorCode(): 0 = no error, nonzero =
error. The errorcode can be resolved to a c string using
GetErrorMessage. Zero is returned in case of error, but note that zero
is also a valid coordinate. One additiona caveat: Labels that are
placed at 0,0 are lables that would have overlapped nucleotides. These
labels should not be drawn.

Parameters:
-----------

i:  is an integer refering to the label to be drawn. This needs to be
a multiple of 10.

Returns:
--------

An int that gives the X coordinate. ";

%feature("docstring")  RNA::GetLabelYCoordinate "int
RNA::GetLabelYCoordinate(const int i)

Get the Y coordinate for placing the nucleotide index label specified
by i.

This function gets the Y coordinate for placing the nucleotide index
label specified by i. The user needs to have determined the
coordinates for a complete structure using DetermineDrawingCoordinates
prior to making this call. This function generates internal error
codes that can be accessed by GetErrorCode(): 0 = no error, nonzero =
error. The errorcode can be resolved to a c string using
GetErrorMessage. Zero is returned in case of error, but note that zero
is also a valid coordinate. One additiona caveat: Labels that are
placed at 0,0 are lables that would have overlapped nucleotides. These
labels should not be drawn.

Parameters:
-----------

i:  is an integer refering to the label to be drawn. This needs to be
a multiple of 10.

Returns:
--------

An int that gives the Y coordinate. ";

%feature("docstring")  RNA::GetNucleotide "char
RNA::GetNucleotide(const int i)

param An integer specifying the nucleotide index (starting at 1 and
ending at GetSequenceLength()). This function generates internal error
codes that can be accessed by GetErrorCode(): 0 = no error, nonzero =
error. The errorcode can be resolved to a c string using
GetErrorMessage. Note that nucleotides are numbered starting at an
index of 1. return The char representing the nucleotide at index i or
'-' if an error occured. ";

%feature("docstring")  RNA::GetSequenceLength "int
RNA::GetSequenceLength()

Get the total length of the sequence.

Returns:
--------

An integer that specifies the total length of the sequence. ";

%feature("docstring")  RNA::GetBackboneType "bool
RNA::GetBackboneType()

Get the backbone type.

This function returns whether the backbone is RNA or DNA. Note that
backbone type is set when calling the constructor.

Returns:
--------

A bool that indicates the backbone (true = RNA, false = DNA). ";

%feature("docstring")  RNA::GetStructure "structure *
RNA::GetStructure()

Access the underlying structure class. This is provided for use with
two sequence methods. Generally, there is no need for end users to use
this function because the RNA class provides an convenient wrapper for
accessing the information in an RNA class.

Returns:
--------

A pointer to structure. ";

%feature("docstring")  RNA::SetProgress "void
RNA::SetProgress(TProgressDialog &Progress)

Provide a TProgressDialog for following calculation progress. A
TProgressDialog class has a public function void update(int percent)
that indicates the progress of a long calculation.

Parameters:
-----------

Progress:  is a TProgressDialog class. ";

%feature("docstring")  RNA::StopProgress "void RNA::StopProgress()

Provide a means to stop using a TProgressDialog. StopProgress tells
the RNA class to no longer follow progress. This should be called if
the TProgressDialog is deleted, so that this class does not make
reference to it. ";

%feature("docstring")  RNA::GetProgress "TProgressDialog *
RNA::GetProgress()

Return the current pointer to TProgressDialog. This is used during
inheritance to provide access to the underlying TProgressDialog. ";

%feature("docstring")  RNA::~RNA "RNA::~RNA()

Destructor.

The destructor automatically cleans up all allocated memory for
predicted or specified structures. ";


// File: class_thermodynamics.xml
%feature("docstring") Thermodynamics "

Thermodynamics Class.

The RNA class provides an encapsulation of the functions and struct
for reading and storing thermodynamic parameters. This includes
methods for changing folding temperatures This class is intended for
use in inheritance for classes that provide functionality.

C++ includes: thermodynamics.h ";

%feature("docstring")  Thermodynamics::Thermodynamics "Thermodynamics::Thermodynamics(const bool ISRNA=true) ";

%feature("docstring")  Thermodynamics::SetTemperature "int
Thermodynamics::SetTemperature(double temperature)

Set the temperature of folding in K.

This function allows the user to specify folding temperatures other
than 310.15 K (37 degrees C). This changes folding free energy changes
that would be returned for existing structures and would alter the set
of structures predicted. When this function is called, the
thermodynamic parameter files are immediately read from disk. These
include both enthalpy parameters (.dh files) and free energy changes
at 310.15 (.dat files). The files must either be at a location
indicated by the $DATAPATH environment variable or in pwd. Changing
the temperature only alters subsequent calculations. For example, if a
structure prediction method has been called, the set of predicted
structures are not changed at the time of a call to SetTemperature.
Likewise, SetTemperature must be called before calling a structure
prediction method if a temperature other than the 310.15 K default is
desired. The function returns an error code where 0 is no error and
non-zero errors can be parsed by by GetErrorMessage() or
GetErrorMessageString() in an inheriting class. ";

%feature("docstring")  Thermodynamics::GetTemperature "double
Thermodynamics::GetTemperature()

Get the current folding temperature in K.

Returns:
--------

A double that indicates the folding temperature in K. ";

%feature("docstring")  Thermodynamics::ReadThermodynamic "int
Thermodynamics::ReadThermodynamic(const char *pathname=NULL)

Function to read the thermodynamic parameters.

This function depends on temp, the current temperature, to determine
in the folding free energies need to be set to other than those read
in files for 310.15 K. Return of zero => no error and a return of non-
zero indicates error. Public functions that need the thermodynamic
parameters call this function automatically. By default, the path to
the thermodynamic paramaters is fetched from the $DATAPATH environment
variable. If a specific path is needed, $DATAPATH is overridden by
specifying the pathname explicitly here as a parameter.

Returns:
--------

An int that indicates whether an error occured.

Parameters:
-----------

pathname:  is a pointer to cstring that indicates the pathname to the
thermodynamnic parameters. By default, this is NULL and the
environment variable $DATAPATH is consulted to get this path. ";

%feature("docstring")  Thermodynamics::GetDatatable "datatable *
Thermodynamics::GetDatatable()

This function is used during inheritance o provide access to the free
energy change parameters. This function generates no error codes.
(Error checking was done for this during construction).

Returns:
--------

A pointer to datatable with free energy change parameters. ";

%feature("docstring")  Thermodynamics::GetEnthalpyTable "datatable *
Thermodynamics::GetEnthalpyTable()

This function is used to provide an enthalpy table. This function will
return a NULL pointer if there is an error reading the tables from
disk. It is important that programs check the status of the pointer
before using it, i.e. make sure it is not NULL.

Returns:
--------

A pointer to datatable with the enthalpy change parameters. ";

%feature("docstring")  Thermodynamics::CopyThermodynamic "void
Thermodynamics::CopyThermodynamic(Thermodynamics *thermo)

Copy thermodynamic parameters from an instance of Thermodynamics
class.

This is generally not needed because functions automatically populate
the parameters from disk. It is helpful, however, when a large number
of calculations with be performed because the parameters can then be
read from disk only once. Note that the source Thermodynamics class
must have been initialized with the \"correct\" ISRNA value and
correct temperature. Return 0 if no error and non-zero errors can be
parsed by GetErrorMessage() or GetErrorMessageString().

Parameters:
-----------

thermo:  is a pointer to Thermodynamics class. That must have already
called the ReadThermodynamics() function. ";

%feature("docstring")  Thermodynamics::GetEnergyRead "bool
Thermodynamics::GetEnergyRead()

Return whether this instance of Thermodynamics has the paremters
populated (either from disk or from another Thermodynamics class).

Returns:
--------

A bool yjay indicates whether the parameters are populated (true =
yes). ";

%feature("docstring")  Thermodynamics::~Thermodynamics "Thermodynamics::~Thermodynamics() ";


// File: class_two_r_n_a.xml
%feature("docstring") TwoRNA "

TwoRNA Class.

The TwoRNA class provides an entry point for all the two sequence
prediction routines of RNAstructure. This contains two instances of
the RNA class to provide the functionality of RNA.

C++ includes: TwoRNA.h ";

%feature("docstring")  TwoRNA::TwoRNA "TwoRNA::TwoRNA(const char
sequence1[], const char sequence2[], bool IsRNA=true)

Constructor - user provides a sequences as c strings.

Input sequences should contain A,C,G,T,U,a,c,g,t,u,x,X. Capitalization
makes no difference. T=t=u=U. If IsRNA is true, the backbone is RNA,
so U is assumed. If IsRNA is false, the backbone is DNA, so T is
assumed. x=X= nucleotide that neither stacks nor pairs. For now, any
unknown nuc is considered 'X'. Both sequences are passed to underlying
RNA classes for each sequence.

Parameters:
-----------

sequence1:  is a NULL terminated c string for sequence 1.

sequence2:  is a NULL terminated c string for sequence 2.

IsRNA:  is a bool that indicates whether these sequences are RNA or
DNA. true= RNA. false=DNA. Default is true. Both sequences must have
the same backbone. ";

%feature("docstring")  TwoRNA::TwoRNA "TwoRNA::TwoRNA(const char
filename1[], const int type1, const char filename2[], const int type2,
bool IsRNA=true)

Constructor - user provides a filenames for existing files as a c
string.

The existing files, specified by filenames, can either be a ct file, a
sequence, or an RNAstructure save file. Therefore, the user provides a
flag for the file type: type = 1 => .ct file, type = 2 => .seq file,
type = 3 => partition function save (.pfs) file, type = 4 => folding
save file (.sav). The file opening is performed by the constructors
for the RNA classes that underlie each sequence. This constructor
generates internal error codes that can be accessed by GetErrorCode()
after the constructor is called. 0 = no error. The errorcode can be
resolved to a c string using GetErrorMessage. Note that the contructor
needs to be explicitly told, via IsRNA, what the backbone is because
files do not store this information. Note also that save files
explicitly store the thermodynamic parameters, therefore changing the
backbone type as compaared to the original calculation will not change
structure predictions.

Parameters:
-----------

filename1:  is a null terminated c string and refers to sequence 1.

filename2:  is a null terminated c string and refers to sequence 2.

type1:  is an integer that indicates the file type for sequence 1.

type2:  is an integer that indicates the file type for sequence 2.

IsRNA:  is a bool that indicates whether these sequences are RNA or
DNA. true= RNA. false=DNA. Default is true. Only one backbone is
allowed for both sequences. ";

%feature("docstring")  TwoRNA::TwoRNA "TwoRNA::TwoRNA()

Constructor Default constructor that requires no parameters. ";

%feature("docstring")  TwoRNA::SetTemperature "int
TwoRNA::SetTemperature(double temperature)

Set the temperature at which the calculation will be performed in K.

This function allows the user to specify folding temperatures other
than 310.15 K (37 degrees C). This changes folding free energy changes
that would be returned for existing structures and would alter the set
of structures predicted. When this function is called, the
thermodynamic parameter files are immediately read from disk. These
include both enthalpy parameters (.dh files) and free energy changes
at 310.15 (.dat files). The files must either be at a location
indicated by the $DATAPATH environment variable or in pwd. Changing
the temperature only alters subsequent calculations. For example, if a
structure prediction method has been called, the set of predicted
structures are not changed at the time of a call to SetTemperature.
Likewise, SetTemperature must be called before calling a structure
prediction method if a temperature other than the 310.15 K default is
desired. The function returns an error code where 0 is no error and
non-zero errors can be parsed by by GetErrorMessage() or
GetErrorMessageString() in an inheriting class. ";

%feature("docstring")  TwoRNA::GetTemperature "double
TwoRNA::GetTemperature()

Get the current folding temperature in K.

Returns:
--------

A double that indicates the folding temperature in K. ";

%feature("docstring")  TwoRNA::GetErrorCode "int
TwoRNA::GetErrorCode()

Return an error code, where a return of zero is no error.

This function returns and error flag that is generated during
construction by RNA(const char &filename, const int type, const bool
IsRNA=true) or from CalculateFreeEnergy(). An error of zero is always
no error. Other codes are errors and a c-string can be fetched for the
error with GetErrorMessage().

Returns:
--------

An integer that provides the error code. ";

%feature("docstring")  TwoRNA::GetErrorMessage "char *
TwoRNA::GetErrorMessage(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

0 = no error 1000 = Error associated with sequence 1 or with a
procedure, function will get message from sequence 1 (the inherited
RNA class). 2000 = Error associated with sequence 2, function will get
message from sequence 2 (the RNA2 class). 3000 = Errors with each
sequence, function will get messages from each.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode().

Returns:
--------

A pointer to a c string that provides an error message or from other
functions that return integer error codes. ";

%feature("docstring")  TwoRNA::GetErrorMessageString "std::string
TwoRNA::GetErrorMessageString(const int error)

Return error messages based on code from GetErrorCode and other error
codes.

Although RNA generally uses c strings, this member function returns a
string that is suitable for interfacing with JAVA, etc. See the error
list in the GetErrorMessage() entry.

Parameters:
-----------

error:  is the integer error code provided by GetErrorCode() or from
other functions that return integer error codes.

Returns:
--------

A string that provides an error message. ";

%feature("docstring")  TwoRNA::GetRNA1 "RNA * TwoRNA::GetRNA1()

return A pointer to the underlying structure class for sequence 1.

Access the underlying RNA class. This is provided for use with two
sequence methods. Generally, there is no need for end users to use
this function.

Returns:
--------

A pointer to the underlying RNA class for sequence 1. ";

%feature("docstring")  TwoRNA::GetRNA2 "RNA * TwoRNA::GetRNA2()

return A pointer to the underlying structure class for sequence 2.

Access the underlying RNA class. This is provided for use with two
sequence methods. Generally, there is no need for end users to use
this function.

Returns:
--------

A pointer to the underlying RNA class for sequence 2. ";

%feature("docstring")  TwoRNA::~TwoRNA "TwoRNA::~TwoRNA()

Destructor.

The destructor cleans up all memory allocation. ";


// File: namespacestd.xml


// File: _dynalign__class_8cpp.xml
%feature("docstring")  std::main "int main(int argc, char *argv[]) ";


// File: _dynalign__object_8cpp.xml


// File: _dynalign__object_8h.xml


// File: _dynalign_s_m_p__class_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: _hybrid_r_n_a_8cpp.xml


// File: _hybrid_r_n_a_8h.xml


// File: _hybrid_r_n_a__class_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: _multifind__object_8cpp.xml


// File: _multifind__object_8h.xml


// File: _multilign__object_8cpp.xml


// File: _multilign__object_8h.xml


// File: _oligo_walk__class_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: _oligowalk__object_8cpp.xml


// File: _oligowalk__object_8h.xml


// File: _prob_scan_8cpp.xml
%feature("docstring")  basestack "basestack_t basestack(double p, int
i, int j, int k, int l) ";

%feature("docstring")  prev_val "PFPRECISION prev_val(int index, int
offset, vector< vector< PFPRECISION > > &arr) ";

%feature("docstring")  hairpin "hairpin_t hairpin(double p, int i,
int j) ";

%feature("docstring")  internal_loop "internal_loop_t
internal_loop(double p, int i, int j, int k, int l) ";

%feature("docstring")  multibranch_loop "multibranch_loop_t
multibranch_loop(int i, int j) ";

%feature("docstring")  add_branch "void add_branch(multibranch_loop_t
&mb, int k, int l) ";

%feature("docstring")  show_hairpins "void show_hairpins(vector<
hairpin_t > hairpins) ";

%feature("docstring")  show_internal_loops "void
show_internal_loops(vector< internal_loop_t > internals) ";

%feature("docstring")  show_bulge_loops "void
show_bulge_loops(vector< internal_loop_t > internals) ";

%feature("docstring")  show_stacks "void show_stacks(vector<
basestack_t > stacks) ";

%feature("docstring")  show_mb_element_array "void
show_mb_element_array(vector< mb_element > e) ";

%feature("docstring")  show_mbl "void show_mbl(multibranch_loop_t mb)
";


// File: _prob_scan_8h.xml
%feature("docstring")  hairpin "hairpin_t hairpin(double p, int i,
int j) ";

%feature("docstring")  internal_loop "internal_loop_t
internal_loop(double p, int i, int j, int k, int l) ";

%feature("docstring")  basestack "basestack_t basestack(double p, int
i, int j, int k, int l) ";

%feature("docstring")  multibranch_loop "multibranch_loop_t
multibranch_loop(int i, int j) ";

%feature("docstring")  add_branch "void add_branch(multibranch_loop_t
&mb, int k, int l) ";

%feature("docstring")  show_hairpins "void show_hairpins(vector<
hairpin_t >) ";

%feature("docstring")  show_stacks "void show_stacks(vector<
basestack_t >) ";

%feature("docstring")  show_internal_loops "void
show_internal_loops(vector< internal_loop_t >) ";

%feature("docstring")  show_bulge_loops "void
show_bulge_loops(vector< internal_loop_t >) ";

%feature("docstring")  show_mbl "void show_mbl(multibranch_loop_t
mbl) ";

%feature("docstring")  show_mb_element_array "void
show_mb_element_array(vector< mb_element >) ";


// File: _r_n_a_8cpp.xml
%feature("docstring")  errmsg "void errmsg(int err, int erri) ";


// File: _r_n_a_8h.xml


// File: _r_n_a__class_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: thermodynamics_8cpp.xml


// File: thermodynamics_8h.xml


// File: _two_r_n_a_8cpp.xml


// File: _two_r_n_a_8h.xml


// File: dir_75b82e7e4a5feb05200b9ad7adf06257.xml


// File: dir_b34ba2ae0e237b9beb7ee2f473e8ea85.xml


// File: dir_2fe3b2dcd0a0c532f0bbc8656c8806ff.xml


// File: dir_3730d91212ae4bedca2b67b4bc3dbf04.xml


// File: dir_eb9ec46fc31d8dcecc8867a4feab15e1.xml

