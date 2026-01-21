 %module RNAstructure_wrap
 %{
 /* Includes the header in the wrapper code */
 #include "../RNA_class/thermodynamics.h"
 #include "../RNA_class/RNA.h"
 #include "../RNA_class/HybridRNA.h"
 #include "../RNA_class/TwoRNA.h"
 #include "../RNA_class/Dynalign_object.h"
 #include "../RNA_class/Multilign_object.h"
 #include "../RNA_class/Oligowalk_object.h"
 #include "../RNA_class/ProbScan.h"
 %}
 /*template instantiations to use STL vectors as python lists*/
 %include <std_string.i>
 %include <std_vector.i>
 %include <std_pair.i>

 %template(IntPair) std::pair<int,int>;
 %template(PairVector) std::vector<std::pair<int,int> >;
 %template(StringVector) std::vector<std::string>;
 %template(StringVectorVector) std::vector<std::vector<std::string> >;

 /* file containing docstrings for C++ methods */
 %include "docstrings.i"

 /* Parse the header file to generate wrappers */
 %include "../RNA_class/thermodynamics.h"
 %include "../RNA_class/RNA.h"
 %include "../RNA_class/HybridRNA.h"
 %include "../RNA_class/TwoRNA.h"
 %include "../RNA_class/Dynalign_object.h"
 %include "../RNA_class/Multilign_object.h"
 %include "../RNA_class/Oligowalk_object.h"
 %include "../RNA_class/ProbScan.h"

%template(HairpinLoopVector) std::vector<hairpin_t>;
%template(InternalLoopVector) std::vector<internal_loop_t>;
%template(MultibranchLoopVector) std::vector<multibranch_loop_t>;
%template(BaseStackVector) std::vector<basestack_t>;
 

/*%extend RNA{
        void reset(){
            $self->ErrorCode = 0;
    }
}*/
