%include "enums.swg"

%typemap(javain) enum SWIGTYPE "$javainput.ordinal()"
%typemap(javaout) enum SWIGTYPE {
    return $javaclassname.values()[$jnicall];
  }
%typemap(javabody) enum SWIGTYPE ""
