%module RNAstructureBackendCalculatorProxy;

%{
#include "RNAstructureBackendCalculator.h"
%}

%include "std_string.i"

//  ****************************************************************************************************
//  Include c++ enums as Java enums (instead of classes).
//  ****************************************************************************************************
	%include "enums.swg"

typedef int RNAInputType;
%include "RNAstructureBackendCalculator.h"