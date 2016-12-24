%module GaussTestMultiD
%{
#define SWIG 1
#include "GaussTestMultiD.h"
#undef SWIG
%}

%include "CDFGenerator.i"
%include "LikelihoodFunction.i"
%include "GaussTestMultiD.h"
