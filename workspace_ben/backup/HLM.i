%module HLM
%{
#define SWIG 1
#include "HLM.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "LikelihoodFunction.i"
%include "HLM.h"
