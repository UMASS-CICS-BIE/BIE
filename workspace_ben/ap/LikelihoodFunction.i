%module LikelihoodFunction

%{
#define SWIG 1
#include "LikelihoodFunction.h"	
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "LikelihoodFunction.h"
