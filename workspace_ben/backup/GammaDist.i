%module GammaDist
%{
#define SWIG 1
#include "GammaDist.h"
#undef SWIG
%}

%include "Distribution.i"
%include "GammaDist.h"
