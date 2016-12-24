%module MixturePrior

%{
#define SWIG 1
#include "MixturePrior.h"	
#undef SWIG
%}
%include "Prior.i"
%include "MixturePrior.h"
