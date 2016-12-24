%module CDFGenerator

%{
#define SWIG 1
#include "CDFGenerator.h"	
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "CDFGenerator.h"
