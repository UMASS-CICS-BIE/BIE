%module Distribution

%{
#define SWIG 1
#include "Distribution.h"	
#undef SWIG
%}


%include "../persistence/Serializable.h"
%include "Distribution.h"
