%module SmplStat

%{
#define SWIG 1
#include "SmplStat.h"	
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "SmplStat.h"