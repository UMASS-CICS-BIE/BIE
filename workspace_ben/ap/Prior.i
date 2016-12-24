%module Prior

%{
#define SWIG 1
#include "Prior.h"	
#undef SWIG
%}

%include "Distribution.i"
%include "Prior.h"
