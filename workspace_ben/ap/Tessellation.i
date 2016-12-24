%module Tesselation

%{
#define SWIG 1
#include "Tessellation.h"	
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Tessellation.h"
