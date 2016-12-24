%module Integration
%{
#define SWIG 1
#include "Integration.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Integration.h"
