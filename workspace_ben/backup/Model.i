%module Model

%{
#define SWIG 1
#include "Model.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Model.h"
