%module Statement

%{
#define SWIG 1
#include "Statement.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Statement.h"
