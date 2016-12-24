%module BSPTree

%{
#define SWIG 1
#include "BSPTree.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "BSPTree.h"
