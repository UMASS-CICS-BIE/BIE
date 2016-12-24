%module KDTree
%{
#define SWIG 1
#include "KDTree.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "BSPTree.i"
%include "KDTree.h"
