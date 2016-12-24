%module KdTessellation
%{
#define SWIG 1
#include "KdTessellation.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Tessellation.h"
%include "KdTessellation.h"
