%module Ensemble
%{
#define SWIG 1
#include "Ensemble.h"
#undef SWIG
%}
%include "../persistence/Serializable.h"
%include "Distribution.i"
%include "Ensemble.h"
