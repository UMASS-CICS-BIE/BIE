%module EnsembleDisc
%{
#define SWIG 1
#include "EnsembleDisc.h"
#undef SWIG
%}

%include "Ensemble.i"
%include "EnsembleDisc.h"
