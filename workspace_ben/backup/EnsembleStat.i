%module EnsembleStat
%{
#define SWIG 1
#include "EnsembleStat.h"
#undef SWIG
%}

%include "Ensemble.i"
%include "EnsembleStat.h"
