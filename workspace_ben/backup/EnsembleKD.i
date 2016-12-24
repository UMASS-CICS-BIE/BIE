%module EnsembleKD
%{
#define SWIG 1
#include "EnsembleKD.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Ensemble.i"
%include "EnsembleKD.h"
