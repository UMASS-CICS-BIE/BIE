%module Histogram1D
%{
#define SWIG 1
#include "Histogram1D.h"
#undef SWIG
%}

%include "../libsrnd/SmplHist.i"
%include "Distribution.i"
%include "Histogram1D.h"
