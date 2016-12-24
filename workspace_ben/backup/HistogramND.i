%module HistogramND
%{
#define SWIG 1
#include "HistogramND.h"
#undef SWIG
%}

%include "Distribution.i"
%include "HistogramND.h"
