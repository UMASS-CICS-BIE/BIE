%module FivePointIntegration
%{
#define SWIG 1
#include "FivePointIntegration.h"
#undef SWIG
%}

%include "Integration.i"
%include "FivePointIntegration.h"
