%module FunnelTest
%{
#define SWIG 1
#include "FunnelTest.h"
#undef SWIG
%}

%include "LikelihoodFunction.i"
%include "FunnelTest.h"
