%module TestMe
%{
#define SWIG 1
#include "TestMe.h"
#undef
%}
%include "TestMe.h"
