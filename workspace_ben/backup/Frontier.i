%module Frontier
%{
#define SWIG 1
#include "Frontier.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Frontier.h"
