%module cliDistribution
%{
#define SWIG 1
#include "cliDistribution.h"
#undef
%}

%include "../persistence/Serializable.h"
%include "Distribution.h"
%include "cliDistribution.h"
