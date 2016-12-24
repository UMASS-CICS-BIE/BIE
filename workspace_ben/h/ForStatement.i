%module ForStatement
%{
#define SWIG 1
#include "ForStatement.h"
#undef SWIG
%}

%include "../persistence/Serializable.h"
%include "Statement.i"
%include "ForStatement.h"
