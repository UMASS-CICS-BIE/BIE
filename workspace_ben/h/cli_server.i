%module cli_server
%{
#define SWIG 1
#include ost
#include "cli_server.h"
#undef SWIG
%}
%include "cli_server.h"
