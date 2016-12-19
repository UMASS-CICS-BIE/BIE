# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_TessToolSender')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_TessToolSender')
    _TessToolSender = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TessToolSender', [dirname(__file__)])
        except ImportError:
            import _TessToolSender
            return _TessToolSender
        if fp is not None:
            try:
                _mod = imp.load_module('_TessToolSender', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TessToolSender = swig_import_helper()
    del swig_import_helper
else:
    import _TessToolSender
del _swig_python_version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class TessToolSender(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TessToolSender, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TessToolSender, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _TessToolSender.new_TessToolSender(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _TessToolSender.delete_TessToolSender
    __del__ = lambda self: None

    def Synchronize(self):
        return _TessToolSender.TessToolSender_Synchronize(self)

    def Detach(self):
        return _TessToolSender.TessToolSender_Detach(self)

    def SetSessionId(self, level, step):
        return _TessToolSender.TessToolSender_SetSessionId(self, level, step)

    def SetLikelihoodComputation(self, arg2):
        return _TessToolSender.TessToolSender_SetLikelihoodComputation(self, arg2)

    def GetMPIComm(self):
        return _TessToolSender.TessToolSender_GetMPIComm(self)

    def SetRecordType(self, rt):
        return _TessToolSender.TessToolSender_SetRecordType(self, rt)

    def SetMPIStream(self, mpistrm):
        return _TessToolSender.TessToolSender_SetMPIStream(self, mpistrm)

    def SetMPIFilter(self, filter):
        return _TessToolSender.TessToolSender_SetMPIFilter(self, filter)

    def SetMPISession(self, session):
        return _TessToolSender.TessToolSender_SetMPISession(self, session)

    def GetRecordType(self):
        return _TessToolSender.TessToolSender_GetRecordType(self)

    def GetMPIStream(self):
        return _TessToolSender.TessToolSender_GetMPIStream(self)

    def GetMPIFilter(self):
        return _TessToolSender.TessToolSender_GetMPIFilter(self)

    def GetMPISession(self):
        return _TessToolSender.TessToolSender_GetMPISession(self)

    def persistTessellation(self):
        return _TessToolSender.TessToolSender_persistTessellation(self)

    def attachFilter(self, *args):
        return _TessToolSender.TessToolSender_attachFilter(self, *args)

    def attachSelectionFilter(self, *args):
        return _TessToolSender.TessToolSender_attachSelectionFilter(self, *args)

    def printFilterChain(self):
        return _TessToolSender.TessToolSender_printFilterChain(self)

    def reinitializeFilterChain(self):
        return _TessToolSender.TessToolSender_reinitializeFilterChain(self)

    def attachMPIFilter(self):
        return _TessToolSender.TessToolSender_attachMPIFilter(self)

    def replaceFilter(self, fil, stackindex):
        return _TessToolSender.TessToolSender_replaceFilter(self, fil, stackindex)

    def attachDefaultFunctionSelectionFilter(self):
        return _TessToolSender.TessToolSender_attachDefaultFunctionSelectionFilter(self)
TessToolSender_swigregister = _TessToolSender.TessToolSender_swigregister
TessToolSender_swigregister(TessToolSender)

# This file is compatible with both classic and new-style classes.

