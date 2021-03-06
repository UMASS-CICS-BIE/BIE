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
        mname = '.'.join((pkg, '_SplatModel3dv')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_SplatModel3dv')
    _SplatModel3dv = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_SplatModel3dv', [dirname(__file__)])
        except ImportError:
            import _SplatModel3dv
            return _SplatModel3dv
        if fp is not None:
            try:
                _mod = imp.load_module('_SplatModel3dv', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _SplatModel3dv = swig_import_helper()
    del swig_import_helper
else:
    import _SplatModel3dv
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

class SplatModel3dv(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SplatModel3dv, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SplatModel3dv, name)
    __repr__ = _swig_repr
    __swig_setmethods__["SIGX"] = _SplatModel3dv.SplatModel3dv_SIGX_set
    __swig_getmethods__["SIGX"] = _SplatModel3dv.SplatModel3dv_SIGX_get
    if _newclass:
        SIGX = _swig_property(_SplatModel3dv.SplatModel3dv_SIGX_get, _SplatModel3dv.SplatModel3dv_SIGX_set)
    __swig_setmethods__["SIGY"] = _SplatModel3dv.SplatModel3dv_SIGY_set
    __swig_getmethods__["SIGY"] = _SplatModel3dv.SplatModel3dv_SIGY_get
    if _newclass:
        SIGY = _swig_property(_SplatModel3dv.SplatModel3dv_SIGY_get, _SplatModel3dv.SplatModel3dv_SIGY_set)
    __swig_setmethods__["V1"] = _SplatModel3dv.SplatModel3dv_V1_set
    __swig_getmethods__["V1"] = _SplatModel3dv.SplatModel3dv_V1_get
    if _newclass:
        V1 = _swig_property(_SplatModel3dv.SplatModel3dv_V1_get, _SplatModel3dv.SplatModel3dv_V1_set)
    __swig_setmethods__["V2"] = _SplatModel3dv.SplatModel3dv_V2_set
    __swig_getmethods__["V2"] = _SplatModel3dv.SplatModel3dv_V2_get
    if _newclass:
        V2 = _swig_property(_SplatModel3dv.SplatModel3dv_V2_get, _SplatModel3dv.SplatModel3dv_V2_set)
    __swig_setmethods__["V3"] = _SplatModel3dv.SplatModel3dv_V3_set
    __swig_getmethods__["V3"] = _SplatModel3dv.SplatModel3dv_V3_get
    if _newclass:
        V3 = _swig_property(_SplatModel3dv.SplatModel3dv_V3_get, _SplatModel3dv.SplatModel3dv_V3_set)

    def __init__(self, mdim, fluxmin, fluxmax):
        this = _SplatModel3dv.new_SplatModel3dv(mdim, fluxmin, fluxmax)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def Initialize(self, *args):
        return _SplatModel3dv.SplatModel3dv_Initialize(self, *args)

    def NormEval(self, *args):
        return _SplatModel3dv.SplatModel3dv_NormEval(self, *args)

    def EvaluateBinned(self, x, y, d):
        return _SplatModel3dv.SplatModel3dv_EvaluateBinned(self, x, y, d)

    def EvaluatePoint(self, x, y, d):
        return _SplatModel3dv.SplatModel3dv_EvaluatePoint(self, x, y, d)

    def ParameterLabels(self):
        return _SplatModel3dv.SplatModel3dv_ParameterLabels(self)
    __swig_destroy__ = _SplatModel3dv.delete_SplatModel3dv
    __del__ = lambda self: None
SplatModel3dv_swigregister = _SplatModel3dv.SplatModel3dv_swigregister
SplatModel3dv_swigregister(SplatModel3dv)
cvar = _SplatModel3dv.cvar

# This file is compatible with both classic and new-style classes.


