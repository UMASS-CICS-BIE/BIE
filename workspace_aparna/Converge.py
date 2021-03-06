# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Converge', [dirname(__file__)])
        except ImportError:
            import _Converge
            return _Converge
        if fp is not None:
            try:
                _mod = imp.load_module('_Converge', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Converge = swig_import_helper()
    del swig_import_helper
else:
    import _Converge
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


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


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0


class Converge(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Converge, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Converge, name)
    __repr__ = _swig_repr
    __swig_setmethods__["count"] = _Converge.Converge_count_set
    __swig_getmethods__["count"] = _Converge.Converge_count_get
    if _newclass:
        count = _swig_property(_Converge.Converge_count_get, _Converge.Converge_count_set)
    __swig_setmethods__["nburn"] = _Converge.Converge_nburn_set
    __swig_getmethods__["nburn"] = _Converge.Converge_nburn_get
    if _newclass:
        nburn = _swig_property(_Converge.Converge_nburn_get, _Converge.Converge_nburn_set)

    def __init__(self, *args):
        this = _Converge.new_Converge(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def Converged(self):
        return _Converge.Converge_Converged(self)

    def New(self, m, d, id):
        return _Converge.Converge_New(self, m, d, id)

    def AccumData(self, *args):
        return _Converge.Converge_AccumData(self, *args)

    def GetLast(self, *args):
        return _Converge.Converge_GetLast(self, *args)

    def IsParallel(self):
        return _Converge.Converge_IsParallel(self)

    def ConvergedIndex(self):
        return _Converge.Converge_ConvergedIndex(self)

    def AccumulateData(self, v, datapoint):
        return _Converge.Converge_AccumulateData(self, v, datapoint)

    def DumpConverged(self, filename):
        return _Converge.Converge_DumpConverged(self, filename)

    def PDF(self, x):
        return _Converge.Converge_PDF(self, x)

    def logPDF(self, x):
        return _Converge.Converge_logPDF(self, x)

    def lower(self):
        return _Converge.Converge_lower(self)

    def upper(self):
        return _Converge.Converge_upper(self)

    def Mean(self):
        return _Converge.Converge_Mean(self)

    def StdDev(self):
        return _Converge.Converge_StdDev(self)

    def Moments(self, m):
        return _Converge.Converge_Moments(self, m)

    def Sample(self):
        return _Converge.Converge_Sample(self)
    __swig_destroy__ = _Converge.delete_Converge
    __del__ = lambda self: None
Converge_swigregister = _Converge.Converge_swigregister
Converge_swigregister(Converge)

# This file is compatible with both classic and new-style classes.


