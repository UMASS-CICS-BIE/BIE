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
            fp, pathname, description = imp.find_module('_GalaxyModelOneDCached', [dirname(__file__)])
        except ImportError:
            import _GalaxyModelOneDCached
            return _GalaxyModelOneDCached
        if fp is not None:
            try:
                _mod = imp.load_module('_GalaxyModelOneDCached', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _GalaxyModelOneDCached = swig_import_helper()
    del swig_import_helper
else:
    import _GalaxyModelOneDCached
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


class GalaxyModelOneDCached(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GalaxyModelOneDCached, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GalaxyModelOneDCached, name)
    __repr__ = _swig_repr

    def __init__(self, ndim, mdim, _dist):
        this = _GalaxyModelOneDCached.new_GalaxyModelOneDCached(ndim, mdim, _dist)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GalaxyModelOneDCached.delete_GalaxyModelOneDCached
    __del__ = lambda self: None

    def NormEval(self, x, y, d):
        return _GalaxyModelOneDCached.GalaxyModelOneDCached_NormEval(self, x, y, d)

    def Evaluate(self, x, y, d):
        return _GalaxyModelOneDCached.GalaxyModelOneDCached_Evaluate(self, x, y, d)
    __swig_setmethods__["numA"] = _GalaxyModelOneDCached.GalaxyModelOneDCached_numA_set
    __swig_getmethods__["numA"] = _GalaxyModelOneDCached.GalaxyModelOneDCached_numA_get
    if _newclass:
        numA = _swig_property(_GalaxyModelOneDCached.GalaxyModelOneDCached_numA_get, _GalaxyModelOneDCached.GalaxyModelOneDCached_numA_set)
    __swig_setmethods__["numH"] = _GalaxyModelOneDCached.GalaxyModelOneDCached_numH_set
    __swig_getmethods__["numH"] = _GalaxyModelOneDCached.GalaxyModelOneDCached_numH_get
    if _newclass:
        numH = _swig_property(_GalaxyModelOneDCached.GalaxyModelOneDCached_numH_get, _GalaxyModelOneDCached.GalaxyModelOneDCached_numH_set)
GalaxyModelOneDCached_swigregister = _GalaxyModelOneDCached.GalaxyModelOneDCached_swigregister
GalaxyModelOneDCached_swigregister(GalaxyModelOneDCached)
cvar = _GalaxyModelOneDCached.cvar

# This file is compatible with both classic and new-style classes.


