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
            fp, pathname, description = imp.find_module('_GalaxyModelTwoDCached', [dirname(__file__)])
        except ImportError:
            import _GalaxyModelTwoDCached
            return _GalaxyModelTwoDCached
        if fp is not None:
            try:
                _mod = imp.load_module('_GalaxyModelTwoDCached', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _GalaxyModelTwoDCached = swig_import_helper()
    del swig_import_helper
else:
    import _GalaxyModelTwoDCached
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


class GalaxyModelTwoDCached(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GalaxyModelTwoDCached, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GalaxyModelTwoDCached, name)
    __repr__ = _swig_repr

    def __init__(self, ndim, mdim, _dist):
        this = _GalaxyModelTwoDCached.new_GalaxyModelTwoDCached(ndim, mdim, _dist)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GalaxyModelTwoDCached.delete_GalaxyModelTwoDCached
    __del__ = lambda self: None

    def NormEval(self, x, y, d):
        return _GalaxyModelTwoDCached.GalaxyModelTwoDCached_NormEval(self, x, y, d)

    def Evaluate(self, x, y, d):
        return _GalaxyModelTwoDCached.GalaxyModelTwoDCached_Evaluate(self, x, y, d)
    __swig_setmethods__["numA"] = _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numA_set
    __swig_getmethods__["numA"] = _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numA_get
    if _newclass:
        numA = _swig_property(_GalaxyModelTwoDCached.GalaxyModelTwoDCached_numA_get, _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numA_set)
    __swig_setmethods__["numH"] = _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numH_set
    __swig_getmethods__["numH"] = _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numH_get
    if _newclass:
        numH = _swig_property(_GalaxyModelTwoDCached.GalaxyModelTwoDCached_numH_get, _GalaxyModelTwoDCached.GalaxyModelTwoDCached_numH_set)
GalaxyModelTwoDCached_swigregister = _GalaxyModelTwoDCached.GalaxyModelTwoDCached_swigregister
GalaxyModelTwoDCached_swigregister(GalaxyModelTwoDCached)
cvar = _GalaxyModelTwoDCached.cvar

# This file is compatible with both classic and new-style classes.

