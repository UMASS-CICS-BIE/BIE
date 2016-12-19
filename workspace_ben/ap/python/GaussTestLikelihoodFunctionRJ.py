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
            fp, pathname, description = imp.find_module('_GaussTestLikelihoodFunctionRJ', [dirname(__file__)])
        except ImportError:
            import _GaussTestLikelihoodFunctionRJ
            return _GaussTestLikelihoodFunctionRJ
        if fp is not None:
            try:
                _mod = imp.load_module('_GaussTestLikelihoodFunctionRJ', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _GaussTestLikelihoodFunctionRJ = swig_import_helper()
    del swig_import_helper
else:
    import _GaussTestLikelihoodFunctionRJ
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


class LikelihoodComputation(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, LikelihoodComputation, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, LikelihoodComputation, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _GaussTestLikelihoodFunctionRJ.new_LikelihoodComputation()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GaussTestLikelihoodFunctionRJ.delete_LikelihoodComputation
    __del__ = lambda self: None

    def Likelihood(self, state, indx):
        return _GaussTestLikelihoodFunctionRJ.LikelihoodComputation_Likelihood(self, state, indx)
LikelihoodComputation_swigregister = _GaussTestLikelihoodFunctionRJ.LikelihoodComputation_swigregister
LikelihoodComputation_swigregister(LikelihoodComputation)

class GaussTestLikelihoodFunction(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GaussTestLikelihoodFunction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GaussTestLikelihoodFunction, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nint"] = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_nint_set
    __swig_getmethods__["nint"] = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_nint_get
    if _newclass:
        nint = _swig_property(_GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_nint_get, _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_nint_set)

    def __init__(self, *args):
        this = _GaussTestLikelihoodFunctionRJ.new_GaussTestLikelihoodFunction(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GaussTestLikelihoodFunctionRJ.delete_GaussTestLikelihoodFunction
    __del__ = lambda self: None

    def useCauchyModel(self, b):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_useCauchyModel(self, b)

    def useCauchyData(self, b):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_useCauchyData(self, b)

    def useNcomp(self):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_useNcomp(self)

    def newModel(self, *args):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_newModel(self, *args)

    def SetDim(self, n):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_SetDim(self, n)

    def LikeProb(self, z, sd, norm, t, s, indx):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_LikeProb(self, z, sd, norm, t, s, indx)

    def CumuProb(self, z, sd, norm, t, s, indx, f):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_CumuProb(self, z, sd, norm, t, s, indx, f)

    def ParameterDescription(self, i):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_ParameterDescription(self, i)
GaussTestLikelihoodFunction_swigregister = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunction_swigregister
GaussTestLikelihoodFunction_swigregister(GaussTestLikelihoodFunction)
cvar = _GaussTestLikelihoodFunctionRJ.cvar

class GaussTestLikelihoodFunctionRJ(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GaussTestLikelihoodFunctionRJ, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GaussTestLikelihoodFunctionRJ, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nint"] = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_nint_set
    __swig_getmethods__["nint"] = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_nint_get
    if _newclass:
        nint = _swig_property(_GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_nint_get, _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_nint_set)

    def __init__(self, *args):
        this = _GaussTestLikelihoodFunctionRJ.new_GaussTestLikelihoodFunctionRJ(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GaussTestLikelihoodFunctionRJ.delete_GaussTestLikelihoodFunctionRJ
    __del__ = lambda self: None

    def useCauchyModel(self, b):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_useCauchyModel(self, b)

    def useCauchyData(self, b):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_useCauchyData(self, b)

    def newModel(self, file):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_newModel(self, file)

    def SetDim(self, n):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_SetDim(self, n)

    def LikeProb(self, z, sd, norm, t, s, indx):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_LikeProb(self, z, sd, norm, t, s, indx)

    def CumuProb(self, z, sd, norm, t, s, indx, f):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_CumuProb(self, z, sd, norm, t, s, indx, f)

    def ParameterDescription(self, i):
        return _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_ParameterDescription(self, i)
GaussTestLikelihoodFunctionRJ_swigregister = _GaussTestLikelihoodFunctionRJ.GaussTestLikelihoodFunctionRJ_swigregister
GaussTestLikelihoodFunctionRJ_swigregister(GaussTestLikelihoodFunctionRJ)

# This file is compatible with both classic and new-style classes.

