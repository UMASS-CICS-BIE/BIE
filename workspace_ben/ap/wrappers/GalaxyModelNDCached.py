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
            fp, pathname, description = imp.find_module('_GalaxyModelNDCached', [dirname(__file__)])
        except ImportError:
            import _GalaxyModelNDCached
            return _GalaxyModelNDCached
        if fp is not None:
            try:
                _mod = imp.load_module('_GalaxyModelNDCached', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _GalaxyModelNDCached = swig_import_helper()
    del swig_import_helper
else:
    import _GalaxyModelNDCached
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


class Serializable(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Serializable, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Serializable, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _GalaxyModelNDCached.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GalaxyModelNDCached.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _GalaxyModelNDCached.Serializable_swigregister
Serializable_swigregister(Serializable)

class Model(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Model, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Model, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_setmethods__["quadtree"] = _GalaxyModelNDCached.Model_quadtree_set
    __swig_getmethods__["quadtree"] = _GalaxyModelNDCached.Model_quadtree_get
    if _newclass:
        quadtree = _swig_property(_GalaxyModelNDCached.Model_quadtree_get, _GalaxyModelNDCached.Model_quadtree_set)
    __swig_setmethods__["maxlevels"] = _GalaxyModelNDCached.Model_maxlevels_set
    __swig_getmethods__["maxlevels"] = _GalaxyModelNDCached.Model_maxlevels_get
    if _newclass:
        maxlevels = _swig_property(_GalaxyModelNDCached.Model_maxlevels_get, _GalaxyModelNDCached.Model_maxlevels_set)
    __swig_setmethods__["qeps"] = _GalaxyModelNDCached.Model_qeps_set
    __swig_getmethods__["qeps"] = _GalaxyModelNDCached.Model_qeps_get
    if _newclass:
        qeps = _swig_property(_GalaxyModelNDCached.Model_qeps_get, _GalaxyModelNDCached.Model_qeps_set)
    __swig_setmethods__["dX0"] = _GalaxyModelNDCached.Model_dX0_set
    __swig_getmethods__["dX0"] = _GalaxyModelNDCached.Model_dX0_get
    if _newclass:
        dX0 = _swig_property(_GalaxyModelNDCached.Model_dX0_get, _GalaxyModelNDCached.Model_dX0_set)
    __swig_setmethods__["dY0"] = _GalaxyModelNDCached.Model_dY0_set
    __swig_getmethods__["dY0"] = _GalaxyModelNDCached.Model_dY0_get
    if _newclass:
        dY0 = _swig_property(_GalaxyModelNDCached.Model_dY0_get, _GalaxyModelNDCached.Model_dY0_set)
    __swig_setmethods__["dX1"] = _GalaxyModelNDCached.Model_dX1_set
    __swig_getmethods__["dX1"] = _GalaxyModelNDCached.Model_dX1_get
    if _newclass:
        dX1 = _swig_property(_GalaxyModelNDCached.Model_dX1_get, _GalaxyModelNDCached.Model_dX1_set)
    __swig_setmethods__["dY1"] = _GalaxyModelNDCached.Model_dY1_set
    __swig_getmethods__["dY1"] = _GalaxyModelNDCached.Model_dY1_get
    if _newclass:
        dY1 = _swig_property(_GalaxyModelNDCached.Model_dY1_get, _GalaxyModelNDCached.Model_dY1_set)
    __swig_setmethods__["numX"] = _GalaxyModelNDCached.Model_numX_set
    __swig_getmethods__["numX"] = _GalaxyModelNDCached.Model_numX_get
    if _newclass:
        numX = _swig_property(_GalaxyModelNDCached.Model_numX_get, _GalaxyModelNDCached.Model_numX_set)
    __swig_setmethods__["numY"] = _GalaxyModelNDCached.Model_numY_set
    __swig_getmethods__["numY"] = _GalaxyModelNDCached.Model_numY_get
    if _newclass:
        numY = _swig_property(_GalaxyModelNDCached.Model_numY_get, _GalaxyModelNDCached.Model_numY_set)
    __swig_destroy__ = _GalaxyModelNDCached.delete_Model
    __del__ = lambda self: None

    def getParameterType(self):
        return _GalaxyModelNDCached.Model_getParameterType(self)

    def setNormKnots(self, nx, ny):
        return _GalaxyModelNDCached.Model_setNormKnots(self, nx, ny)

    def setQuadTreeParams(self, dx0, dy0, dx1, dy1, eps, mlev):
        return _GalaxyModelNDCached.Model_setQuadTreeParams(self, dx0, dy0, dx1, dy1, eps, mlev)

    def Initialize(self, arg2):
        return _GalaxyModelNDCached.Model_Initialize(self, arg2)

    def ResetCache(self):
        return _GalaxyModelNDCached.Model_ResetCache(self)

    def NormEval(self, *args):
        return _GalaxyModelNDCached.Model_NormEval(self, *args)

    def NormEvalMeasure(self, X, Y):
        return _GalaxyModelNDCached.Model_NormEvalMeasure(self, X, Y)

    def Evaluate(self, x, y, d):
        return _GalaxyModelNDCached.Model_Evaluate(self, x, y, d)

    def EvaluateBinned(self, x, y, d):
        return _GalaxyModelNDCached.Model_EvaluateBinned(self, x, y, d)

    def EvaluatePoint(self, x, y, d):
        return _GalaxyModelNDCached.Model_EvaluatePoint(self, x, y, d)

    def ParameterLabels(self):
        return _GalaxyModelNDCached.Model_ParameterLabels(self)

    def ExtendedLabels(self):
        return _GalaxyModelNDCached.Model_ExtendedLabels(self)
    binned = _GalaxyModelNDCached.Model_binned
    point = _GalaxyModelNDCached.Model_point
Model_swigregister = _GalaxyModelNDCached.Model_swigregister
Model_swigregister(Model)
cvar = _GalaxyModelNDCached.cvar

class GalaxyModelND(Model):
    __swig_setmethods__ = {}
    for _s in [Model]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, GalaxyModelND, name, value)
    __swig_getmethods__ = {}
    for _s in [Model]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, GalaxyModelND, name)
    __repr__ = _swig_repr

    def __init__(self, ndim, mdim, histo):
        this = _GalaxyModelNDCached.new_GalaxyModelND(ndim, mdim, histo)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GalaxyModelNDCached.delete_GalaxyModelND
    __del__ = lambda self: None

    def SetKnots(self, num):
        return _GalaxyModelNDCached.GalaxyModelND_SetKnots(self, num)

    def CacheLimit(self, num):
        return _GalaxyModelNDCached.GalaxyModelND_CacheLimit(self, num)

    def SetExtinction(self, A, Z):
        return _GalaxyModelNDCached.GalaxyModelND_SetExtinction(self, A, Z)

    def Initialize(self, *args):
        return _GalaxyModelNDCached.GalaxyModelND_Initialize(self, *args)

    def ResetCache(self):
        return _GalaxyModelNDCached.GalaxyModelND_ResetCache(self)

    def NormEval(self, *args):
        return _GalaxyModelNDCached.GalaxyModelND_NormEval(self, *args)

    def NormEvalMeasure(self, x, y):
        return _GalaxyModelNDCached.GalaxyModelND_NormEvalMeasure(self, x, y)

    def Evaluate(self, x, y, d):
        return _GalaxyModelNDCached.GalaxyModelND_Evaluate(self, x, y, d)

    def ParameterLabels(self):
        return _GalaxyModelNDCached.GalaxyModelND_ParameterLabels(self)
    __swig_setmethods__["LMAG"] = _GalaxyModelNDCached.GalaxyModelND_LMAG_set
    __swig_getmethods__["LMAG"] = _GalaxyModelNDCached.GalaxyModelND_LMAG_get
    if _newclass:
        LMAG = _swig_property(_GalaxyModelNDCached.GalaxyModelND_LMAG_get, _GalaxyModelNDCached.GalaxyModelND_LMAG_set)
    __swig_setmethods__["HMAG"] = _GalaxyModelNDCached.GalaxyModelND_HMAG_set
    __swig_getmethods__["HMAG"] = _GalaxyModelNDCached.GalaxyModelND_HMAG_get
    if _newclass:
        HMAG = _swig_property(_GalaxyModelNDCached.GalaxyModelND_HMAG_get, _GalaxyModelNDCached.GalaxyModelND_HMAG_set)
    __swig_setmethods__["A1"] = _GalaxyModelNDCached.GalaxyModelND_A1_set
    __swig_getmethods__["A1"] = _GalaxyModelNDCached.GalaxyModelND_A1_get
    if _newclass:
        A1 = _swig_property(_GalaxyModelNDCached.GalaxyModelND_A1_get, _GalaxyModelNDCached.GalaxyModelND_A1_set)
    __swig_setmethods__["Z1"] = _GalaxyModelNDCached.GalaxyModelND_Z1_set
    __swig_getmethods__["Z1"] = _GalaxyModelNDCached.GalaxyModelND_Z1_get
    if _newclass:
        Z1 = _swig_property(_GalaxyModelNDCached.GalaxyModelND_Z1_get, _GalaxyModelNDCached.GalaxyModelND_Z1_set)
    __swig_setmethods__["K0"] = _GalaxyModelNDCached.GalaxyModelND_K0_set
    __swig_getmethods__["K0"] = _GalaxyModelNDCached.GalaxyModelND_K0_get
    if _newclass:
        K0 = _swig_property(_GalaxyModelNDCached.GalaxyModelND_K0_get, _GalaxyModelNDCached.GalaxyModelND_K0_set)
    __swig_setmethods__["NUM"] = _GalaxyModelNDCached.GalaxyModelND_NUM_set
    __swig_getmethods__["NUM"] = _GalaxyModelNDCached.GalaxyModelND_NUM_get
    if _newclass:
        NUM = _swig_property(_GalaxyModelNDCached.GalaxyModelND_NUM_get, _GalaxyModelNDCached.GalaxyModelND_NUM_set)
    __swig_setmethods__["ALPHA"] = _GalaxyModelNDCached.GalaxyModelND_ALPHA_set
    __swig_getmethods__["ALPHA"] = _GalaxyModelNDCached.GalaxyModelND_ALPHA_get
    if _newclass:
        ALPHA = _swig_property(_GalaxyModelNDCached.GalaxyModelND_ALPHA_get, _GalaxyModelNDCached.GalaxyModelND_ALPHA_set)
    __swig_setmethods__["BETA"] = _GalaxyModelNDCached.GalaxyModelND_BETA_set
    __swig_getmethods__["BETA"] = _GalaxyModelNDCached.GalaxyModelND_BETA_get
    if _newclass:
        BETA = _swig_property(_GalaxyModelNDCached.GalaxyModelND_BETA_get, _GalaxyModelNDCached.GalaxyModelND_BETA_set)
    __swig_setmethods__["R0"] = _GalaxyModelNDCached.GalaxyModelND_R0_set
    __swig_getmethods__["R0"] = _GalaxyModelNDCached.GalaxyModelND_R0_get
    if _newclass:
        R0 = _swig_property(_GalaxyModelNDCached.GalaxyModelND_R0_get, _GalaxyModelNDCached.GalaxyModelND_R0_set)
    __swig_setmethods__["RMAX"] = _GalaxyModelNDCached.GalaxyModelND_RMAX_set
    __swig_getmethods__["RMAX"] = _GalaxyModelNDCached.GalaxyModelND_RMAX_get
    if _newclass:
        RMAX = _swig_property(_GalaxyModelNDCached.GalaxyModelND_RMAX_get, _GalaxyModelNDCached.GalaxyModelND_RMAX_set)
    __swig_setmethods__["AK"] = _GalaxyModelNDCached.GalaxyModelND_AK_set
    __swig_getmethods__["AK"] = _GalaxyModelNDCached.GalaxyModelND_AK_get
    if _newclass:
        AK = _swig_property(_GalaxyModelNDCached.GalaxyModelND_AK_get, _GalaxyModelNDCached.GalaxyModelND_AK_set)
    __swig_setmethods__["AMIN"] = _GalaxyModelNDCached.GalaxyModelND_AMIN_set
    __swig_getmethods__["AMIN"] = _GalaxyModelNDCached.GalaxyModelND_AMIN_get
    if _newclass:
        AMIN = _swig_property(_GalaxyModelNDCached.GalaxyModelND_AMIN_get, _GalaxyModelNDCached.GalaxyModelND_AMIN_set)
    __swig_setmethods__["AMAX"] = _GalaxyModelNDCached.GalaxyModelND_AMAX_set
    __swig_getmethods__["AMAX"] = _GalaxyModelNDCached.GalaxyModelND_AMAX_get
    if _newclass:
        AMAX = _swig_property(_GalaxyModelNDCached.GalaxyModelND_AMAX_get, _GalaxyModelNDCached.GalaxyModelND_AMAX_set)
    __swig_setmethods__["HMIN"] = _GalaxyModelNDCached.GalaxyModelND_HMIN_set
    __swig_getmethods__["HMIN"] = _GalaxyModelNDCached.GalaxyModelND_HMIN_get
    if _newclass:
        HMIN = _swig_property(_GalaxyModelNDCached.GalaxyModelND_HMIN_get, _GalaxyModelNDCached.GalaxyModelND_HMIN_set)
    __swig_setmethods__["HMAX"] = _GalaxyModelNDCached.GalaxyModelND_HMAX_set
    __swig_getmethods__["HMAX"] = _GalaxyModelNDCached.GalaxyModelND_HMAX_get
    if _newclass:
        HMAX = _swig_property(_GalaxyModelNDCached.GalaxyModelND_HMAX_get, _GalaxyModelNDCached.GalaxyModelND_HMAX_set)
    __swig_setmethods__["BASISDATA"] = _GalaxyModelNDCached.GalaxyModelND_BASISDATA_set
    __swig_getmethods__["BASISDATA"] = _GalaxyModelNDCached.GalaxyModelND_BASISDATA_get
    if _newclass:
        BASISDATA = _swig_property(_GalaxyModelNDCached.GalaxyModelND_BASISDATA_get, _GalaxyModelNDCached.GalaxyModelND_BASISDATA_set)
GalaxyModelND_swigregister = _GalaxyModelNDCached.GalaxyModelND_swigregister
GalaxyModelND_swigregister(GalaxyModelND)

class GalaxyModelNDCached(GalaxyModelND):
    __swig_setmethods__ = {}
    for _s in [GalaxyModelND]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, GalaxyModelNDCached, name, value)
    __swig_getmethods__ = {}
    for _s in [GalaxyModelND]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, GalaxyModelNDCached, name)
    __repr__ = _swig_repr

    def __init__(self, ndim, mdim, _dist):
        this = _GalaxyModelNDCached.new_GalaxyModelNDCached(ndim, mdim, _dist)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GalaxyModelNDCached.delete_GalaxyModelNDCached
    __del__ = lambda self: None

    def NormEval(self, x, y, d):
        return _GalaxyModelNDCached.GalaxyModelNDCached_NormEval(self, x, y, d)

    def EvaluateBinned(self, x, y, d):
        return _GalaxyModelNDCached.GalaxyModelNDCached_EvaluateBinned(self, x, y, d)
    __swig_setmethods__["numA"] = _GalaxyModelNDCached.GalaxyModelNDCached_numA_set
    __swig_getmethods__["numA"] = _GalaxyModelNDCached.GalaxyModelNDCached_numA_get
    if _newclass:
        numA = _swig_property(_GalaxyModelNDCached.GalaxyModelNDCached_numA_get, _GalaxyModelNDCached.GalaxyModelNDCached_numA_set)
    __swig_setmethods__["numH"] = _GalaxyModelNDCached.GalaxyModelNDCached_numH_set
    __swig_getmethods__["numH"] = _GalaxyModelNDCached.GalaxyModelNDCached_numH_get
    if _newclass:
        numH = _swig_property(_GalaxyModelNDCached.GalaxyModelNDCached_numH_get, _GalaxyModelNDCached.GalaxyModelNDCached_numH_set)
GalaxyModelNDCached_swigregister = _GalaxyModelNDCached.GalaxyModelNDCached_swigregister
GalaxyModelNDCached_swigregister(GalaxyModelNDCached)

# This file is compatible with both classic and new-style classes.


