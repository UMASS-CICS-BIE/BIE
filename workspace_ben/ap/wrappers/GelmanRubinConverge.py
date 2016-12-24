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
            fp, pathname, description = imp.find_module('_GelmanRubinConverge', [dirname(__file__)])
        except ImportError:
            import _GelmanRubinConverge
            return _GelmanRubinConverge
        if fp is not None:
            try:
                _mod = imp.load_module('_GelmanRubinConverge', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _GelmanRubinConverge = swig_import_helper()
    del swig_import_helper
else:
    import _GelmanRubinConverge
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
        this = _GelmanRubinConverge.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GelmanRubinConverge.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _GelmanRubinConverge.Serializable_swigregister
Serializable_swigregister(Serializable)

class Distribution(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Distribution, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Distribution, name)
    __repr__ = _swig_repr
    UniformAdd = _GelmanRubinConverge.Distribution_UniformAdd
    UniformMult = _GelmanRubinConverge.Distribution_UniformMult
    NormalAdd = _GelmanRubinConverge.Distribution_NormalAdd
    NormalMult = _GelmanRubinConverge.Distribution_NormalMult
    Undefined = _GelmanRubinConverge.Distribution_Undefined

    def __init__(self):
        this = _GelmanRubinConverge.new_Distribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GelmanRubinConverge.delete_Distribution
    __del__ = lambda self: None

    def New(self):
        return _GelmanRubinConverge.Distribution_New(self)

    def PDF(self, arg2):
        return _GelmanRubinConverge.Distribution_PDF(self, arg2)

    def logPDF(self, arg2):
        return _GelmanRubinConverge.Distribution_logPDF(self, arg2)

    def CDF(self, arg2):
        return _GelmanRubinConverge.Distribution_CDF(self, arg2)

    def lower(self):
        return _GelmanRubinConverge.Distribution_lower(self)

    def upper(self):
        return _GelmanRubinConverge.Distribution_upper(self)

    def Mean(self):
        return _GelmanRubinConverge.Distribution_Mean(self)

    def StdDev(self):
        return _GelmanRubinConverge.Distribution_StdDev(self)

    def Moments(self, arg2):
        return _GelmanRubinConverge.Distribution_Moments(self, arg2)

    def Sample(self):
        return _GelmanRubinConverge.Distribution_Sample(self)

    def setWidth(self, x):
        return _GelmanRubinConverge.Distribution_setWidth(self, x)

    def Type(self):
        return _GelmanRubinConverge.Distribution_Type(self)

    def Dim(self):
        return _GelmanRubinConverge.Distribution_Dim(self)
Distribution_swigregister = _GelmanRubinConverge.Distribution_swigregister
Distribution_swigregister(Distribution)

class SampleDistribution(Distribution):
    __swig_setmethods__ = {}
    for _s in [Distribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, SampleDistribution, name, value)
    __swig_getmethods__ = {}
    for _s in [Distribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, SampleDistribution, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_destroy__ = _GelmanRubinConverge.delete_SampleDistribution
    __del__ = lambda self: None

    def CDF(self, arg2):
        return _GelmanRubinConverge.SampleDistribution_CDF(self, arg2)

    def New(self):
        return _GelmanRubinConverge.SampleDistribution_New(self)

    def AccumData(self, *args):
        return _GelmanRubinConverge.SampleDistribution_AccumData(self, *args)

    def ComputeDistribution(self):
        return _GelmanRubinConverge.SampleDistribution_ComputeDistribution(self)

    def numberData(self):
        return _GelmanRubinConverge.SampleDistribution_numberData(self)

    def getValue(self, i):
        return _GelmanRubinConverge.SampleDistribution_getValue(self, i)

    def getdim(self, i):
        return _GelmanRubinConverge.SampleDistribution_getdim(self, i)

    def getRecordType(self):
        return _GelmanRubinConverge.SampleDistribution_getRecordType(self)

    def getDataSetSize(self):
        return _GelmanRubinConverge.SampleDistribution_getDataSetSize(self)
SampleDistribution_swigregister = _GelmanRubinConverge.SampleDistribution_swigregister
SampleDistribution_swigregister(SampleDistribution)

class BinnedDistribution(SampleDistribution):
    __swig_setmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, BinnedDistribution, name, value)
    __swig_getmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, BinnedDistribution, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_destroy__ = _GelmanRubinConverge.delete_BinnedDistribution
    __del__ = lambda self: None

    def New(self):
        return _GelmanRubinConverge.BinnedDistribution_New(self)

    def getLow(self, i):
        return _GelmanRubinConverge.BinnedDistribution_getLow(self, i)

    def getHigh(self, i):
        return _GelmanRubinConverge.BinnedDistribution_getHigh(self, i)

    def CDF(self, arg2):
        return _GelmanRubinConverge.BinnedDistribution_CDF(self, arg2)
BinnedDistribution_swigregister = _GelmanRubinConverge.BinnedDistribution_swigregister
BinnedDistribution_swigregister(BinnedDistribution)

class PointDistribution(SampleDistribution):
    __swig_setmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, PointDistribution, name, value)
    __swig_getmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, PointDistribution, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _GelmanRubinConverge.new_PointDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GelmanRubinConverge.delete_PointDistribution
    __del__ = lambda self: None

    def New(self):
        return _GelmanRubinConverge.PointDistribution_New(self)

    def getRecordType(self):
        return _GelmanRubinConverge.PointDistribution_getRecordType(self)

    def AccumData(self, v, datapoint):
        return _GelmanRubinConverge.PointDistribution_AccumData(self, v, datapoint)

    def getdim(self, i):
        return _GelmanRubinConverge.PointDistribution_getdim(self, i)

    def numberData(self):
        return _GelmanRubinConverge.PointDistribution_numberData(self)

    def getValue(self, i):
        return _GelmanRubinConverge.PointDistribution_getValue(self, i)

    def Point(self):
        return _GelmanRubinConverge.PointDistribution_Point(self)

    def CDF(self, v):
        return _GelmanRubinConverge.PointDistribution_CDF(self, v)
PointDistribution_swigregister = _GelmanRubinConverge.PointDistribution_swigregister
PointDistribution_swigregister(PointDistribution)

class NullDistribution(SampleDistribution):
    __swig_setmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, NullDistribution, name, value)
    __swig_getmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, NullDistribution, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _GelmanRubinConverge.new_NullDistribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _GelmanRubinConverge.delete_NullDistribution
    __del__ = lambda self: None

    def New(self):
        return _GelmanRubinConverge.NullDistribution_New(self)
NullDistribution_swigregister = _GelmanRubinConverge.NullDistribution_swigregister
NullDistribution_swigregister(NullDistribution)

class Converge(SampleDistribution):
    __swig_setmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Converge, name, value)
    __swig_getmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Converge, name)
    __repr__ = _swig_repr
    __swig_setmethods__["count"] = _GelmanRubinConverge.Converge_count_set
    __swig_getmethods__["count"] = _GelmanRubinConverge.Converge_count_get
    if _newclass:
        count = _swig_property(_GelmanRubinConverge.Converge_count_get, _GelmanRubinConverge.Converge_count_set)
    __swig_setmethods__["nburn"] = _GelmanRubinConverge.Converge_nburn_set
    __swig_getmethods__["nburn"] = _GelmanRubinConverge.Converge_nburn_get
    if _newclass:
        nburn = _swig_property(_GelmanRubinConverge.Converge_nburn_get, _GelmanRubinConverge.Converge_nburn_set)

    def __init__(self, *args):
        this = _GelmanRubinConverge.new_Converge(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def Converged(self):
        return _GelmanRubinConverge.Converge_Converged(self)

    def New(self, m, d, id):
        return _GelmanRubinConverge.Converge_New(self, m, d, id)

    def AccumData(self, *args):
        return _GelmanRubinConverge.Converge_AccumData(self, *args)

    def GetLast(self, *args):
        return _GelmanRubinConverge.Converge_GetLast(self, *args)

    def IsParallel(self):
        return _GelmanRubinConverge.Converge_IsParallel(self)

    def ConvergedIndex(self):
        return _GelmanRubinConverge.Converge_ConvergedIndex(self)

    def AccumulateData(self, v, datapoint):
        return _GelmanRubinConverge.Converge_AccumulateData(self, v, datapoint)

    def DumpConverged(self, filename):
        return _GelmanRubinConverge.Converge_DumpConverged(self, filename)

    def PDF(self, x):
        return _GelmanRubinConverge.Converge_PDF(self, x)

    def logPDF(self, x):
        return _GelmanRubinConverge.Converge_logPDF(self, x)

    def lower(self):
        return _GelmanRubinConverge.Converge_lower(self)

    def upper(self):
        return _GelmanRubinConverge.Converge_upper(self)

    def Mean(self):
        return _GelmanRubinConverge.Converge_Mean(self)

    def StdDev(self):
        return _GelmanRubinConverge.Converge_StdDev(self)

    def Moments(self, m):
        return _GelmanRubinConverge.Converge_Moments(self, m)

    def Sample(self):
        return _GelmanRubinConverge.Converge_Sample(self)
    __swig_destroy__ = _GelmanRubinConverge.delete_Converge
    __del__ = lambda self: None
Converge_swigregister = _GelmanRubinConverge.Converge_swigregister
Converge_swigregister(Converge)


def inv_student_t_1sided(alpha, df):
    return _GelmanRubinConverge.inv_student_t_1sided(alpha, df)
inv_student_t_1sided = _GelmanRubinConverge.inv_student_t_1sided
class GelmanRubinConverge(Converge):
    __swig_setmethods__ = {}
    for _s in [Converge]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, GelmanRubinConverge, name, value)
    __swig_getmethods__ = {}
    for _s in [Converge]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, GelmanRubinConverge, name)
    __repr__ = _swig_repr
    __swig_setmethods__["ngood"] = _GelmanRubinConverge.GelmanRubinConverge_ngood_set
    __swig_getmethods__["ngood"] = _GelmanRubinConverge.GelmanRubinConverge_ngood_get
    if _newclass:
        ngood = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_ngood_get, _GelmanRubinConverge.GelmanRubinConverge_ngood_set)
    __swig_setmethods__["nskip"] = _GelmanRubinConverge.GelmanRubinConverge_nskip_set
    __swig_getmethods__["nskip"] = _GelmanRubinConverge.GelmanRubinConverge_nskip_get
    if _newclass:
        nskip = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_nskip_get, _GelmanRubinConverge.GelmanRubinConverge_nskip_set)
    __swig_setmethods__["alpha"] = _GelmanRubinConverge.GelmanRubinConverge_alpha_set
    __swig_getmethods__["alpha"] = _GelmanRubinConverge.GelmanRubinConverge_alpha_get
    if _newclass:
        alpha = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_alpha_get, _GelmanRubinConverge.GelmanRubinConverge_alpha_set)
    __swig_setmethods__["maxoutlier"] = _GelmanRubinConverge.GelmanRubinConverge_maxoutlier_set
    __swig_getmethods__["maxoutlier"] = _GelmanRubinConverge.GelmanRubinConverge_maxoutlier_get
    if _newclass:
        maxoutlier = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_maxoutlier_get, _GelmanRubinConverge.GelmanRubinConverge_maxoutlier_set)
    __swig_setmethods__["noutlier"] = _GelmanRubinConverge.GelmanRubinConverge_noutlier_set
    __swig_getmethods__["noutlier"] = _GelmanRubinConverge.GelmanRubinConverge_noutlier_get
    if _newclass:
        noutlier = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_noutlier_get, _GelmanRubinConverge.GelmanRubinConverge_noutlier_set)
    __swig_setmethods__["rtol"] = _GelmanRubinConverge.GelmanRubinConverge_rtol_set
    __swig_getmethods__["rtol"] = _GelmanRubinConverge.GelmanRubinConverge_rtol_get
    if _newclass:
        rtol = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_rtol_get, _GelmanRubinConverge.GelmanRubinConverge_rtol_set)
    __swig_setmethods__["verbose"] = _GelmanRubinConverge.GelmanRubinConverge_verbose_set
    __swig_getmethods__["verbose"] = _GelmanRubinConverge.GelmanRubinConverge_verbose_get
    if _newclass:
        verbose = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_verbose_get, _GelmanRubinConverge.GelmanRubinConverge_verbose_set)
    __swig_setmethods__["poffset"] = _GelmanRubinConverge.GelmanRubinConverge_poffset_set
    __swig_getmethods__["poffset"] = _GelmanRubinConverge.GelmanRubinConverge_poffset_get
    if _newclass:
        poffset = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_poffset_get, _GelmanRubinConverge.GelmanRubinConverge_poffset_set)
    __swig_setmethods__["debug"] = _GelmanRubinConverge.GelmanRubinConverge_debug_set
    __swig_getmethods__["debug"] = _GelmanRubinConverge.GelmanRubinConverge_debug_get
    if _newclass:
        debug = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_debug_get, _GelmanRubinConverge.GelmanRubinConverge_debug_set)
    __swig_setmethods__["maxkp"] = _GelmanRubinConverge.GelmanRubinConverge_maxkp_set
    __swig_getmethods__["maxkp"] = _GelmanRubinConverge.GelmanRubinConverge_maxkp_get
    if _newclass:
        maxkp = _swig_property(_GelmanRubinConverge.GelmanRubinConverge_maxkp_get, _GelmanRubinConverge.GelmanRubinConverge_maxkp_set)

    def __init__(self, m, d, id):
        this = _GelmanRubinConverge.new_GelmanRubinConverge(m, d, id)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def Converged(self):
        return _GelmanRubinConverge.GelmanRubinConverge_Converged(self)

    def IsParallel(self):
        return _GelmanRubinConverge.GelmanRubinConverge_IsParallel(self)

    def setMax(self, n):
        return _GelmanRubinConverge.GelmanRubinConverge_setMax(self, n)

    def setNskip(self, n):
        return _GelmanRubinConverge.GelmanRubinConverge_setNskip(self, n)

    def setNoutlier(self, n):
        return _GelmanRubinConverge.GelmanRubinConverge_setNoutlier(self, n)

    def setMaxout(self, n):
        return _GelmanRubinConverge.GelmanRubinConverge_setMaxout(self, n)

    def setNgood(self, n):
        return _GelmanRubinConverge.GelmanRubinConverge_setNgood(self, n)

    def Quiet(self):
        return _GelmanRubinConverge.GelmanRubinConverge_Quiet(self)

    def setPoffset(self, z):
        return _GelmanRubinConverge.GelmanRubinConverge_setPoffset(self, z)

    def ConvergedIndex(self):
        return _GelmanRubinConverge.GelmanRubinConverge_ConvergedIndex(self)

    def setRhatMax(self, r):
        return _GelmanRubinConverge.GelmanRubinConverge_setRhatMax(self, r)

    def setAlpha(self, a):
        return _GelmanRubinConverge.GelmanRubinConverge_setAlpha(self, a)
    __swig_getmethods__["setMaxK"] = lambda x: _GelmanRubinConverge.GelmanRubinConverge_setMaxK
    if _newclass:
        setMaxK = staticmethod(_GelmanRubinConverge.GelmanRubinConverge_setMaxK)

    def AccumData(self, values, states):
        return _GelmanRubinConverge.GelmanRubinConverge_AccumData(self, values, states)

    def GetLast(self, values, states):
        return _GelmanRubinConverge.GelmanRubinConverge_GetLast(self, values, states)

    def ComputeDistribution(self):
        return _GelmanRubinConverge.GelmanRubinConverge_ComputeDistribution(self)

    def New(self, *args):
        return _GelmanRubinConverge.GelmanRubinConverge_New(self, *args)
    __swig_destroy__ = _GelmanRubinConverge.delete_GelmanRubinConverge
    __del__ = lambda self: None
GelmanRubinConverge_swigregister = _GelmanRubinConverge.GelmanRubinConverge_swigregister
GelmanRubinConverge_swigregister(GelmanRubinConverge)
cvar = _GelmanRubinConverge.cvar

def GelmanRubinConverge_setMaxK(n):
    return _GelmanRubinConverge.GelmanRubinConverge_setMaxK(n)
GelmanRubinConverge_setMaxK = _GelmanRubinConverge.GelmanRubinConverge_setMaxK

# This file is compatible with both classic and new-style classes.

