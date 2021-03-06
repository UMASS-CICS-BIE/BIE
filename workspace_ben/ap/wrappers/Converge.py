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


class Serializable(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Serializable, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Serializable, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _Converge.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Converge.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _Converge.Serializable_swigregister
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
    UniformAdd = _Converge.Distribution_UniformAdd
    UniformMult = _Converge.Distribution_UniformMult
    NormalAdd = _Converge.Distribution_NormalAdd
    NormalMult = _Converge.Distribution_NormalMult
    Undefined = _Converge.Distribution_Undefined

    def __init__(self):
        this = _Converge.new_Distribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Converge.delete_Distribution
    __del__ = lambda self: None

    def New(self):
        return _Converge.Distribution_New(self)

    def PDF(self, arg2):
        return _Converge.Distribution_PDF(self, arg2)

    def logPDF(self, arg2):
        return _Converge.Distribution_logPDF(self, arg2)

    def CDF(self, arg2):
        return _Converge.Distribution_CDF(self, arg2)

    def lower(self):
        return _Converge.Distribution_lower(self)

    def upper(self):
        return _Converge.Distribution_upper(self)

    def Mean(self):
        return _Converge.Distribution_Mean(self)

    def StdDev(self):
        return _Converge.Distribution_StdDev(self)

    def Moments(self, arg2):
        return _Converge.Distribution_Moments(self, arg2)

    def Sample(self):
        return _Converge.Distribution_Sample(self)

    def setWidth(self, x):
        return _Converge.Distribution_setWidth(self, x)

    def Type(self):
        return _Converge.Distribution_Type(self)

    def Dim(self):
        return _Converge.Distribution_Dim(self)
Distribution_swigregister = _Converge.Distribution_swigregister
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
    __swig_destroy__ = _Converge.delete_SampleDistribution
    __del__ = lambda self: None

    def CDF(self, arg2):
        return _Converge.SampleDistribution_CDF(self, arg2)

    def New(self):
        return _Converge.SampleDistribution_New(self)

    def AccumData(self, *args):
        return _Converge.SampleDistribution_AccumData(self, *args)

    def ComputeDistribution(self):
        return _Converge.SampleDistribution_ComputeDistribution(self)

    def numberData(self):
        return _Converge.SampleDistribution_numberData(self)

    def getValue(self, i):
        return _Converge.SampleDistribution_getValue(self, i)

    def getdim(self, i):
        return _Converge.SampleDistribution_getdim(self, i)

    def getRecordType(self):
        return _Converge.SampleDistribution_getRecordType(self)

    def getDataSetSize(self):
        return _Converge.SampleDistribution_getDataSetSize(self)
SampleDistribution_swigregister = _Converge.SampleDistribution_swigregister
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
    __swig_destroy__ = _Converge.delete_BinnedDistribution
    __del__ = lambda self: None

    def New(self):
        return _Converge.BinnedDistribution_New(self)

    def getLow(self, i):
        return _Converge.BinnedDistribution_getLow(self, i)

    def getHigh(self, i):
        return _Converge.BinnedDistribution_getHigh(self, i)

    def CDF(self, arg2):
        return _Converge.BinnedDistribution_CDF(self, arg2)
BinnedDistribution_swigregister = _Converge.BinnedDistribution_swigregister
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
        this = _Converge.new_PointDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Converge.delete_PointDistribution
    __del__ = lambda self: None

    def New(self):
        return _Converge.PointDistribution_New(self)

    def getRecordType(self):
        return _Converge.PointDistribution_getRecordType(self)

    def AccumData(self, v, datapoint):
        return _Converge.PointDistribution_AccumData(self, v, datapoint)

    def getdim(self, i):
        return _Converge.PointDistribution_getdim(self, i)

    def numberData(self):
        return _Converge.PointDistribution_numberData(self)

    def getValue(self, i):
        return _Converge.PointDistribution_getValue(self, i)

    def Point(self):
        return _Converge.PointDistribution_Point(self)

    def CDF(self, v):
        return _Converge.PointDistribution_CDF(self, v)
PointDistribution_swigregister = _Converge.PointDistribution_swigregister
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
        this = _Converge.new_NullDistribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Converge.delete_NullDistribution
    __del__ = lambda self: None

    def New(self):
        return _Converge.NullDistribution_New(self)
NullDistribution_swigregister = _Converge.NullDistribution_swigregister
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


