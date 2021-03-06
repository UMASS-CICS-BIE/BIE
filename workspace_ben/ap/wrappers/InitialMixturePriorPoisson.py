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
            fp, pathname, description = imp.find_module('_InitialMixturePriorPoisson', [dirname(__file__)])
        except ImportError:
            import _InitialMixturePriorPoisson
            return _InitialMixturePriorPoisson
        if fp is not None:
            try:
                _mod = imp.load_module('_InitialMixturePriorPoisson', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _InitialMixturePriorPoisson = swig_import_helper()
    del swig_import_helper
else:
    import _InitialMixturePriorPoisson
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
        this = _InitialMixturePriorPoisson.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _InitialMixturePriorPoisson.Serializable_swigregister
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
    UniformAdd = _InitialMixturePriorPoisson.Distribution_UniformAdd
    UniformMult = _InitialMixturePriorPoisson.Distribution_UniformMult
    NormalAdd = _InitialMixturePriorPoisson.Distribution_NormalAdd
    NormalMult = _InitialMixturePriorPoisson.Distribution_NormalMult
    Undefined = _InitialMixturePriorPoisson.Distribution_Undefined

    def __init__(self):
        this = _InitialMixturePriorPoisson.new_Distribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_Distribution
    __del__ = lambda self: None

    def New(self):
        return _InitialMixturePriorPoisson.Distribution_New(self)

    def PDF(self, arg2):
        return _InitialMixturePriorPoisson.Distribution_PDF(self, arg2)

    def logPDF(self, arg2):
        return _InitialMixturePriorPoisson.Distribution_logPDF(self, arg2)

    def CDF(self, arg2):
        return _InitialMixturePriorPoisson.Distribution_CDF(self, arg2)

    def lower(self):
        return _InitialMixturePriorPoisson.Distribution_lower(self)

    def upper(self):
        return _InitialMixturePriorPoisson.Distribution_upper(self)

    def Mean(self):
        return _InitialMixturePriorPoisson.Distribution_Mean(self)

    def StdDev(self):
        return _InitialMixturePriorPoisson.Distribution_StdDev(self)

    def Moments(self, arg2):
        return _InitialMixturePriorPoisson.Distribution_Moments(self, arg2)

    def Sample(self):
        return _InitialMixturePriorPoisson.Distribution_Sample(self)

    def setWidth(self, x):
        return _InitialMixturePriorPoisson.Distribution_setWidth(self, x)

    def Type(self):
        return _InitialMixturePriorPoisson.Distribution_Type(self)

    def Dim(self):
        return _InitialMixturePriorPoisson.Distribution_Dim(self)
Distribution_swigregister = _InitialMixturePriorPoisson.Distribution_swigregister
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
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_SampleDistribution
    __del__ = lambda self: None

    def CDF(self, arg2):
        return _InitialMixturePriorPoisson.SampleDistribution_CDF(self, arg2)

    def New(self):
        return _InitialMixturePriorPoisson.SampleDistribution_New(self)

    def AccumData(self, *args):
        return _InitialMixturePriorPoisson.SampleDistribution_AccumData(self, *args)

    def ComputeDistribution(self):
        return _InitialMixturePriorPoisson.SampleDistribution_ComputeDistribution(self)

    def numberData(self):
        return _InitialMixturePriorPoisson.SampleDistribution_numberData(self)

    def getValue(self, i):
        return _InitialMixturePriorPoisson.SampleDistribution_getValue(self, i)

    def getdim(self, i):
        return _InitialMixturePriorPoisson.SampleDistribution_getdim(self, i)

    def getRecordType(self):
        return _InitialMixturePriorPoisson.SampleDistribution_getRecordType(self)

    def getDataSetSize(self):
        return _InitialMixturePriorPoisson.SampleDistribution_getDataSetSize(self)
SampleDistribution_swigregister = _InitialMixturePriorPoisson.SampleDistribution_swigregister
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
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_BinnedDistribution
    __del__ = lambda self: None

    def New(self):
        return _InitialMixturePriorPoisson.BinnedDistribution_New(self)

    def getLow(self, i):
        return _InitialMixturePriorPoisson.BinnedDistribution_getLow(self, i)

    def getHigh(self, i):
        return _InitialMixturePriorPoisson.BinnedDistribution_getHigh(self, i)

    def CDF(self, arg2):
        return _InitialMixturePriorPoisson.BinnedDistribution_CDF(self, arg2)
BinnedDistribution_swigregister = _InitialMixturePriorPoisson.BinnedDistribution_swigregister
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
        this = _InitialMixturePriorPoisson.new_PointDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_PointDistribution
    __del__ = lambda self: None

    def New(self):
        return _InitialMixturePriorPoisson.PointDistribution_New(self)

    def getRecordType(self):
        return _InitialMixturePriorPoisson.PointDistribution_getRecordType(self)

    def AccumData(self, v, datapoint):
        return _InitialMixturePriorPoisson.PointDistribution_AccumData(self, v, datapoint)

    def getdim(self, i):
        return _InitialMixturePriorPoisson.PointDistribution_getdim(self, i)

    def numberData(self):
        return _InitialMixturePriorPoisson.PointDistribution_numberData(self)

    def getValue(self, i):
        return _InitialMixturePriorPoisson.PointDistribution_getValue(self, i)

    def Point(self):
        return _InitialMixturePriorPoisson.PointDistribution_Point(self)

    def CDF(self, v):
        return _InitialMixturePriorPoisson.PointDistribution_CDF(self, v)
PointDistribution_swigregister = _InitialMixturePriorPoisson.PointDistribution_swigregister
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
        this = _InitialMixturePriorPoisson.new_NullDistribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_NullDistribution
    __del__ = lambda self: None

    def New(self):
        return _InitialMixturePriorPoisson.NullDistribution_New(self)
NullDistribution_swigregister = _InitialMixturePriorPoisson.NullDistribution_swigregister
NullDistribution_swigregister(NullDistribution)

class Prior(Distribution):
    __swig_setmethods__ = {}
    for _s in [Distribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Prior, name, value)
    __swig_getmethods__ = {}
    for _s in [Distribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Prior, name)
    __repr__ = _swig_repr
    __swig_setmethods__["width_factor"] = _InitialMixturePriorPoisson.Prior_width_factor_set
    __swig_getmethods__["width_factor"] = _InitialMixturePriorPoisson.Prior_width_factor_get
    if _newclass:
        width_factor = _swig_property(_InitialMixturePriorPoisson.Prior_width_factor_get, _InitialMixturePriorPoisson.Prior_width_factor_set)
    __swig_setmethods__["max_iter"] = _InitialMixturePriorPoisson.Prior_max_iter_set
    __swig_getmethods__["max_iter"] = _InitialMixturePriorPoisson.Prior_max_iter_get
    if _newclass:
        max_iter = _swig_property(_InitialMixturePriorPoisson.Prior_max_iter_get, _InitialMixturePriorPoisson.Prior_max_iter_set)

    def __init__(self, *args):
        this = _InitialMixturePriorPoisson.new_Prior(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_Prior
    __del__ = lambda self: None

    def InitialStates(self, e):
        return _InitialMixturePriorPoisson.Prior_InitialStates(self, e)

    def SI(self):
        return _InitialMixturePriorPoisson.Prior_SI(self)

    def Dim(self):
        return _InitialMixturePriorPoisson.Prior_Dim(self)

    def New(self):
        return _InitialMixturePriorPoisson.Prior_New(self)

    def PDF(self, arg2):
        return _InitialMixturePriorPoisson.Prior_PDF(self, arg2)

    def logPDF(self, arg2):
        return _InitialMixturePriorPoisson.Prior_logPDF(self, arg2)

    def lower(self):
        return _InitialMixturePriorPoisson.Prior_lower(self)

    def upper(self):
        return _InitialMixturePriorPoisson.Prior_upper(self)

    def Mean(self):
        return _InitialMixturePriorPoisson.Prior_Mean(self)

    def StdDev(self):
        return _InitialMixturePriorPoisson.Prior_StdDev(self)

    def Moments(self, arg2):
        return _InitialMixturePriorPoisson.Prior_Moments(self, arg2)

    def Sample(self, *args):
        return _InitialMixturePriorPoisson.Prior_Sample(self, *args)

    def SamplePrior(self, arg2):
        return _InitialMixturePriorPoisson.Prior_SamplePrior(self, arg2)

    def SampleProposal(self, ch, width):
        return _InitialMixturePriorPoisson.Prior_SampleProposal(self, ch, width)

    def EnforceBounds(self, ch):
        return _InitialMixturePriorPoisson.Prior_EnforceBounds(self, ch)

    def EnforceBoundsState(self, s):
        return _InitialMixturePriorPoisson.Prior_EnforceBoundsState(self, s)
Prior_swigregister = _InitialMixturePriorPoisson.Prior_swigregister
Prior_swigregister(Prior)
cvar = _InitialMixturePriorPoisson.cvar

class MixturePrior(Prior):
    __swig_setmethods__ = {}
    for _s in [Prior]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, MixturePrior, name, value)
    __swig_getmethods__ = {}
    for _s in [Prior]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, MixturePrior, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_setmethods__["ITMAX"] = _InitialMixturePriorPoisson.MixturePrior_ITMAX_set
    __swig_getmethods__["ITMAX"] = _InitialMixturePriorPoisson.MixturePrior_ITMAX_get
    if _newclass:
        ITMAX = _swig_property(_InitialMixturePriorPoisson.MixturePrior_ITMAX_get, _InitialMixturePriorPoisson.MixturePrior_ITMAX_set)
    __swig_setmethods__["min_prob"] = _InitialMixturePriorPoisson.MixturePrior_min_prob_set
    __swig_getmethods__["min_prob"] = _InitialMixturePriorPoisson.MixturePrior_min_prob_get
    if _newclass:
        min_prob = _swig_property(_InitialMixturePriorPoisson.MixturePrior_min_prob_get, _InitialMixturePriorPoisson.MixturePrior_min_prob_set)
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_MixturePrior
    __del__ = lambda self: None

    def setMaxIter(self, i):
        return _InitialMixturePriorPoisson.MixturePrior_setMaxIter(self, i)

    def New(self):
        return _InitialMixturePriorPoisson.MixturePrior_New(self)

    def PDF(self, arg2):
        return _InitialMixturePriorPoisson.MixturePrior_PDF(self, arg2)

    def logPDF(self, arg2):
        return _InitialMixturePriorPoisson.MixturePrior_logPDF(self, arg2)

    def logPDFMarginal(self, M, n, V):
        return _InitialMixturePriorPoisson.MixturePrior_logPDFMarginal(self, M, n, V)

    def lower(self):
        return _InitialMixturePriorPoisson.MixturePrior_lower(self)

    def upper(self):
        return _InitialMixturePriorPoisson.MixturePrior_upper(self)

    def SampleM(self):
        return _InitialMixturePriorPoisson.MixturePrior_SampleM(self)

    def Sample(self, *args):
        return _InitialMixturePriorPoisson.MixturePrior_Sample(self, *args)

    def SamplePrior(self, *args):
        return _InitialMixturePriorPoisson.MixturePrior_SamplePrior(self, *args)

    def SampleProposal(self, ch, width):
        return _InitialMixturePriorPoisson.MixturePrior_SampleProposal(self, ch, width)

    def BirthWeight(self, M):
        return _InitialMixturePriorPoisson.MixturePrior_BirthWeight(self, M)

    def BirthWeightPDF(self, M, p):
        return _InitialMixturePriorPoisson.MixturePrior_BirthWeightPDF(self, M, p)

    def Exists(self, m):
        return _InitialMixturePriorPoisson.MixturePrior_Exists(self, m)
MixturePrior_swigregister = _InitialMixturePriorPoisson.MixturePrior_swigregister
MixturePrior_swigregister(MixturePrior)

class InitialMixturePrior(MixturePrior):
    __swig_setmethods__ = {}
    for _s in [MixturePrior]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, InitialMixturePrior, name, value)
    __swig_getmethods__ = {}
    for _s in [MixturePrior]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, InitialMixturePrior, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _InitialMixturePriorPoisson.new_InitialMixturePrior(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def New(self):
        return _InitialMixturePriorPoisson.InitialMixturePrior_New(self)

    def PDF(self, arg2):
        return _InitialMixturePriorPoisson.InitialMixturePrior_PDF(self, arg2)

    def logPDF(self, arg2):
        return _InitialMixturePriorPoisson.InitialMixturePrior_logPDF(self, arg2)

    def logPDFMarginal(self, M, n, V):
        return _InitialMixturePriorPoisson.InitialMixturePrior_logPDFMarginal(self, M, n, V)

    def lower(self):
        return _InitialMixturePriorPoisson.InitialMixturePrior_lower(self)

    def upper(self):
        return _InitialMixturePriorPoisson.InitialMixturePrior_upper(self)

    def Mean(self, m):
        return _InitialMixturePriorPoisson.InitialMixturePrior_Mean(self, m)

    def StdDev(self, m):
        return _InitialMixturePriorPoisson.InitialMixturePrior_StdDev(self, m)

    def Moments(self, m, n):
        return _InitialMixturePriorPoisson.InitialMixturePrior_Moments(self, m, n)

    def Sample(self, m):
        return _InitialMixturePriorPoisson.InitialMixturePrior_Sample(self, m)

    def SamplePrior(self, *args):
        return _InitialMixturePriorPoisson.InitialMixturePrior_SamplePrior(self, *args)
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_InitialMixturePrior
    __del__ = lambda self: None
InitialMixturePrior_swigregister = _InitialMixturePriorPoisson.InitialMixturePrior_swigregister
InitialMixturePrior_swigregister(InitialMixturePrior)

class InitialMixturePriorPoisson(InitialMixturePrior):
    __swig_setmethods__ = {}
    for _s in [InitialMixturePrior]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, InitialMixturePriorPoisson, name, value)
    __swig_getmethods__ = {}
    for _s in [InitialMixturePrior]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, InitialMixturePriorPoisson, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _InitialMixturePriorPoisson.new_InitialMixturePriorPoisson(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _InitialMixturePriorPoisson.delete_InitialMixturePriorPoisson
    __del__ = lambda self: None

    def New(self):
        return _InitialMixturePriorPoisson.InitialMixturePriorPoisson_New(self)

    def SamplePrior(self, *args):
        return _InitialMixturePriorPoisson.InitialMixturePriorPoisson_SamplePrior(self, *args)
InitialMixturePriorPoisson_swigregister = _InitialMixturePriorPoisson.InitialMixturePriorPoisson_swigregister
InitialMixturePriorPoisson_swigregister(InitialMixturePriorPoisson)

# This file is compatible with both classic and new-style classes.


