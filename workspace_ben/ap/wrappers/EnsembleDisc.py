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
            fp, pathname, description = imp.find_module('_EnsembleDisc', [dirname(__file__)])
        except ImportError:
            import _EnsembleDisc
            return _EnsembleDisc
        if fp is not None:
            try:
                _mod = imp.load_module('_EnsembleDisc', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _EnsembleDisc = swig_import_helper()
    del swig_import_helper
else:
    import _EnsembleDisc
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
        this = _EnsembleDisc.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _EnsembleDisc.Serializable_swigregister
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
    UniformAdd = _EnsembleDisc.Distribution_UniformAdd
    UniformMult = _EnsembleDisc.Distribution_UniformMult
    NormalAdd = _EnsembleDisc.Distribution_NormalAdd
    NormalMult = _EnsembleDisc.Distribution_NormalMult
    Undefined = _EnsembleDisc.Distribution_Undefined

    def __init__(self):
        this = _EnsembleDisc.new_Distribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_Distribution
    __del__ = lambda self: None

    def New(self):
        return _EnsembleDisc.Distribution_New(self)

    def PDF(self, arg2):
        return _EnsembleDisc.Distribution_PDF(self, arg2)

    def logPDF(self, arg2):
        return _EnsembleDisc.Distribution_logPDF(self, arg2)

    def CDF(self, arg2):
        return _EnsembleDisc.Distribution_CDF(self, arg2)

    def lower(self):
        return _EnsembleDisc.Distribution_lower(self)

    def upper(self):
        return _EnsembleDisc.Distribution_upper(self)

    def Mean(self):
        return _EnsembleDisc.Distribution_Mean(self)

    def StdDev(self):
        return _EnsembleDisc.Distribution_StdDev(self)

    def Moments(self, arg2):
        return _EnsembleDisc.Distribution_Moments(self, arg2)

    def Sample(self):
        return _EnsembleDisc.Distribution_Sample(self)

    def setWidth(self, x):
        return _EnsembleDisc.Distribution_setWidth(self, x)

    def Type(self):
        return _EnsembleDisc.Distribution_Type(self)

    def Dim(self):
        return _EnsembleDisc.Distribution_Dim(self)
Distribution_swigregister = _EnsembleDisc.Distribution_swigregister
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
    __swig_destroy__ = _EnsembleDisc.delete_SampleDistribution
    __del__ = lambda self: None

    def CDF(self, arg2):
        return _EnsembleDisc.SampleDistribution_CDF(self, arg2)

    def New(self):
        return _EnsembleDisc.SampleDistribution_New(self)

    def AccumData(self, *args):
        return _EnsembleDisc.SampleDistribution_AccumData(self, *args)

    def ComputeDistribution(self):
        return _EnsembleDisc.SampleDistribution_ComputeDistribution(self)

    def numberData(self):
        return _EnsembleDisc.SampleDistribution_numberData(self)

    def getValue(self, i):
        return _EnsembleDisc.SampleDistribution_getValue(self, i)

    def getdim(self, i):
        return _EnsembleDisc.SampleDistribution_getdim(self, i)

    def getRecordType(self):
        return _EnsembleDisc.SampleDistribution_getRecordType(self)

    def getDataSetSize(self):
        return _EnsembleDisc.SampleDistribution_getDataSetSize(self)
SampleDistribution_swigregister = _EnsembleDisc.SampleDistribution_swigregister
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
    __swig_destroy__ = _EnsembleDisc.delete_BinnedDistribution
    __del__ = lambda self: None

    def New(self):
        return _EnsembleDisc.BinnedDistribution_New(self)

    def getLow(self, i):
        return _EnsembleDisc.BinnedDistribution_getLow(self, i)

    def getHigh(self, i):
        return _EnsembleDisc.BinnedDistribution_getHigh(self, i)

    def CDF(self, arg2):
        return _EnsembleDisc.BinnedDistribution_CDF(self, arg2)
BinnedDistribution_swigregister = _EnsembleDisc.BinnedDistribution_swigregister
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
        this = _EnsembleDisc.new_PointDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_PointDistribution
    __del__ = lambda self: None

    def New(self):
        return _EnsembleDisc.PointDistribution_New(self)

    def getRecordType(self):
        return _EnsembleDisc.PointDistribution_getRecordType(self)

    def AccumData(self, v, datapoint):
        return _EnsembleDisc.PointDistribution_AccumData(self, v, datapoint)

    def getdim(self, i):
        return _EnsembleDisc.PointDistribution_getdim(self, i)

    def numberData(self):
        return _EnsembleDisc.PointDistribution_numberData(self)

    def getValue(self, i):
        return _EnsembleDisc.PointDistribution_getValue(self, i)

    def Point(self):
        return _EnsembleDisc.PointDistribution_Point(self)

    def CDF(self, v):
        return _EnsembleDisc.PointDistribution_CDF(self, v)
PointDistribution_swigregister = _EnsembleDisc.PointDistribution_swigregister
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
        this = _EnsembleDisc.new_NullDistribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_NullDistribution
    __del__ = lambda self: None

    def New(self):
        return _EnsembleDisc.NullDistribution_New(self)
NullDistribution_swigregister = _EnsembleDisc.NullDistribution_swigregister
NullDistribution_swigregister(NullDistribution)

class StateData(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, StateData, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, StateData, name)
    __repr__ = _swig_repr
    __swig_setmethods__["prob"] = _EnsembleDisc.StateData_prob_set
    __swig_getmethods__["prob"] = _EnsembleDisc.StateData_prob_get
    if _newclass:
        prob = _swig_property(_EnsembleDisc.StateData_prob_get, _EnsembleDisc.StateData_prob_set)
    __swig_setmethods__["p"] = _EnsembleDisc.StateData_p_set
    __swig_getmethods__["p"] = _EnsembleDisc.StateData_p_get
    if _newclass:
        p = _swig_property(_EnsembleDisc.StateData_p_get, _EnsembleDisc.StateData_p_set)

    def __init__(self, *args):
        this = _EnsembleDisc.new_StateData(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def Broadcast(self):
        return _EnsembleDisc.StateData_Broadcast(self)
    __swig_destroy__ = _EnsembleDisc.delete_StateData
    __del__ = lambda self: None
StateData_swigregister = _EnsembleDisc.StateData_swigregister
StateData_swigregister(StateData)

class PairIndex(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, PairIndex, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, PairIndex, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _EnsembleDisc.new_PairIndex()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_setmethods__["index"] = _EnsembleDisc.PairIndex_index_set
    __swig_getmethods__["index"] = _EnsembleDisc.PairIndex_index_get
    if _newclass:
        index = _swig_property(_EnsembleDisc.PairIndex_index_get, _EnsembleDisc.PairIndex_index_set)
    __swig_setmethods__["value"] = _EnsembleDisc.PairIndex_value_set
    __swig_getmethods__["value"] = _EnsembleDisc.PairIndex_value_get
    if _newclass:
        value = _swig_property(_EnsembleDisc.PairIndex_value_get, _EnsembleDisc.PairIndex_value_set)

    def __eq__(self, t):
        return _EnsembleDisc.PairIndex___eq__(self, t)
    __swig_destroy__ = _EnsembleDisc.delete_PairIndex
    __del__ = lambda self: None
PairIndex_swigregister = _EnsembleDisc.PairIndex_swigregister
PairIndex_swigregister(PairIndex)

class StateCache(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, StateCache, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, StateCache, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _EnsembleDisc.new_StateCache(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_StateCache
    __del__ = lambda self: None

    def Mean(self, m):
        return _EnsembleDisc.StateCache_Mean(self, m)

    def Sample(self, m):
        return _EnsembleDisc.StateCache_Sample(self, m)

    def logPDF(self, P):
        return _EnsembleDisc.StateCache_logPDF(self, P)

    def size(self):
        return _EnsembleDisc.StateCache_size(self)

    def stateSize(self):
        return _EnsembleDisc.StateCache_stateSize(self)
StateCache_swigregister = _EnsembleDisc.StateCache_swigregister
StateCache_swigregister(StateCache)

class Ensemble(SampleDistribution):
    __swig_setmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Ensemble, name, value)
    __swig_getmethods__ = {}
    for _s in [SampleDistribution]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Ensemble, name)
    __repr__ = _swig_repr
    __swig_setmethods__["verbose_debug"] = _EnsembleDisc.Ensemble_verbose_debug_set
    __swig_getmethods__["verbose_debug"] = _EnsembleDisc.Ensemble_verbose_debug_get
    if _newclass:
        verbose_debug = _swig_property(_EnsembleDisc.Ensemble_verbose_debug_get, _EnsembleDisc.Ensemble_verbose_debug_set)
    __swig_setmethods__["verbose_default"] = _EnsembleDisc.Ensemble_verbose_default_set
    __swig_getmethods__["verbose_default"] = _EnsembleDisc.Ensemble_verbose_default_get
    if _newclass:
        verbose_default = _swig_property(_EnsembleDisc.Ensemble_verbose_default_get, _EnsembleDisc.Ensemble_verbose_default_set)
    __swig_setmethods__["keep"] = _EnsembleDisc.Ensemble_keep_set
    __swig_getmethods__["keep"] = _EnsembleDisc.Ensemble_keep_get
    if _newclass:
        keep = _swig_property(_EnsembleDisc.Ensemble_keep_get, _EnsembleDisc.Ensemble_keep_set)
    __swig_setmethods__["minsub"] = _EnsembleDisc.Ensemble_minsub_set
    __swig_getmethods__["minsub"] = _EnsembleDisc.Ensemble_minsub_get
    if _newclass:
        minsub = _swig_property(_EnsembleDisc.Ensemble_minsub_get, _EnsembleDisc.Ensemble_minsub_set)
    __swig_setmethods__["thresh"] = _EnsembleDisc.Ensemble_thresh_set
    __swig_getmethods__["thresh"] = _EnsembleDisc.Ensemble_thresh_get
    if _newclass:
        thresh = _swig_property(_EnsembleDisc.Ensemble_thresh_get, _EnsembleDisc.Ensemble_thresh_set)
    __swig_setmethods__["continuous"] = _EnsembleDisc.Ensemble_continuous_set
    __swig_getmethods__["continuous"] = _EnsembleDisc.Ensemble_continuous_get
    if _newclass:
        continuous = _swig_property(_EnsembleDisc.Ensemble_continuous_get, _EnsembleDisc.Ensemble_continuous_set)
    __swig_setmethods__["key_pos"] = _EnsembleDisc.Ensemble_key_pos_set
    __swig_getmethods__["key_pos"] = _EnsembleDisc.Ensemble_key_pos_get
    if _newclass:
        key_pos = _swig_property(_EnsembleDisc.Ensemble_key_pos_get, _EnsembleDisc.Ensemble_key_pos_set)

    def __init__(self, *args):
        this = _EnsembleDisc.new_Ensemble(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_Ensemble
    __del__ = lambda self: None

    def Reset(self, *args):
        return _EnsembleDisc.Ensemble_Reset(self, *args)

    def setDimensions(self, si):
        return _EnsembleDisc.Ensemble_setDimensions(self, si)

    def setContinuous(self):
        return _EnsembleDisc.Ensemble_setContinuous(self)

    def setMaxRange(self):
        return _EnsembleDisc.Ensemble_setMaxRange(self)

    def setNKeep(self, n):
        return _EnsembleDisc.Ensemble_setNKeep(self, n)

    def Order(self, n):
        return _EnsembleDisc.Ensemble_Order(self, n)

    def setVerboseOn(self):
        return _EnsembleDisc.Ensemble_setVerboseOn(self)

    def setVerboseOff(self):
        return _EnsembleDisc.Ensemble_setVerboseOff(self)

    def setVerboseDebugOn(self):
        return _EnsembleDisc.Ensemble_setVerboseDebugOn(self)

    def setVerboseDebugOff(self):
        return _EnsembleDisc.Ensemble_setVerboseDebugOff(self)

    def AccumData(self, *args):
        return _EnsembleDisc.Ensemble_AccumData(self, *args)

    def AccumulateData(self, v, datapoint):
        return _EnsembleDisc.Ensemble_AccumulateData(self, v, datapoint)

    def ComputeDistribution(self, *args):
        return _EnsembleDisc.Ensemble_ComputeDistribution(self, *args)

    def stats(self, n, m, s):
        return _EnsembleDisc.Ensemble_stats(self, n, m, s)

    def Nstates(self):
        return _EnsembleDisc.Ensemble_Nstates(self)

    def Npopped(self):
        return _EnsembleDisc.Ensemble_Npopped(self)

    def PDFSubspace(self, m):
        return _EnsembleDisc.Ensemble_PDFSubspace(self, m)

    def getValidSubspaceMixtures(self):
        return _EnsembleDisc.Ensemble_getValidSubspaceMixtures(self)

    def getStates(self):
        return _EnsembleDisc.Ensemble_getStates(self)

    def getStateInfo(self):
        return _EnsembleDisc.Ensemble_getStateInfo(self)

    def logPDF(self, p):
        return _EnsembleDisc.Ensemble_logPDF(self, p)

    def logPDFMarginal(self, m, n, V):
        return _EnsembleDisc.Ensemble_logPDFMarginal(self, m, n, V)

    def New(self):
        return _EnsembleDisc.Ensemble_New(self)

    def StdDev(self, *args):
        return _EnsembleDisc.Ensemble_StdDev(self, *args)

    def Mean(self, *args):
        return _EnsembleDisc.Ensemble_Mean(self, *args)

    def Sample(self, *args):
        return _EnsembleDisc.Ensemble_Sample(self, *args)

    def sampleSubspace(self):
        return _EnsembleDisc.Ensemble_sampleSubspace(self)

    def Width(self, m):
        return _EnsembleDisc.Ensemble_Width(self, m)

    def WidthVector(self, m):
        return _EnsembleDisc.Ensemble_WidthVector(self, m)

    def printWidth(self, *args):
        return _EnsembleDisc.Ensemble_printWidth(self, *args)

    def StdDevMarginal(self, m, n):
        return _EnsembleDisc.Ensemble_StdDevMarginal(self, m, n)

    def MeanMarginal(self, m, n):
        return _EnsembleDisc.Ensemble_MeanMarginal(self, m, n)

    def SampleMarginal(self, m, n):
        return _EnsembleDisc.Ensemble_SampleMarginal(self, m, n)

    def WidthMarginal(self, m, n):
        return _EnsembleDisc.Ensemble_WidthMarginal(self, m, n)

    def PrintDiag(self, *args):
        return _EnsembleDisc.Ensemble_PrintDiag(self, *args)

    def Exists(self, m):
        return _EnsembleDisc.Ensemble_Exists(self, m)

    def Broadcast(self):
        return _EnsembleDisc.Ensemble_Broadcast(self)

    def enableDump(self):
        return _EnsembleDisc.Ensemble_enableDump(self)

    def disableDump(self):
        return _EnsembleDisc.Ensemble_disableDump(self)

    def binaryDump(self):
        return _EnsembleDisc.Ensemble_binaryDump(self)

    def asciiDump(self):
        return _EnsembleDisc.Ensemble_asciiDump(self)

    def statesToFile(self, filename, binary):
        return _EnsembleDisc.Ensemble_statesToFile(self, filename, binary)

    def dumpStates(self, *args):
        return _EnsembleDisc.Ensemble_dumpStates(self, *args)

    def PrintDensity(self, *args):
        return _EnsembleDisc.Ensemble_PrintDensity(self, *args)

    def getMaxLikeState(self):
        return _EnsembleDisc.Ensemble_getMaxLikeState(self)

    def getMaxProbState(self):
        return _EnsembleDisc.Ensemble_getMaxProbState(self)
Ensemble_swigregister = _EnsembleDisc.Ensemble_swigregister
Ensemble_swigregister(Ensemble)
cvar = _EnsembleDisc.cvar

class EnsembleDisc(Ensemble):
    __swig_setmethods__ = {}
    for _s in [Ensemble]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, EnsembleDisc, name, value)
    __swig_getmethods__ = {}
    for _s in [Ensemble]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, EnsembleDisc, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nnear"] = _EnsembleDisc.EnsembleDisc_nnear_set
    __swig_getmethods__["nnear"] = _EnsembleDisc.EnsembleDisc_nnear_get
    if _newclass:
        nnear = _swig_property(_EnsembleDisc.EnsembleDisc_nnear_get, _EnsembleDisc.EnsembleDisc_nnear_set)
    __swig_setmethods__["errortol"] = _EnsembleDisc.EnsembleDisc_errortol_set
    __swig_getmethods__["errortol"] = _EnsembleDisc.EnsembleDisc_errortol_get
    if _newclass:
        errortol = _swig_property(_EnsembleDisc.EnsembleDisc_errortol_get, _EnsembleDisc.EnsembleDisc_errortol_set)
    __swig_setmethods__["mindist"] = _EnsembleDisc.EnsembleDisc_mindist_set
    __swig_getmethods__["mindist"] = _EnsembleDisc.EnsembleDisc_mindist_get
    if _newclass:
        mindist = _swig_property(_EnsembleDisc.EnsembleDisc_mindist_get, _EnsembleDisc.EnsembleDisc_mindist_set)
    __swig_setmethods__["minwidth"] = _EnsembleDisc.EnsembleDisc_minwidth_set
    __swig_getmethods__["minwidth"] = _EnsembleDisc.EnsembleDisc_minwidth_get
    if _newclass:
        minwidth = _swig_property(_EnsembleDisc.EnsembleDisc_minwidth_get, _EnsembleDisc.EnsembleDisc_minwidth_set)
    __swig_setmethods__["scaled"] = _EnsembleDisc.EnsembleDisc_scaled_set
    __swig_getmethods__["scaled"] = _EnsembleDisc.EnsembleDisc_scaled_get
    if _newclass:
        scaled = _swig_property(_EnsembleDisc.EnsembleDisc_scaled_get, _EnsembleDisc.EnsembleDisc_scaled_set)
    __swig_setmethods__["mscale"] = _EnsembleDisc.EnsembleDisc_mscale_set
    __swig_getmethods__["mscale"] = _EnsembleDisc.EnsembleDisc_mscale_get
    if _newclass:
        mscale = _swig_property(_EnsembleDisc.EnsembleDisc_mscale_get, _EnsembleDisc.EnsembleDisc_mscale_set)

    def __init__(self, *args):
        this = _EnsembleDisc.new_EnsembleDisc(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleDisc.delete_EnsembleDisc
    __del__ = lambda self: None

    def Reset(self, *args):
        return _EnsembleDisc.EnsembleDisc_Reset(self, *args)

    def setDimensions(self, si):
        return _EnsembleDisc.EnsembleDisc_setDimensions(self, si)

    def setTarget(self, n):
        return _EnsembleDisc.EnsembleDisc_setTarget(self, n)

    def setBucket(self, n):
        return _EnsembleDisc.EnsembleDisc_setBucket(self, n)

    def setNnear(self, n):
        return _EnsembleDisc.EnsembleDisc_setNnear(self, n)

    def setRelative(self, b):
        return _EnsembleDisc.EnsembleDisc_setRelative(self, b)

    def setDensityScale(self, f):
        return _EnsembleDisc.EnsembleDisc_setDensityScale(self, f)

    def setKernelScale(self, s):
        return _EnsembleDisc.EnsembleDisc_setKernelScale(self, s)

    def ComputeDistribution(self, *args):
        return _EnsembleDisc.EnsembleDisc_ComputeDistribution(self, *args)

    def ComputeDistributionMarginal(self, out):
        return _EnsembleDisc.EnsembleDisc_ComputeDistributionMarginal(self, out)

    def preComputeDensity(self, *args):
        return _EnsembleDisc.EnsembleDisc_preComputeDensity(self, *args)

    def New(self):
        return _EnsembleDisc.EnsembleDisc_New(self)

    def PDF(self, arg2):
        return _EnsembleDisc.EnsembleDisc_PDF(self, arg2)

    def logPDF(self, arg2):
        return _EnsembleDisc.EnsembleDisc_logPDF(self, arg2)

    def logPDFMarginal(self, m, n, V):
        return _EnsembleDisc.EnsembleDisc_logPDFMarginal(self, m, n, V)

    def lower(self):
        return _EnsembleDisc.EnsembleDisc_lower(self)

    def upper(self):
        return _EnsembleDisc.EnsembleDisc_upper(self)

    def Mean(self, m):
        return _EnsembleDisc.EnsembleDisc_Mean(self, m)

    def StdDev(self, m):
        return _EnsembleDisc.EnsembleDisc_StdDev(self, m)

    def Moments(self, m, k):
        return _EnsembleDisc.EnsembleDisc_Moments(self, m, k)

    def Sample(self, m):
        return _EnsembleDisc.EnsembleDisc_Sample(self, m)

    def SampleMarginal(self, m, n):
        return _EnsembleDisc.EnsembleDisc_SampleMarginal(self, m, n)

    def SampleOne(self, m, j):
        return _EnsembleDisc.EnsembleDisc_SampleOne(self, m, j)

    def PrintDiag(self, *args):
        return _EnsembleDisc.EnsembleDisc_PrintDiag(self, *args)

    def PrintDensity(self, *args):
        return _EnsembleDisc.EnsembleDisc_PrintDensity(self, *args)
EnsembleDisc_swigregister = _EnsembleDisc.EnsembleDisc_swigregister
EnsembleDisc_swigregister(EnsembleDisc)

# This file is compatible with both classic and new-style classes.

