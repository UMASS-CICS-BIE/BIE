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
            fp, pathname, description = imp.find_module('_cliDistribution', [dirname(__file__)])
        except ImportError:
            import _cliDistribution
            return _cliDistribution
        if fp is not None:
            try:
                _mod = imp.load_module('_cliDistribution', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _cliDistribution = swig_import_helper()
    del swig_import_helper
else:
    import _cliDistribution
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
        this = _cliDistribution.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _cliDistribution.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _cliDistribution.Serializable_swigregister
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
    UniformAdd = _cliDistribution.Distribution_UniformAdd
    UniformMult = _cliDistribution.Distribution_UniformMult
    NormalAdd = _cliDistribution.Distribution_NormalAdd
    NormalMult = _cliDistribution.Distribution_NormalMult
    Undefined = _cliDistribution.Distribution_Undefined

    def __init__(self):
        this = _cliDistribution.new_Distribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _cliDistribution.delete_Distribution
    __del__ = lambda self: None

    def New(self):
        return _cliDistribution.Distribution_New(self)

    def PDF(self, arg2):
        return _cliDistribution.Distribution_PDF(self, arg2)

    def logPDF(self, arg2):
        return _cliDistribution.Distribution_logPDF(self, arg2)

    def CDF(self, arg2):
        return _cliDistribution.Distribution_CDF(self, arg2)

    def lower(self):
        return _cliDistribution.Distribution_lower(self)

    def upper(self):
        return _cliDistribution.Distribution_upper(self)

    def Mean(self):
        return _cliDistribution.Distribution_Mean(self)

    def StdDev(self):
        return _cliDistribution.Distribution_StdDev(self)

    def Moments(self, arg2):
        return _cliDistribution.Distribution_Moments(self, arg2)

    def Sample(self):
        return _cliDistribution.Distribution_Sample(self)

    def setWidth(self, x):
        return _cliDistribution.Distribution_setWidth(self, x)

    def Type(self):
        return _cliDistribution.Distribution_Type(self)

    def Dim(self):
        return _cliDistribution.Distribution_Dim(self)
Distribution_swigregister = _cliDistribution.Distribution_swigregister
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
    __swig_destroy__ = _cliDistribution.delete_SampleDistribution
    __del__ = lambda self: None

    def CDF(self, arg2):
        return _cliDistribution.SampleDistribution_CDF(self, arg2)

    def New(self):
        return _cliDistribution.SampleDistribution_New(self)

    def AccumData(self, *args):
        return _cliDistribution.SampleDistribution_AccumData(self, *args)

    def ComputeDistribution(self):
        return _cliDistribution.SampleDistribution_ComputeDistribution(self)

    def numberData(self):
        return _cliDistribution.SampleDistribution_numberData(self)

    def getValue(self, i):
        return _cliDistribution.SampleDistribution_getValue(self, i)

    def getdim(self, i):
        return _cliDistribution.SampleDistribution_getdim(self, i)

    def getRecordType(self):
        return _cliDistribution.SampleDistribution_getRecordType(self)

    def getDataSetSize(self):
        return _cliDistribution.SampleDistribution_getDataSetSize(self)
SampleDistribution_swigregister = _cliDistribution.SampleDistribution_swigregister
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
    __swig_destroy__ = _cliDistribution.delete_BinnedDistribution
    __del__ = lambda self: None

    def New(self):
        return _cliDistribution.BinnedDistribution_New(self)

    def getLow(self, i):
        return _cliDistribution.BinnedDistribution_getLow(self, i)

    def getHigh(self, i):
        return _cliDistribution.BinnedDistribution_getHigh(self, i)

    def CDF(self, arg2):
        return _cliDistribution.BinnedDistribution_CDF(self, arg2)
BinnedDistribution_swigregister = _cliDistribution.BinnedDistribution_swigregister
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
        this = _cliDistribution.new_PointDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _cliDistribution.delete_PointDistribution
    __del__ = lambda self: None

    def New(self):
        return _cliDistribution.PointDistribution_New(self)

    def getRecordType(self):
        return _cliDistribution.PointDistribution_getRecordType(self)

    def AccumData(self, v, datapoint):
        return _cliDistribution.PointDistribution_AccumData(self, v, datapoint)

    def getdim(self, i):
        return _cliDistribution.PointDistribution_getdim(self, i)

    def numberData(self):
        return _cliDistribution.PointDistribution_numberData(self)

    def getValue(self, i):
        return _cliDistribution.PointDistribution_getValue(self, i)

    def Point(self):
        return _cliDistribution.PointDistribution_Point(self)

    def CDF(self, v):
        return _cliDistribution.PointDistribution_CDF(self, v)
PointDistribution_swigregister = _cliDistribution.PointDistribution_swigregister
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
        this = _cliDistribution.new_NullDistribution()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _cliDistribution.delete_NullDistribution
    __del__ = lambda self: None

    def New(self):
        return _cliDistribution.NullDistribution_New(self)
NullDistribution_swigregister = _cliDistribution.NullDistribution_swigregister
NullDistribution_swigregister(NullDistribution)

class cliDistribution(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, cliDistribution, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, cliDistribution, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _cliDistribution.new_cliDistribution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this

    def getval(self, pos):
        return _cliDistribution.cliDistribution_getval(self, pos)

    def setval(self, i, val):
        return _cliDistribution.cliDistribution_setval(self, i, val)
    __swig_destroy__ = _cliDistribution.delete_cliDistribution
    __del__ = lambda self: None
cliDistribution_swigregister = _cliDistribution.cliDistribution_swigregister
cliDistribution_swigregister(cliDistribution)

# This file is compatible with both classic and new-style classes.


