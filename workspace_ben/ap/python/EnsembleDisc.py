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


class EnsembleDisc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, EnsembleDisc, name, value)
    __swig_getmethods__ = {}
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
cvar = _EnsembleDisc.cvar

# This file is compatible with both classic and new-style classes.

