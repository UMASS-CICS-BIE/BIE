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
            fp, pathname, description = imp.find_module('_EnsembleKD', [dirname(__file__)])
        except ImportError:
            import _EnsembleKD
            return _EnsembleKD
        if fp is not None:
            try:
                _mod = imp.load_module('_EnsembleKD', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _EnsembleKD = swig_import_helper()
    del swig_import_helper
else:
    import _EnsembleKD
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


class EnsembleKD(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, EnsembleKD, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, EnsembleKD, name)
    __repr__ = _swig_repr
    __swig_setmethods__["ncut"] = _EnsembleKD.EnsembleKD_ncut_set
    __swig_getmethods__["ncut"] = _EnsembleKD.EnsembleKD_ncut_get
    if _newclass:
        ncut = _swig_property(_EnsembleKD.EnsembleKD_ncut_get, _EnsembleKD.EnsembleKD_ncut_set)
    __swig_setmethods__["minsub"] = _EnsembleKD.EnsembleKD_minsub_set
    __swig_getmethods__["minsub"] = _EnsembleKD.EnsembleKD_minsub_get
    if _newclass:
        minsub = _swig_property(_EnsembleKD.EnsembleKD_minsub_get, _EnsembleKD.EnsembleKD_minsub_set)

    def __init__(self, *args):
        this = _EnsembleKD.new_EnsembleKD(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _EnsembleKD.delete_EnsembleKD
    __del__ = lambda self: None

    def Reset(self, *args):
        return _EnsembleKD.EnsembleKD_Reset(self, *args)

    def setBucketSize(self, n):
        return _EnsembleKD.EnsembleKD_setBucketSize(self, n)

    def ComputeDistribution(self, *args):
        return _EnsembleKD.EnsembleKD_ComputeDistribution(self, *args)

    def ComputeDistributionMarginal(self, out):
        return _EnsembleKD.EnsembleKD_ComputeDistributionMarginal(self, out)

    def New(self):
        return _EnsembleKD.EnsembleKD_New(self)

    def getDim(self):
        return _EnsembleKD.EnsembleKD_getDim(self)

    def PDF(self, arg2):
        return _EnsembleKD.EnsembleKD_PDF(self, arg2)

    def logPDF(self, arg2):
        return _EnsembleKD.EnsembleKD_logPDF(self, arg2)

    def logPDFMarginal(self, m, n, V):
        return _EnsembleKD.EnsembleKD_logPDFMarginal(self, m, n, V)

    def lower(self):
        return _EnsembleKD.EnsembleKD_lower(self)

    def upper(self):
        return _EnsembleKD.EnsembleKD_upper(self)

    def Mean(self, m):
        return _EnsembleKD.EnsembleKD_Mean(self, m)

    def StdDev(self, m):
        return _EnsembleKD.EnsembleKD_StdDev(self, m)

    def Moments(self, m, k):
        return _EnsembleKD.EnsembleKD_Moments(self, m, k)

    def Sample(self, m):
        return _EnsembleKD.EnsembleKD_Sample(self, m)

    def SampleMarginal(self, m, n):
        return _EnsembleKD.EnsembleKD_SampleMarginal(self, m, n)

    def SampleOne(self, m, j):
        return _EnsembleKD.EnsembleKD_SampleOne(self, m, j)

    def PrintDiag(self, *args):
        return _EnsembleKD.EnsembleKD_PrintDiag(self, *args)

    def PrintDensity(self, dim1, dim2, num1, num2, file):
        return _EnsembleKD.EnsembleKD_PrintDensity(self, dim1, dim2, num1, num2, file)
EnsembleKD_swigregister = _EnsembleKD.EnsembleKD_swigregister
EnsembleKD_swigregister(EnsembleKD)
cvar = _EnsembleKD.cvar

# This file is compatible with both classic and new-style classes.


