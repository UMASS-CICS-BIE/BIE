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
            fp, pathname, description = imp.find_module('_DifferentialEvolution', [dirname(__file__)])
        except ImportError:
            import _DifferentialEvolution
            return _DifferentialEvolution
        if fp is not None:
            try:
                _mod = imp.load_module('_DifferentialEvolution', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _DifferentialEvolution = swig_import_helper()
    del swig_import_helper
else:
    import _DifferentialEvolution
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


class DifferentialEvolution(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DifferentialEvolution, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DifferentialEvolution, name)
    __repr__ = _swig_repr
    parallel = _DifferentialEvolution.DifferentialEvolution_parallel
    serial = _DifferentialEvolution.DifferentialEvolution_serial
    uniform = _DifferentialEvolution.DifferentialEvolution_uniform
    normal = _DifferentialEvolution.DifferentialEvolution_normal
    cauchy = _DifferentialEvolution.DifferentialEvolution_cauchy
    __swig_setmethods__["cntrl"] = _DifferentialEvolution.DifferentialEvolution_cntrl_set
    __swig_getmethods__["cntrl"] = _DifferentialEvolution.DifferentialEvolution_cntrl_get
    if _newclass:
        cntrl = _swig_property(_DifferentialEvolution.DifferentialEvolution_cntrl_get, _DifferentialEvolution.DifferentialEvolution_cntrl_set)
    __swig_setmethods__["minmc"] = _DifferentialEvolution.DifferentialEvolution_minmc_set
    __swig_getmethods__["minmc"] = _DifferentialEvolution.DifferentialEvolution_minmc_get
    if _newclass:
        minmc = _swig_property(_DifferentialEvolution.DifferentialEvolution_minmc_get, _DifferentialEvolution.DifferentialEvolution_minmc_set)
    __swig_setmethods__["state_iter"] = _DifferentialEvolution.DifferentialEvolution_state_iter_set
    __swig_getmethods__["state_iter"] = _DifferentialEvolution.DifferentialEvolution_state_iter_get
    if _newclass:
        state_iter = _swig_property(_DifferentialEvolution.DifferentialEvolution_state_iter_get, _DifferentialEvolution.DifferentialEvolution_state_iter_set)
    __swig_setmethods__["tiny"] = _DifferentialEvolution.DifferentialEvolution_tiny_set
    __swig_getmethods__["tiny"] = _DifferentialEvolution.DifferentialEvolution_tiny_get
    if _newclass:
        tiny = _swig_property(_DifferentialEvolution.DifferentialEvolution_tiny_get, _DifferentialEvolution.DifferentialEvolution_tiny_set)
    __swig_setmethods__["jfreq"] = _DifferentialEvolution.DifferentialEvolution_jfreq_set
    __swig_getmethods__["jfreq"] = _DifferentialEvolution.DifferentialEvolution_jfreq_get
    if _newclass:
        jfreq = _swig_property(_DifferentialEvolution.DifferentialEvolution_jfreq_get, _DifferentialEvolution.DifferentialEvolution_jfreq_set)
    __swig_setmethods__["linear"] = _DifferentialEvolution.DifferentialEvolution_linear_set
    __swig_getmethods__["linear"] = _DifferentialEvolution.DifferentialEvolution_linear_get
    if _newclass:
        linear = _swig_property(_DifferentialEvolution.DifferentialEvolution_linear_get, _DifferentialEvolution.DifferentialEvolution_linear_set)
    __swig_setmethods__["nfifo"] = _DifferentialEvolution.DifferentialEvolution_nfifo_set
    __swig_getmethods__["nfifo"] = _DifferentialEvolution.DifferentialEvolution_nfifo_get
    if _newclass:
        nfifo = _swig_property(_DifferentialEvolution.DifferentialEvolution_nfifo_get, _DifferentialEvolution.DifferentialEvolution_nfifo_set)

    def __init__(self, *args):
        this = _DifferentialEvolution.new_DifferentialEvolution(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _DifferentialEvolution.delete_DifferentialEvolution
    __del__ = lambda self: None

    def NewState(self, *args):
        return _DifferentialEvolution.DifferentialEvolution_NewState(self, *args)

    def Initialize(self):
        return _DifferentialEvolution.DifferentialEvolution_Initialize(self)

    def Reinitialize(self, width, mp):
        return _DifferentialEvolution.DifferentialEvolution_Reinitialize(self, width, mp)

    def SetControl(self, arg2):
        return _DifferentialEvolution.DifferentialEvolution_SetControl(self, arg2)

    def GetValue(self):
        return _DifferentialEvolution.DifferentialEvolution_GetValue(self)

    def GetPrior(self):
        return _DifferentialEvolution.DifferentialEvolution_GetPrior(self)

    def GetLikelihood(self):
        return _DifferentialEvolution.DifferentialEvolution_GetLikelihood(self)

    def GetState(self):
        return _DifferentialEvolution.DifferentialEvolution_GetState(self)

    def ReportState(self):
        return _DifferentialEvolution.DifferentialEvolution_ReportState(self)

    def LogState(self, *args):
        return _DifferentialEvolution.DifferentialEvolution_LogState(self, *args)

    def PrintState(self):
        return _DifferentialEvolution.DifferentialEvolution_PrintState(self)

    def PrintStepDiagnostic(self):
        return _DifferentialEvolution.DifferentialEvolution_PrintStepDiagnostic(self)

    def PrintStateDiagnostic(self):
        return _DifferentialEvolution.DifferentialEvolution_PrintStateDiagnostic(self)

    def GetStat(self):
        return _DifferentialEvolution.DifferentialEvolution_GetStat(self)

    def NewNumber(self, arg2):
        return _DifferentialEvolution.DifferentialEvolution_NewNumber(self, arg2)

    def NewGamma(self, p):
        return _DifferentialEvolution.DifferentialEvolution_NewGamma(self, p)

    def EnableLogging(self):
        return _DifferentialEvolution.DifferentialEvolution_EnableLogging(self)

    def SetLinearMapping(self, b):
        return _DifferentialEvolution.DifferentialEvolution_SetLinearMapping(self, b)

    def SetJumpFreq(self, n):
        return _DifferentialEvolution.DifferentialEvolution_SetJumpFreq(self, n)

    def AdditionalInfo(self):
        return _DifferentialEvolution.DifferentialEvolution_AdditionalInfo(self)
DifferentialEvolution_swigregister = _DifferentialEvolution.DifferentialEvolution_swigregister
DifferentialEvolution_swigregister(DifferentialEvolution)
cvar = _DifferentialEvolution.cvar

# This file is compatible with both classic and new-style classes.


