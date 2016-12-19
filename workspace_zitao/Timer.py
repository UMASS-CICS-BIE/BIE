# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_Timer')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_Timer')
    _Timer = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Timer', [dirname(__file__)])
        except ImportError:
            import _Timer
            return _Timer
        if fp is not None:
            try:
                _mod = imp.load_module('_Timer', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Timer = swig_import_helper()
    del swig_import_helper
else:
    import _Timer
del _swig_python_version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

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


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class TimeElapsed(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TimeElapsed, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TimeElapsed, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Timer.new_TimeElapsed(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def getUserTime(self):
        return _Timer.TimeElapsed_getUserTime(self)

    def getSystemTime(self):
        return _Timer.TimeElapsed_getSystemTime(self)

    def getTotalTime(self):
        return _Timer.TimeElapsed_getTotalTime(self)

    def getRealTime(self):
        return _Timer.TimeElapsed_getRealTime(self)

    def __call__(self):
        return _Timer.TimeElapsed___call__(self)
    __swig_destroy__ = _Timer.delete_TimeElapsed
    __del__ = lambda self: None
TimeElapsed_swigregister = _Timer.TimeElapsed_swigregister
TimeElapsed_swigregister(TimeElapsed)

class Timer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Timer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Timer, name)
    __repr__ = _swig_repr

    def __init__(self, precision=False):
        this = _Timer.new_Timer(precision)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def Seconds(self):
        return _Timer.Timer_Seconds(self)

    def Microseconds(self):
        return _Timer.Timer_Microseconds(self)

    def Precision(self):
        return _Timer.Timer_Precision(self)

    def start(self):
        return _Timer.Timer_start(self)

    def restart(self):
        return _Timer.Timer_restart(self)

    def stop(self):
        return _Timer.Timer_stop(self)

    def reset(self):
        return _Timer.Timer_reset(self)

    def getTime(self):
        return _Timer.Timer_getTime(self)

    def isStarted(self):
        return _Timer.Timer_isStarted(self)
    __swig_destroy__ = _Timer.delete_Timer
    __del__ = lambda self: None
Timer_swigregister = _Timer.Timer_swigregister
Timer_swigregister(Timer)

# This file is compatible with both classic and new-style classes.

