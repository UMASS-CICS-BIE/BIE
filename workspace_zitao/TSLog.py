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
        mname = '.'.join((pkg, '_TSLog')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_TSLog')
    _TSLog = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TSLog', [dirname(__file__)])
        except ImportError:
            import _TSLog
            return _TSLog
        if fp is not None:
            try:
                _mod = imp.load_module('_TSLog', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TSLog = swig_import_helper()
    del swig_import_helper
else:
    import _TSLog
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

class TSLog(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TSLog, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TSLog, name)
    __repr__ = _swig_repr
    circBufNumEntries = _TSLog.TSLog_circBufNumEntries
    __swig_setmethods__["circBuf"] = _TSLog.TSLog_circBuf_set
    __swig_getmethods__["circBuf"] = _TSLog.TSLog_circBuf_get
    if _newclass:
        circBuf = _swig_property(_TSLog.TSLog_circBuf_get, _TSLog.TSLog_circBuf_set)
    __swig_setmethods__["circBufPos"] = _TSLog.TSLog_circBufPos_set
    __swig_getmethods__["circBufPos"] = _TSLog.TSLog_circBufPos_get
    if _newclass:
        circBufPos = _swig_property(_TSLog.TSLog_circBufPos_get, _TSLog.TSLog_circBufPos_set)
    __swig_setmethods__["startTime"] = _TSLog.TSLog_startTime_set
    __swig_getmethods__["startTime"] = _TSLog.TSLog_startTime_get
    if _newclass:
        startTime = _swig_property(_TSLog.TSLog_startTime_get, _TSLog.TSLog_startTime_set)
    __swig_setmethods__["disabled"] = _TSLog.TSLog_disabled_set
    __swig_getmethods__["disabled"] = _TSLog.TSLog_disabled_get
    if _newclass:
        disabled = _swig_property(_TSLog.TSLog_disabled_get, _TSLog.TSLog_disabled_set)

    def __init__(self, *args):
        this = _TSLog.new_TSLog(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def LogIt(self, buf):
        return _TSLog.TSLog_LogIt(self, buf)

    def LogInt(self, buf, i):
        return _TSLog.TSLog_LogInt(self, buf, i)

    def LogDouble(self, buf, d):
        return _TSLog.TSLog_LogDouble(self, buf, d)

    def PrintIt(self):
        return _TSLog.TSLog_PrintIt(self)
    __swig_destroy__ = _TSLog.delete_TSLog
    __del__ = lambda self: None
TSLog_swigregister = _TSLog.TSLog_swigregister
TSLog_swigregister(TSLog)

# This file is compatible with both classic and new-style classes.

