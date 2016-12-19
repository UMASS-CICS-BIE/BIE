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
        mname = '.'.join((pkg, '_SquareTile')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_SquareTile')
    _SquareTile = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_SquareTile', [dirname(__file__)])
        except ImportError:
            import _SquareTile
            return _SquareTile
        if fp is not None:
            try:
                _mod = imp.load_module('_SquareTile', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _SquareTile = swig_import_helper()
    del swig_import_helper
else:
    import _SquareTile
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

class SquareTile(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SquareTile, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SquareTile, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _SquareTile.new_SquareTile(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def New(self):
        return _SquareTile.SquareTile_New(self)

    def setup(self, arg2, arg3, arg4, arg5):
        return _SquareTile.SquareTile_setup(self, arg2, arg3, arg4, arg5)

    def corners(self, arg2, arg3, arg4, arg5):
        return _SquareTile.SquareTile_corners(self, arg2, arg3, arg4, arg5)

    def X(self, arg2, arg3):
        return _SquareTile.SquareTile_X(self, arg2, arg3)

    def Y(self, arg2, arg3):
        return _SquareTile.SquareTile_Y(self, arg2, arg3)

    def measure(self, arg2, arg3):
        return _SquareTile.SquareTile_measure(self, arg2, arg3)

    def InTile(self, arg2, arg3):
        return _SquareTile.SquareTile_InTile(self, arg2, arg3)

    def _print(self, outputstream):
        return _SquareTile.SquareTile__print(self, outputstream)

    def printforplot(self, outputstream):
        return _SquareTile.SquareTile_printforplot(self, outputstream)

    def overlapsWith(self, arg2):
        return _SquareTile.SquareTile_overlapsWith(self, arg2)

    def toString(self):
        return _SquareTile.SquareTile_toString(self)
    __swig_destroy__ = _SquareTile.delete_SquareTile
    __del__ = lambda self: None
SquareTile_swigregister = _SquareTile.SquareTile_swigregister
SquareTile_swigregister(SquareTile)

# This file is compatible with both classic and new-style classes.

