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
        mname = '.'.join((pkg, '_Tessellation')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_Tessellation')
    _Tessellation = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Tessellation', [dirname(__file__)])
        except ImportError:
            import _Tessellation
            return _Tessellation
        if fp is not None:
            try:
                _mod = imp.load_module('_Tessellation', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Tessellation = swig_import_helper()
    del swig_import_helper
else:
    import _Tessellation
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

class twodcoords(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, twodcoords, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, twodcoords, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Tessellation.new_twodcoords(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_setmethods__["x"] = _Tessellation.twodcoords_x_set
    __swig_getmethods__["x"] = _Tessellation.twodcoords_x_get
    if _newclass:
        x = _swig_property(_Tessellation.twodcoords_x_get, _Tessellation.twodcoords_x_set)
    __swig_setmethods__["y"] = _Tessellation.twodcoords_y_set
    __swig_getmethods__["y"] = _Tessellation.twodcoords_y_get
    if _newclass:
        y = _swig_property(_Tessellation.twodcoords_y_get, _Tessellation.twodcoords_y_set)
    __swig_destroy__ = _Tessellation.delete_twodcoords
    __del__ = lambda self: None
twodcoords_swigregister = _Tessellation.twodcoords_swigregister
twodcoords_swigregister(twodcoords)

class Tessellation(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Tessellation, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Tessellation, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _Tessellation.delete_Tessellation
    __del__ = lambda self: None

    def GetTile(self, tileid):
        return _Tessellation.Tessellation_GetTile(self, tileid)

    def FindAll(self, x, y, found):
        return _Tessellation.Tessellation_FindAll(self, x, y, found)

    def GetRootTiles(self):
        return _Tessellation.Tessellation_GetRootTiles(self)

    def GetRootNodes(self):
        return _Tessellation.Tessellation_GetRootNodes(self)

    def IsValidTileID(self, tileid):
        return _Tessellation.Tessellation_IsValidTileID(self, tileid)

    def MinID(self):
        return _Tessellation.Tessellation_MinID(self)

    def MaxID(self):
        return _Tessellation.Tessellation_MaxID(self)

    def NumberTiles(self):
        return _Tessellation.Tessellation_NumberTiles(self)

    def PrintPreOrder(self, *args):
        return _Tessellation.Tessellation_PrintPreOrder(self, *args)

    def GetLimits(self, x1, x2, y1, y2):
        return _Tessellation.Tessellation_GetLimits(self, x1, x2, y1, y2)
Tessellation_swigregister = _Tessellation.Tessellation_swigregister
Tessellation_swigregister(Tessellation)

# This file is compatible with both classic and new-style classes.


