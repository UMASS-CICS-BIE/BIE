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
            fp, pathname, description = imp.find_module('_Tesselation', [dirname(__file__)])
        except ImportError:
            import _Tesselation
            return _Tesselation
        if fp is not None:
            try:
                _mod = imp.load_module('_Tesselation', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Tesselation = swig_import_helper()
    del swig_import_helper
else:
    import _Tesselation
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
        this = _Tesselation.new_Serializable()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Tesselation.delete_Serializable
    __del__ = lambda self: None
Serializable_swigregister = _Tesselation.Serializable_swigregister
Serializable_swigregister(Serializable)

class compx(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, compx, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, compx, name)
    __repr__ = _swig_repr

    def __call__(self, a, b):
        return _Tesselation.compx___call__(self, a, b)

    def __init__(self):
        this = _Tesselation.new_compx()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Tesselation.delete_compx
    __del__ = lambda self: None
compx_swigregister = _Tesselation.compx_swigregister
compx_swigregister(compx)

class compy(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, compy, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, compy, name)
    __repr__ = _swig_repr

    def __call__(self, a, b):
        return _Tesselation.compy___call__(self, a, b)

    def __init__(self):
        this = _Tesselation.new_compy()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Tesselation.delete_compy
    __del__ = lambda self: None
compy_swigregister = _Tesselation.compy_swigregister
compy_swigregister(compy)
cvar = _Tesselation.cvar

class compxy(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, compxy, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, compxy, name)
    __repr__ = _swig_repr

    def __call__(self, a, b):
        return _Tesselation.compxy___call__(self, a, b)

    def __init__(self):
        this = _Tesselation.new_compxy()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Tesselation.delete_compxy
    __del__ = lambda self: None
compxy_swigregister = _Tesselation.compxy_swigregister
compxy_swigregister(compxy)

class hashCoords(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, hashCoords, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, hashCoords, name)
    __repr__ = _swig_repr

    def __call__(self, a):
        return _Tesselation.hashCoords___call__(self, a)

    def __init__(self):
        this = _Tesselation.new_hashCoords()
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_destroy__ = _Tesselation.delete_hashCoords
    __del__ = lambda self: None
hashCoords_swigregister = _Tesselation.hashCoords_swigregister
hashCoords_swigregister(hashCoords)

class twodcoords(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, twodcoords, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, twodcoords, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Tesselation.new_twodcoords(*args)
        try:
            self.this.append(this)
        except Exception:
            self.this = this
    __swig_setmethods__["x"] = _Tesselation.twodcoords_x_set
    __swig_getmethods__["x"] = _Tesselation.twodcoords_x_get
    if _newclass:
        x = _swig_property(_Tesselation.twodcoords_x_get, _Tesselation.twodcoords_x_set)
    __swig_setmethods__["y"] = _Tesselation.twodcoords_y_set
    __swig_getmethods__["y"] = _Tesselation.twodcoords_y_get
    if _newclass:
        y = _swig_property(_Tesselation.twodcoords_y_get, _Tesselation.twodcoords_y_set)
    __swig_destroy__ = _Tesselation.delete_twodcoords
    __del__ = lambda self: None
twodcoords_swigregister = _Tesselation.twodcoords_swigregister
twodcoords_swigregister(twodcoords)

class Tessellation(Serializable):
    __swig_setmethods__ = {}
    for _s in [Serializable]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Tessellation, name, value)
    __swig_getmethods__ = {}
    for _s in [Serializable]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, Tessellation, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _Tesselation.delete_Tessellation
    __del__ = lambda self: None

    def GetTile(self, tileid):
        return _Tesselation.Tessellation_GetTile(self, tileid)

    def FindAll(self, x, y, found):
        return _Tesselation.Tessellation_FindAll(self, x, y, found)

    def GetRootTiles(self):
        return _Tesselation.Tessellation_GetRootTiles(self)

    def GetRootNodes(self):
        return _Tesselation.Tessellation_GetRootNodes(self)

    def IsValidTileID(self, tileid):
        return _Tesselation.Tessellation_IsValidTileID(self, tileid)

    def MinID(self):
        return _Tesselation.Tessellation_MinID(self)

    def MaxID(self):
        return _Tesselation.Tessellation_MaxID(self)

    def NumberTiles(self):
        return _Tesselation.Tessellation_NumberTiles(self)

    def PrintPreOrder(self, *args):
        return _Tesselation.Tessellation_PrintPreOrder(self, *args)

    def GetLimits(self, x1, x2, y1, y2):
        return _Tesselation.Tessellation_GetLimits(self, x1, x2, y1, y2)
Tessellation_swigregister = _Tesselation.Tessellation_swigregister
Tessellation_swigregister(Tessellation)

# This file is compatible with both classic and new-style classes.

