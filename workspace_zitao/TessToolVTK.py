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
        mname = '.'.join((pkg, '_TessToolVTK')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_TessToolVTK')
    _TessToolVTK = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TessToolVTK', [dirname(__file__)])
        except ImportError:
            import _TessToolVTK
            return _TessToolVTK
        if fp is not None:
            try:
                _mod = imp.load_module('_TessToolVTK', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TessToolVTK = swig_import_helper()
    del swig_import_helper
else:
    import _TessToolVTK
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

class vtkGtkCallback(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vtkGtkCallback, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vtkGtkCallback, name)
    __repr__ = _swig_repr
    if _newclass:
        New = staticmethod(_TessToolVTK.vtkGtkCallback_New)
    else:
        New = _TessToolVTK.vtkGtkCallback_New

    def Execute(self, caller, eid, cdata):
        return _TessToolVTK.vtkGtkCallback_Execute(self, caller, eid, cdata)
    __swig_setmethods__["textMapper"] = _TessToolVTK.vtkGtkCallback_textMapper_set
    __swig_getmethods__["textMapper"] = _TessToolVTK.vtkGtkCallback_textMapper_get
    if _newclass:
        textMapper = _swig_property(_TessToolVTK.vtkGtkCallback_textMapper_get, _TessToolVTK.vtkGtkCallback_textMapper_set)
    __swig_setmethods__["textActor"] = _TessToolVTK.vtkGtkCallback_textActor_set
    __swig_getmethods__["textActor"] = _TessToolVTK.vtkGtkCallback_textActor_get
    if _newclass:
        textActor = _swig_property(_TessToolVTK.vtkGtkCallback_textActor_get, _TessToolVTK.vtkGtkCallback_textActor_set)
    __swig_setmethods__["renWin"] = _TessToolVTK.vtkGtkCallback_renWin_set
    __swig_getmethods__["renWin"] = _TessToolVTK.vtkGtkCallback_renWin_get
    if _newclass:
        renWin = _swig_property(_TessToolVTK.vtkGtkCallback_renWin_get, _TessToolVTK.vtkGtkCallback_renWin_set)

    def __init__(self):
        this = _TessToolVTK.new_vtkGtkCallback()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _TessToolVTK.delete_vtkGtkCallback
    __del__ = lambda self: None
vtkGtkCallback_swigregister = _TessToolVTK.vtkGtkCallback_swigregister
vtkGtkCallback_swigregister(vtkGtkCallback)
cvar = _TessToolVTK.cvar

def vtkGtkCallback_New():
    return _TessToolVTK.vtkGtkCallback_New()
vtkGtkCallback_New = _TessToolVTK.vtkGtkCallback_New

class TessToolVTK(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TessToolVTK, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TessToolVTK, name)
    __repr__ = _swig_repr
    __swig_setmethods__["default_height"] = _TessToolVTK.TessToolVTK_default_height_set
    __swig_getmethods__["default_height"] = _TessToolVTK.TessToolVTK_default_height_get
    if _newclass:
        default_height = _swig_property(_TessToolVTK.TessToolVTK_default_height_get, _TessToolVTK.TessToolVTK_default_height_set)
    __swig_setmethods__["default_width"] = _TessToolVTK.TessToolVTK_default_width_set
    __swig_getmethods__["default_width"] = _TessToolVTK.TessToolVTK_default_width_get
    if _newclass:
        default_width = _swig_property(_TessToolVTK.TessToolVTK_default_width_get, _TessToolVTK.TessToolVTK_default_width_set)

    def __init__(self, *args):
        this = _TessToolVTK.new_TessToolVTK(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _TessToolVTK.delete_TessToolVTK
    __del__ = lambda self: None

    def run(self):
        return _TessToolVTK.TessToolVTK_run(self)

    def SetSourceStream(self, mpistrm):
        return _TessToolVTK.TessToolVTK_SetSourceStream(self, mpistrm)

    def Render(self):
        return _TessToolVTK.TessToolVTK_Render(self)

    def getSemaphore(self):
        return _TessToolVTK.TessToolVTK_getSemaphore(self)

    def setProducerSemaphore(self, sem):
        return _TessToolVTK.TessToolVTK_setProducerSemaphore(self, sem)

    def isTessellationInitialized(self):
        return _TessToolVTK.TessToolVTK_isTessellationInitialized(self)

    def setHueLow(self, l):
        return _TessToolVTK.TessToolVTK_setHueLow(self, l)

    def setHueHigh(self, l):
        return _TessToolVTK.TessToolVTK_setHueHigh(self, l)

    def setSaturationLow(self, l):
        return _TessToolVTK.TessToolVTK_setSaturationLow(self, l)

    def setSaturationHigh(self, l):
        return _TessToolVTK.TessToolVTK_setSaturationHigh(self, l)

    def setValueLow(self, l):
        return _TessToolVTK.TessToolVTK_setValueLow(self, l)

    def setValueHigh(self, l):
        return _TessToolVTK.TessToolVTK_setValueHigh(self, l)

    def setAlphaLow(self, l):
        return _TessToolVTK.TessToolVTK_setAlphaLow(self, l)

    def setAlphaHigh(self, l):
        return _TessToolVTK.TessToolVTK_setAlphaHigh(self, l)

    def setPersistenceDir(self, dir):
        return _TessToolVTK.TessToolVTK_setPersistenceDir(self, dir)

    def initVTK(self):
        return _TessToolVTK.TessToolVTK_initVTK(self)

    def getVTKCellPicker(self):
        return _TessToolVTK.TessToolVTK_getVTKCellPicker(self)

    def getVTKRenderer(self):
        return _TessToolVTK.TessToolVTK_getVTKRenderer(self)

    def getVTKRenderWindow(self):
        return _TessToolVTK.TessToolVTK_getVTKRenderWindow(self)

    def getVTKRenderWindowInteractor(self):
        return _TessToolVTK.TessToolVTK_getVTKRenderWindowInteractor(self)

    def setScalarName(self, *args):
        return _TessToolVTK.TessToolVTK_setScalarName(self, *args)

    def setScalarInfoName(self, *args):
        return _TessToolVTK.TessToolVTK_setScalarInfoName(self, *args)
TessToolVTK_swigregister = _TessToolVTK.TessToolVTK_swigregister
TessToolVTK_swigregister(TessToolVTK)

class vtkBIECallback(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vtkBIECallback, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vtkBIECallback, name)
    __repr__ = _swig_repr
    if _newclass:
        New = staticmethod(_TessToolVTK.vtkBIECallback_New)
    else:
        New = _TessToolVTK.vtkBIECallback_New

    def Execute(self, caller, eid, cdata):
        return _TessToolVTK.vtkBIECallback_Execute(self, caller, eid, cdata)
    __swig_setmethods__["textMapper"] = _TessToolVTK.vtkBIECallback_textMapper_set
    __swig_getmethods__["textMapper"] = _TessToolVTK.vtkBIECallback_textMapper_get
    if _newclass:
        textMapper = _swig_property(_TessToolVTK.vtkBIECallback_textMapper_get, _TessToolVTK.vtkBIECallback_textMapper_set)
    __swig_setmethods__["textActor"] = _TessToolVTK.vtkBIECallback_textActor_set
    __swig_getmethods__["textActor"] = _TessToolVTK.vtkBIECallback_textActor_get
    if _newclass:
        textActor = _swig_property(_TessToolVTK.vtkBIECallback_textActor_get, _TessToolVTK.vtkBIECallback_textActor_set)
    __swig_setmethods__["renWin"] = _TessToolVTK.vtkBIECallback_renWin_set
    __swig_getmethods__["renWin"] = _TessToolVTK.vtkBIECallback_renWin_get
    if _newclass:
        renWin = _swig_property(_TessToolVTK.vtkBIECallback_renWin_get, _TessToolVTK.vtkBIECallback_renWin_set)

    def __init__(self):
        this = _TessToolVTK.new_vtkBIECallback()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _TessToolVTK.delete_vtkBIECallback
    __del__ = lambda self: None
vtkBIECallback_swigregister = _TessToolVTK.vtkBIECallback_swigregister
vtkBIECallback_swigregister(vtkBIECallback)

def vtkBIECallback_New():
    return _TessToolVTK.vtkBIECallback_New()
vtkBIECallback_New = _TessToolVTK.vtkBIECallback_New

# This file is compatible with both classic and new-style classes.


