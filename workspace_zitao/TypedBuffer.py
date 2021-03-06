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
        mname = '.'.join((pkg, '_TypedBuffer')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_TypedBuffer')
    _TypedBuffer = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TypedBuffer', [dirname(__file__)])
        except ImportError:
            import _TypedBuffer
            return _TypedBuffer
        if fp is not None:
            try:
                _mod = imp.load_module('_TypedBuffer', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TypedBuffer = swig_import_helper()
    del swig_import_helper
else:
    import _TypedBuffer
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

class TypedBuffer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    if _newclass:
        createNew = staticmethod(_TypedBuffer.TypedBuffer_createNew)
    else:
        createNew = _TypedBuffer.TypedBuffer_createNew
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer
    __del__ = lambda self: None

    def setStringValue(self, arg2):
        return _TypedBuffer.TypedBuffer_setStringValue(self, arg2)

    def setIntValue(self, arg2):
        return _TypedBuffer.TypedBuffer_setIntValue(self, arg2)

    def setRealValue(self, arg2):
        return _TypedBuffer.TypedBuffer_setRealValue(self, arg2)

    def setBoolValue(self, arg2):
        return _TypedBuffer.TypedBuffer_setBoolValue(self, arg2)

    def setIntArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_setIntArrayValue(self, *args)

    def setRealArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_setRealArrayValue(self, *args)

    def setBoolArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_setBoolArrayValue(self, *args)

    def getStringValue(self):
        return _TypedBuffer.TypedBuffer_getStringValue(self)

    def getIntValue(self):
        return _TypedBuffer.TypedBuffer_getIntValue(self)

    def getRealValue(self):
        return _TypedBuffer.TypedBuffer_getRealValue(self)

    def getBoolValue(self):
        return _TypedBuffer.TypedBuffer_getBoolValue(self)

    def getIntArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_getIntArrayValue(self, *args)

    def getRealArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_getRealArrayValue(self, *args)

    def getBoolArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_getBoolArrayValue(self, *args)

    def getType(self):
        return _TypedBuffer.TypedBuffer_getType(self)

    def hasValue(self):
        return _TypedBuffer.TypedBuffer_hasValue(self)

    def reset(self):
        return _TypedBuffer.TypedBuffer_reset(self)

    def toString(self):
        return _TypedBuffer.TypedBuffer_toString(self)
TypedBuffer_swigregister = _TypedBuffer.TypedBuffer_swigregister
TypedBuffer_swigregister(TypedBuffer)

def TypedBuffer_createNew(arg2):
    return _TypedBuffer.TypedBuffer_createNew(arg2)
TypedBuffer_createNew = _TypedBuffer.TypedBuffer_createNew

class TypedBuffer_String(TypedBuffer):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_String, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_String, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _TypedBuffer.new_TypedBuffer_String()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def setStringValue(self, arg2):
        return _TypedBuffer.TypedBuffer_String_setStringValue(self, arg2)

    def getStringValue(self):
        return _TypedBuffer.TypedBuffer_String_getStringValue(self)

    def toString(self):
        return _TypedBuffer.TypedBuffer_String_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_String
    __del__ = lambda self: None
TypedBuffer_String_swigregister = _TypedBuffer.TypedBuffer_String_swigregister
TypedBuffer_String_swigregister(TypedBuffer_String)

class TypedBuffer_Int(TypedBuffer):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_Int, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_Int, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _TypedBuffer.new_TypedBuffer_Int()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def setIntValue(self, arg2):
        return _TypedBuffer.TypedBuffer_Int_setIntValue(self, arg2)

    def getIntValue(self):
        return _TypedBuffer.TypedBuffer_Int_getIntValue(self)

    def toString(self):
        return _TypedBuffer.TypedBuffer_Int_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_Int
    __del__ = lambda self: None
TypedBuffer_Int_swigregister = _TypedBuffer.TypedBuffer_Int_swigregister
TypedBuffer_Int_swigregister(TypedBuffer_Int)

class TypedBuffer_Real(TypedBuffer):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_Real, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_Real, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _TypedBuffer.new_TypedBuffer_Real()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def setRealValue(self, arg2):
        return _TypedBuffer.TypedBuffer_Real_setRealValue(self, arg2)

    def getRealValue(self):
        return _TypedBuffer.TypedBuffer_Real_getRealValue(self)

    def toString(self):
        return _TypedBuffer.TypedBuffer_Real_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_Real
    __del__ = lambda self: None
TypedBuffer_Real_swigregister = _TypedBuffer.TypedBuffer_Real_swigregister
TypedBuffer_Real_swigregister(TypedBuffer_Real)

class TypedBuffer_Bool(TypedBuffer):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_Bool, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_Bool, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _TypedBuffer.new_TypedBuffer_Bool()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def setBoolValue(self, arg2):
        return _TypedBuffer.TypedBuffer_Bool_setBoolValue(self, arg2)

    def getBoolValue(self):
        return _TypedBuffer.TypedBuffer_Bool_getBoolValue(self)

    def toString(self):
        return _TypedBuffer.TypedBuffer_Bool_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_Bool
    __del__ = lambda self: None
TypedBuffer_Bool_swigregister = _TypedBuffer.TypedBuffer_Bool_swigregister
TypedBuffer_Bool_swigregister(TypedBuffer_Bool)

class TypedBuffer_IntArray(TypedBuffer_Int):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer_Int]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_IntArray, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer_Int]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_IntArray, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _TypedBuffer.new_TypedBuffer_IntArray(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def getIntArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_IntArray_getIntArrayValue(self, *args)

    def setIntArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_IntArray_setIntArrayValue(self, *args)

    def toString(self):
        return _TypedBuffer.TypedBuffer_IntArray_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_IntArray
    __del__ = lambda self: None
TypedBuffer_IntArray_swigregister = _TypedBuffer.TypedBuffer_IntArray_swigregister
TypedBuffer_IntArray_swigregister(TypedBuffer_IntArray)

class TypedBuffer_RealArray(TypedBuffer_Real):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer_Real]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_RealArray, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer_Real]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_RealArray, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _TypedBuffer.new_TypedBuffer_RealArray(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def getRealArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_RealArray_getRealArrayValue(self, *args)

    def setRealArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_RealArray_setRealArrayValue(self, *args)

    def toString(self):
        return _TypedBuffer.TypedBuffer_RealArray_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_RealArray
    __del__ = lambda self: None
TypedBuffer_RealArray_swigregister = _TypedBuffer.TypedBuffer_RealArray_swigregister
TypedBuffer_RealArray_swigregister(TypedBuffer_RealArray)

class TypedBuffer_BoolArray(TypedBuffer_Bool):
    __swig_setmethods__ = {}
    for _s in [TypedBuffer_Bool]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, TypedBuffer_BoolArray, name, value)
    __swig_getmethods__ = {}
    for _s in [TypedBuffer_Bool]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, TypedBuffer_BoolArray, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _TypedBuffer.new_TypedBuffer_BoolArray(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def getBoolArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_BoolArray_getBoolArrayValue(self, *args)

    def setBoolArrayValue(self, *args):
        return _TypedBuffer.TypedBuffer_BoolArray_setBoolArrayValue(self, *args)

    def toString(self):
        return _TypedBuffer.TypedBuffer_BoolArray_toString(self)
    __swig_destroy__ = _TypedBuffer.delete_TypedBuffer_BoolArray
    __del__ = lambda self: None
TypedBuffer_BoolArray_swigregister = _TypedBuffer.TypedBuffer_BoolArray_swigregister
TypedBuffer_BoolArray_swigregister(TypedBuffer_BoolArray)

# This file is compatible with both classic and new-style classes.


