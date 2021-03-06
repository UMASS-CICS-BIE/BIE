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
            fp, pathname, description = imp.find_module('_bieTags', [dirname(__file__)])
        except ImportError:
            import _bieTags
            return _bieTags
        if fp is not None:
            try:
                _mod = imp.load_module('_bieTags', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _bieTags = swig_import_helper()
    del swig_import_helper
else:
    import _bieTags
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



_bieTags.BIE_ENDOFDATA_MARKER_swigconstant(_bieTags)
BIE_ENDOFDATA_MARKER = _bieTags.BIE_ENDOFDATA_MARKER

_bieTags.BIE_SYNCHRONIZE_TAG_swigconstant(_bieTags)
BIE_SYNCHRONIZE_TAG = _bieTags.BIE_SYNCHRONIZE_TAG

_bieTags.BIE_SYNCHRONIZE_REQUEST_swigconstant(_bieTags)
BIE_SYNCHRONIZE_REQUEST = _bieTags.BIE_SYNCHRONIZE_REQUEST

_bieTags.BIE_INITSOCKET_TAG_swigconstant(_bieTags)
BIE_INITSOCKET_TAG = _bieTags.BIE_INITSOCKET_TAG

_bieTags.BIE_INITSOCKET_REQUEST_swigconstant(_bieTags)
BIE_INITSOCKET_REQUEST = _bieTags.BIE_INITSOCKET_REQUEST

_bieTags.BIE_DATATYPE_TAG_swigconstant(_bieTags)
BIE_DATATYPE_TAG = _bieTags.BIE_DATATYPE_TAG

_bieTags.BIE_DATATYPE_REQUEST_swigconstant(_bieTags)
BIE_DATATYPE_REQUEST = _bieTags.BIE_DATATYPE_REQUEST

_bieTags.BIE_DATA_TAG_swigconstant(_bieTags)
BIE_DATA_TAG = _bieTags.BIE_DATA_TAG

_bieTags.BIE_DATA_REQUEST_swigconstant(_bieTags)
BIE_DATA_REQUEST = _bieTags.BIE_DATA_REQUEST

_bieTags.TESS_CMND_GET_TESSALATION_swigconstant(_bieTags)
TESS_CMND_GET_TESSALATION = _bieTags.TESS_CMND_GET_TESSALATION

_bieTags.TESS_CMND_ADD_FILTER_swigconstant(_bieTags)
TESS_CMND_ADD_FILTER = _bieTags.TESS_CMND_ADD_FILTER

_bieTags.TESS_CMND_CLR_FILTER_swigconstant(_bieTags)
TESS_CMND_CLR_FILTER = _bieTags.TESS_CMND_CLR_FILTER

_bieTags.TESS_CMND_GET_NODEINFO_swigconstant(_bieTags)
TESS_CMND_GET_NODEINFO = _bieTags.TESS_CMND_GET_NODEINFO

_bieTags.TESS_CMND_GET_NODETILES_swigconstant(_bieTags)
TESS_CMND_GET_NODETILES = _bieTags.TESS_CMND_GET_NODETILES

_bieTags.TESS_CMND_GET_TILEINFO_swigconstant(_bieTags)
TESS_CMND_GET_TILEINFO = _bieTags.TESS_CMND_GET_TILEINFO

_bieTags.PERSISTENCE_DIR_swigconstant(_bieTags)
PERSISTENCE_DIR = _bieTags.PERSISTENCE_DIR

_bieTags.TESSELLATION_STORE_swigconstant(_bieTags)
TESSELLATION_STORE = _bieTags.TESSELLATION_STORE

_bieTags.COMMAND_MASK_swigconstant(_bieTags)
COMMAND_MASK = _bieTags.COMMAND_MASK

_bieTags.SWITCH_CLION_MASK_swigconstant(_bieTags)
SWITCH_CLION_MASK = _bieTags.SWITCH_CLION_MASK

_bieTags.SUSPEND_SIMULATION_MASK_swigconstant(_bieTags)
SUSPEND_SIMULATION_MASK = _bieTags.SUSPEND_SIMULATION_MASK

_bieTags.SAMPLE_NEXT_MASK_swigconstant(_bieTags)
SAMPLE_NEXT_MASK = _bieTags.SAMPLE_NEXT_MASK

_bieTags.DETACH_TESSTOOL_MASK_swigconstant(_bieTags)
DETACH_TESSTOOL_MASK = _bieTags.DETACH_TESSTOOL_MASK
# This file is compatible with both classic and new-style classes.


