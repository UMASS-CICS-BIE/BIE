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
        mname = '.'.join((pkg, '_TileRequestHandlerThread')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_TileRequestHandlerThread')
    _TileRequestHandlerThread = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TileRequestHandlerThread', [dirname(__file__)])
        except ImportError:
            import _TileRequestHandlerThread
            return _TileRequestHandlerThread
        if fp is not None:
            try:
                _mod = imp.load_module('_TileRequestHandlerThread', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TileRequestHandlerThread = swig_import_helper()
    del swig_import_helper
else:
    import _TileRequestHandlerThread
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

class TileRequestHandlerThread(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TileRequestHandlerThread, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TileRequestHandlerThread, name)
    __repr__ = _swig_repr

    def __init__(self, numWorkersInTileGroup, p_GroupIdToSystemId, p_SystemIdToGroupId, p_comms, numberSlavesReady):
        this = _TileRequestHandlerThread.new_TileRequestHandlerThread(numWorkersInTileGroup, p_GroupIdToSystemId, p_SystemIdToGroupId, p_comms, numberSlavesReady)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def reInitialize(self, tileGroupSize, groupIdToSystemId, systemIdToGroupId, comms, numberSlavesReady):
        return _TileRequestHandlerThread.TileRequestHandlerThread_reInitialize(self, tileGroupSize, groupIdToSystemId, systemIdToGroupId, comms, numberSlavesReady)

    def setTileMasterWorkerThread(self, tileMasterWorkerThread):
        return _TileRequestHandlerThread.TileRequestHandlerThread_setTileMasterWorkerThread(self, tileMasterWorkerThread)

    def masterWorkerRegister(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_masterWorkerRegister(self)

    def masterWorkerAvailable(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_masterWorkerAvailable(self)

    def resetPause(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_resetPause(self)

    def unPause(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_unPause(self)

    def stopPausedThread(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_stopPausedThread(self)
    __swig_destroy__ = _TileRequestHandlerThread.delete_TileRequestHandlerThread
    __del__ = lambda self: None

    def run(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_run(self)

    def enqueueRequest(self, request):
        return _TileRequestHandlerThread.TileRequestHandlerThread_enqueueRequest(self, request)

    def readAndProcessMessage(self, status):
        return _TileRequestHandlerThread.TileRequestHandlerThread_readAndProcessMessage(self, status)

    def allWorking(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_allWorking(self)

    def noneWorking(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_noneWorking(self)

    def noMPISlaves(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_noMPISlaves(self)

    def workAvailable(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_workAvailable(self)

    def final(self):
        return _TileRequestHandlerThread.TileRequestHandlerThread_final(self)
    __swig_setmethods__["MPIThreadLoopEvent"] = _TileRequestHandlerThread.TileRequestHandlerThread_MPIThreadLoopEvent_set
    __swig_getmethods__["MPIThreadLoopEvent"] = _TileRequestHandlerThread.TileRequestHandlerThread_MPIThreadLoopEvent_get
    if _newclass:
        MPIThreadLoopEvent = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_MPIThreadLoopEvent_get, _TileRequestHandlerThread.TileRequestHandlerThread_MPIThreadLoopEvent_set)
    __swig_setmethods__["startTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_startTime_set
    __swig_getmethods__["startTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_startTime_get
    if _newclass:
        startTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_startTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_startTime_set)
    __swig_setmethods__["endTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_endTime_set
    __swig_getmethods__["endTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_endTime_get
    if _newclass:
        endTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_endTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_endTime_set)
    __swig_setmethods__["histo"] = _TileRequestHandlerThread.TileRequestHandlerThread_histo_set
    __swig_getmethods__["histo"] = _TileRequestHandlerThread.TileRequestHandlerThread_histo_get
    if _newclass:
        histo = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_histo_get, _TileRequestHandlerThread.TileRequestHandlerThread_histo_set)
    __swig_setmethods__["workerWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_workerWaitTime_set
    __swig_getmethods__["workerWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_workerWaitTime_get
    if _newclass:
        workerWaitTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_workerWaitTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_workerWaitTime_set)
    __swig_setmethods__["workWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_workWaitTime_set
    __swig_getmethods__["workWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_workWaitTime_get
    if _newclass:
        workWaitTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_workWaitTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_workWaitTime_set)
    __swig_setmethods__["mutexWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_mutexWaitTime_set
    __swig_getmethods__["mutexWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_mutexWaitTime_get
    if _newclass:
        mutexWaitTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_mutexWaitTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_mutexWaitTime_set)
    __swig_setmethods__["probeTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_probeTime_set
    __swig_getmethods__["probeTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_probeTime_get
    if _newclass:
        probeTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_probeTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_probeTime_set)
    __swig_setmethods__["IprobeTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_IprobeTime_set
    __swig_getmethods__["IprobeTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_IprobeTime_get
    if _newclass:
        IprobeTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_IprobeTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_IprobeTime_set)
    __swig_setmethods__["loopEventWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_loopEventWaitTime_set
    __swig_getmethods__["loopEventWaitTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_loopEventWaitTime_get
    if _newclass:
        loopEventWaitTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_loopEventWaitTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_loopEventWaitTime_set)
    __swig_setmethods__["thrqMutexHoldTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_thrqMutexHoldTime_set
    __swig_getmethods__["thrqMutexHoldTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_thrqMutexHoldTime_get
    if _newclass:
        thrqMutexHoldTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_thrqMutexHoldTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_thrqMutexHoldTime_set)
    __swig_setmethods__["recvTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_recvTime_set
    __swig_getmethods__["recvTime"] = _TileRequestHandlerThread.TileRequestHandlerThread_recvTime_get
    if _newclass:
        recvTime = _swig_property(_TileRequestHandlerThread.TileRequestHandlerThread_recvTime_get, _TileRequestHandlerThread.TileRequestHandlerThread_recvTime_set)
TileRequestHandlerThread_swigregister = _TileRequestHandlerThread.TileRequestHandlerThread_swigregister
TileRequestHandlerThread_swigregister(TileRequestHandlerThread)

# This file is compatible with both classic and new-style classes.


