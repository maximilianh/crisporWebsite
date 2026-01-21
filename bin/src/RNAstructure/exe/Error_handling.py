#automated error checking for RNAstructure python interface
from __future__ import print_function
import inspect
from functools import wraps
from collections import defaultdict
debug = False
class StructureError(Exception): pass
class RNAstructureInternalError(Exception):pass
lookup_exceptions = defaultdict(lambda:RuntimeError,
        {   1:IOError,
            2:IOError,
            3:IndexError,
            4:IndexError,
            5:EnvironmentError,
            6:StructureError,
            7:StructureError,
            8:StructureError,
            9:StructureError,
            10:ValueError,
            11:ValueError,
            12:ValueError,
            13:IOError,
            14:RNAstructureInternalError,
            15:ValueError,
            16:ValueError,
            17:ValueError,
            18:ValueError,
            19:ValueError,
            20:ValueError,
            21:RNAstructureInternalError,
            22:RNAstructureInternalError,
            23:ValueError,
            24:ValueError,
            25:ValueError,
            26:ValueError
            })

def check_for_errors(method):
    @wraps(method)
    def RNAstructure_error_checker(self,*args,**kwargs):
        if debug: print ("checking for errors in %s" % method.__name__)
        ret = method(self,*args,**kwargs)
        error = self.GetErrorCode()
        self.ResetError()
        if error != 0:
            raise lookup_exceptions[error]("Error in %s: " % method.__name__ +
                                self.GetErrorMessage(error))
        return ret
    return RNAstructure_error_checker

def check_for_init_errors(method):
    @wraps(method)
    def RNAstructure_error_checker(self,*args):
        if debug: print ("checking for errors in %s" % method.__name__)
        ret = method(self,*args)
        error = self.GetErrorCode()
        if error != 0:
            raise RuntimeError("Error in call to %s.%s: " % (self.__name__,method.__name__) +
                                self.GetErrorMessage(error))
        return ret
    return RNAstructure_error_checker

def is_init(method):
    result =  inspect.ismethod(method) and method.__name__=="__init__"
    if inspect.ismethod(method):
        pass
    return result

def not_excluded(method):
    excluded = ["__repr__","__setattr__","__getattr__","__str__","__init__","<lambda>","swig_repr",
                "GetErrorCode","GetErrorMessage","GetErrorMessageString","ResetError","fromFile","fromString"]
    result = (inspect.ismethod(method) or inspect.isfunction(method)) and method.__name__ not in excluded
    if inspect.ismethod(method):
        if debug: print ("checking if", method.__name__ , "should be excluded: ",result)
    return result

def decorate_methods(decorator,methodtype):
    def decorate(cls):
        for attr in inspect.getmembers(cls, methodtype):
            if debug: print ("decorating %s!" % attr[0])
            setattr(cls, attr[0], decorator(getattr(cls, attr[0])))
        return cls
    return decorate
