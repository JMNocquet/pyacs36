"""
Exception class for pyacs.sol
"""
from colors import red

###############################################################################
class PyacsSolError(Exception):
###############################################################################
    pass

###############################################################################
class PyacsSol_BadLsnx(PyacsSolError):
###############################################################################

    def __init__(self,lsnx):
        self.lsnx = lsnx
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] empty or bad lsnx file for -lsinex option: %s" % 
                          ( self.lsnx ) ) 
        return( return_str )

###############################################################################
class PyacsSol_HelmertError(PyacsSolError):
###############################################################################

    def __init__(self, method_name,lib,msg):
        self.msg = msg
        self.method_name = method_name    
        self.lib = lib
        #self.sinex = sinex
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: %s " % 
                          ( self.method_name , self.lib , self.msg ) ) 
        return( return_str )
