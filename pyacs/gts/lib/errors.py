"""
Exception class for Gts
"""
from colors import red

###############################################################################
class GtsError(Exception):
###############################################################################
    pass

###############################################################################
class GtsInputDataNone(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: input gts has no data for site: %s" % 
                          ( self.method_name, self.lib , self.gts.code ) ) 
        return( return_str )
    
###############################################################################
class GtsInputData_xyzNone(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: input gts has no data_xyz for site: %s" % 
                          ( self.method_name, self.lib , self.gts.code ) ) 
        return( return_str )
        
###############################################################################
class GtsInputData_allNone(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: input gts has neither data nor data_xyz for site: %s" % 
                          ( self.method_name, self.lib , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputDataTypeError(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data is not a 2D numpy array. gts.data.type: %s for site: %s " % 
                          ( self.method_name, self.lib , type(self.gts.data).__name__ , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputDataBadDim(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data is not a 2D numpy array. gts.data.ndim: %d for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data.ndim , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputDataBadNcolumns(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data has wrong number of columns (must be in [7,10]). gts.data.shape[1]: %d for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data.shape[1] , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputDataBadNrows(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data has wrong number of rows (must be >0). gts.data.shape[0]: %d for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data.shape[0] , self.gts.code) ) 
        return( return_str )

###############################################################################
class GtsInputData_xyzTypeError(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data_xyz is not a 2D numpy array. gts.data.type: %s for site: %s " % 
                          ( self.method_name, self.lib , type(self.gts.data_xyz).__name__ , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputData_xyzBadDim(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data_xyz is not a 2D numpy array. gts.data_xyz.ndim: %d for site: %s " % 
                          ( self.method_name, self.lib , self.gts.data_xyz.ndim , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputData_xyzBadNcolumns(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data_xyz has wrong number of columns (must be in [7,10]). gts.data_xyz.shape[1]: %d for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data_xyz.shape[1] , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputData_xyzBadNrows(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data_xyz has wrong number of rows (must be >0). gts.data_xyz.shape[0]: %d for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data_xyz.shape[0] , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsInputDataDiffShape(GtsError):
###############################################################################
    
    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data and gts.data_xyz have different shapes. gts.data.shape=(%d,%d) and gts.data_xyz.shape=(%d,%d) for site: %s" % 
                          ( self.method_name, self.lib , self.gts.data.shape[0] , self.gts.data.shape[1], self.gts.data_xyz.shape[0],self.gts.data_xyz.shape[1] , self.gts.code) ) 
        return( return_str )

###############################################################################
class GtsInputDataDiffDate(GtsError):
###############################################################################

    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS WARNING] method %s from %s: gts.data and gts.data_xyz have different dates. Check gts.data[:,0] and gts.data_xyz[:,0] for site: %s" % 
                          ( self.method_name, self.lib , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsCDataError(GtsError):
###############################################################################
    def __init__(self,method_name,lib,gts):
        self.method_name = method_name    
        self.gts = gts
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] method %s from %s: initial check with cdata failed. Check .data and or .data_xyz for site: %s" % 
                          ( self.method_name, self.lib , self.gts.code ) ) 
        return( return_str )

###############################################################################
class GtsReadFileError(GtsError):
###############################################################################
    def __init__(self,method_name,lib,file_name):
        self.method_name = method_name    
        self.file_name = file_name
        self.lib = lib
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] method %s from %s: Could not read file: %s" % 
                          ( self.method_name, self.lib , self.file_name   ) ) 
        return( return_str )


###############################################################################
class GtsMethodError(GtsError):
###############################################################################
    def __init__(self, method_name, lib, file_name):
        self.method_name = method_name
        self.file_name = file_name
        self.lib = lib

    def __str__(self):
        return_str = red("[PYACS ERROR] running method %s from %s" %
                         (self.method_name, self.lib))
        return (return_str)
