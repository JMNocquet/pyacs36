"""Base class for exceptions in PYACS lib core module"""

from colors import red

###############################################################################
class PyacsError(BaseException):
###############################################################################
    """Base class for exceptions in PYACS lib core module"""
    pass

###############################################################################
class BadOption(Exception):
###############################################################################
    """Exception raised for an option not allowed."""
    pass

###############################################################################
class ReadingFile(PyacsError):
###############################################################################

    def __init__(self,file_name):
        self.file_name = file_name
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] Could not read: %s" %
                          ( self.file_name ) ) 
        return( return_str )
