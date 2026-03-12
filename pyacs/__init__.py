__version__ = '0.66.20'
import sys
import importlib
_this = sys.modules[__name__]
_verbose_mod = importlib.import_module('pyacs.verbose')
_this.verbose = _verbose_mod.verbose
_debug_mod = importlib.import_module('pyacs.debug')
_this.debug = _debug_mod.debug
import pyacs.lib
