###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['message',\
                     'warning',\
                     'error',
                     'debug_message',
                     'verbose_message']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

