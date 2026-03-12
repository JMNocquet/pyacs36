###################################################################
def copy(self):
###################################################################
    """Return a deep copy of this Sgts.

    Returns
    -------
    Sgts
        New Sgts instance.
    """
    
    import copy

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    new_Sgts=copy.deepcopy( self )

    return( new_Sgts )
