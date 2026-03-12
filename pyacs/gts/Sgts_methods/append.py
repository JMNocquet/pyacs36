###################################################################
def append(self, gts):
###################################################################
    """Append a Gts to this Sgts (in-place).

    Parameters
    ----------
    gts : Gts
        Gts instance to append (key = gts.code).

    Returns
    -------
    None
    """

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


#    import inspect
#    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    self.__dict__[gts.code]= gts
