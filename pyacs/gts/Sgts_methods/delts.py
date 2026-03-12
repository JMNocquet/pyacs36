###################################################################
def delts(self,code):
###################################################################
    """Remove one time series from this Sgts (in-place).

    Parameters
    ----------
    code : str
        Site code to remove.

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

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    delattr(self, code)
