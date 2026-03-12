###################################################################
def delnone(self):
###################################################################
    """Remove all series with .data is None.

    Returns
    -------
    Sgts
        New Sgts without None-data series.
    """

    # import
    from pyacs.gts.Sgts import Sgts

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    ngts = self.copy()

    VERBOSE("Sgts has %d gts" % self.n())
    for code in self.lcode():
        if self.__dict__[code].data is None:
            MESSAGE("No data for gts: %s. Removing it." % code)
            ngts.delts(code)
    VERBOSE("Sgts now has %d gts" % self.n())

    return ngts
