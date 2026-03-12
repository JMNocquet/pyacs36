def dates(self,unit='decyear'):
    """Return start and end dates of this Sgts.

    Parameters
    ----------
    unit : str, optional
        'decyear' or 'datetime'. Default is 'decyear'.

    Returns
    -------
    tuple
        (start_date, end_date) in the chosen unit.
    """

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    start_date = 9999.
    end_date   = 0.

    for wts in self.lGts():
        start_date = min(start_date,wts.data[0,0])
        end_date = max(end_date,wts.data[-1,0])

    return start_date,end_date
