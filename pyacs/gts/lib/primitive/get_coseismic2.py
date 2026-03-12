###################################################################
def get_coseismic2(self,
                   eq_date,
                   window_days=5,
                   sample_after=1,
                   method='median',
                   exclude_eq_day=True,
                   verbose=False):
    ###################################################################
    """
    Get coseismic displacement at a given date.

    Coseismic displacement is the position difference between the median (or mean) of
    window_days before the earthquake and the median (or mean) of sample_after samples
    after the earthquake date.

    Parameters
    ----------
    eq_date : datetime or float
        Earthquake date as datetime instance or decimal year.
    window_days : int, optional
        Number of days before the earthquake used to compute position before.
    sample_after : int, optional
        Number of data after the earthquake used to compute position after.
    method : str, optional
        Method to compute positions: 'median' or 'mean'.
    exclude_eq_day : bool, optional
        If True, exclude the day of the earthquake from the data.
    verbose : bool, optional
        Verbose mode.
    """
    # debug
    debug = False

    # import
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # verbose

    logging.getLogger("my_logger").setLevel(logging.WARNING)

    if not(verbose):
        # WARNING LEVEL: ONLY MAJOR STEPS AND WARNINGS WILL BE PRINTED
        logging.getLogger("my_logger").setLevel(logging.WARNING)
    else:
        if verbose:
            # VERBOSE MODE
            logging.getLogger("my_logger").setLevel(logging.INFO)

    if debug:
        # DEBUG MODE
        logging.getLogger("my_logger").setLevel(logging.DEBUG)

    if logging.getLogger("my_logger").getEffectiveLevel() == logging.WARNING:
        MESSAGE("verbose level: WARNING ONLY")

    if logging.getLogger("my_logger").getEffectiveLevel() == logging.INFO:
        MESSAGE("verbose level: VERBOSE")

    if logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG:
        MESSAGE("verbose level: DEBUG")

    # convert time series to pandas
    pd_ts = self.to_pandas_ts()

    # extract the desired period
