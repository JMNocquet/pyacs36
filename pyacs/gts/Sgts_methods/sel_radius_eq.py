###################################################################
def sel_radius_eq(self, eq, time=None, verbose=True):
    ###################################################################
    """Select time series potentially impacted by an earthquake.

    Parameters
    ----------
    eq : list
        [lon, lat, magnitude]; epicenter in decimal degrees.
    time : float, optional
        Earthquake time in decimal year. Default is None.
    verbose : bool, optional
        If True, print progress. Default is True.

    Returns
    -------
    Sgts
        New Sgts with offset_dates set at eq time.

    Notes
    -----
    Radius from NGL formula: radius = 10**(magnitude/2 - 0.79).
    """

    # import
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import numpy as np
    from pyacs.lib.gmtpoint import GMT_Point

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # decipher center
    if not isinstance( eq, list ):
        ERROR("eq should be a list of [lon, lat, magnitude]", exit=True)

    # compute radius
    radius = np.power(10, (eq[2]/2 - 0.79))

    # select time series
    new_ts = self.sel_radius(eq[:2], [0,radius])

    # add offset_dates at eq time
    if time is not None:
        new_ts = new_ts.gts('add_offsets_dates', [time])

    return new_ts