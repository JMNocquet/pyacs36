###################################################################
def sel_radius(self, center, range, verbose=True):
    ###################################################################
    """Select time series for sites within a radius range around center.

    Parameters
    ----------
    center : list or str
        [lon, lat] in decimal degrees, or site code.
    range : list or float
        [min_radius_km, max_radius_km], or single max_radius_km, or site code.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    Sgts
        New Sgts instance.
    """

    # import
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    from pyacs.lib.gmtpoint import GMT_Point

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # decipher center
    if isinstance( center, str ):
        [lon_center,lat_center] = [self.__dict__[center].lon, self.__dict__[center].lat ]
    else:
        [lon_center, lat_center] = center
    # decipher range
    if isinstance( range , float ):
        range = [0,range]

    # start loop

    CNTR = GMT_Point(code='XXXX', lon=lon_center, lat=lat_center, he=0.)
    lsel = []

    for code in self.lcode():
        XXXX = GMT_Point(code='XXXX', lon=self.__dict__[code].lon, lat=self.__dict__[code].lat, he=0.)
        sdistance = CNTR.spherical_distance(XXXX) / 1.E3
        if sdistance >= range[0] and sdistance < range[1]:
            # print("%s %s : %10.1lf " % (EQ__.code,code,EQ__.spherical_distance(XXXX)/1.E3))
            lsel.append(code)

    # ensure the provided site is included if range[0] == 0
    if int(range[0]) == 0 and not(center in lsel):
        lsel.append(center)

    return self.sub(linclude=lsel)
