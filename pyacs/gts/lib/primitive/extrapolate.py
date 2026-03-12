def extrapolate(self,idates,kind='linear'):
    """

    Extrapolate the time series at dates outside the observation span.

    Parameters
    ----------
    idates : list or ndarray
        Dates in decimal year where extrapolation will be performed.
    kind : str, optional
        Interpolation kind for scipy.interpolate.interp1d (e.g. 'linear', 'nearest', 'cubic').

    Returns
    -------
    Gts
        New Gts with extrapolated values at the given dates.
    """

    # import
    import numpy as np
    from scipy import interpolate
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.lib.astrotime as at

    data = np.copy(self.data)

    # it is better to work in seconds
    np_date_s = at.decyear2seconds(self.data[:,0],rounding='day')
    np_idate_s = at.decyear2seconds(np.array(idates),rounding='day')

    # check edates
    sdate = np_date_s[0]
    edate = np_date_s[-1]

    iidate=[]
    for date in np_idate_s:
        if date >= sdate and date <= edate:
            as_datetime = at.seconds2datetime(date)
            WARNING("Provided date %.4lf (%s) is within the observation period of time series %s" %
                    (at.datetime2decyear(as_datetime),as_datetime.isoformat(),self.code))
            WARNING("Removing it from interpolated dates")
        else:
            iidate.append(date)

    tnew = np.sort(np.array(iidate))

    # interpolate
    # north
    f = interpolate.interp1d( np_date_s, self.data[:,1], kind=kind,fill_value='extrapolate')
    ynewn = f(tnew)
    # east
    f = interpolate.interp1d( np_date_s, self.data[:,2], kind=kind,fill_value='extrapolate')
    ynewe = f(tnew)
    # up
    f = interpolate.interp1d( np_date_s, self.data[:,3], kind=kind,fill_value='extrapolate')
    ynewu = f(tnew)

    # adds interpolated values
    # get uncertainty
    sigma = np.median( self.data[:,4:7], axis=0 )
    # fills value to be added
    NEW_OBS = np.zeros((tnew.shape[0],10))
    NEW_OBS[:,0] = at.datetime2decyear( at.seconds2datetime( tnew ) )
    NEW_OBS[:, 1] = ynewn
    NEW_OBS[:, 2] = ynewe
    NEW_OBS[:, 3] = ynewu
    NEW_OBS[:, 4:7] = sigma.reshape(-1,3)

    return self.add_obs( NEW_OBS[:,0] , NEW_OBS[:,1:] )