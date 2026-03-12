###################################################################
def displacement(self,sdate=None,edate=None,window=None,method='median',speriod=[],eperiod=[], rounding='day' , verbose=True):
###################################################################
    """
    Calculate displacement between two dates or two periods.

    Parameters
    ----------
    sdate : float, optional
        Start date in decimal year.
    edate : float, optional
        End date in decimal year.
    window : float, optional
        Time window in days for searching available dates.
    method : str, optional
        Method to calculate position: 'median' or 'mean'. Default is 'median'.
    speriod : list, optional
        Period [start, end] for calculating the start position.
    eperiod : list, optional
        Period [start, end] for calculating the end position.
    rounding : str, optional
        Rounding for dates: 'day', 'hour', 'minute', or 'second'. Default is 'day'.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    ndarray
        Displacement as np.array([dn, de, du, sdn, sde, sdu]).
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.lib.astrotime as at
    from colors import red

    # check consistency of passed arguments
    if ( sdate is not None ) and ( speriod != [] ):
        print(red("[ERROR] provide sdate or speriod not both"))
    if ( edate is not None ) and ( eperiod != [] ):
        print(red("[ERROR] provide edate or eperiod not both"))

    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( )

    # convert dates to seconds
    np_seconds = at.decyear2seconds( self.data[:,0] , rounding=rounding )

    # convert window to seconds
    if window is not None:
        window_s = int( window * 86400. )
    else:
        window_s = 1

    # fills speriod and eperiod

    if speriod == []:
        speriod = [ at.decyear2seconds( sdate , rounding=rounding ) - window_s , at.decyear2seconds( sdate , rounding=rounding ) + window_s ]
    else:
        speriod = [ at.decyear2seconds( speriod[0] , rounding=rounding ) - window_s , at.decyear2seconds( speriod[1] , rounding=rounding ) + window_s ]

    if eperiod == []:
        eperiod = [at.decyear2seconds(edate, rounding=rounding) - window_s,
                   at.decyear2seconds(edate, rounding=rounding) + window_s]
    else:
        eperiod = [at.decyear2seconds(eperiod[0], rounding=rounding) - window_s,
                   at.decyear2seconds(eperiod[1], rounding=rounding) + window_s]


    # get values
    lindex = np.argwhere( ( np_seconds >= speriod[0] ) & ( np_seconds <= speriod[1] ) )
    np_start = self.data[lindex,1:7].reshape(-1,6)

    lindex = np.argwhere((np_seconds >= eperiod[0]) & (np_seconds <= eperiod[1]))
    np_end = self.data[lindex, 1:7].reshape(-1,6)

    if method=='median':
        pos_sdate=np.median( np_start , axis=0)
        pos_edate=np.median( np_end   , axis=0)
    else:
        pos_sdate=np.mean( np_start , axis=0)
        pos_edate=np.mean( np_end   , axis=0)


    disp = pos_edate[:3]-pos_sdate[:3]
    std_disp = np.sqrt(pos_sdate**2+pos_edate**2)[3:7]
    
    displacement=[disp[0],disp[1],disp[2],std_disp[0],std_disp[1],std_disp[2]]
    
    return(np.array(displacement))
