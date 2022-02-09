###################################################################
def displacement(self,sdate=None,edate=None,window=None,method='median',speriod=[],eperiod=[], rounding='day' , verbose=True):
###################################################################
    """
    Calculates displacements between two dates or two periods

    :param sdate: start date in decimal year
    :param edate: start date in decimal year
    :param window: time window in days for searching available dates
    :param method: method to calculate the position. 'median' or 'mean'. default is 'median'.
    :param speriod: period for calculating the start position
    :param eperiod: period for calculating the end position
    :param rounding: rounding for dates. Choose among 'day','hour','minute' or 'second'. default is 'day'.
    :param verbose: verbose mode
    :return: displacement as np.array([dn,de,du,sdn,sde,sdu])
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
