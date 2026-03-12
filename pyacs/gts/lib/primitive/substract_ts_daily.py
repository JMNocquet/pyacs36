###################################################################
def substract_ts_daily(self,ts,verbose=True):
###################################################################
    """Subtract the given time series from the current one (sample-by-sample).

    Parameters
    ----------
    ts : Gts
        Time series to subtract (Gts instance).
    verbose : bool, optional
        If True, print progress. Default is True.

    Returns
    -------
    Gts
        New Gts instance (current minus ts).

    Notes
    -----
    Assumes daily time series; dates are matched.
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.lib.astrotime as at

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    if self.data is None:
        ERROR("%s time series to be processed has not data" % self.code)
        return( self )

    # check data is not None
    
    if ts.data is None:
        ERROR("%s user provided time series has no data" % ts.code)
        new_data = None

    # find common dates

    np_mjd_1 = at.decyear2mjd( self.data[:,0] ).astype( int )
    np_mjd_2 = at.decyear2mjd( ts.data[:,0] ).astype( int )
    
    np_mjd_common_dates = np.intersect1d( np_mjd_1 ,  np_mjd_2 )

    VERBOSE("%d common dates found for %s " % (np_mjd_common_dates.shape[0] , self.code ))

    # test whether there are common dates
    
    if np_mjd_common_dates.shape[0] == 0:
        warning = ("No common dates between for %s " % self.code)
        WARNING(warning)
        new_data = None

    else:
        # extract

        mask1 = np.isin( np_mjd_1, np_mjd_common_dates )
        data1 = self.data[ mask1 ]
        
        mask2 = np.isin( np_mjd_2, np_mjd_common_dates )
        data2 = ts.data[ mask2 ]

        
        
        if data1.shape != data2.shape:
            ERROR("Extracted data have different shapes")
            new_data = None
        
        # ENU
        new_data = data1 - data2
        # dates
        new_data[:,0] = data1[:,0]
        # uncertainties
        new_data[:,4:8] = np.sqrt( data1[:,4:8]**2 + data2[:,4:8]**2 )
        new_data[:,8:] = 0.
                
    #new_code=self.code+'_'+ts.code

    new_Gts=self.copy(data=new_data)
    #new_Gts.code=new_code

    # .data_xyz set to None
    new_Gts.data_xyz = None


    return(new_Gts)
