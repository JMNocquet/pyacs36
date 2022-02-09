###################################################################
def substract_ts_daily(self,ts,verbose=True):
###################################################################
    """
    substract the ts provided as argument to the current time series
    
    :param ts: time series to be substracted as a Gts instance
    :param verbose: verbose mode
    :return : new Gts
    
    :note: this method assumes daily time series
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.lib.astrotime as at
    
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )

    # check data is not None
    
    try:
        if ts.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,ts)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )


    # find common dates

    np_mjd_1 = at.decyear2mjd( self.data[:,0] ).astype( int )
    np_mjd_2 = at.decyear2mjd( ts.data[:,0] ).astype( int )
    
    np_mjd_common_dates = np.intersect1d( np_mjd_1 ,  np_mjd_2 )

    if verbose:
        print('-- ', np_mjd_common_dates.shape[0] , ' common dates found.')

    # test whether there are common dates
    
    if np_mjd_common_dates.shape[0] == 0:
        print('! WARNING No common dates between ',self.code,' from ',self.ifile)
        print('!!! and ',ts.code, ' from ',ts.ifile)
        
        new_data = None

    else:
        # extract

        mask1 = np.isin( np_mjd_1, np_mjd_common_dates )
        data1 = self.data[ mask1 ]
        
        mask2 = np.isin( np_mjd_2, np_mjd_common_dates )
        data2 = ts.data[ mask2 ]

        
        if data1.shape != data2.shape:
            print("!!! ERROR. Extracted data have different shapes")
            new_data = None
        
        # ENU
        new_data = data1 - data2
        # dates
        new_data[:,0] = data1[:,0]
        # uncertainties
        new_data[:,4:8] = np.sqrt( data1[:,4:8]**2 + data2[:,4:8]**2 )
        new_data[:,8:] = 0.
                
    new_code=self.code+'_'+ts.code

    new_Gts=self.copy(data=new_data)
    new_Gts.code=new_code

    # .data_xyz set to None
    new_Gts.data_xyz = None


    return(new_Gts)
