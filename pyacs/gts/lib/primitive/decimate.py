###################################################################
def decimate(self,time_step=30.,dates=[],method='median',verbose=False):
###################################################################
    """
    decimate a time series
    
    :param time_step: time step in days
    :param dates: list of dates where point are forced to be written regardless time_step
    :param method: method used to be used to calculated the position. choose among ['median','mean','exact']
    :param verbose: verbose mode
    
    :return : new Gts  
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

    # start and end date
    start_decyear = self.data[0,0]
    end_decyear = self.data[-1,0]

    start_mjd = at.decyear2mjd(start_decyear)
    end_mjd = at.decyear2mjd(end_decyear)
    
    # new_gts to be filled
    new_gts=self.copy()
    
    new_gts.data=np.zeros((1,self.data.shape[1]))
    
    # loop on dates by time_step
    
    for date in np.arange(start_mjd ,end_mjd + time_step / 2. , time_step ):
        
        sub_ts=self.\
                    extract_periods( [at.mjd2decyear(date-time_step/2.),at.mjd2decyear(date+time_step/2.)] ,\
                                      verbose = verbose)
                
        if sub_ts.data is None:
            continue
        
        if method == 'median':
            new_obs = np.median(sub_ts.data,axis=0)
        if method == 'mean':
            new_obs = np.mean(sub_ts.data,axis=0)
        if method == 'exact':
            new_obs = sub_ts.data[int(sub_ts.data.shape[0])-1,:]

    
        new_gts.data = np.vstack( ( new_gts.data, new_obs.reshape(1,-1) ) )

    # case dates are provided
    
    for date in dates:
        sub_ts=self.extract_periods( [date,end_decyear] , verbose = verbose)
        
        if sub_ts.data is not None:
            new_obs = sub_ts.data[0,:]
            if verbose:
                print("-- adding observation at user requested date: %lf = %s = doy %d" % 
                      ( date, '-'.join(map(str, at.decyear2cal(date))), at.decyear2dayno(date)))
                
            if new_gts.data is not None:
                new_gts.data = np.vstack( ( new_gts.data, new_obs.reshape(1,-1) ) )
            else:
                new_gts.data = new_obs.new_obs.reshape(1,-1)

    # remove first obs
    
    if new_gts.data.shape[0] == 1:
        if verbose:
            print('!!! decimated Gts has no date')
        new_gts.data = None
        new_gts.data_xyz = None
        return(new_gts)
    
    else:
        
        new_gts.data = np.delete(new_gts.data,0,axis=0)
        new_gts.neu2xyz()
        
        # reorder
        new_gts.reorder(verbose=verbose)
        if verbose:
            print('-- decimated Gts has ',new_gts.data.shape[0],' entries')

        return(new_gts)
