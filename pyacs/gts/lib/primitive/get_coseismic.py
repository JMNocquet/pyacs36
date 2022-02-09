###################################################################
def get_coseismic(self,eq_date,window_days=5,sample_after=1,method='median',in_place=False):
###################################################################
    """
    Get coseismic displacement at a given date.
    Coseismic displacement is estimated as the position difference between the median of window_days before 
    the earthquake date and the median of sample_after samples after the earthquake date. 
    
    note: only median method implemented
    
    """
    
    # import 
    import inspect
    import numpy as np
    import pyacs.lib.astrotime
    from datetime import timedelta



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

    # eq date as datetime
    eq_datetime = pyacs.lib.astrotime.decyear2datetime( eq_date )
    # get the eq day at 00:00:00 and then go backward by window_days
    wdatetime = eq_datetime.replace(hour=0,minute=0,second=0,microsecond=0) - timedelta(days =  window_days )
    # convert back to decimal year
    wdecyear = pyacs.lib.astrotime.datetime2decyear(wdatetime)
    # make extraction between eq - window_days and eq_date
    lindex = np.where( ( self.data[:,0] > wdecyear ) & ( self.data[:,0] < eq_date ) )[0]
    data_before_eq = self.data[lindex,1:4]
    # take sample_after values after eq_date
    lindex = np.where( self.data[:,0] > eq_date )[0][0]
    data_after_eq = self.data[lindex:lindex + sample_after,1:4]
    # get offset values
    disp = np.median(data_after_eq,axis=0) - np.median(data_before_eq,axis=0)
    # calculates std before eq
    std_sdate=np.std( data_before_eq , axis=0 )
    # calculates std after eq
    std_edate=np.std( data_after_eq , axis=0 )
    
    # if std cannot be calculated or is 0.
    if np.max(std_sdate)==0.0:
        lindex = np.where( ( self.data[:,0] > wdecyear ) & ( self.data[:,0] < eq_date ) )[0]
        std_sdate=np.mean( self.data[lindex,4:7],axis=0)

    if np.max(std_edate)==0.0:
        lindex = np.where( self.data[:,0] > eq_date )[0][0]
        std_sdate=np.mean( self.data[lindex,4:7],axis=0)

    std_disp=np.sqrt(std_edate**2+std_sdate**2)

    new_Gts=self.copy( )
    
    # modify the output time series
    new_Gts.data_xyz = None
    new_Gts.data[:,1:4] = new_Gts.data[:,1:4] - np.median( data_before_eq , axis=0 )
    
    if in_place:
        self.offsets_values=np.array([eq_date]+disp.tolist()+std_disp.tolist()).reshape(1,7)
        return(self)
        del new_Gts
    else:
        new_Gts.offsets_values=np.array([eq_date]+disp.tolist()+std_disp.tolist()).reshape(1,7)
        return(new_Gts)
