###################################################################
def substract_ts(self,ts,tol=0.05,verbose=True):
###################################################################
    """
    substract the ts provided as argument to the current time series
    
    :param ts: time series to be substracted as a Gts instance
    :param tol: date tolerance to decide whether two dates are identical in both time series. default = 1/4 day
    :param verbose: verbose mode
    :return : new Gts
    
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.gts.Gts

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

    tol = tol /365.25 # tolerance in decimal year

    # find common dates
    
    #common_dates=np.intersect1d(np.round(self.data[:,0],decimals=7),np.round(ts.data[:,0],decimals=7))

    lindex=pyacs.gts.Gts.get_index_from_dates(self.data[:,0], ts.data, tol)

    if verbose:
        print('-- ',len(lindex),' common dates found.')
         
    common_dates = self.data[lindex,0]

    # test whether there are common dates
    
    if common_dates.shape[0] == 0:
        print('! WARNING No common dates between ',self.code,' from ',self.ifile)
        print('!!! and ',ts.code, ' from ',ts.ifile)
        
        new_data = None

    else:
        # extract
        extracted_target=self.extract_dates(common_dates,tol=tol)
        extracted_ref=ts.extract_dates(common_dates,tol=tol)
    
        new_data=extracted_target.data - extracted_ref.data
        
        new_data[:,0]=common_dates 
        
        new_data[:,4]=np.sqrt(extracted_target.data[:,4]**2+extracted_ref.data[:,4]**2) 
        new_data[:,5]=np.sqrt(extracted_target.data[:,5]**2+extracted_ref.data[:,5]**2) 
        new_data[:,6]=np.sqrt(extracted_target.data[:,6]**2+extracted_ref.data[:,6]**2) 
                
    new_code=self.code+'_'+ts.code

    new_Gts=self.copy(data=new_data)
    new_Gts.code=new_code
    # .data_xyz set to None
    new_Gts.data_xyz = None

    return(new_Gts)
