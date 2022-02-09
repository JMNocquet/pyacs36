###################################################################
def extract_ndates_after_date(self,date,n,verbose=False):
###################################################################
    """
    Extract n values after a given date
    If n values are not available, returns all available values after date
    .data is set to None if no value at all is available

    :param date: date in decimal year
    :param n: number of observations to be extracted
    
    :return: a new Gts

    """

    # import 
    import inspect
    import numpy as np

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

    # copy
    new_Gts=self.copy( data_xyz=None )
    new_Gts.data=None
    
    try:
        sel_data=np.copy(self.data[np.where((self.data[:,0]>date))])
    except:
        print("-- time series ",self.code," does not have data after ",date)
        return(new_Gts)
    
    # extract data
    if sel_data.shape[0]>n:
        sel_data=sel_data[:n,:]
    elif sel_data.shape[0]>0:
        sel_data=sel_data[:,:]
    else:
        sel_data = None
        if verbose:
            print("-- time series ",self.code," does not have data after ",date)
    
    new_Gts.data=sel_data
    
    return(new_Gts)    
