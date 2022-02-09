
###################################################################
def extract_ndates_around_date(self,date,n):
###################################################################
    """
    Extract n values before and n values after a given date
    If n values are not available, returns all available values
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
        sel_data_after=np.copy(self.data[np.where((self.data[:,0]>date))])
    except:
        print("-- time series ",self.code," does not have data after ",date)
        return(new_Gts)
    
    
    try:
        sel_data_before=np.copy(self.data[np.where((self.data[:,0]<date))])
    except:
        print("-- time series ",self.code," does not have data before ",date)
        return(new_Gts)
    
    # extract data
    if sel_data_after.shape[0]>n:
        sel_data_after=sel_data_after[:n,:]
    
    if sel_data_before.shape[0]>n:
        sel_data_before=sel_data_before[-n:,:]

    
    sel_data=np.vstack(( sel_data_before, sel_data_after))
        
    new_Gts.data=sel_data
    
    return(new_Gts)    
