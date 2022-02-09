###################################################################
def extract_ndates_before_date(self,date,n,verbose=False):
###################################################################
    """
    Extract n values before a given date
    If n values are not available, returns all available values before date
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
    
    # extract data
    sel_data     = np.copy(self.data[ np.where( self.data[:,0] < date ) ] )
    if self.data_xyz is not None:
        sel_data_xyz = np.copy(self.data_xyz[ np.where( self.data[:,0] < date ) ] )
        
    # extract data
    if sel_data.shape[0]>n:
        sel_data=sel_data[-n:,:]
        if self.data_xyz is not None:
            sel_data_xyz=sel_data_xyz[-n:,:]
    elif sel_data.shape[0]>0:
        sel_data=sel_data[:,:]
        if self.data_xyz is not None:
            sel_data_xyz=sel_data_xyz[:,:]
    else:
        sel_data = None
        if verbose:
            print("-- time series ",self.code," does not have data before ",date)

    new_Gts.data=sel_data
    if self.data_xyz is not None:
        new_Gts.data_xyz = sel_data_xyz
        
    return(new_Gts)    
