###################################################################
def differentiate(self):
###################################################################
    """
    differentiate the current time series
    :return: the differentiated time series as a new Gts object
    :note : differentiation is made on .data. .data_xyz is set to None.
    """
    # import 
    import inspect
    import numpy as np
    import pyacs.lib.astrotime

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

    new_Gts=self.copy()
    data=self.data
    data_diff=np.diff(data,axis=0)
    
    # for dates
    data_mjd = pyacs.lib.astrotime.decyear2mjd( data[:,0] )
    new_data_mjd = ( data_mjd[0:-1] + data_mjd[1:] ) / 2. 
    
#    data_diff[:,0]=data[0:len(data[:,0])-1,0]+data_diff[:,0]

    data_diff[:,0] = pyacs.lib.astrotime.mjd2decyear( new_data_mjd )

    new_Gts.data=data_diff.copy()
    
    # set .data_xyz to None
    
    new_Gts.data_xyz = None
    
    return(new_Gts)
