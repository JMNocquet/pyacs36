###################################################################
def add_offsets_dates(self,offsets_dates,in_place=False):
###################################################################
    """
    add_offsets_dates to a time series
    if in_place = True then replace the current time series
    """
    
    ### import
    
    import inspect
    import numpy as np
    
    # check offsets_dates type
    try:
        if type( offsets_dates ) not in [np.ndarray, list]:
            # raise exception
            from pyacs.gts.lib.errors import GtsInputDataTypeError
            raise GtsInputDataTypeError(inspect.stack()[0][3],__name__,self)
    except GtsInputDataTypeError as error:
        # print PYACS WARNING
        print( error )
        return( self )


    if in_place:
        self.offsets_dates = offsets_dates
        return( self )

    else:
        new_ts = self.copy()
        new_ts.offsets_dates = offsets_dates
        return( new_ts )
