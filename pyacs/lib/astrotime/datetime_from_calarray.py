"""
Build a datetime array from an array of calendar dates.
"""

def datetime_from_calarray(calarray, hour=12, minute=0, sec=0,ms=0):
    """
    Build a datetime array from an array of calendar dates.

    Parameters
    ----------
    calarray : array-like, 1D or 2D
        Columns order: year, month, mday, hour, minute, second, microsecond.
    hour : int, optional
        Default hour if not in calarray (default 12).
    minute : int, optional
        Default minute (default 0).
    sec : int, optional
        Default second (default 0).
    ms : int, optional
        Default millisecond (default 0).

    Returns
    -------
    ndarray
        1D array of datetime instances.

    Notes
    -----
    year, month, mday are mandatory. If hour/minute/sec/ms not provided,
    defaults are hour=12, minute=0, sec=0, ms=0 (middle of day).
    """

    #import
    import numpy as np
    from datetime import datetime

    # convert to numpy array if a list of list was provided
    calarray = np.array( calarray )
    

    # case array is 1D
    if calarray.ndim == 1:
        calarray = calarray.reshape(-1,calarray.size)
        
    # prepare a formatted array
    
    fmt_array = np.zeros( (calarray.shape[0] , 7) , dtype = int )
    
    # default value for hour
    
    fmt_array[:,3] = 12
    
    # fill the year, month, mday columns
    
    fmt_array[:,:3] = calarray[:,:3].astype(int)
    
    # if calarray is not fully integer, then we assume that microseconds are not provided, but seconds are decimal
    if  calarray.dtype is not np.dtype('int'):
        # if there are 7 columns, then microseconds are provided and dtype was accidentally of float type
        if calarray.shape[1] == 7:
            fmt_array = calarray.astype(int)
        # else, the decimal part of seconds are converted to microseconds
        else:
            seconds = calarray[:,5].astype(int)
            microseconds = ( np.modf(calarray[:,5])[0] * 1.E6).astype(int) 
            
            fmt_array[:,3] = calarray[:,3].astype(int)
            fmt_array[:,4] = calarray[:,4].astype(int)
            fmt_array[:,5] = seconds
            fmt_array[:,6] = microseconds
            
    return( np.array(list(map(datetime,fmt_array[:,0],fmt_array[:,1],fmt_array[:,2],fmt_array[:,3],fmt_array[:,4],fmt_array[:,5],fmt_array[:,6]))) )
