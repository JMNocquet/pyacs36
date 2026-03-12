"""
Convert seconds since 1980-01-01 00:00:00 to datetime object(s).
"""

def seconds2datetime( seconds ):
    """
    Convert seconds since 1980-01-01 00:00:00 to datetime object(s).

    Parameters
    ----------
    seconds : int or array-like
        Seconds since 1980-01-01 00:00:00.

    Returns
    -------
    datetime or ndarray
        Datetime instance(s).
    """
    # import  
    import numpy as np
    from datetime import datetime , timedelta
    
    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)

    if isinstance( seconds , np.ndarray):
        np_datetime =   np.array(list(map(lambda x: ref_date_time + timedelta(seconds=x) ,  seconds.astype(int).tolist() ))) 
    else:
        np_datetime =   ref_date_time + timedelta(seconds=int(seconds)) 
    return np_datetime
