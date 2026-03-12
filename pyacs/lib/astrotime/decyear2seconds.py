"""
Convert decimal year to seconds since 1980-01-01 00:00:00.
"""

def decyear2seconds( np_decyear , rounding='day' ):
    """
    Convert decimal year to seconds since 1980-01-01 00:00:00.

    Parameters
    ----------
    np_decyear : float or array-like
        Decimal year(s).
    rounding : {'day', 'hour', 'minute', 'second'}, optional
        Rounding: 'day' at 12:00:00, 'hour' at 00:00, 'minute' at current
        minute with 00 seconds, 'second' at integer second.

    Returns
    -------
    int or np.ndarray
        Seconds (scalar) for scalar input, else 1D integer (np.int64) array.
    """

    # import  
    from .decyear2datetime import decyear2datetime
    import numpy as np
    from datetime import datetime, timedelta
    
    input_is_scalar = np.isscalar(np_decyear) or (
        isinstance(np_decyear, np.ndarray) and np_decyear.ndim == 0
    )
    
    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)
    
    np_datetime = decyear2datetime( np_decyear )

    if not( isinstance(np_datetime, np.ndarray) ):
        np_datetime = np.array([ np_datetime ])

    # handle round
    
    for i in np.arange( np_datetime.shape[0] ):
        
        if rounding == 'day':
            np_datetime[i] = np_datetime[i].replace(hour= 12, minute = 0, second = 0, microsecond = 0)
        if rounding == 'hour':
            if np_datetime[i].minute >= 30:
                np_datetime[i] = np_datetime[i] + timedelta( minutes=30 ) 
            np_datetime[i] = np_datetime[i].replace( minute = 0, second = 0, microsecond = 0)
        if rounding == 'minute':
            if np_datetime[i].second >= 30:
                np_datetime[i] = np_datetime[i] + timedelta( seconds=30 ) 
            np_datetime[i] = np_datetime[i].replace( second = 0, microsecond = 0)
        if rounding == 'second':
            np_datetime[i] = np_datetime[i].replace( microsecond = 0)

    
    # convert to timedelta and then to seconds

    np_dates_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)

    # return scalar for scalar input
    if input_is_scalar:
        return int(np_dates_s[0])
    return np_dates_s
