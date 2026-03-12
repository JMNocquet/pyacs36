"""
Convert Python datetime to modified Julian day.
"""

def datetime2mjd(datetime):
    """
    Convert Python datetime to modified Julian day.

    Parameters
    ----------
    datetime : datetime or array-like
        Datetime instance(s).

    Returns
    -------
    float or ndarray
        Modified Julian day.
    """

    from .hmsmicros2uts import hmsmicros2uts
    from .uts2ut import uts2ut
    from .cal2mjd import cal2mjd
    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)
    
    if isinstance(datetime, np.ndarray):
        mjd=np.array(list(map(datetime2mjd,datetime)))
    
    else:

        
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        h=datetime.hour
        mn=datetime.minute
        s=datetime.second
        microsecond=datetime.microsecond
        uts=hmsmicros2uts(h,mn,s=s,microsecond=microsecond)
        ut=uts2ut(uts)
        mjd=cal2mjd(mday, month, year, ut=ut)
        
    return mjd
