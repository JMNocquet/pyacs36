"""
Convert Python datetime to calendar date.
"""

def datetime2cal(datetime):
    """
    Convert Python datetime to calendar date.

    Parameters
    ----------
    datetime : datetime or array-like
        Datetime instance(s).

    Returns
    -------
    year : int or ndarray
        Year.
    month : int or ndarray
        Month (1-12).
    mday : int or ndarray
        Day of month.
    ut : float or ndarray
        Day fraction.
    """
    
    from .hmsmicros2uts import hmsmicros2uts
    from .uts2ut import uts2ut
    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)

    if isinstance(datetime, np.ndarray):
        [year,month,mday,ut]=np.array(list(map(datetime2cal,datetime))).T

    else:
    
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        uts=hmsmicros2uts(datetime.hour,datetime.minute,datetime.second,datetime.microsecond)
        ut=uts2ut(uts)

    return( np.array(year, dtype=int)+0,\
            np.array(month, dtype=int)+0,\
            np.array(mday, dtype=int)+0,\
            ut )
