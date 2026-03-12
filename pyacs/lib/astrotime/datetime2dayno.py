"""
Convert Python datetime to day of year, year, and day fraction.
"""

def datetime2dayno(datetime):
    """
    Convert Python datetime to day of year, year, and day fraction.

    Parameters
    ----------
    datetime : datetime or array-like
        Datetime instance(s).

    Returns
    -------
    year : int or ndarray
        Year.
    dayno : int or ndarray
        Day of year.
    ut : float or ndarray
        Day fraction.
    """

    from .hmsmicros2uts import hmsmicros2uts
    from .uts2ut import uts2ut
    from .cal2dayno import cal2dayno
    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)
    
    if isinstance(datetime, np.ndarray):
        [year,dayno,ut]=np.array(list(map(datetime2dayno,datetime))).T
    
    else:
    
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        uts=hmsmicros2uts(datetime.hour,datetime.minute,datetime.second,datetime.microsecond)
        ut=uts2ut(uts)
        dayno=cal2dayno(mday, month, year)

    return( np.array(year , dtype=int)+0,\
            np.array(dayno, dtype=int)+0,\
            ut)
