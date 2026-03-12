"""
Convert a calendar date to a Python datetime instance.
"""

def cal2datetime(day,month,year,ut=0.5):
    """
    Convert a calendar date to a Python datetime instance.

    Parameters
    ----------
    day : int or array-like
        Day of month.
    month : int or array-like
        Month (1-12).
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction in [0., 1.[. Default 0.5.

    Returns
    -------
    datetime or ndarray of datetime
        Datetime instance(s). At 12:00 (ut=0.5) by default.
    """
    
    from ._common import __monthOK, __dayOK
    from .cal2dayno import cal2dayno
    from .dayno2datetime import dayno2datetime
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=day*0.0+ut
        return_datetime=np.array(list(map(cal2datetime,day,month,year,ut)))

    else:
    
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year) 

        dayno=cal2dayno(day,month,year)
        return_datetime=dayno2datetime(dayno,year,ut=ut)
        
    return return_datetime
