"""
Return the day of year (doy).
"""

def cal2dayno (day, month, year):
    """
    Return the day of year (doy).

    Parameters
    ----------
    day : int or array-like
        Day of month.
    month : int or array-like
        Month (1-12).
    year : int or array-like
        Year.

    Returns
    -------
    int or ndarray
        Day of year.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.cal2dayno(1,1,2000)
    >>> 1
    >>> pyacs.cal2dayno(29,2,2000)
    >>> 60
    >>> pyacs.cal2dayno(29.,2,2000)
    >>> 60.
    >>> pyacs.cal2dayno(29,2,2001)
    >>> !!!  29  out of range
    >>> !!! bad day
    """
    
    from ._common import __monthOK, __dayOK, days
    from .leap_year import leap_year
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        dayno=np.array(list(map(cal2dayno,day, month, year)))

    else:
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year) 
    
        month=month-1 # For array indexing
    
        if (leap_year(year)):days[1] = 29
        else:days[1] = 28
    
    
        dayno = day
        mon=0
        while mon<month:
            dayno = dayno + days[mon]
            mon=mon+1

    return( np.array(dayno, dtype=int)+0 )
