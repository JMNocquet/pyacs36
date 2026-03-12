"""
Convert a calendar date (universal time) to modified Julian day (JD-2400000.5).
"""

def cal2mjd(day, month, year, ut=0.5):
    """
    Convert a calendar date (universal time) to modified Julian day (JD-2400000.5).

    Parameters
    ----------
    day : int or array-like
        Day of month.
    month : int or array-like
        Month (1-12).
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction in [0., 1.[. Default 0.5 (middle of day).

    Returns
    -------
    float or ndarray
        Modified Julian day.

    Notes
    -----
    When ut is not provided, the middle of the day is used.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.cal2mjd(29,2,2000)
    >>> 51603.5
    >>> pyacs.cal2mjd(29,2,2000,ut=0.0)
    >>> 51603.0
    """

    from ._common import __monthOK, __dayOK, __utOK
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=day*0.0+ut
        mjd=np.array(list(map(cal2mjd,day, month, year, ut)))

    else:
        
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year)
        __utOK(ut)
    
        if (month <= 2):
            m = int(month+9)
            y = int(year-1)
        else:
            m = int(month-3)
            y = int(year)
    
        c = int(y/100)
        y = y-c*100
        x1 = int(146097.0*c/4.0)
        x2 = int(1461.0*y/4.0)
        x3 = int((153.0*m+2.0)/5.0)

        mjd=int(x1+x2+x3+day-678882)
        mjd+=ut

    return(mjd)
