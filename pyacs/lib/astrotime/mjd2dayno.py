"""
Convert modified Julian day to year and day of year (universal time).
"""

def mjd2dayno(mjd):
    """
    Convert modified Julian day to year and day of year (universal time).

    Parameters
    ----------
    mjd : float or array-like
        Modified Julian day.

    Returns
    -------
    dayno : int or ndarray
        Day of year.
    year : int or ndarray
        Year.
    ut : float or ndarray
        Day fraction.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.mjd2dayno(51603.0)
    >>> (60, 2000, 0.0)
    >>> pyacs.mjd2dayno(51603.5)
    >>> (60, 2000, 0.5)
    """
    
    from .mjd2cal import mjd2cal
    from .cal2dayno import cal2dayno
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
   
    if isinstance(mjd, np.ndarray):
        [dayno,year,ut]=np.array(list(map(mjd2dayno,mjd))).T
    
    else:
        
        (day, month, year, ut) = mjd2cal(mjd)
        dayno=cal2dayno(day,month,year)
    return  dayno, year, ut
