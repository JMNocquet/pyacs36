"""
Convert day of year and year to modified Julian day (JD - 2400000.5).
"""

def dayno2mjd(dayno, year, ut=0.5):
    """
    Convert day of year and year to modified Julian day (JD - 2400000.5).

    Parameters
    ----------
    dayno : int or array-like
        Day of year.
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction. Default 0.5 (middle of day).

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
    >>> pyacs.dayno2mjd(60,2000)
    >>> 51603.5
    >>> pyacs.dayno2mjd(60.,2000,ut=0)
    >>> 51603.0
    """

    from ._common import __daynoOK, __utOK
    from .dayno2cal import dayno2cal
    from .cal2mjd import cal2mjd
    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)
   
    if isinstance(dayno, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=dayno*0.0+ut
        mjd=np.array(list(map(dayno2mjd,dayno, year, ut)))

    else:

        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
        __utOK(ut)
    
        (day, month) = dayno2cal(dayno, year)
        mjd=cal2mjd(day, month, year, ut)
    
    return mjd
