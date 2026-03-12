"""
Convert decimal year to day of year and day fraction.
"""

def decyear2dayno(decyear):
    """
    Convert decimal year to day of year and day fraction.

    Parameters
    ----------
    decyear : float or array-like
        Decimal year.

    Returns
    -------
    noday : int or ndarray
        Day of year.
    ut : float or ndarray
        Day fraction.

    Notes
    -----
    The input decimal year is assumed to account for leap years; the day of
    year is the decimal part of decyear times the number of days in the year.
    """

    from .cal2dayno import cal2dayno
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)
    
    if isinstance(decyear, np.ndarray):
        [noday,ut]=np.array(list(map(decyear2dayno,decyear))).T

    else:
    
        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)

    return( np.array(noday, dtype=int)+0,ut)
