"""
Convert decimal year to calendar date.
"""

def decyear2cal(decyear):
    """
    Convert decimal year to calendar date.

    Parameters
    ----------
    decyear : float or array-like
        Decimal year.

    Returns
    -------
    mday : int or ndarray
        Day of month.
    month : int or ndarray
        Month (1-12).
    ut : float or ndarray
        Day fraction.

    Notes
    -----
    The input decimal year is assumed to account for leap years; the day of
    year is the decimal part of decyear times the number of days in the year.
    """

    from .cal2dayno import cal2dayno
    from .dayno2cal import dayno2cal
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)

    if isinstance(decyear, np.ndarray):
        [mday,month,ut]=np.array(list(map(decyear2cal,decyear))).T
    
    else:

        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)
        (mday,month)=dayno2cal(noday,year)

    # this trick ensures that a single decyear provided as float
    # returns int,int,float for monthday, month and ut
    # while a list or np.array of decyear will return
    # int 1D np array for monthday and month and float 1D np array for ut
    
    return( np.array(mday,dtype=int)+0,
            np.array(month,dtype=int)+0,
            ut)
