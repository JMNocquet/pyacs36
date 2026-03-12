"""
Convert day of year and year to decimal year.
"""

def dayno2decyear(dayno,year,ut=0.5):
    """
    Convert day of year and year to decimal year.

    Parameters
    ----------
    dayno : int or array-like
        Day of year.
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction. Default 0.5 (12:00).

    Returns
    -------
    float or ndarray
        Decimal year.

    Notes
    -----
    Accounts for leap years (366 days). Unless ut is specified, the result
    is at 12:00 (ut=0.5).

    Examples
    --------
    >>> import pyacs
    >>> pyacs.dayno2decyear(60,1999)
    >>> 1999.16301369863
    >>> pyacs.dayno2decyear(60,2000)
    >>> 2000.1625683060108
    >>> pyacs.dayno2decyear(60,2001)
    >>> 2001.16301369863
    >>> pyacs.dayno2decyear(60,2001,ut=0.)
    >>> 2001.1616438356164
    """
    
    from ._common import __daynoOK, __utOK
    from .cal2dayno import cal2dayno
    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)

    
    if isinstance(dayno, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=dayno*0.0+ut
        decyear=np.array(list(map(dayno2decyear,dayno, year, ut)))

    else:
    
        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
        __utOK(ut) 
    
        nday=cal2dayno(31,12,year)
        dayno=dayno-1
        
        dayno=dayno+ut

        decyear = float(year)+float(dayno)/float(nday)
    
    return decyear
