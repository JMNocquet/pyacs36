"""
Convert decimal year to modified Julian day.
"""

def decyear2mjd(decyear):
    """
    Convert decimal year to modified Julian day.

    Parameters
    ----------
    decyear : float or array-like
        Decimal year.

    Returns
    -------
    float or ndarray
        Modified Julian day.
    """

    from .cal2dayno import cal2dayno
    from .dayno2mjd import dayno2mjd
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)
    
    if isinstance(decyear, np.ndarray):
        mjd=np.array(list(map(decyear2mjd,decyear)))

    else:

        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)
        mjd=dayno2mjd(noday,year,ut)
        
    return mjd
