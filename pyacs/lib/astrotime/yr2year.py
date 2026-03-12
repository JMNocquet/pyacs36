"""
Convert two-digit year to four-digit year.
"""

def yr2year(yr):
    """
    Convert two-digit year to four-digit year.

    Parameters
    ----------
    yr : int or array-like
        2-digit year (e.g. 01).

    Returns
    -------
    int or ndarray
        4-digit year (e.g. 2001). yr+2000 if yr<80, else yr+1900.

    Notes
    -----
    Heritage from geodesy convention (GPS era). See examples.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.yr2year(96)
    >>> 1996
    >>> pyacs.yr2year(11)
    >>> 2011
    """

    import numpy as np

    if isinstance(yr, list):
        yr=np.array(yr)

    if isinstance(yr, np.ndarray):
        year=np.array(list(map(yr2year,yr)))

    else:
    
        if yr<80:year=yr+2000
        else:year=yr+1900

    return(year)
