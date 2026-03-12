"""
Convert 4-digit year to 2-digit year.
"""

def year2yr(year):
    """
    Convert 4-digit year to 2-digit year.

    Parameters
    ----------
    year : int or array-like
        4-digit year (e.g. 2001).

    Returns
    -------
    int or ndarray
        2-digit year (e.g. 01). year-2000 if year>1999, else year-1900.

    Notes
    -----
    Heritage from geodesy convention (GPS in 1980). See examples.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.year2yr(1996)
    >>> 96
    >>> pyacs.year2yr(2011)
    >>> 11
    """

    import numpy as np

    if isinstance(year, list):
        year=np.array(year)

    if isinstance(year, np.ndarray):
        yr=np.array(list(map(year2yr,year)))

    else:
    
        
        if year > 1999:yr=year-2000
        else:yr=year-1900

    return(yr)
