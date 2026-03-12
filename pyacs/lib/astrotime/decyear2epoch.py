"""
Convert decimal year to formatted epoch string 'yy:doy:sec' (%02d:%03d:%05d).
"""

def decyear2epoch(decyear):
    """
    Convert decimal year to formatted epoch string 'yy:doy:sec' (%02d:%03d:%05d).

    Parameters
    ----------
    decyear : float or array-like
        Decimal year.

    Returns
    -------
    str or ndarray
        Epoch string(s) 'yy:doy:sec'.

    Notes
    -----
    The input decimal year is assumed to account for leap years; the day of
    year is the decimal part of decyear times the number of days in the year.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.cal2decyear(1,3,2001)
    >>> 2001.16301369863
    >>> pyacs.decyear2epoch(2001.16301369863)
    >>> '01:060:43200'
    >>> pyacs.decyear2epoch(2000.16301369863)
    >>> '00:060:57284'
    >>> pyacs.decyear2epoch(pyacs.cal2decyear(1,3,2000))
    >>> '00:061:43200'
    """
    
    from .decyear2dayno import decyear2dayno
    from .year2yr import year2yr
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array( decyear )

    if isinstance(decyear, np.ndarray):
        epoch=np.array(list(map(decyear2epoch,decyear)))

    else:
    
        (doy,ut)=decyear2dayno(decyear)
        seconds=ut*(60*60*24)
        yr=year2yr(decyear)
        epoch=("%02d:%03d:%05d" %(yr,doy,round(seconds)))

    return(epoch)
