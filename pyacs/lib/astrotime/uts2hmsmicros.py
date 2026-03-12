"""
Convert seconds since 00:00:00.0 (0-86400) to hours, minutes, seconds, microseconds.
"""

def uts2hmsmicros(uts):
    """
    Convert seconds since 00:00:00.0 (0-86400) to hours, minutes, seconds, microseconds.

    Parameters
    ----------
    uts : float or array-like
        Seconds since 00:00:00.0, in [0., 86400.[.

    Returns
    -------
    h : int or ndarray
        Hours.
    mn : int or ndarray
        Minutes.
    s : int or ndarray
        Seconds.
    microsecond : float or ndarray
        Microseconds.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.uts2hmsmicros(86399.999999)
    >>> (23, 59, 59, 999999.0000069374)
    >>> pyacs.uts2hmsmicros(86400)
    >>> !!!  86400 uts out of range [0- 86400.0 [
    """

    from ._common import __utsOK
    import numpy as np
    
    if isinstance(uts, list):
        uts=np.array(uts)

    if isinstance(uts, np.ndarray):
        [h,mn,s,microsecond]=np.array(list(map(uts2hmsmicros,uts))).T

    else:

        __utsOK(uts)
        
        nsec_hour= 60 * 60
        nsec_min= 60
        
        int_sec=int(uts)
        
        h = int_sec // nsec_hour
        r = int_sec % nsec_hour

        mn = r // nsec_min
        s =  r % nsec_min

        microsecond= (uts - float(int_sec)) * 1.E6

        return( np.array(h, dtype=int)+0,\
                np.array(mn, dtype=int)+0,\
                np.array(s,dtype=int)+0,\
                microsecond)
