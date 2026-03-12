"""
Convert hours, minutes, seconds, microseconds to seconds since 00:00:00.0.
"""

def hmsmicros2uts(h,mn,s=0,microsecond=0):
    """
    Convert hours, minutes, seconds, microseconds to seconds since 00:00:00.0.

    Parameters
    ----------
    h : int
        Hours.
    mn : int
        Minutes.
    s : int, optional
        Seconds (default 0).
    microsecond : int or float, optional
        Microseconds (default 0).

    Returns
    -------
    float
        Seconds since 00:00:00.0 (uts).
    """

    nsec_hour=60 * 60
    nsec_min=60

    return((h*nsec_hour+mn*nsec_min+s)+(microsecond)/1.E6)
