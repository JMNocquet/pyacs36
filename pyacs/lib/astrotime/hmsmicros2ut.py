"""
Convert hours, minutes, seconds, microseconds to fractional day.
"""

def hmsmicros2ut(h,mn,s=0,microsecond=0):
    """
    Convert hours, minutes, seconds, microseconds to fractional day.

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
        Fractional day (ut) in [0., 1.[.
    """

    from .hmsmicros2uts import hmsmicros2uts
    from .uts2ut import uts2ut

    uts=hmsmicros2uts(h,mn,s=s,microsecond=microsecond)

    return(uts2ut(uts))
