"""
Convert JPL date to decimal year (seconds since 2000-01-01 12:00 TAI).
"""

def jpldate_2_decyear(jpldate):
    """
    Convert JPL date to decimal year.

    JPL dates are seconds since 2000-01-01 12:00 (TAI).

    Parameters
    ----------
    jpldate : float or array-like
        JPL date (seconds since 2000-01-01 12:00 TAI).

    Returns
    -------
    float or ndarray
        Decimal year.
    """

    from .jpldate_2_datetime import jpldate_2_datetime
    from .datetime2decyear import datetime2decyear

    return datetime2decyear(jpldate_2_datetime(jpldate))
