"""
Convert JPL date to datetime object (seconds since 2000-01-01 12:00 TAI).
"""

def jpldate_2_datetime(jpldate):
    """
    Convert JPL date to datetime object.

    JPL dates are seconds since 2000-01-01 12:00 (TAI).

    Parameters
    ----------
    jpldate : float or array-like
        JPL date (seconds since 2000-01-01 12:00 TAI).

    Returns
    -------
    datetime or ndarray
        Datetime instance(s).
    """

    from datetime import datetime
    from .datetime2seconds import datetime2seconds
    from .seconds2datetime import seconds2datetime
    import numpy as np

    # check jpldate type
    if isinstance(jpldate, list):
        jpldate = np.array(jpldate)

    zero_jpl_datetime = datetime(2000, 1, 1, 12, 0)
    zero_jpl_in_pyacs_seconds = datetime2seconds(zero_jpl_datetime)
    jpldate_in_pyacs_seconds = zero_jpl_in_pyacs_seconds + jpldate

    return seconds2datetime(jpldate_in_pyacs_seconds)
