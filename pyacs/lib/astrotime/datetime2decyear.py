"""
Convert Python datetime to decimal year.
"""

def datetime2decyear(datetime):
    """
    Convert Python datetime to decimal year.

    Parameters
    ----------
    datetime : datetime or array-like
        Datetime instance(s).

    Returns
    -------
    float or ndarray
        Decimal year.
    """

    from .datetime2mjd import datetime2mjd
    from .mjd2decyear import mjd2decyear

    return mjd2decyear( datetime2mjd(datetime) )
