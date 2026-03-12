"""
Convert datetime to seconds since 1980-01-01 00:00:00.
"""

def datetime2seconds( np_datetime ):
    """
    Convert datetime to seconds since 1980-01-01 00:00:00.

    Parameters
    ----------
    np_datetime : datetime or array-like
        Datetime instance(s).

    Returns
    -------
    ndarray
        1D array of seconds (int64).
    """


    # import
    import numpy as np
    from datetime import datetime, timedelta

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)


    if not (isinstance(np_datetime, np.ndarray)):
        np_datetime = np.array([np_datetime])


    # convert to timedelta and then to seconds

    np_dates_s = np.array(list(map(int, [x.total_seconds() for x in (np_datetime - ref_date_time)])), dtype=np.int64)

    # return

    return (np_dates_s)
