"""
Convert decimal year to a Python datetime instance.
"""

def decyear2datetime(decyear):
    """
    Convert decimal year to a Python datetime instance.

    Parameters
    ----------
    decyear : float or array-like
        Decimal year.

    Returns
    -------
    datetime or ndarray of datetime
        Datetime instance(s).
    """

    from .decyear2dayno import decyear2dayno
    from .dayno2datetime import dayno2datetime
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array( decyear )

    if isinstance(decyear, np.ndarray):
        return_datetime=np.array(list(map(decyear2datetime,decyear)))

    else:

        year=int(decyear)
        (dayno,ut)=decyear2dayno(decyear)
        return_datetime=dayno2datetime(dayno,year,ut=ut)

    return return_datetime
