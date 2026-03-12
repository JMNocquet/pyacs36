"""
Convert modified Julian day to a Python datetime instance.
"""

def mjd2datetime(mjd):
    """
    Convert modified Julian day to a Python datetime instance.

    Parameters
    ----------
    mjd : float or array-like
        Modified Julian day.

    Returns
    -------
    datetime or ndarray of datetime
        Datetime instance(s).
    """

    from .mjd2dayno import mjd2dayno
    from .dayno2datetime import dayno2datetime
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
  
    if isinstance(mjd, np.ndarray):
        return_datetime=np.array(list(map(mjd2datetime,mjd)))
                                 
    else:
    
        (dayno,year,ut)=mjd2dayno(mjd)
        return_datetime=dayno2datetime(dayno,year,ut=ut)
    
    return return_datetime
