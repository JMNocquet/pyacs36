"""
Convert modified Julian day to decimal year.
"""

def mjd2decyear(mjd):
    """
    Convert modified Julian day to decimal year.

    Parameters
    ----------
    mjd : float or array-like
        Modified Julian day.

    Returns
    -------
    float or ndarray
        Decimal year.
    """
    
    from .mjd2dayno import mjd2dayno
    from .dayno2decyear import dayno2decyear
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
 
    if isinstance(mjd, np.ndarray):
        decyear=np.array(list(map(mjd2decyear,mjd)))
    
    else:
        (dayno,year,ut)=mjd2dayno(mjd)
        decyear = dayno2decyear(dayno,year,ut=ut)
    
    return decyear
