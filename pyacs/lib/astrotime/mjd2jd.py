"""
Convert modified Julian day to Julian day.
"""

def mjd2jd(mjd):
    """
    Convert modified Julian day to Julian day.

    Parameters
    ----------
    mjd : float
        Modified Julian day.

    Returns
    -------
    float
        mjd + 2400000.5.
    """
    
    return mjd+2400000.5
