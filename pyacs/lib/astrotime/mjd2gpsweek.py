"""
Convert modified Julian day to GPS week and GPS day of week.
"""

def mjd2gpsweek(mjd):
    """
    Convert modified Julian day to GPS week and GPS day of week.

    Parameters
    ----------
    mjd : float or array-like
        Modified Julian day.

    Returns
    -------
    gps_week : int or ndarray
        GPS week number.
    day_of_week : int or ndarray
        Day of week (0-6).
    """
    
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
    
    if isinstance(mjd, np.ndarray):
        [gweek,dow]=np.array(list(map(mjd2gpsweek,mjd))).T
    
    else:
    
        mjd010580 = 44243
        mjd_d = ( int(mjd) - mjd010580 ) - 1
        gweek = mjd_d // 7
        dow = mjd_d % 7

    return(gweek,dow)
