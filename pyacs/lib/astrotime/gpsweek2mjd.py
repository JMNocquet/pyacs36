"""
Convert GPS week and day of week to modified Julian day.
"""

def gpsweek2mjd(gpsweek,dow):
    """
    Convert GPS week and day of week to modified Julian day.

    Parameters
    ----------
    gpsweek : int or array-like
        GPS week number.
    dow : int or array-like
        Day of week (0=Sunday, 6=Saturday).

    Returns
    -------
    float or ndarray
        Modified Julian day.
    """

    import numpy as np

    if isinstance(gpsweek, list):
        gpsweek=np.array(gpsweek)
        dow=np.array(dow)


    if isinstance(gpsweek, np.ndarray):
        mjd=np.array(list(map(gpsweek2mjd,gpsweek,dow))).T
    
    else:
    
        mjd010580 = 44243
        mjd = ( gpsweek * 7 ) + dow + mjd010580 + 1

    return( mjd )
