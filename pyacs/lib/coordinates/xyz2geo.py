#from pyacs.lib.errors import OptionError
import numpy as np


###################################################################
def xyz2geo(x, y, z, A=6378137.0, E2=0.006694380022903, unit="radians"):
###################################################################
    """Convert geocentric cartesian (XYZ) to geodetic coordinates (lon, lat, h).

    Parameters
    ----------
    x : float or array_like
        X coordinate(s) in meters.
    y : float or array_like
        Y coordinate(s) in meters.
    z : float or array_like
        Z coordinate(s) in meters.
    A : float, optional
        Semi-major axis (equatorial radius), meters. Default 6378137.0 (GRS80).
    E2 : float, optional
        Squared eccentricity. Default 0.006694380022903 (GRS80).
    unit : {'radians', 'dec_deg'}, optional
        Output units for longitude and latitude. Default is 'radians'.

    Returns
    -------
    long : float or ndarray
        Longitude in requested unit.
    lat : float or ndarray
        Latitude in requested unit.
    he : float or ndarray
        Height above ellipsoid in meters.

    Notes
    -----
    Default ellipsoid is GRS80 (WGS84): A=6378137 m, E2=0.006694380022903,
    flattening F = 1 - sqrt(1-E2). Reference: Bowring (1985), Survey Review.
    """    
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
   
    if unit not in ["radians", "dec_deg"]:
        ERROR("unit option must be in [radians,dec_deg];unit=", unit, exit=True)

    F = 1.0 - np.sqrt(1 - E2)
    
    TP = np.sqrt(x ** 2 + y ** 2)
    R = np.sqrt(TP ** 2 + z ** 2)
 
    TMU = np.arctan2((z / TP) * ((1.0 - F) + (E2 * A) / R), 1)
    RLBDA = np.arctan2(y, x)

    S3 = (np.sin(TMU)) ** 3
    C3 = (np.cos(TMU)) ** 3
    T1 = z * (1 - F) + E2 * A * S3
    T2 = (1 - F) * (TP - E2 * A * C3)


    RPHI = np.arctan2(T1, T2)
    RHE = TP * (np.cos(RPHI)) + z * (np.sin(RPHI))
    RHE = RHE - A * (np.sqrt(1 - E2 * (np.sin(RPHI) ** 2)))

    if unit == "dec_deg":
        RLBDA = np.degrees(RLBDA)
        RPHI = np.degrees(RPHI)

    return (RLBDA, RPHI, RHE)


