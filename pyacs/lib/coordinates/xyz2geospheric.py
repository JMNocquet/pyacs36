from pyacs.lib.errors import OptionError
import numpy as np


###################################################################
def xyz2geospheric(x, y, z, unit="radians"):
###################################################################
    """Convert geocentric XYZ to geo-spherical (longitude, latitude, radius).

    Parameters
    ----------
    x : float or array_like
        X coordinate in meters.
    y : float or array_like
        Y coordinate in meters.
    z : float or array_like
        Z coordinate in meters.
    unit : {'radians', 'dec_deg'}, optional
        Units for longitude and latitude. Default is 'radians'.

    Returns
    -------
    longitude : float or ndarray
        Longitude in requested unit.
    latitude : float or ndarray
        Latitude in requested unit (geocentric, not co-latitude).
    R : float or ndarray
        Radius from Earth center in meters.

    Notes
    -----
    Latitude here is geocentric (angle from equatorial plane), not the usual
    spherical co-latitude.
    """

    if unit not in ["radians", "dec_deg"]:
        raise OptionError("unit option must be in [radians,dec_deg];unit=", unit)
    
    R = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    RLBDA = np.arctan2(y, x)
    
    Req = np.sqrt(x ** 2 + y ** 2)
    if Req != 0.0:
        RPHI = np.arctan(z / Req)
    else:
        RPHI = np.pi / 2.0 * z / np.sqrt(z ** 2)

    if unit == "dec_deg":
        RLBDA = np.degrees(RLBDA)
        RPHI = np.degrees(RPHI)
    
    return (RLBDA, RPHI, R)


