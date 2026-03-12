from pyacs.lib.errors import OptionError
import numpy as np


###################################################################
def geo2xyz(llambda, phi, he, unit="radians", A=6378137.0, E2=0.006694380022903):
###################################################################
    """Convert geodetic coordinates (lon, lat, h) to geocentric XYZ.

    Parameters
    ----------
    llambda : float or array_like
        Longitude.
    phi : float or array_like
        Latitude.
    he : float or array_like
        Ellipsoidal height in meters.
    unit : {'radians', 'dec_deg'}, optional
        Units for longitude and latitude. Default is 'radians'.
    A : float, optional
        Semi-major axis (equatorial radius), meters. Default 6378137.0 (GRS80).
    E2 : float, optional
        Squared eccentricity. Default 0.006694380022903 (GRS80).

    Returns
    -------
    x, y, z : float or ndarray
        Geocentric cartesian coordinates in meters.

    Notes
    -----
    Default ellipsoid is GRS80 (WGS84): A=6378137 m, E2=0.006694380022903,
    flattening F = 1 - sqrt(1-E2).
    """

    from .wnorm import wnorm

    if unit not in ["radians", "dec_deg"]:
        raise OptionError("unit option must be in [radians,dec_deg];unit=", unit)

    if unit == "dec_deg":
        llambda = np.radians(llambda)
        phi = np.radians(phi)

    a = A
    e2 = E2

    wn = wnorm(phi, A=a, E2=e2)
   
    xx = (wn + he) * np.cos(phi) * np.cos(llambda)
    yy = (wn + he) * np.cos(phi) * np.sin(llambda)
    zz = (wn * (1.0 - e2) + he) * np.sin(phi)
   
    return (xx, yy, zz)


