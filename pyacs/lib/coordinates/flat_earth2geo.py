import numpy as np
from pyproj import Transformer


###################################################################
def flat_earth2geo(x, y):
##################################################################
    """Convert web Mercator coordinates (km) to geographical coordinates.

    Uses pyproj with EPSG:3857 (Web Mercator) to EPSG:4326 (WGS84).

    Parameters
    ----------
    x : float or array_like
        Easting in km (Web Mercator).
    y : float or array_like
        Northing in km (Web Mercator).

    Returns
    -------
    lon : ndarray
        Longitude in decimal degrees.
    lat : ndarray
        Latitude in decimal degrees.
    """

    x = np.array(x)
    y = np.array(y)

    TRAN_3857_TO_4326 = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    return TRAN_3857_TO_4326.transform(x * 1.0e3, y * 1.0e3)


