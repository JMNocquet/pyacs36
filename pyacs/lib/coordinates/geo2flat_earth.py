from pyproj import Transformer
import numpy as np


###################################################################
def geo2flat_earth(longitude, latitude):
##################################################################
    """Convert geographical coordinates to Web Mercator (flat earth).

    Uses pyproj with EPSG:4326 (WGS84) to EPSG:3857 (Web Mercator).

    Parameters
    ----------
    longitude : float or array_like
        Longitude in decimal degrees.
    latitude : float or array_like
        Latitude in decimal degrees.

    Returns
    -------
    x : ndarray
        Easting in km.
    y : ndarray
        Northing in km.
    """

    longitude = np.array(longitude)
    latitude = np.array(latitude)

    TRAN_4326_TO_3857 = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    (x, y) = TRAN_4326_TO_3857.transform(longitude, latitude)
    return x * 1.0e-3, y * 1.0e-3


