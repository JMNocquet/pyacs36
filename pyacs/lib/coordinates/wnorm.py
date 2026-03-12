import numpy as np


###################################################################
def wnorm(phi, A=6378137.0, E2=0.006694380022903):
###################################################################
    """
    Calculate the geodetic radius of curvature normal to the ellipsoid.

    Parameters
    ----------
    phi : float or ndarray
        Geodetic latitude in radians.
    A : float, optional
        Semi-major axis of the ellipsoid (default WGS84 in m).
    E2 : float, optional
        First eccentricity squared (default WGS84).

    Returns
    -------
    float or ndarray
        Radius of curvature (m).
    """

    a = A
    e2 = E2
  
    wd = np.sqrt(1.0 - e2 * np.sin(phi) ** 2) 
    result = a / wd
    
    return result


