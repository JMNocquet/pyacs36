import numpy as np


###################################################################
def xyz_spherical_distance(x1, y1, z1, x2, y2, z2, Rt=6.371E6):
###################################################################
    """Spherical distance between two points in geocentric XYZ.

    Parameters
    ----------
    x1, y1, z1 : float
        First point in geocentric cartesian coordinates, meters.
    x2, y2, z2 : float
        Second point in geocentric cartesian coordinates, meters.
    Rt : float, optional
        Mean Earth radius in meters. Default 6.371e6.

    Returns
    -------
    float
        Distance in meters.
    """
 
    v1 = np.array([x1, y1, z1], dtype=float)
    v2 = np.array([x2, y2, z2], dtype=float)
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    sp = np.dot(v1, v2)
    sp = np.around(sp, decimals=13)
    return np.fabs(np.arccos(sp) * Rt)


