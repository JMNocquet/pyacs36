"""Design matrix relating horizontal velocity to rotation rate vector."""

import numpy as np
import pyacs.lib.coordinates


def pole_matrix(coor):
    """Design matrix relating horizontal velocity to a rotation rate vector.

    For n sites with coordinates [lon, lat] (decimal degrees), returns matrix W
    such that ``np.dot(W, w)`` is the flattened [ve1, vn1, ve2, vn2, ...] in m/yr
    for rotation rate vector `w` in rad/yr.

    Parameters
    ----------
    coor : array_like, shape (n, 2)
        Site coordinates as [longitude, latitude] in decimal degrees.

    Returns
    -------
    ndarray, shape (2*n, 3)
        Pole matrix such that (W @ w) gives [ve1, vn1, ve2, vn2, ...] in m/yr.
    """
    coor = coor.reshape(-1, 2)

    pole_matrix_out = None

    for i in np.arange(coor.shape[0]):
        he = 0.0
        (x, y, z) = pyacs.lib.coordinates.geo2xyz(coor[i, 0], coor[i, 1], he, unit='dec_deg')
        (l, p, _) = pyacs.lib.coordinates.xyz2geospheric(x, y, z, unit='radians')
        R = pyacs.lib.coordinates.mat_rot_general_to_local(l, p)

        Ai = np.zeros([3, 3], float)
        Ai[0, 1] = z
        Ai[0, 2] = -y
        Ai[1, 0] = -z
        Ai[1, 2] = x
        Ai[2, 0] = y
        Ai[2, 1] = -x

        RAi = np.dot(R, Ai)[:2, :]

        if pole_matrix_out is None:
            pole_matrix_out = RAi
        else:
            pole_matrix_out = np.vstack((pole_matrix_out, RAi))

    return pole_matrix_out
