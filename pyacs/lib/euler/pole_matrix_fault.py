"""Design matrix relating fault slip components to rotation rate vector."""

import numpy as np
import pyacs.lib.coordinates


def pole_matrix_fault(coor, strike, order=None):
    """Design matrix relating fault slip components to a rotation rate vector.

    For n sites with [lon, lat] (decimal degrees) and strike (counter-clockwise
    from north), returns matrix W such that ``np.dot(W, w)`` is the flattened
    [ss1, ns1, ss2, ns2, ...] in m/yr for rotation rate vector `w` in rad/yr.

    Parameters
    ----------
    coor : array_like, shape (n, 2)
        Site coordinates as [longitude, latitude] in decimal degrees.
    strike : array_like, shape (n,) or (n, 1)
        Fault strike at each site, counter-clockwise from north, decimal degrees.
    order : {None, 'SS_DS'}, optional
        If provided, return matrix ordered as strike-slip then dip-slip rows.

    Returns
    -------
    ndarray, shape (2*n, 3)
        Pole matrix. ``np.dot(W, w).reshape(-1, 2)`` gives [along-strike, normal]
        components in two columns, in m/yr. Values refer to the hanging-wall
        block (right-hand rule polygon).
    """
    coor = coor.reshape(-1, 2)

    rstrike = np.radians(strike)

    if order is not None:
        pmf_SS = np.zeros((strike.shape[0], 3))
        pmf_DS = np.zeros((strike.shape[0], 3))

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

        RF = np.array([
            [np.sin(rstrike[i]), np.cos(rstrike[i])],
            [-np.cos(rstrike[i]), np.sin(rstrike[i])]
        ])

        RW = np.dot(RF, RAi)

        if order is not None:
            pmf_SS[i] = RW[0, :]
            pmf_DS[i] = RW[1, :]

        if pole_matrix_out is None:
            pole_matrix_out = RW
        else:
            pole_matrix_out = np.vstack((pole_matrix_out, RW))

    if order is not None:
        return np.vstack((pmf_SS, pmf_DS))
    else:
        return pole_matrix_out
