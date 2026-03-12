"""Euler pole prediction for GMT_Point."""

import numpy as np
import pyacs.lib.coordinates as Coordinates
import pyacs.lib.euler


def pole(self, W=None, SW=None, type_euler='rot', option='predict'):
    """Apply or remove an Euler pole prediction to this point's velocity.

    Parameters
    ----------
    W : array_like
        Euler pole: [lon, lat, omega] (deg/Myr) or [wx, wy, wz] (rad/yr).
    SW : array_like, optional
        Variance-covariance of pole (3x3). If provided, updates SVe, SVn.
    type_euler : str, optional
        'euler' (geographic) or 'rot' (cartesian rad/yr). Default is 'rot'.
    option : str, optional
        'predict' (velocity = pole only), 'remove' (subtract), 'add'. Default is 'predict'.

    Returns
    -------
    GMT_Point
        New point with modified velocity.
    """
    ltype = ['euler', 'rot']
    if type_euler not in ltype:
        raise ValueError(
            "!!! Error. type_euler should be either euler " +
            "or rot (Euler pole (long. lat in dec. degrees and deg/Myr)" +
            " or cartesian rotation rate vector in rad/yr. type_euler = %s )"
            % type_euler)

    loption = ['remove', 'predict', 'add']
    if option not in loption:
        raise ValueError(
            "!!! Error. option should be either predict, remove or  add. option = %s ",
            option)

    if type_euler == 'euler':
        if len(W.shape) != 2:
            W = W.reshape(3, 1)

        (lon, lat, omega) = (W[0, 0], W[1, 0], W[2, 0])
        (wx, wy, wz) = pyacs.lib.euler.euler2rot(lon, lat, omega)

        W = np.zeros((3, 1))
        W[0, 0] = wx
        W[1, 0] = wy
        W[2, 0] = wz

    (x, y, z) = Coordinates.geo2xyz(self.lon, self.lat, self.he, unit='dec_deg')
    (l, p, _) = Coordinates.xyz2geospheric(x, y, z, unit='radians')
    R = Coordinates.mat_rot_general_to_local(l, p)

    Ai = np.zeros([3, 3], float)
    Ai[0, 1] = z
    Ai[0, 2] = -y
    Ai[1, 0] = -z
    Ai[1, 2] = x
    Ai[2, 0] = y
    Ai[2, 1] = -x

    RAi = np.dot(R, Ai)

    Pi = np.dot(RAi, W)

    pve = Pi[0, 0] * 1.E3
    pvn = Pi[1, 0] * 1.E3

    N = self.copy()

    if option == 'remove':
        N.Ve = self.Ve - pve
        N.Vn = self.Vn - pvn
    if option == 'add':
        N.Ve = self.Ve + pve
        N.Vn = self.Vn + pvn
    if option == 'predict':
        N.Ve = pve
        N.Vn = pvn

    if SW is not None:
        VCV = SW
        VCV_PREDICTION_ENU = np.dot(np.dot(RAi, VCV), RAi.T)
        (N.SVe, N.SVn) = (np.sqrt(VCV_PREDICTION_ENU[0, 0]) * 1000.0, np.sqrt(VCV_PREDICTION_ENU[1, 1]) * 1000.0)

    return N
