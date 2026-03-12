###############################################################################
def v_to_n_ss(ve, vn, strike):
###############################################################################
    """Decompose horizontal velocity into normal and strike-slip components.

    Point is in the left-domain with respect to the fault. Shortening and
    right-lateral are positive; extension and left-lateral are negative.

    Parameters
    ----------
    ve : float or array-like
        East component of motion (any unit).
    vn : float or array-like
        North component of motion (any unit).
    strike : float or array-like
        Fault strike in decimal degrees.

    Returns
    -------
    normal : float or numpy.ndarray
        Normal component, same units as ve, vn.
    strike_slip : float or numpy.ndarray
        Strike-slip component, same units as ve, vn.
    """

    import numpy as np

    azimuth_radian = np.arctan2(ve, vn)
    radian_strike = np.radians(strike)

    angle = azimuth_radian - radian_strike

    v_i = np.sqrt(ve ** 2 + vn ** 2)
    ss = np.cos(angle) * v_i
    normal = np.sin(angle) * v_i

    return (normal, ss)
