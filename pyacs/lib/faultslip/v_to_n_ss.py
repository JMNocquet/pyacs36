###############################################################################
def v_to_n_ss(ve, vn, strike):
###############################################################################
    """

    for a given relative horizontal velocity (ve,vn) between two blocks separated by a fault of a given strike
    , returns normal and strike-slip component of motion.
    The convention is that the point is in the left-domain with respect to the fault.
    Shortening & right-lateral are positive, extension and left-lateral negative

    :param ve,vn: east and north components of motion (any unit)
    :param strike: in decimal degrees
    :returns normal,strike-slip: same units as ve,vn

    """

    import numpy as np

    azimuth_radian = np.arctan2(ve, vn)
    radian_strike = np.radians(strike)

    angle = azimuth_radian - radian_strike

    v_i = np.sqrt(ve ** 2 + vn ** 2)
    ss = np.cos(angle) * v_i
    normal = np.sin(angle) * v_i

    return (normal, ss)
