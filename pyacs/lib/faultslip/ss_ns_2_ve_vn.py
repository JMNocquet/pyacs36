###################################################################
def ss_ns_2_ve_vn(ss, ns, strike):
###################################################################
    """
    Converts strike-slip and normal slip components of a fault to ve, vn

    :param ss: vector of strike-slip
    :param ns: vector of normal-slip
    :param strike: vector of fault strike, counter-clockwise from north in degrees

    :return: ve, vn, same unit as ss and ds

    :note:
    ss is assume to be positive for left-lateral motion
    ds will be positive for reverse motion
    """

    import numpy as np

    # strike in radians
    rstrike = np.radians(strike)
    # rcos and rsin
    rcos = np.cos(rstrike)
    rsin = np.sin(rstrike)

    return ss * rsin - ns * rcos, ss * rcos + ns * rsin

