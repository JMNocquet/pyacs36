###################################################################
def ss_ns_2_ve_vn(ss, ns, strike):
###################################################################
    """Convert strike-slip and normal-slip components to east/north velocity.

    Parameters
    ----------
    ss : array-like
        Strike-slip component (positive for left-lateral).
    ns : array-like
        Normal-slip component (positive for reverse motion).
    strike : array-like
        Fault strike, counter-clockwise from north in degrees.

    Returns
    -------
    ve : numpy.ndarray
        East component, same unit as ss and ns.
    vn : numpy.ndarray
        North component, same unit as ss and ns.

    Notes
    -----
    ss is assumed positive for left-lateral motion; ns positive for reverse.
    """

    import numpy as np

    # strike in radians
    rstrike = np.radians(strike)
    # rcos and rsin
    rcos = np.cos(rstrike)
    rsin = np.sin(rstrike)

    return ss * rsin - ns * rcos, ss * rcos + ns * rsin

