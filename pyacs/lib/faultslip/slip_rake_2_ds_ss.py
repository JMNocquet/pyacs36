###################################################################
def slip_rake_2_ds_ss(slip, rake):
###################################################################
    """Convert slip magnitude and rake to dip-slip and strike-slip components.

    Parameters
    ----------
    slip : float or array-like
        Slip magnitude.
    rake : float or array-like
        Rake in degrees.

    Returns
    -------
    ds : float or numpy.ndarray
        Dip-slip component (positive for reverse, rake in [0, 180]).
    ss : float or numpy.ndarray
        Strike-slip component (positive for left-lateral).

    Notes
    -----
    ds positive for reverse motion; ss positive for left-lateral motion.
    """

    import numpy as np

    return (slip * np.sin(np.radians(rake)), slip * np.cos(np.radians(rake)))
