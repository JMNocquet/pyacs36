###################################################################
def rake_from_slip_az(strike, dip, slipdir , motion_type ):
###################################################################
    """Predict rake from horizontal slip direction and motion style.

    Input parameters can be single float or 1D numpy array.

    Parameters
    ----------
    strike : float or numpy.ndarray
        Fault strike from north in decimal degrees.
    dip : float or numpy.ndarray
        Fault dip from north in decimal degrees.
    slipdir : float or numpy.ndarray
        Slip direction (azimuth) in decimal degrees.
    motion_type : str
        Motion style (e.g. 'inverse', 'normal', 'leftlateral', 'rightlateral').

    Returns
    -------
    float or numpy.ndarray
        Rake in decimal degrees.

    Notes
    -----
    Style needs to be provided to ensure the correct sense of slip.
    """

    import pyacs.lib.gmtpoint
    import numpy as np
    import pyacs.lib.faultslip


    (ve, vn) = (np.sin(np.radians(slipdir)) , np.cos(np.radians(slipdir)) )

    rake = pyacs.lib.faultslip.v_to_rake(ve, vn, strike, dip, motion_type)

    return (rake)
