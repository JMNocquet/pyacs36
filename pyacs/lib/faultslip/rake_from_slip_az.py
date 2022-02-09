###################################################################
def rake_from_slip_az(strike, dip, slipdir , motion_type ):
###################################################################
    """
    predicts rake for a given fault from the horizontal slip direction and motion style
    input parameters can be single float or 1D numpy array

    :param longitude,latitude: in decimal degrees
    :param strike: fault strike from north in decimal degrees
    :param dip: fault dip from north in decimal degrees
    :param euler: Euler pole as a string '/long/lat/w/style' (style among 'inverse', 'normal', 'leftlateral','rightlateral')

    :return rake: in decimal degrees
    :note: style needs to be provided to ensure the correct sense of slip.
    """

    import pyacs.lib.gmtpoint
    import numpy as np
    import pyacs.lib.faultslip


    (ve, vn) = (np.sin(np.radians(slipdir)) , np.cos(np.radians(slipdir)) )

    rake = pyacs.lib.faultslip.v_to_rake(ve, vn, strike, dip, motion_type)

    return (rake)
