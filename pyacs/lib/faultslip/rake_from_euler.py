###################################################################
def rake_from_euler(longitude, latitude, strike, dip, euler):
###################################################################
    """Predict rake for a fault from Euler pole, position and fault strike/dip.

    Parameters
    ----------
    longitude : float or numpy.ndarray
        Longitude in decimal degrees.
    latitude : float or numpy.ndarray
        Latitude in decimal degrees.
    strike : float or numpy.ndarray
        Fault strike from north in decimal degrees.
    dip : float or numpy.ndarray
        Fault dip from north in decimal degrees.
    euler : str
        Euler pole as '/long/lat/w/style' (style: 'inverse', 'normal', 'leftlateral', 'rightlateral').

    Returns
    -------
    float or numpy.ndarray
        Rake in decimal degrees.

    Notes
    -----
    Style needs to be provided to ensure the correct sense of slip.
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyacs.lib.gmtpoint
    import numpy as np
    import pyacs.lib.faultslip


    ###############################################################################
    # INPUT AS 1D NUMPY ARRAY
    ###############################################################################

    if isinstance(longitude, np.ndarray) and longitude.ndim == 1:
        rake = np.zeros(longitude.shape[0])
        for i in np.arange(longitude.shape[0]):
            rake[i] = pyacs.lib.faultslip.rake_from_euler(longitude[i], latitude[i], strike[i], dip[i], euler)

    else:

    ###############################################################################
    # INPUT AS FLOAT
    ###############################################################################

        M = pyacs.lib.gmtpoint.GMT_Point(code='XXXX', lon=longitude, lat=latitude, he=0., Ve=0., Vn=0., Vu=0., \
                                         SVe=0, SVn=0, SVu=0, SVen=0, \
                                         Cv_xyz=None, Cv_enu=None, index=None)

        _tmp, elon, elat, ew, motion_type = euler.split('/')

        W = np.array(list(map(float, [elon, elat, ew])))

        N = M.pole(W=W, type_euler='euler', option='predict')

        (ve, vn) = (N.Ve, N.Vn)

        rake = pyacs.lib.faultslip.v_to_rake(ve, vn, strike, dip, motion_type)

    return (rake)
