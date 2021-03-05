###################################################################
def rake_from_euler(longitude, latitude, strike, dip, euler):
###################################################################
    """
    predicts rake for a given fault according to an Euler pole, position and strike of the fault

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

    M = pyacs.lib.gmtpoint.GMT_Point(code='XXXX', lon=longitude, lat=latitude, he=0., Ve=0., Vn=0., Vu=0., \
                                     SVe=0, SVn=0, SVu=0, SVen=0, \
                                     Cv_xyz=None, Cv_enu=None, index=None)

    _tmp, elon, elat, ew, motion_type = euler.split('/')

    W = np.array(list(map(float, [elon, elat, ew])))

    N = M.pole(W=W, type_euler='euler', option='predict')

    (ve, vn) = (N.Ve, N.Vn)

    rake = pyacs.lib.faultslip.v_to_rake(ve, vn, strike, dip, motion_type)

    return (rake)
