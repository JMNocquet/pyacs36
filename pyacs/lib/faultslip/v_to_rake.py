###############################################################################
def v_to_rake(ve, vn, strike, dip, style='inverse'):
###############################################################################

    """
    for a given relative horizontal velocity (ve,vn) between two blocks separated by a fault of a given strike & dip
    , returns the expected rake for the eq. Input can be single float or 1D numpy array

    :param ve,vn: east and north components of motion (any unit)
    :param strike,dip: in decimal degrees
    :param style: string among {'inverse', 'normal', 'leftlateral','rightlateral'}
    :returns rake: in decimal degrees

    :note: Because only providing ve,vn is ambiguous (we don't know whether if ve,vn is hanging wall motion wrt the footwall or the opposite) \
    an additional information (style) in the form of one of {'inverse', 'normal', 'leftlateral','rightlateral'} must be provided.
    """

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np

    ###############################################################################
    # INPUT AS 1D NUMPY ARRAY
    ###############################################################################

    if isinstance(ve,np.ndarray) and ve.ndim ==1:
        rake = np.zeros( ve.shape[0] )
        for i in np.arange( ve.shape[0] ):
            rake[i] = v_to_rake(ve[i], vn[i], strike[i], dip[i], style=style[i])

    else:

    ###############################################################################
    # INPUT AS FLOAT
    ###############################################################################

        if dip == 90.:
            import sys
            print("!!! ERROR: dip is 90 degrees. Cannot calculate rake")
            sys.exit()

        r_dip = np.radians(dip)
        ralpha = np.radians(90. - strike)
        v = np.array([ve, vn])
        v_norm = v / np.sqrt(v[0] ** 2 + v[1] ** 2)

        R = np.array([[np.cos(ralpha), -np.sin(ralpha)], [np.sin(ralpha), np.cos(ralpha)]])

        VRAKE = np.dot(R.T, v_norm)
        VRAKE[1] = VRAKE[1] / np.cos(r_dip)

        rake_1 = np.degrees(np.arctan2(VRAKE[1], VRAKE[0]))

        rake_2 = rake_1 + 180.0
        if rake_2 > 180: rake_2 = rake_2 - 360.

        rakes = np.array([rake_1, rake_2])

        if style is not None:

            (rake_negative, rake_positive) = (np.min(rakes), np.max(rakes))
            if np.abs(rake_1) <= 90:
                rake_leftlateral = rake_1;
                rake_rightlateral = rake_2
            else:
                rake_leftlateral = rake_2;
                rake_rightlateral = rake_1

            if style == 'inverse': rake = rake_positive
            if style == 'normal': rake = rake_negative
            if style == 'leftlateral': rake = rake_leftlateral
            if style == 'rightlateral': rake = rake_rightlateral

        else:
            rake = rake_1

    return (rake)
