###############################################################################
def strike_dip_rake_to_dir(strike,dip,rake):
###############################################################################
    """
    for a given (strike,dip,rake)  returns the azimuth of the horizontal motion
    :param strike,dip,rake: in decimal degrees
    :returns: direction: slip direction (modulo 180 degrees)
    """


    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyacs.lib.faultslip.unit_slip

    ###############################################################################
    # INPUT AS 1D NUMPY ARRAY
    ###############################################################################

    if isinstance(strike, np.ndarray) and strike.ndim == 1:
        dir_slip = np.zeros(strike.shape[0])
        for i in np.arange(strike.shape[0]):
            dir_slip[i] = strike_dip_rake_to_dir(strike[i],dip[i],rake[i])

    else:

    ###############################################################################
    # INPUT AS FLOAT
    ###############################################################################

        # OLD & INCORRECT STUFF - CHANGE ON 03/04/2021
        # r_rake=np.radians(rake)
        # r_dip=np.radians(dip)
        # dir_slip=strike-np.degrees(np.arctan(np.tan(r_rake)*np.cos(r_dip)))
        #
        # if dir_slip > 180.:dir_slip=dir_slip-180.0
        # if dir_slip > 180.:dir_slip=dir_slip-180.0
        # if dir_slip<0:dir_slip=dir_slip+180.0

        (ue,un,uu) = pyacs.lib.faultslip.unit_slip( strike , dip, rake )
        dir_slip = np.degrees(np.arctan2( ue, un ))

    return(dir_slip)
