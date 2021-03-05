###############################################################################
def strike_dip_rake_to_dir(strike,dip,rake):
###############################################################################
    """
    for a given (strike,dip,rake)  returns the azimuth of the horizontal motion
    :param strike,dip,rake: in decimal degrees
    :returns: direction: slip direction (modulo 180 degrees)
    """


    import numpy as np
    r_rake=np.radians(rake)
    r_dip=np.radians(dip)
    dir_slip=strike-np.degrees(np.arctan(np.tan(r_rake)*np.cos(r_dip)))

    if dir_slip > 180.:dir_slip=dir_slip-180.0
    if dir_slip > 180.:dir_slip=dir_slip-180.0
    if dir_slip<0:dir_slip=dir_slip+180.0

    return(dir_slip)
