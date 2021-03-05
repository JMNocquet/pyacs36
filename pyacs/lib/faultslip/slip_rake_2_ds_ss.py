###################################################################
def slip_rake_2_ds_ss(slip, rake):
###################################################################
    """
    Converts a slip vector provided as slip and rake to dip slip and strike-slip components

    :param slip: slip
    :param rake: rake in degrees

    :return ( slip*np.sin( np.radians(rake) ) ,  slip*np.cos( np.radians(rake) ))

    :note:
    ds will be positive for rake in [0,180] that is reverse motion
    ss will be positive for left-lateral motion

    """

    import numpy as np

    return (slip * np.sin(np.radians(rake)), slip * np.cos(np.radians(rake)))
