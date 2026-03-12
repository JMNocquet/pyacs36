import numpy as np


###############################################################################
def azimuth_to_en(azimuth):
###############################################################################
    """Convert azimuth (clockwise from north) to east and north unit vector.

    Parameters
    ----------
    azimuth : float
        Azimuth in decimal degrees.

    Returns
    -------
    east : float
        East component of unit vector.
    north : float
        North component of unit vector.

    Notes
    -----
    Calculation uses a spherical approximation.
    """

    r_az = np.radians(azimuth)
    return (np.sin(r_az), np.cos(r_az))


