"""Convert rotation rate vector to Euler pole in geographical coordinates."""

import numpy as np


def rot2euler(Wx, Wy, Wz):
    """Convert a rotation rate vector to an Euler pole in geographical coordinates.

    Converts Wx, Wy, Wz (radians/yr) into longitude, latitude (degrees) and
    angular velocity (deg/Myr).

    Parameters
    ----------
    Wx : float
        X component of rotation rate vector (geocentric cartesian), radians/yr.
    Wy : float
        Y component of rotation rate vector (geocentric cartesian), radians/yr.
    Wz : float
        Z component of rotation rate vector (geocentric cartesian), radians/yr.

    Returns
    -------
    W_long : float
        Longitude of the Euler pole in decimal degrees.
    W_lat : float
        Latitude of the Euler pole in decimal degrees.
    W_omega : float
        Angular velocity in decimal degrees per Myr.

    Notes
    -----
    Longitude and latitude are relative to the sphere, not the ellipsoid.
    Euler pole and rigid rotations only have sense on a sphere.
    """
    W = np.sqrt(Wx**2 + Wy**2 + Wz**2)
    W_lat = 90.0 - np.arccos(Wz / W) * 180.0 / np.pi

    if Wx > 0:
        W_long = np.arctan(Wy / Wx) * 180.0 / np.pi
    elif Wx < 0:
        if Wy > 0:
            W_long = np.arctan(Wy / Wx) * 180.0 / np.pi + 180.0
        else:
            W_long = np.arctan(Wy / Wx) * 180.0 / np.pi - 180.0
    else:
        if Wy > 0:
            W_long = 90.0
        else:
            W_long = -90.0

    W_omega = W / np.pi * 180 * 1.0E6

    return (W_long, W_lat, W_omega)
