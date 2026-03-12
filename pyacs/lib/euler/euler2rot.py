"""Convert Euler pole to cartesian rotation rate vector."""

import numpy as np


def euler2rot(lon, lat, omega):
    """Convert an Euler pole to a cartesian rotation rate vector.

    Converts longitude, latitude (degrees) and angular velocity (deg/Myr)
    into Wx, Wy, Wz in radians/yr.

    Parameters
    ----------
    lon : float
        Longitude of the Euler pole in decimal degrees.
    lat : float
        Latitude of the Euler pole in decimal degrees.
    omega : float
        Angular velocity in decimal degrees per Myr.

    Returns
    -------
    wx : float
        X component of rotation rate, radians/yr.
    wy : float
        Y component of rotation rate, radians/yr.
    wz : float
        Z component of rotation rate, radians/yr.

    Notes
    -----
    Longitude and latitude are relative to the sphere, not the ellipsoid.
    """
    wx = np.cos(np.radians(lat)) * np.cos(np.radians(lon)) * np.radians(omega) * 1.0E-6
    wy = np.cos(np.radians(lat)) * np.sin(np.radians(lon)) * np.radians(omega) * 1.0E-6
    wz = np.sin(np.radians(lat)) * np.radians(omega) * 1.0E-6

    return (wx, wy, wz)
