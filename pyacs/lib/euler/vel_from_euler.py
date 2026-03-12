"""Horizontal velocity at a point predicted from an Euler pole."""

import numpy as np
import pyacs.lib.coordinates

from .euler2rot import euler2rot


def vel_from_euler(lonp, latp, lon_euler, lat_euler, omega_euler, use_spherical=False):
    """Return horizontal velocity predicted at a point from an Euler pole.

    Parameters
    ----------
    lonp : float
        Longitude of the point where velocity is predicted, decimal degrees.
    latp : float
        Latitude of the point where velocity is predicted, decimal degrees.
    lon_euler : float
        Longitude of the Euler pole, decimal degrees.
    lat_euler : float
        Latitude of the Euler pole, decimal degrees.
    omega_euler : float
        Angular velocity of the Euler pole, decimal degrees per Myr.
    use_spherical : bool, optional
        If True, compute site position with a spherical Earth (radius 6371 km)
        so that pole and site use the same geometry. Default False (ellipsoidal).

    Returns
    -------
    ve : float
        East velocity at (lonp, latp), mm/yr.
    vn : float
        North velocity at (lonp, latp), mm/yr.
    """
    (wx, wy, wz) = euler2rot(lon_euler, lat_euler, omega_euler)

    W = np.zeros((3, 1))
    W[0, 0] = wx
    W[1, 0] = wy
    W[2, 0] = wz

    if use_spherical:
        R_earth = 6371000.0  # m, mean radius
        lat_rad = np.radians(latp)
        lon_rad = np.radians(lonp)
        x = R_earth * np.cos(lat_rad) * np.cos(lon_rad)
        y = R_earth * np.cos(lat_rad) * np.sin(lon_rad)
        z = R_earth * np.sin(lat_rad)
    else:
        he = 0.0
        (x, y, z) = pyacs.lib.coordinates.geo2xyz(lonp, latp, he, unit='dec_deg')
    (l, p, _) = pyacs.lib.coordinates.xyz2geospheric(x, y, z, unit='radians')
    R = pyacs.lib.coordinates.mat_rot_general_to_local(l, p)

    Ai = np.zeros([3, 3], float)
    Ai[0, 1] = z
    Ai[0, 2] = -y
    Ai[1, 0] = -z
    Ai[1, 2] = x
    Ai[2, 0] = y
    Ai[2, 1] = -x

    RAi = np.dot(R, Ai)

    V = np.dot(RAi, W) * 1.E3  # to get mm/yr

    return (V[0, 0], V[1, 0])
