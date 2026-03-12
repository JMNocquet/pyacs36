"""
Unit tests for pyacs.lib.coordinates.

Reference values from GAMIT convertc/tform where noted.
"""

import numpy as np
import pytest

from pyacs.lib import coordinates as COO
from pyacs.lib import units as UNITS


# -----------------------------------------------------------------------------
# geo2xyz / xyz2geo (geodetic <-> geocentric)
# -----------------------------------------------------------------------------


def test_geo2xyz_dec_deg():
    """geo2xyz with decimal degrees; GAMIT convertc reference."""
    x, y, z = COO.geo2xyz(-112.567, 43.36, 27.2, unit="dec_deg")
    # GAMIT convertc: -1782430.03490 -4288974.18902 4356684.78382
    ref = (-1782430.03490, -4288974.18902, 4356684.78382)
    assert np.allclose((x, y, z), ref, atol=0.00002)


def test_xyz2geo_radians():
    """xyz2geo returns radians by default; GAMIT tform reference."""
    lon, lat, h = COO.xyz2geo(918129.451, -4346071.255, 4561977.839)
    # GAMIT tform (deg): 45.95580023 281.92863291 200.9022
    assert np.allclose(np.degrees(lon) % 360.0, 281.92863291, atol=1e-6)
    assert np.allclose(np.degrees(lat), 45.95580023, atol=1e-6)
    assert np.allclose(h, 200.9022, atol=0.0001)


def test_xyz2geo_dec_deg():
    """xyz2geo with unit='dec_deg' returns degrees (lon in [-180, 180])."""
    lon, lat, h = COO.xyz2geo(
        918129.451, -4346071.255, 4561977.839, unit="dec_deg"
    )
    # GAMIT tform 281.93° ≡ -78.07° (pyacs uses [-180, 180])
    assert np.allclose(lon, -78.07136709, atol=1e-6)
    assert np.allclose(lat, 45.95580023, atol=1e-6)
    assert np.allclose(h, 200.9022, atol=0.0001)


def test_geo2xyz_xyz2geo_round_trip():
    """Round-trip geo -> xyz -> geo in radians and dec_deg."""
    lon_deg, lat_deg, he = -112.567, 43.36, 27.2
    x, y, z = COO.geo2xyz(lon_deg, lat_deg, he, unit="dec_deg")
    lon_out, lat_out, he_out = COO.xyz2geo(x, y, z, unit="dec_deg")
    assert np.allclose((lon_out, lat_out, he_out), (lon_deg, lat_deg, he), atol=1e-8)


# -----------------------------------------------------------------------------
# xyz2geospheric (XYZ -> longitude, latitude, radius)
# -----------------------------------------------------------------------------


def test_xyz2geospheric_radians():
    """xyz2geospheric; GAMIT tform reference for same XYZ."""
    rlon, rlat, R = COO.xyz2geospheric(918129.451, -4346071.255, 4561977.839)
    # Geocentric lat/lon/radius; tform N45 45 48.48702 W78 4 16.92154 6367333.7313
    assert np.allclose(R, 6367333.7313, atol=0.0001)
    deg_lon = np.degrees(rlon)
    deg_lat = np.degrees(rlat)
    if deg_lon > 180:
        deg_lon -= 360
    assert np.allclose(deg_lon, -78.07137, atol=1e-4)
    assert np.allclose(deg_lat, 45.76347, atol=1e-4)


def test_xyz2geospheric_deg_min_sec():
    """xyz2geospheric plus radians2deg_mn_sec; tform N45 45 48.48702 W78 4 16.92154."""
    rlon, rlat, h = COO.xyz2geospheric(918129.451, -4346071.255, 4561977.839)
    deg, mn, sec = UNITS.radians2deg_mn_sec(rlon, angle_range="-180-180")
    assert deg == -78
    assert mn == 4
    assert np.allclose(sec, 16.92154, atol=0.00001)
    deg, mn, sec = UNITS.radians2deg_mn_sec(rlat, angle_range="-180-180")
    assert deg == 45
    assert mn == 45
    assert np.allclose(sec, 48.48702, atol=0.00001)
    assert np.allclose(h, 6367333.7313, atol=0.0001)


# -----------------------------------------------------------------------------
# denu_at_x0y0z0_to_xyz and delta_xyz -> delta_NEU (round-trip via rotation)
# -----------------------------------------------------------------------------


def test_denu_at_x0y0z0_to_xyz_and_back_to_neu():
    """Local NEU offset -> XYZ offset -> back to NEU using rotation matrix."""
    x0, y0, z0 = COO.geo2xyz(-112.567, 43.36, 27.2, unit="dec_deg")
    dn, de, du = 10.0, 5.0, -2.0

    xn, yn, zn = COO.denu_at_x0y0z0_to_xyz(de, dn, du, x0, y0, z0)
    dx, dy, dz = xn - x0, yn - y0, zn - z0

    # Recover NEU from delta_xyz using rotation (general-to-local)
    lam, phi, _ = COO.xyz2geo(x0, y0, z0)
    R = COO.mat_rot_general_to_local(lam, phi)
    enu = np.dot(R, np.array([dx, dy, dz]))
    # R rows are E, N, U
    rn, re, ru = enu[1], enu[0], enu[2]

    assert np.allclose([rn, re, ru], [dn, de, du], atol=1e-6)


# -----------------------------------------------------------------------------
# Rotation matrices
# -----------------------------------------------------------------------------


def test_mat_rot_general_to_local_orthogonal():
    """Rotation matrix R is orthogonal: R @ R.T = I."""
    lam, phi = np.radians(10.0), np.radians(45.0)
    R = COO.mat_rot_general_to_local(lam, phi)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)


def test_mat_rot_local_to_general_is_transpose():
    """mat_rot_local_to_general is transpose of mat_rot_general_to_local."""
    lam, phi = np.radians(10.0), np.radians(45.0)
    R_gl = COO.mat_rot_general_to_local(lam, phi)
    R_lg = COO.mat_rot_local_to_general(lam, phi)
    assert np.allclose(R_lg, R_gl.T, atol=1e-10)


# -----------------------------------------------------------------------------
# azimuth_to_en
# -----------------------------------------------------------------------------


def test_azimuth_to_en():
    """Azimuth (deg) to east/north unit vector; azimuth 0 = north."""
    east, north = COO.azimuth_to_en(0.0)
    assert np.allclose((east, north), (0.0, 1.0), atol=1e-10)
    east, north = COO.azimuth_to_en(90.0)
    assert np.allclose((east, north), (1.0, 0.0), atol=1e-10)


# -----------------------------------------------------------------------------
# wnorm
# -----------------------------------------------------------------------------


def test_wnorm_pole_and_equator():
    """wnorm (prime vertical radius N) at pole and equator."""
    # At equator, N = a / sqrt(1 - e² sin²(0)) = a
    phi_equator = 0.0
    w_eq = COO.wnorm(phi_equator)
    a = 6378137.0
    assert np.allclose(w_eq, a, atol=1.0)
    # At pole, N = a / sqrt(1 - e²) > a
    phi_pole = np.pi / 2
    w_pole = COO.wnorm(phi_pole)
    assert w_pole > w_eq


# -----------------------------------------------------------------------------
# Spherical distances
# -----------------------------------------------------------------------------


def test_xyz_spherical_distance_same_point():
    """Distance between same point is zero."""
    x, y, z = COO.geo2xyz(-112.567, 43.36, 0.0, unit="dec_deg")
    d = COO.xyz_spherical_distance(x, y, z, x, y, z)
    assert np.allclose(d, 0.0, atol=1e-6)


def test_geo_spherical_distance_same_point():
    """Geo spherical distance between same point is zero."""
    d = COO.geo_spherical_distance(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, unit="radians")
    assert np.allclose(d, 0.0, atol=1e-6)


# -----------------------------------------------------------------------------
# geo2flat_earth / flat_earth2geo (Web Mercator)
# -----------------------------------------------------------------------------


def test_geo2flat_earth_flat_earth2geo_round_trip():
    """Round-trip WGS84 <-> Web Mercator (km)."""
    lon, lat = 7.0, 43.5
    x_km, y_km = COO.geo2flat_earth(lon, lat)
    lon_out, lat_out = COO.flat_earth2geo(x_km, y_km)
    assert np.allclose(lon_out, lon, atol=1e-6)
    assert np.allclose(lat_out, lat, atol=1e-6)
