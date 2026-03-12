"""
Unit tests for pyacs.lib.gmtpoint (GMT_Point: position and velocity).
"""

import numpy as np

from pyacs.lib.gmtpoint import GMT_Point


# -----------------------------------------------------------------------------
# GMT_Point creation and attributes
# -----------------------------------------------------------------------------


def test_gmt_point_creation():
    """GMT_Point stores code, lon, lat, he and optional velocity."""
    p = GMT_Point(code="TEST", lon=2.0, lat=45.0, he=100.0)
    assert p.code == "TEST"
    assert p.lon == 2.0
    assert p.lat == 45.0
    assert p.he == 100.0
    assert p.Ve is None
    assert p.Vn is None


def test_gmt_point_with_velocity():
    """GMT_Point with Ve, Vn, Vu and sigmas."""
    p = GMT_Point(
        code="SITE",
        lon=-78.0,
        lat=46.0,
        Ve=1.0,
        Vn=2.0,
        Vu=0.5,
        SVe=0.1,
        SVn=0.2,
        SVu=0.05,
    )
    assert p.Ve == 1.0
    assert p.Vn == 2.0
    assert p.Vu == 0.5
    assert p.SVe == 0.1
    assert p.SVn == 0.2


# -----------------------------------------------------------------------------
# copy
# -----------------------------------------------------------------------------


def test_copy():
    """copy() returns a new point with same attributes."""
    p = GMT_Point(code="A", lon=0.0, lat=0.0, Ve=1.0, Vn=0.0, index=5)
    q = p.copy()
    assert q is not p
    assert q.code == p.code
    assert q.lon == p.lon
    assert q.lat == p.lat
    assert q.Ve == p.Ve
    assert q.Vn == p.Vn
    assert q.index == p.index


# -----------------------------------------------------------------------------
# get_info
# -----------------------------------------------------------------------------


def test_get_info_returns_string():
    """get_info(display=False) returns a non-empty string."""
    p = GMT_Point(code="X", lon=5.0, lat=43.0, Ve=0.0, Vn=0.0)
    info = p.get_info(display=False)
    assert isinstance(info, str)
    assert "X" in info
    assert "5" in info
    assert "43" in info


# -----------------------------------------------------------------------------
# magaz
# -----------------------------------------------------------------------------


def test_magaz_magnitude_azimuth():
    """magaz() returns magnitude and azimuth (clockwise from North)."""
    p = GMT_Point(code="M", lon=0.0, lat=0.0, Ve=0.0, Vn=10.0)
    mag, az = p.magaz()
    assert np.isclose(mag, 10.0, atol=1e-10)
    # North: Ve=0, Vn>0 -> atan2(0, 10) = 0°
    assert np.isclose(az, 0.0, atol=1e-6)


def test_magaz_east():
    """Velocity due east: azimuth 90°."""
    p = GMT_Point(code="E", lon=0.0, lat=0.0, Ve=5.0, Vn=0.0)
    mag, az = p.magaz()
    assert np.isclose(mag, 5.0, atol=1e-10)
    assert np.isclose(az, 90.0, atol=1e-6)


def test_magaz_magnitude():
    """Magnitude is sqrt(Ve^2 + Vn^2)."""
    p = GMT_Point(code="P", lon=0.0, lat=0.0, Ve=3.0, Vn=4.0)
    mag, _ = p.magaz()
    assert np.isclose(mag, 5.0, atol=1e-10)


# -----------------------------------------------------------------------------
# spherical_distance
# -----------------------------------------------------------------------------


def test_spherical_distance_same_point():
    """Spherical distance from a point to itself is zero."""
    p = GMT_Point(code="A", lon=2.0, lat=45.0)
    d = p.spherical_distance(p)
    assert np.isclose(d, 0.0, atol=1e-6)


def test_spherical_distance_antipodal():
    """Spherical distance between antipodal points is pi * R_earth."""
    p = GMT_Point(code="A", lon=0.0, lat=0.0)
    q = GMT_Point(code="B", lon=180.0, lat=0.0)
    d = p.spherical_distance(q)
    R = 6.371e6
    expected = np.pi * R
    assert np.isclose(d, expected, rtol=1e-4)


def test_spherical_distance_small():
    """Spherical distance between two close points is positive and ~ deg-to-m."""
    p = GMT_Point(code="A", lon=0.0, lat=45.0)
    q = GMT_Point(code="B", lon=0.001, lat=45.0)
    d = p.spherical_distance(q)
    # At 45° lat, 1 deg lon ~ 78 km; 0.001 deg ~ 78 m
    assert d > 0
    assert 50 < d < 150  # meters


# -----------------------------------------------------------------------------
# midpoint
# -----------------------------------------------------------------------------


def test_midpoint_symmetric():
    """Midpoint of A and B is same as midpoint of B and A (lon/lat)."""
    a = GMT_Point(code="A", lon=0.0, lat=0.0)
    b = GMT_Point(code="B", lon=2.0, lat=4.0)
    m_ab = a.midpoint(b, code="M1")
    m_ba = b.midpoint(a, code="M2")
    assert np.isclose(m_ab.lon, m_ba.lon, atol=1e-8)
    assert np.isclose(m_ab.lat, m_ba.lat, atol=1e-8)


def test_midpoint_between_same():
    """Midpoint of a point with itself is that point."""
    p = GMT_Point(code="P", lon=10.0, lat=20.0)
    m = p.midpoint(p, code="M")
    assert np.isclose(m.lon, 10.0, atol=1e-8)
    assert np.isclose(m.lat, 20.0, atol=1e-8)


def test_midpoint_has_code():
    """Midpoint accepts custom code."""
    a = GMT_Point(code="A", lon=0.0, lat=0.0)
    b = GMT_Point(code="B", lon=1.0, lat=1.0)
    m = a.midpoint(b, code="MID")
    assert m.code == "MID"


# -----------------------------------------------------------------------------
# rotate_vel
# -----------------------------------------------------------------------------


def test_rotate_vel_magnitude_preserved():
    """Rotating velocity preserves magnitude."""
    p = GMT_Point(code="R", lon=0.0, lat=0.0, Ve=1.0, Vn=0.0)
    q = p.rotate_vel(np.radians(30.0), unit="radians")
    mag_q, _ = q.magaz()
    mag_p, _ = p.magaz()
    assert np.isclose(mag_q, mag_p, atol=1e-10)


def test_rotate_vel_round_trip():
    """Rotate by theta then by -theta recovers original Ve, Vn."""
    p = GMT_Point(code="R", lon=0.0, lat=0.0, Ve=2.0, Vn=3.0)
    theta_deg = 40.0
    q = p.rotate_vel(theta_deg, unit="degrees")
    r = q.rotate_vel(-theta_deg, unit="degrees")
    assert np.isclose(r.Ve, p.Ve, atol=1e-10)
    assert np.isclose(r.Vn, p.Vn, atol=1e-10)


def test_rotate_vel_degrees():
    """rotate_vel with unit='degrees' accepts angle in degrees."""
    p = GMT_Point(code="R", lon=0.0, lat=0.0, Ve=1.0, Vn=0.0)
    q = p.rotate_vel(90.0, unit="degrees")
    # Initial: (Ve, Vn) = (1, 0). Rotate 90° CCW -> (-0, 1) = (0, 1) in (E,N)
    # Implementation uses R = [[cos,-sin],[sin,cos]] so R@[1,0]=[-sin,cos]=[-1,0] for 90°? No: cos(90)=0, sin(90)=1, so [0,-1] and [1,0] -> (0*1 + (-1)*0, 1*1+0*0) = (0, 1). So (0, 1) = North. So Ve=0, Vn=1.
    assert np.isclose(q.Ve, 0.0, atol=1e-10)
    assert np.isclose(q.Vn, 1.0, atol=1e-10)


# -----------------------------------------------------------------------------
# assign_index / get_index
# -----------------------------------------------------------------------------


def test_assign_index_get_index():
    """assign_index sets index; get_index returns it."""
    p = GMT_Point(code="I", lon=0.0, lat=0.0)
    p.assign_index(42)
    assert p.get_index() == 42


def test_assign_index_returns_self():
    """assign_index returns self for chaining."""
    p = GMT_Point(code="I", lon=0.0, lat=0.0)
    out = p.assign_index(7)
    assert out is p


# -----------------------------------------------------------------------------
# pole (Euler pole prediction)
# -----------------------------------------------------------------------------


def test_pole_predict_rot_at_pole():
    """At Euler pole (rotation axis), predicted velocity is zero (or very small)."""
    # Euler pole at North Pole: lon=0, lat=90. Rotation vector W = (0, 0, omega) in rad/yr.
    # Point at North Pole: velocity = omega x r = 0.
    p = GMT_Point(code="NP", lon=0.0, lat=90.0, Ve=0.0, Vn=0.0)
    W = np.array([[0.0], [0.0], [1e-6]])  # (3,1) rad/yr, rotation around z
    q = p.pole(W=W, type_euler="rot", option="predict")
    assert np.isclose(q.Ve, 0.0, atol=1e-6)
    assert np.isclose(q.Vn, 0.0, atol=1e-6)


def test_pole_predict_rot_equator():
    """At equator, rotation around z gives eastward velocity."""
    # Point on equator at lon=0, lat=0. W = (0, 0, omega) -> v = omega x r, eastward.
    p = GMT_Point(code="EQ", lon=0.0, lat=0.0, Ve=0.0, Vn=0.0)
    omega_rad_yr = 1e-6
    W = np.array([[0.0], [0.0], [omega_rad_yr]])
    q = p.pole(W=W, type_euler="rot", option="predict")
    # Velocity magnitude ~ omega * R_earth in m/s, then *1e3 for mm/yr: omega*R*1e3 mm/yr
    R_m = 6.371e6
    expected_mm_yr = omega_rad_yr * R_m * 1e3  # conversion in code: pve = Pi[0,0]*1e3
    assert q.Ve is not None
    assert q.Vn is not None
    # Eastward at (0,0): Ve > 0, Vn ≈ 0
    assert q.Ve > 0
    assert np.isclose(q.Vn, 0.0, atol=0.1)


def test_pole_remove():
    """option='remove' subtracts pole prediction from velocity."""
    p = GMT_Point(code="S", lon=10.0, lat=45.0, Ve=2.0, Vn=1.0)
    W = np.array([[0.0], [0.0], [0.0]])  # zero pole (3,1)
    q = p.pole(W=W, type_euler="rot", option="remove")
    assert np.isclose(q.Ve, 2.0, atol=1e-10)
    assert np.isclose(q.Vn, 1.0, atol=1e-10)


def test_pole_add():
    """option='add' adds pole prediction to velocity."""
    p = GMT_Point(code="S", lon=10.0, lat=45.0, Ve=0.0, Vn=0.0)
    W = np.array([[0.0], [0.0], [1e-6]])
    q_pred = p.pole(W=W, type_euler="rot", option="predict")
    q_add = p.pole(W=W, type_euler="rot", option="add")
    assert np.isclose(q_add.Ve, q_pred.Ve, atol=1e-10)
    assert np.isclose(q_add.Vn, q_pred.Vn, atol=1e-10)
