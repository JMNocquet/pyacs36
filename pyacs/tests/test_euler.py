"""
Unit tests for pyacs.lib.euler (Euler pole, rotation vectors, velocities).
"""

import numpy as np

from pyacs.lib import euler as EUL
from pyacs.lib import coordinates as COO


# -----------------------------------------------------------------------------
# euler2rot / rot2euler
# -----------------------------------------------------------------------------


def test_euler2rot_rot2euler_round_trip():
    """Euler pole (lon, lat, omega deg/Myr) <-> rotation vector (rad/yr) round-trip."""
    lon, lat, omega = 10.0, 45.0, 0.5
    wx, wy, wz = EUL.euler2rot(lon, lat, omega)
    lon_out, lat_out, omega_out = EUL.rot2euler(wx, wy, wz)
    assert np.isclose(lon_out, lon, atol=1e-6)
    assert np.isclose(lat_out, lat, atol=1e-6)
    assert np.isclose(omega_out, omega, atol=1e-6)


def test_euler2rot_pole_at_north():
    """Euler pole at North pole (0, 90) gives rotation vector along Z."""
    wx, wy, wz = EUL.euler2rot(0.0, 90.0, 1.0)
    assert np.isclose(wx, 0.0, atol=1e-15)
    assert np.isclose(wy, 0.0, atol=1e-15)
    assert wz > 0
    # omega 1 deg/Myr = 1e-6 deg/yr -> rad/yr = 1e-6 * pi/180
    assert np.isclose(wz, np.radians(1e-6), atol=1e-15)


def test_rot2euler_zero_omega():
    """Zero rotation vector: rot2euler(0,0,0) can be singular; use small non-zero."""
    wx, wy, wz = EUL.euler2rot(30.0, 20.0, 0.1)
    lon, lat, omega = EUL.rot2euler(wx, wy, wz)
    assert np.isclose(lon, 30.0, atol=1e-5)
    assert np.isclose(lat, 20.0, atol=1e-5)
    assert np.isclose(omega, 0.1, atol=1e-5)


# -----------------------------------------------------------------------------
# vel_from_euler
# -----------------------------------------------------------------------------


def test_vel_from_euler_at_pole():
    """Velocity at the Euler pole is zero when using spherical coordinates."""
    lon_pole, lat_pole = 10.0, 45.0
    omega = 0.5
    ve, vn = EUL.vel_from_euler(
        lon_pole, lat_pole, lon_pole, lat_pole, omega, use_spherical=True
    )
    assert np.isclose(ve, 0.0, atol=1e-10)
    assert np.isclose(vn, 0.0, atol=1e-10)


def test_vel_from_euler_vs_pole_matrix():
    """vel_from_euler matches pole_matrix @ w (with unit conversion)."""
    lon_site, lat_site = 7.0, 44.0
    lon_euler, lat_euler, omega = 10.0, 45.0, 0.3
    ve_mm, vn_mm = EUL.vel_from_euler(
        lon_site, lat_site, lon_euler, lat_euler, omega
    )
    wx, wy, wz = EUL.euler2rot(lon_euler, lat_euler, omega)
    w = np.array([wx, wy, wz])
    W = EUL.pole_matrix(np.array([[lon_site, lat_site]]))
    v_m_yr = W @ w
    # vel_from_euler returns mm/yr; pole_matrix returns m/yr
    assert np.isclose(ve_mm, v_m_yr[0] * 1e3, atol=1e-3)
    assert np.isclose(vn_mm, v_m_yr[1] * 1e3, atol=1e-3)


# -----------------------------------------------------------------------------
# pole_matrix
# -----------------------------------------------------------------------------


def test_pole_matrix_shape():
    """pole_matrix has shape (2*n, 3) for n sites."""
    coor = np.array([[0.0, 45.0], [10.0, 44.0]])
    W = EUL.pole_matrix(coor)
    assert W.shape == (4, 3)


def test_pole_matrix_single_site():
    """Single site: pole_matrix has shape (2, 3)."""
    coor = np.array([[5.0, 43.0]])
    W = EUL.pole_matrix(coor)
    assert W.shape == (2, 3)


# -----------------------------------------------------------------------------
# pole_matrix_fault
# -----------------------------------------------------------------------------


def test_pole_matrix_fault_shape():
    """pole_matrix_fault has shape (2*n, 3) for n sites."""
    coor = np.array([[0.0, 45.0], [10.0, 44.0]])
    strike = np.array([0.0, 90.0])
    W = EUL.pole_matrix_fault(coor, strike)
    assert W.shape == (4, 3)


def test_pole_matrix_fault_vs_pole_matrix_strike_zero():
    """pole_matrix_fault with strike=0 matches pole_matrix (E,N = E,N)."""
    coor = np.array([[7.0, 44.0]])
    strike = np.array([0.0])
    W_fault = EUL.pole_matrix_fault(coor, strike)
    W_enu = EUL.pole_matrix(coor)
    # Strike 0: first row is N (north), second is -E (east). So [sin(0), cos(0)] = [0,1] -> N; [-cos(0), sin(0)] = [-1, 0] -> -E. So row0 = N, row1 = -E. So W_fault @ w = [vn, -ve] while W_enu @ w = [ve, vn]. So W_fault[0] = W_enu[1], W_fault[1] = -W_enu[0]. Let me check the code again. RF = [[sin(strike), cos(strike)], [-cos(strike), sin(strike)]]. For strike=0: RF = [[0, 1], [-1, 0]]. So RW = RF @ RAi. RAi rows are [E_component, N_component] of the design. So RW row0 = 0*E + 1*N = N, RW row1 = -1*E + 0*N = -E. So W_fault gives [N, -E] in the two rows. W_enu gives [E, N]. So W_fault[0] = W_enu[1], W_fault[1] = -W_enu[0]. So we can test that.
    np.testing.assert_allclose(W_fault[0], W_enu[1], atol=1e-10)
    np.testing.assert_allclose(W_fault[1], -W_enu[0], atol=1e-10)
