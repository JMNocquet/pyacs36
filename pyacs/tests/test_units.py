"""
Unit tests for pyacs.lib.units (angle and seismic moment conversions).
"""

import math

import numpy as np

import pyacs.lib.units as U


# -----------------------------------------------------------------------------
# mas2rad / rad2mas
# -----------------------------------------------------------------------------


def test_mas2rad_one():
    """1 mas in radians (docstring value)."""
    rad = U.mas2rad(1.0)
    expected = 2.0 * math.pi / 360.0 / 60.0 / 60.0 / 1000.0
    assert np.isclose(rad, expected, atol=1e-15)
    assert np.isclose(rad, 4.84813681109536e-09, atol=1e-15)


def test_rad2mas_one():
    """1 rad in milliarcseconds."""
    mas = U.rad2mas(1.0)
    expected = 360.0 * 60.0 * 60.0 * 1000.0 / (2.0 * math.pi)
    assert np.isclose(mas, expected, atol=1e-6)


def test_mas2rad_rad2mas_round_trip():
    """mas2rad then rad2mas recovers original value."""
    mas_in = 1234.5
    rad = U.mas2rad(mas_in)
    mas_out = U.rad2mas(rad)
    assert np.isclose(mas_out, mas_in, atol=1e-10)


def test_rad2mas_mas2rad_round_trip():
    """rad2mas then mas2rad recovers original value."""
    rad_in = 1e-6
    mas = U.rad2mas(rad_in)
    rad_out = U.mas2rad(mas)
    assert np.isclose(rad_out, rad_in, atol=1e-18)


def test_mas2rad_array():
    """mas2rad accepts array and returns array."""
    mas = np.array([1.0, 2.0, 0.0])
    rad = U.mas2rad(mas)
    assert isinstance(rad, np.ndarray)
    assert rad.shape == (3,)
    assert np.isclose(rad[0], U.mas2rad(1.0), atol=1e-15)
    assert np.isclose(rad[1], U.mas2rad(2.0), atol=1e-15)
    assert rad[2] == 0.0


def test_rad2mas_scalar():
    """rad2mas with scalar returns scalar (float)."""
    mas = U.rad2mas(0.0)
    assert np.isscalar(mas) or (isinstance(mas, np.ndarray) and mas.ndim == 0)
    assert np.isclose(mas, 0.0, atol=1e-10)


# -----------------------------------------------------------------------------
# radians2deg_mn_sec
# -----------------------------------------------------------------------------


def test_radians2deg_mn_sec_zero():
    """Zero radians -> 0 deg, 0 min, 0 sec."""
    deg, mn, sec = U.radians2deg_mn_sec(0.0)
    assert deg == 0
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_90():
    """pi/2 radians -> 90 deg, 0 min, 0 sec."""
    deg, mn, sec = U.radians2deg_mn_sec(math.pi / 2.0)
    assert deg == 90
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_180():
    """pi radians -> 180 deg (in -180-180 range)."""
    deg, mn, sec = U.radians2deg_mn_sec(math.pi, angle_range='-180-180')
    assert deg == 180
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_negative():
    """-pi/2 -> -90 deg, 0 min, 0 sec."""
    deg, mn, sec = U.radians2deg_mn_sec(-math.pi / 2.0)
    assert deg == -90
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_0_360():
    """angle_range 0-360: 270° stays 270."""
    deg, mn, sec = U.radians2deg_mn_sec(3.0 * math.pi / 2.0, angle_range='0-360')
    assert deg == 270
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_270_minus180_180():
    """270° in -180-180 range -> -90 deg."""
    deg, mn, sec = U.radians2deg_mn_sec(3.0 * math.pi / 2.0, angle_range='-180-180')
    assert deg == -90
    assert mn == 0
    assert np.isclose(sec, 0.0, atol=1e-10)


def test_radians2deg_mn_sec_with_minutes_seconds():
    """Angle that has non-zero minutes and seconds."""
    # 1.5 deg = 1 deg 30 min 0 sec
    rad = np.radians(1.5)
    deg, mn, sec = U.radians2deg_mn_sec(rad)
    assert deg == 1
    assert mn == 30
    assert np.isclose(sec, 0.0, atol=1e-10)


# -----------------------------------------------------------------------------
# moment_to_magnitude / magnitude_to_moment
# -----------------------------------------------------------------------------


def test_moment_to_magnitude_zero():
    """M0 = 10^9.05 N.m -> Mw = 0 (standard definition)."""
    moment = 10.0 ** 9.05
    mw = U.moment_to_magnitude(moment)
    assert np.isclose(mw, 0.0, atol=1e-10)


def test_magnitude_to_moment_zero():
    """Mw = 0 -> M0 = 10^9.05 N.m."""
    moment = U.magnitude_to_moment(0.0)
    assert np.isclose(moment, 10.0 ** 9.05, atol=1e-6)


def test_moment_magnitude_round_trip():
    """moment_to_magnitude then magnitude_to_moment recovers M0."""
    m0 = 1e20  # large event
    mw = U.moment_to_magnitude(m0)
    m0_back = U.magnitude_to_moment(mw)
    assert np.isclose(m0_back, m0, rtol=1e-10)


def test_magnitude_moment_round_trip():
    """magnitude_to_moment then moment_to_magnitude recovers Mw."""
    mw = 7.0
    m0 = U.magnitude_to_moment(mw)
    mw_back = U.moment_to_magnitude(m0)
    assert np.isclose(mw_back, mw, atol=1e-10)


def test_moment_to_magnitude_formula():
    """Mw = 2/3 * (log10(M0) - 9.05)."""
    m0 = 1e21
    mw = U.moment_to_magnitude(m0)
    expected = 2.0 / 3.0 * (math.log10(m0) - 9.05)
    assert np.isclose(mw, expected, atol=1e-10)


def test_magnitude_to_moment_formula():
    """M0 = 10^(1.5*Mw + 9.05)."""
    mw = 6.0
    m0 = U.magnitude_to_moment(mw)
    expected = 10.0 ** (1.5 * mw + 9.05)
    assert np.isclose(m0, expected, atol=1e-6)


def test_moment_to_magnitude_array():
    """moment_to_magnitude accepts array."""
    moments = np.array([1e18, 1e20, 1e22])
    mw = U.moment_to_magnitude(moments)
    assert mw.shape == (3,)
    assert mw[1] > mw[0]
    assert mw[2] > mw[1]


def test_magnitude_to_moment_array():
    """magnitude_to_moment accepts array."""
    mw = np.array([5.0, 6.0, 7.0])
    moments = U.magnitude_to_moment(mw)
    assert moments.shape == (3,)
    assert moments[1] > moments[0]
    assert moments[2] > moments[1]
