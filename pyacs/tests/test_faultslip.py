"""
Unit tests for pyacs.lib.faultslip (fault geometry, slip, rake, strike).
"""

import numpy as np

import pyacs.lib.faultslip as FS


# -----------------------------------------------------------------------------
# unit_slip
# -----------------------------------------------------------------------------


def test_unit_slip_unit_norm():
    """Unit slip vector has magnitude 1."""
    strike, dip, rake = 45.0, 60.0, -90.0
    ue, un, uu = FS.unit_slip(strike, dip, rake)
    assert np.isclose(np.sqrt(ue**2 + un**2 + uu**2), 1.0, atol=1e-10)


def test_unit_slip_pure_strike_slip_rake_zero():
    """Rake 0: slip along strike (horizontal)."""
    strike, dip, rake = 90.0, 90.0, 0.0
    ue, un, uu = FS.unit_slip(strike, dip, rake)
    # Strike 90 = fault strikes E-W, so slip along strike is E or W
    assert np.isclose(uu, 0.0, atol=1e-10)
    assert np.isclose(np.abs(ue), 1.0, atol=1e-10)
    assert np.isclose(un, 0.0, atol=1e-10)


def test_unit_slip_pure_dip_slip_rake_90():
    """Rake 90: reverse (up-dip); vertical fault has slip in horizontal plane."""
    strike, dip, rake = 0.0, 45.0, 90.0
    ue, un, uu = FS.unit_slip(strike, dip, rake)
    # Up-dip component; should have non-zero U
    assert not np.isclose(uu, 0.0, atol=1e-6)


# -----------------------------------------------------------------------------
# unit_normal
# -----------------------------------------------------------------------------


def test_unit_normal_unit_norm():
    """Unit normal has magnitude 1."""
    strike, dip = 30.0, 50.0
    ue, un, uu = FS.unit_normal(strike, dip)
    assert np.isclose(np.sqrt(ue**2 + un**2 + uu**2), 1.0, atol=1e-10)


def test_unit_normal_perpendicular_to_unit_slip_rake_zero():
    """For rake 0 or 180, slip is in the fault plane so normal · slip = 0."""
    strike, dip = 20.0, 40.0
    ne, nn, nu = FS.unit_normal(strike, dip)
    for rake in (0.0, 180.0):
        se, sn, su = FS.unit_slip(strike, dip, rake)
        dot = ne * se + nn * sn + nu * su
        assert np.isclose(dot, 0.0, atol=1e-10)


# -----------------------------------------------------------------------------
# slip_rake_2_ds_ss / convention
# -----------------------------------------------------------------------------


def test_slip_rake_2_ds_ss_pure_strike_slip():
    """Rake 0: all slip is strike-slip; dip-slip zero."""
    slip, rake = 1.0, 0.0
    ds, ss = FS.slip_rake_2_ds_ss(slip, rake)
    assert np.isclose(ds, 0.0, atol=1e-10)
    assert np.isclose(ss, 1.0, atol=1e-10)


def test_slip_rake_2_ds_ss_pure_dip_slip():
    """Rake 90: all slip is dip-slip; strike-slip zero."""
    slip, rake = 1.0, 90.0
    ds, ss = FS.slip_rake_2_ds_ss(slip, rake)
    assert np.isclose(ds, 1.0, atol=1e-10)
    assert np.isclose(ss, 0.0, atol=1e-10)


def test_slip_rake_2_ds_ss_magnitude():
    """Slip magnitude squared = ds^2 + ss^2."""
    slip, rake = 2.0, 35.0
    ds, ss = FS.slip_rake_2_ds_ss(slip, rake)
    assert np.isclose(np.sqrt(ds**2 + ss**2), slip, atol=1e-10)


# -----------------------------------------------------------------------------
# strike_dip_rake_to_dir
# -----------------------------------------------------------------------------


def test_strike_dip_rake_to_dir_consistent_with_unit_slip():
    """Slip direction from strike_dip_rake_to_dir matches unit_slip horizontal."""
    strike, dip, rake = 100.0, 50.0, -45.0
    slip_dir = FS.strike_dip_rake_to_dir(strike, dip, rake)
    ue, un, uu = FS.unit_slip(strike, dip, rake)
    expected_az = np.degrees(np.arctan2(ue, un))
    assert np.isclose(np.mod(slip_dir - expected_az, 360.0), 0.0, atol=1e-6) or \
           np.isclose(np.mod(slip_dir - expected_az, 360.0), 360.0, atol=1e-6)


# -----------------------------------------------------------------------------
# v_to_n_ss / ss_ns_2_ve_vn
# -----------------------------------------------------------------------------


def test_v_to_n_ss_ss_ns_2_ve_vn_round_trip():
    """Decompose ve,vn into normal/strike-slip and back.
    v_to_n_ss uses normal axis (cos(strike), -sin(strike)); ss_ns_2_ve_vn uses
    ns axis (-cos(strike), sin(strike)), so pass -normal for correct round-trip.
    """
    ve, vn = 1.0, 2.0
    strike = 30.0
    normal, ss = FS.v_to_n_ss(ve, vn, strike)
    ve_back, vn_back = FS.ss_ns_2_ve_vn(ss, -normal, strike)
    assert np.isclose(ve_back, ve, atol=1e-10)
    assert np.isclose(vn_back, vn, atol=1e-10)


def test_v_to_n_ss_along_strike():
    """Motion along strike: normal component zero."""
    strike = 45.0
    # Velocity perpendicular to normal = along strike: e.g. (sin(strike), cos(strike))
    ve = np.sin(np.radians(strike))
    vn = np.cos(np.radians(strike))
    normal, ss_comp = FS.v_to_n_ss(ve, vn, strike)
    assert np.isclose(normal, 0.0, atol=1e-10)
    assert np.isclose(np.sqrt(ve**2 + vn**2), np.abs(ss_comp), atol=1e-10)


# -----------------------------------------------------------------------------
# v_to_rake
# -----------------------------------------------------------------------------


def test_v_to_rake_pure_strike_slip():
    """Horizontal motion along strike (strike 90 = E-W, motion E) gives rake 0 or 180."""
    strike, dip = 90.0, 80.0
    ve, vn = 1.0, 0.0  # motion due east = along strike
    rake = FS.v_to_rake(ve, vn, strike, dip, style="leftlateral")
    assert np.isclose(np.abs(rake), 0.0, atol=1e-6) or np.isclose(np.abs(rake), 180.0, atol=1e-6)


def test_v_to_rake_round_trip_via_unit_slip():
    """Given strike, dip, rake: unit_slip gives direction; v_to_rake recovers rake."""
    strike, dip, rake = 120.0, 60.0, 45.0
    ue, un, uu = FS.unit_slip(strike, dip, rake)
    # Horizontal part (normalize)
    h = np.sqrt(ue**2 + un**2)
    if h < 1e-10:
        return  # skip if purely vertical
    ve, vn = ue / h, un / h
    rake_back = FS.v_to_rake(ve, vn, strike, dip, style="inverse")
    # Rake is defined modulo 180 for direction
    assert np.isclose(np.mod(rake_back - rake + 180, 360) - 180, 0.0, atol=2.0) or \
           np.isclose(np.mod(rake_back + rake + 180, 360) - 180, 0.0, atol=2.0)


# -----------------------------------------------------------------------------
# rake_from_slip_az
# -----------------------------------------------------------------------------


def test_rake_from_slip_az_along_strike():
    """Slip direction = strike (along fault): rake 0 or 180."""
    strike, dip = 45.0, 70.0
    slipdir = strike
    rake_ll = FS.rake_from_slip_az(strike, dip, slipdir, "leftlateral")
    rake_rl = FS.rake_from_slip_az(strike, dip, slipdir, "rightlateral")
    assert np.isclose(np.abs(rake_ll), 0.0, atol=1e-5) or np.isclose(np.abs(rake_ll), 180.0, atol=1e-5)
    assert np.isclose(np.abs(rake_rl), 0.0, atol=1e-5) or np.isclose(np.abs(rake_rl), 180.0, atol=1e-5)


# -----------------------------------------------------------------------------
# geo_to_strike
# -----------------------------------------------------------------------------


def test_geo_to_strike_north_south_segment():
    """Segment due N-S has strike 0 (or 180)."""
    # Start at (0, 0), end at (0, 1) -> north
    strike = FS.geo_to_strike(0.0, 0.0, 0.0, 1.0)
    assert np.isclose(strike, 0.0, atol=1e-5)


def test_geo_to_strike_east_west_segment():
    """Segment due E-W has strike 90 (or -90)."""
    strike = FS.geo_to_strike(0.0, 0.0, 1.0, 0.0)
    assert np.isclose(np.abs(strike), 90.0, atol=1e-5)


# -----------------------------------------------------------------------------
# fault_to_corners
# -----------------------------------------------------------------------------


def test_fault_to_corners_shape():
    """fault_to_corners returns (4, 3) array of [lon, lat, depth]."""
    corners = FS.fault_to_corners(
        lon=7.0, lat=44.0, depth=-10.0,
        length=20.0, width=10.0, strike=90.0, dip=45.0
    )
    assert corners.shape == (4, 3)
    assert np.all(corners[:, 2] <= 0.0)  # depth non-positive
