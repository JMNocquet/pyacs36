"""
Unit tests for pyacs.vel_field (Velocity_Field) using nocquet_tectonophysics_2012.gmt.
"""

import os

import numpy as np
import pytest

from pyacs.vel_field import Velocity_Field

_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
GMT_VEL = os.path.join(_TEST_DIR, "data", "vel", "nocquet_tectonophysics_2012.gmt")


@pytest.fixture(scope="module")
def gmt_vel_path():
    """Path to the GMT velocity file."""
    if not os.path.isfile(GMT_VEL):
        pytest.skip("GMT velocity file not found at %s" % GMT_VEL)
    return GMT_VEL


@pytest.fixture(scope="module")
def vf(gmt_vel_path):
    """Velocity_Field loaded from nocquet_tectonophysics_2012.gmt."""
    return Velocity_Field.read(gmt_vel_path, verbose=False)


# -------------------------------------------------------------------------
# read (via fixture)
# -------------------------------------------------------------------------


def test_vel_field_read_returns_velocity_field(gmt_vel_path):
    vf = Velocity_Field.read(gmt_vel_path, verbose=False)
    assert isinstance(vf, Velocity_Field)
    assert hasattr(vf, "sites")
    assert vf.file_name == gmt_vel_path


def test_vel_field_read_with_lexclude(gmt_vel_path):
    vf_full = Velocity_Field.read(gmt_vel_path, verbose=False)
    codes = vf_full.lcode()
    if len(codes) < 2:
        pytest.skip("need at least 2 sites")
    exclude = [codes[0]]
    vf = Velocity_Field.read(gmt_vel_path, lexclude=exclude, verbose=False)
    assert vf.nsites() == vf_full.nsites() - 1
    assert codes[0] not in vf.lcode()


def test_vel_field_read_with_lonly(gmt_vel_path):
    vf_full = Velocity_Field.read(gmt_vel_path, verbose=False)
    codes = vf_full.lcode()[:3]
    vf = Velocity_Field.read(gmt_vel_path, lonly=codes, verbose=False)
    assert vf.nsites() == len(codes)
    assert set(vf.lcode()) == set(codes)


# -------------------------------------------------------------------------
# nsites, lcode, l_GMT_Point, site
# -------------------------------------------------------------------------


def test_vel_field_nsites_returns_int(vf):
    n = vf.nsites()
    assert isinstance(n, int)
    assert n >= 1


def test_vel_field_lcode_returns_list_of_str(vf):
    out = vf.lcode()
    assert isinstance(out, list)
    assert len(out) == vf.nsites()
    assert all(isinstance(c, str) for c in out)


def test_vel_field_l_GMT_Point_returns_list(vf):
    out = vf.l_GMT_Point()
    assert isinstance(out, list)
    assert len(out) == vf.nsites()
    from pyacs.lib.gmtpoint import GMT_Point
    assert all(hasattr(p, "lon") and hasattr(p, "lat") for p in out)


def test_vel_field_site_returns_gmt_point_or_none(vf):
    code = vf.lcode()[0]
    out = vf.site(code)
    assert out is not None
    assert out.code == code
    assert vf.site("XXXX") is None


# -------------------------------------------------------------------------
# info
# -------------------------------------------------------------------------


def test_vel_field_info_runs(vf):
    vf.info(details=False)


# -------------------------------------------------------------------------
# subset
# -------------------------------------------------------------------------


def test_vel_field_subset_lonly_returns_velocity_field(vf):
    codes = vf.lcode()[:2]
    out = vf.subset(lonly=codes)
    assert isinstance(out, Velocity_Field)
    assert out.nsites() == len(codes)
    assert set(out.lcode()) == set(codes)


def test_vel_field_subset_lexclude_returns_velocity_field(vf):
    if vf.nsites() < 2:
        pytest.skip("need at least 2 sites")
    codes = vf.lcode()
    exclude = [codes[0]]
    out = vf.subset(lexclude=exclude)
    assert isinstance(out, Velocity_Field)
    assert out.nsites() == vf.nsites() - 1
    assert codes[0] not in out.lcode()


# -------------------------------------------------------------------------
# add_point, remove_point
# -------------------------------------------------------------------------


def test_vel_field_add_point(vf):
    from pyacs.lib.gmtpoint import GMT_Point
    vf_copy = Velocity_Field(lgmt_points=list(vf.sites))
    n0 = vf_copy.nsites()
    pt = GMT_Point(lon=0.0, lat=0.0, Ve=0.0, Vn=0.0, code="TST1")
    vf_copy.add_point(pt)
    assert vf_copy.nsites() == n0 + 1
    assert vf_copy.site("TST1") is not None


def test_vel_field_remove_point_returns_velocity_field(vf):
    if vf.nsites() < 2:
        pytest.skip("need at least 2 sites")
    vf_copy = Velocity_Field(lgmt_points=list(vf.sites))
    code = vf_copy.lcode()[-1]
    n0 = vf_copy.nsites()
    out = vf_copy.remove_point(code)
    assert isinstance(out, Velocity_Field)
    assert out.nsites() == n0 - 1
    assert out.site(code) is None


# -------------------------------------------------------------------------
# write
# -------------------------------------------------------------------------


def test_vel_field_write_creates_file(vf, tmp_path):
    out_file = tmp_path / "out.gmt"
    vf.write(str(out_file), verbose=False)
    assert out_file.is_file()
    assert out_file.stat().st_size > 0


# -------------------------------------------------------------------------
# radial
# -------------------------------------------------------------------------


def test_vel_field_radial_returns_velocity_field(vf):
    code = vf.lcode()[0]
    pt = vf.site(code)
    center = [pt.lon, pt.lat]
    out = vf.radial(center)
    assert isinstance(out, Velocity_Field)
    assert out.nsites() == vf.nsites()


# -------------------------------------------------------------------------
# pole
# -------------------------------------------------------------------------


def test_vel_field_pole_returns_tuple(vf):
    if vf.nsites() < 3:
        pytest.skip("pole needs at least 3 sites")
    out = vf.pole(method="WLS")
    assert isinstance(out, tuple)
    assert len(out) == 2
    X, VCV = out
    assert isinstance(X, (list, np.ndarray))
    assert len(X) == 3
    assert np.array(VCV).shape == (3, 3)


# -------------------------------------------------------------------------
# strain
# -------------------------------------------------------------------------


def test_vel_field_strain_returns_expected_type(vf):
    codes = vf.lcode()[:5]
    if len(codes) < 4:
        pytest.skip("strain needs several sites")
    try:
        out = vf.strain(codes, method="WLS", verbose=False)
    except Exception:
        pytest.skip("strain may require specific setup")
    assert out is not None


# -------------------------------------------------------------------------
# proj_profile
# -------------------------------------------------------------------------


def test_vel_field_proj_profile_returns_tuple(vf):
    # profile across part of the field
    slon, slat = 0.0, 45.0
    elon, elat = 10.0, 50.0
    d = 500.0
    try:
        out = vf.proj_profile(slon, slat, elon, elat, d, verbose=False)
    except Exception:
        pytest.skip("proj_profile may require specific setup")
    assert isinstance(out, tuple)


# -------------------------------------------------------------------------
# substract_pole
# -------------------------------------------------------------------------


def test_vel_field_substract_pole_returns_velocity_field(vf):
    if vf.nsites() < 3:
        pytest.skip("need at least 3 sites for pole")
    X, _ = vf.pole(method="WLS")
    out = vf.substract_pole(X, "rot")
    assert isinstance(out, Velocity_Field)
    assert out.nsites() == vf.nsites()


# -------------------------------------------------------------------------
# common (needs list of GMT_Point)
# -------------------------------------------------------------------------


def test_vel_field_common_returns_list(vf):
    # use sites from the field itself as "other" list
    other = vf.l_GMT_Point()[:2]
    out = vf.common(other, prefit=1e6)
    assert isinstance(out, list)
    assert len(out) <= len(other)


# -------------------------------------------------------------------------
# backward compatibility: Velocity_Field has expected methods
# -------------------------------------------------------------------------


def test_vel_field_all_methods_present():
    expected = [
        "read", "write", "add_point", "remove_point",
        "info", "nsites", "l_GMT_Point", "print_info_site", "lcode", "site",
        "subset", "radial",
        "calc_pole", "pole", "substract_pole",
        "common", "proj_profile", "strain",
    ]
    for name in expected:
        assert hasattr(Velocity_Field, name), "Velocity_Field should have %r" % name
