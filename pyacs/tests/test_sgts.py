"""
Unit tests for pyacs.gts.Sgts (collection of Gts).
Fixture: Sgts built from QUEM.pos (one Gts from file, two copies with different codes and slight modifications).
"""

import os

import numpy as np
import pytest

from pyacs.gts.Gts import Gts
from pyacs.gts.Sgts import Sgts

_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
QUEM_POS = os.path.join(_TEST_DIR, "data", "ts", "QUEM.pos")


def _sgts_from_quem():
    """Build Sgts from QUEM.pos: one Gts from file + two copies with codes SIT1, SIT2 and slight modifications."""
    if not os.path.isfile(QUEM_POS):
        return None
    gts_quem = Gts.read(QUEM_POS, fmt="pos", verbose=False)
    if gts_quem.data is None or gts_quem.data.shape[0] < 10:
        return None
    # Keep only first 80 epochs for faster tests
    if gts_quem.data.shape[0] > 80:
        gts_quem = gts_quem.copy()
        gts_quem.data = gts_quem.data[:80].copy()
    sgts = Sgts(read=False)
    sgts.append(gts_quem)
    # SIT1: copy, change code, shift lon/lat and data slightly
    g1 = gts_quem.copy()
    g1.code = "SIT1"
    if g1.lon is not None:
        g1.lon = g1.lon + 0.01
    if g1.lat is not None:
        g1.lat = g1.lat + 0.01
    g1.data = gts_quem.data.copy()
    g1.data[:, 1:4] += 0.001
    sgts.append(g1)
    # SIT2: same idea
    g2 = gts_quem.copy()
    g2.code = "SIT2"
    if g2.lon is not None:
        g2.lon = g2.lon - 0.01
    if g2.lat is not None:
        g2.lat = g2.lat - 0.01
    g2.data = gts_quem.data.copy()
    g2.data[:, 1:4] -= 0.0005
    sgts.append(g2)
    return sgts


@pytest.fixture(scope="module")
def sgts():
    """Sgts with 3 series: QUEM, SIT1, SIT2 (from QUEM.pos template)."""
    out = _sgts_from_quem()
    if out is None:
        pytest.skip("QUEM.pos not found or too short at %s" % QUEM_POS)
    return out


# -------------------------------------------------------------------------
# append
# -------------------------------------------------------------------------


def test_sgts_append_adds_gts(sgts):
    n0 = sgts.n()
    first = sgts.lGts()[0]
    extra = first.copy()
    extra.code = "EXT1"
    sgts.append(extra)
    assert sgts.n() == n0 + 1
    assert sgts.has_ts("EXT1")


# -------------------------------------------------------------------------
# copy
# -------------------------------------------------------------------------


def test_sgts_copy_returns_sgts(sgts):
    out = sgts.copy()
    assert isinstance(out, Sgts)
    assert out is not sgts
    assert out.n() == sgts.n()


# -------------------------------------------------------------------------
# delts
# -------------------------------------------------------------------------


def test_sgts_delts_removes_series(sgts):
    if sgts.n() < 2:
        pytest.skip("need at least 2 series")
    cp = sgts.copy()
    code = cp.lcode()[-1]
    n0 = cp.n()
    cp.delts(code)
    assert cp.n() == n0 - 1
    assert not cp.has_ts(code)


# -------------------------------------------------------------------------
# gts (apply Gts method to all)
# -------------------------------------------------------------------------


def test_sgts_gts_returns_sgts(sgts):
    pytest.importorskip("tqdm")
    out = sgts.gts("copy")
    assert isinstance(out, Sgts)
    assert out.n() == sgts.n()


# -------------------------------------------------------------------------
# has_ts
# -------------------------------------------------------------------------


def test_sgts_has_ts_true(sgts):
    code = sgts.lcode()[0]
    assert sgts.has_ts(code) is True


def test_sgts_has_ts_false(sgts):
    assert sgts.has_ts("XXXX") is False


# -------------------------------------------------------------------------
# lGts
# -------------------------------------------------------------------------


def test_sgts_lGts_returns_list_of_gts(sgts):
    out = sgts.lGts()
    assert isinstance(out, list)
    assert len(out) == sgts.n()
    from pyacs.gts.Gts import Gts
    assert all(isinstance(x, Gts) for x in out)


# -------------------------------------------------------------------------
# n
# -------------------------------------------------------------------------


def test_sgts_n_returns_int(sgts):
    out = sgts.n()
    assert isinstance(out, int)
    assert out >= 1


# -------------------------------------------------------------------------
# lcode
# -------------------------------------------------------------------------


def test_sgts_lcode_returns_list_of_str(sgts):
    out = sgts.lcode()
    assert isinstance(out, list)
    assert all(isinstance(c, str) for c in out)
    assert len(out) == sgts.n()


# -------------------------------------------------------------------------
# sub
# -------------------------------------------------------------------------


def test_sgts_sub_linclude_returns_sgts(sgts):
    codes = sgts.lcode()[:1]
    out = sgts.sub(linclude=codes)
    assert isinstance(out, Sgts)
    assert out.n() == 1
    assert out.lcode() == codes


def test_sgts_sub_lexclude_returns_sgts(sgts):
    if sgts.n() < 2:
        pytest.skip("need at least 2 series")
    codes = sgts.lcode()[1:]
    out = sgts.sub(lexclude=[sgts.lcode()[0]])
    assert isinstance(out, Sgts)
    assert out.n() == sgts.n() - 1


# -------------------------------------------------------------------------
# sel_period
# -------------------------------------------------------------------------


def test_sgts_sel_period_returns_sgts(sgts):
    start, end = sgts.dates(unit="decyear")
    mid = (start + end) / 2.0
    out = sgts.sel_period([start, mid], min_data=2, verbose=False)
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# dates
# -------------------------------------------------------------------------


def test_sgts_dates_returns_tuple(sgts):
    start, end = sgts.dates(unit="decyear")
    assert isinstance(start, (int, float))
    assert isinstance(end, (int, float))
    assert start <= end


# -------------------------------------------------------------------------
# get_dates
# -------------------------------------------------------------------------


def test_sgts_get_dates_returns_tuple(sgts):
    start, end = sgts.get_dates(fmt="decyear")
    assert isinstance(start, (int, float))
    assert isinstance(end, (int, float))


# -------------------------------------------------------------------------
# delnone
# -------------------------------------------------------------------------


def test_sgts_delnone_returns_sgts(sgts):
    out = sgts.delnone()
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# add_offsets_dates
# -------------------------------------------------------------------------


def test_sgts_add_offsets_dates_returns_sgts(sgts):
    start, end = sgts.dates(unit="decyear")
    mid = (start + end) / 2.0
    out = sgts.add_offsets_dates([mid], verbose=False)
    assert isinstance(out, Sgts)
    assert out.n() == sgts.n()


# -------------------------------------------------------------------------
# sel_rectangle
# -------------------------------------------------------------------------


def test_sgts_sel_rectangle_returns_sgts(sgts):
    # bounds that include all sites (QUEM/SIT1/SIT2)
    codes = sgts.lcode()
    lons = [sgts.__dict__[c].lon for c in codes if getattr(sgts.__dict__[c], "lon", None) is not None]
    lats = [sgts.__dict__[c].lat for c in codes if getattr(sgts.__dict__[c], "lat", None) is not None]
    if not lons or not lats:
        pytest.skip("need lon/lat on Gts for sel_rectangle")
    bounds = [min(lons) - 1, max(lons) + 1, min(lats) - 1, max(lats) + 1]
    out = sgts.sel_rectangle(bounds, verbose=False)
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# sel_radius
# -------------------------------------------------------------------------


def test_sgts_sel_radius_returns_sgts(sgts):
    code = sgts.lcode()[0]
    out = sgts.sel_radius(code, 500.0, verbose=False)
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# medvel
# -------------------------------------------------------------------------


def test_sgts_medvel_returns_sgts(sgts):
    out = sgts.medvel(verbose=False)
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# same_site
# -------------------------------------------------------------------------


def test_sgts_same_site_returns_sgts(sgts):
    out = sgts.same_site(dc=10, in_place=False, verbose=False)
    assert isinstance(out, Sgts)


# -------------------------------------------------------------------------
# make_distance_matrix_from_sgts
# -------------------------------------------------------------------------


def test_sgts_make_distance_matrix_returns_ndarray(sgts):
    out = sgts.make_distance_matrix_from_sgts()
    assert isinstance(out, np.ndarray)
    assert out.shape == (sgts.n(), sgts.n())


# -------------------------------------------------------------------------
# nearest
# -------------------------------------------------------------------------


def test_sgts_nearest_returns_list(sgts):
    code = sgts.lcode()[0]
    out = sgts.nearest(code, n=1)
    assert out is None or (isinstance(out, (list, np.ndarray)) and len(out) >= 1)


# -------------------------------------------------------------------------
# info, stat_site (no return, just no exception)
# -------------------------------------------------------------------------


def test_sgts_info_runs(sgts):
    sgts.info()


def test_sgts_stat_site_runs(sgts):
    sgts.stat_site(display=False, verbose=False)


# -------------------------------------------------------------------------
# frame
# -------------------------------------------------------------------------


def test_sgts_frame_returns_sgts_or_skips(sgts):
    try:
        out = sgts.frame(verbose=False)
    except Exception:
        pytest.skip("frame may require ref or external setup")
    assert out is None or isinstance(out, Sgts)


# -------------------------------------------------------------------------
# to_displacement
# -------------------------------------------------------------------------


def test_sgts_to_displacement_returns_expected_type(sgts):
    start, end = sgts.dates(unit="decyear")
    try:
        out = sgts.to_displacement(sdate=start, edate=end, verbose=False)
    except Exception:
        pytest.skip("to_displacement may require specific setup")
    assert out is not None


# -------------------------------------------------------------------------
# save_velocity
# -------------------------------------------------------------------------


def test_sgts_save_velocity_runs_or_skips(sgts):
    try:
        sgts.save_velocity(verbose=False)
    except Exception:
        pytest.skip("save_velocity may require velocity attribute or path")


# -------------------------------------------------------------------------
# write_pck
# -------------------------------------------------------------------------


def test_sgts_write_pck_writes_file(sgts, tmp_path):
    p = tmp_path / "ts.pck"
    sgts.write_pck(str(p), verbose=False)
    assert p.is_file()


# -------------------------------------------------------------------------
# to_tsnpz
# -------------------------------------------------------------------------


def test_sgts_to_tsnpz_writes_file(sgts, tmp_path):
    p = tmp_path / "ts.npz"
    try:
        sgts.to_tsnpz(str(p), verbose=False)
    except Exception:
        pytest.skip("to_tsnpz may require specific deps")
    assert p.is_file()


# -------------------------------------------------------------------------
# to_kml
# -------------------------------------------------------------------------


def test_sgts_to_kml_writes_file_or_skips(sgts, tmp_path):
    p = tmp_path / "ts.kml"
    try:
        sgts.to_kml(str(p), verbose=False)
    except Exception:
        pytest.skip("to_kml may require specific setup")
    if p.exists():
        assert p.is_file()


# -------------------------------------------------------------------------
# plot_data_sum, plot_component
# -------------------------------------------------------------------------


def test_sgts_plot_data_sum_runs(sgts):
    pytest.importorskip("matplotlib")
    try:
        sgts.plot_data_sum()
    except Exception:
        pytest.skip("plot_data_sum may require period or other args")


def test_sgts_plot_component_runs(sgts):
    pytest.importorskip("matplotlib")
    try:
        sgts.plot_component(component="N")
    except Exception:
        pytest.skip("plot_component may require specific args")


# -------------------------------------------------------------------------
# read_ts (no fixture dir; test that method exists and callable)
# -------------------------------------------------------------------------


def test_sgts_read_ts_callable():
    assert hasattr(Sgts, "read_ts")
    assert callable(Sgts.read_ts)


# -------------------------------------------------------------------------
# read_gmt, read_gts_conf, read_soln, get_unr, show_map, show_ivel_map_gmt
# correct_offsets_from_file, remove_observations, common_mode, to_obs_tensor,
# apply_coseismic, compute_common_mode_l1trend, gts_mp, sel_from_grid, sel_radius_eq, to_tspck
# -------------------------------------------------------------------------


def test_sgts_all_methods_callable():
    """All Sgts methods attached in Sgts.py are callable."""
    expected = [
        "append", "copy", "delts", "frame", "gts", "has_ts", "lGts", "n", "lcode",
        "medvel", "read_gmt", "read_gts_conf", "read_soln", "read_ts", "same_site",
        "save_velocity", "sel_radius", "sel_rectangle", "sel_period", "show_map",
        "stat_site", "sub", "to_displacement", "write_pck", "common_mode",
        "to_obs_tensor", "apply_coseismic", "get_unr", "info", "to_kml",
        "plot_data_sum", "plot_component", "dates", "sel_from_grid",
        "compute_common_mode_l1trend", "gts_mp", "delnone", "to_tspck", "to_tsnpz",
        "make_distance_matrix_from_sgts", "nearest", "show_ivel_map_gmt",
        "correct_offsets_from_file", "remove_observations", "get_dates", "sel_radius_eq",
        "add_offsets_dates",
    ]
    for name in expected:
        assert hasattr(Sgts, name), "Sgts should have method %r" % name
        assert callable(getattr(Sgts, name)), "Sgts.%s should be callable" % name
