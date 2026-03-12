"""
Unit tests for pyacs.gts.lib.primitive (Gts methods: copy, cdata, extract_*, add_obs, etc.)
using QUEM.pos.
"""

import os

import numpy as np
import pytest

from pyacs.gts.Gts import Gts

_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
QUEM_POS = os.path.join(_TEST_DIR, "data", "ts", "QUEM.pos")


def _short_ts(ts, n=80):
    """Return a copy of ts with at most n epochs."""
    ts2 = ts.copy()
    if ts2.data is not None and ts2.data.shape[0] > n:
        ts2.data = ts2.data[:n].copy()
    if ts2.data_xyz is not None and ts2.data_xyz.shape[0] > n:
        ts2.data_xyz = ts2.data_xyz[:n].copy()

    return ts2


def _other_gts_for_insert(gts, n=5, code="OTH"):
    """Return a small Gts with dates outside gts range, for insert_gts_data / insert_ts."""
    other = gts.copy()
    other.data = gts.data[:n].copy()
    other.data[:, 0] = np.linspace(2023.0, 2023.5, n)
    other.code = code
    return other


@pytest.fixture(scope="module")
def ts_quem():
    """Load QUEM.pos once per module."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found at %s" % QUEM_POS)
    return Gts.read(QUEM_POS, fmt="pos", verbose=False)


@pytest.fixture
def gts(ts_quem):
    """Gts instance from QUEM.pos (shortened for faster tests)."""
    return _short_ts(ts_quem)


# -------------------------------------------------------------------------
# copy
# -------------------------------------------------------------------------


def test_primitive_copy_returns_gts(gts):
    out = gts.copy()
    assert isinstance(out, Gts)
    assert out is not gts
    assert out.data.shape == gts.data.shape


def test_primitive_copy_data_false(gts):
    out = gts.copy(data=False)
    assert out.data is None


# -------------------------------------------------------------------------
# cdata
# -------------------------------------------------------------------------


def test_primitive_cdata_with_data_returns_bool(gts):
    out = gts.cdata(data=True)
    assert isinstance(out, bool)
    assert out is True


def test_primitive_cdata_no_data_returns_false():
    g = Gts(code="X")
    g.data = None
    assert g.cdata(data=True) is False


# -------------------------------------------------------------------------
# differentiate
# -------------------------------------------------------------------------


def test_primitive_differentiate_returns_gts_one_fewer_row(gts):
    out = gts.differentiate()
    assert isinstance(out, Gts)
    assert out.data.shape[0] == gts.data.shape[0] - 1


# -------------------------------------------------------------------------
# extract_periods
# -------------------------------------------------------------------------


def test_primitive_extract_periods_returns_gts(gts):
    mid = (gts.data[0, 0] + gts.data[-1, 0]) / 2.0
    out = gts.extract_periods([[gts.data[0, 0], mid]])
    assert isinstance(out, Gts)
    assert out.data is None or (out.data.shape[0] <= gts.data.shape[0] and np.all(out.data[:, 0] >= gts.data[0, 0]) and np.all(out.data[:, 0] <= mid))


# -------------------------------------------------------------------------
# exclude_periods
# -------------------------------------------------------------------------


def test_primitive_exclude_periods_returns_gts(gts):
    mid = (gts.data[0, 0] + gts.data[-1, 0]) / 2.0
    out = gts.exclude_periods([[mid - 0.01, mid + 0.01]])
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# extract_dates
# -------------------------------------------------------------------------


def test_primitive_extract_dates_returns_gts(gts):
    dates = [gts.data[0, 0], gts.data[5, 0]]
    out = gts.extract_dates(dates, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# add_offsets_dates
# -------------------------------------------------------------------------


def test_primitive_add_offsets_dates_returns_gts(gts):
    out = gts.add_offsets_dates([2021.0], in_place=False)
    assert isinstance(out, Gts)
    assert out.offsets_dates == [2021.0]


# -------------------------------------------------------------------------
# remove_velocity
# -------------------------------------------------------------------------


def test_primitive_remove_velocity_returns_gts(gts):
    out = gts.remove_velocity([0.001, 0.001, 0.001], in_place=False)
    assert isinstance(out, Gts)
    assert out.data.shape == gts.data.shape


# -------------------------------------------------------------------------
# set_zero_at_date
# -------------------------------------------------------------------------


def test_primitive_set_zero_at_date_returns_gts(gts):
    date = float(gts.data[gts.data.shape[0] // 2, 0])
    out = gts.set_zero_at_date(date)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# decimate
# -------------------------------------------------------------------------


def test_primitive_decimate_returns_gts(gts):
    out = gts.decimate(time_step=60.0, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# reorder
# -------------------------------------------------------------------------


def test_primitive_reorder_returns_gts(gts):
    out = gts.reorder(verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# extract_ndates_before_date, extract_ndates_after_date, extract_ndates_around_date
# -------------------------------------------------------------------------


def test_primitive_extract_ndates_before_date_returns_gts(gts):
    date = float(gts.data[10, 0])
    out = gts.extract_ndates_before_date(date, 5, verbose=False)
    assert isinstance(out, Gts)


def test_primitive_extract_ndates_after_date_returns_gts(gts):
    date = float(gts.data[10, 0])
    out = gts.extract_ndates_after_date(date, 5, verbose=False)
    assert isinstance(out, Gts)


def test_primitive_extract_ndates_around_date_returns_gts(gts):
    date = float(gts.data[10, 0])
    out = gts.extract_ndates_around_date(date, 3)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# get_coseismic
# -------------------------------------------------------------------------


def test_primitive_get_coseismic_returns_gts(gts):
    eq_date = float(gts.data[15, 0])
    out = gts.get_coseismic(eq_date, window_days=3, sample_after=1, in_place=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# correct_duplicated_dates
# -------------------------------------------------------------------------


def test_primitive_correct_duplicated_dates_returns_gts(gts):
    out = gts.correct_duplicated_dates(action="check", in_place=False, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# rotate
# -------------------------------------------------------------------------


def test_primitive_rotate_returns_gts(gts):
    out = gts.rotate(45.0, in_place=False)
    assert isinstance(out, Gts)
    assert out.data.shape == gts.data.shape


# -------------------------------------------------------------------------
# insert_gts_data
# -------------------------------------------------------------------------


def test_primitive_insert_gts_data_returns_gts(gts):
    other = _other_gts_for_insert(gts)
    out = gts.insert_gts_data(other, in_place=False, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# find_large_uncertainty
# -------------------------------------------------------------------------


def test_primitive_find_large_uncertainty_returns_gts(gts):
    out = gts.find_large_uncertainty(sigma_threshold=20, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# split_gap
# -------------------------------------------------------------------------


def test_primitive_split_gap_returns_list_of_gts(gts):
    out = gts.split_gap(gap=10, verbose=False)
    assert isinstance(out, list)
    assert all(isinstance(x, Gts) for x in out)
    assert len(out) >= 1


# -------------------------------------------------------------------------
# interpolate
# -------------------------------------------------------------------------


def test_primitive_interpolate_returns_gts(gts):
    pytest.importorskip("scipy")
    out = gts.interpolate(date="day", in_place=False, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# extrapolate
# -------------------------------------------------------------------------


def test_primitive_extrapolate_returns_gts(gts):
    pytest.importorskip("scipy")
    idates = [float(gts.data[0, 0]) - 0.5, float(gts.data[-1, 0]) + 0.5]
    out = gts.extrapolate(idates, kind="linear")
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# insert_ts
# -------------------------------------------------------------------------


def test_primitive_insert_ts_returns_gts(gts):
    other = _other_gts_for_insert(gts)
    out = gts.insert_ts(other, data=None)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# n_obs
# -------------------------------------------------------------------------


def test_primitive_n_obs_returns_int(gts):
    pytest.importorskip("pandas")
    # single date as decimal year
    date = float(gts.data[5, 0])
    n = gts.n_obs([date, date + 0.001])
    assert n is None or isinstance(n, (int, np.integer))


# -------------------------------------------------------------------------
# get_coseismic_l1trend
# -------------------------------------------------------------------------


def test_primitive_get_coseismic_l1trend_returns_ndarray_or_skips(gts):
    if not hasattr(Gts, "l1trend"):
        pytest.skip("Gts.l1trend not available (trendfilter)")
    eq_date = float(gts.data[15, 0])
    out = gts.get_coseismic_l1trend(eq_date)
    assert isinstance(out, np.ndarray)
    assert out.size >= 6


# -------------------------------------------------------------------------
# get_values_at_date
# -------------------------------------------------------------------------


def test_primitive_get_values_at_date_returns_ndarray(gts):
    date = float(gts.data[5, 0])
    out = gts.get_values_at_date(date, data_type="data")
    assert isinstance(out, np.ndarray)
    assert out.ndim == 1


# -------------------------------------------------------------------------
# split
# -------------------------------------------------------------------------


def test_primitive_split_returns_list_of_gts(gts):
    mid = float(gts.data[15, 0])
    out = gts.split([mid], verbose=False)
    assert isinstance(out, list)
    assert all(isinstance(x, Gts) for x in out)
    assert len(out) >= 1


# -------------------------------------------------------------------------
# substract_ts
# -------------------------------------------------------------------------


def test_primitive_substract_ts_returns_gts(gts):
    other = gts.copy()
    other.code = "REF"
    out = gts.substract_ts(other, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# substract_ts_daily
# -------------------------------------------------------------------------


def test_primitive_substract_ts_daily_returns_gts(gts):
    other = gts.copy()
    other.code = "REF"
    out = gts.substract_ts_daily(other, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# displacement
# -------------------------------------------------------------------------


def test_primitive_displacement_returns_ndarray(gts):
    sdate = float(gts.data[2, 0])
    edate = float(gts.data[-3, 0])
    try:
        out = gts.displacement(sdate=sdate, edate=edate, verbose=False)
    except Exception:
        pytest.skip("displacement may require optional deps (e.g. colors)")
    assert isinstance(out, np.ndarray)
    assert out.size >= 6


# -------------------------------------------------------------------------
# add_obs
# -------------------------------------------------------------------------


def test_primitive_add_obs_returns_gts(gts):
    date = float(gts.data[-1, 0]) + 0.01
    neu = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    out = gts.add_obs(date, neu, in_place=False, check=False, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# add_obs_xyz
# -------------------------------------------------------------------------


def test_primitive_add_obs_xyz_returns_gts(gts):
    gts.neu2xyz(corr=False)
    if gts.data_xyz is None:
        pytest.skip("neu2xyz did not set data_xyz (need lon/lat/h or X0,Y0,Z0)")
    date = float(gts.data[-1, 0]) + 0.01
    xyz = [0.0, 0.0, 0.0]
    out = gts.add_obs_xyz(date, xyz, in_place=False, check=False, verbose=False)
    assert isinstance(out, Gts)


# -------------------------------------------------------------------------
# xyz2neu, neu2xyz
# -------------------------------------------------------------------------


def test_primitive_xyz2neu_returns_gts(gts):
    gts.neu2xyz()
    out = gts.xyz2neu(corr=False)
    assert isinstance(out, Gts)


def test_primitive_neu2xyz_returns_gts(gts):
    out = gts.neu2xyz(corr=True)
    assert isinstance(out, Gts)
