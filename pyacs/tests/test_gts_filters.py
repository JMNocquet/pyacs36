"""
Unit tests for pyacs.gts.lib.filters functions using QUEM.pos as toy data.
"""

import os

import numpy as np
import pytest

from pyacs.gts.Gts import Gts

_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
QUEM_POS = os.path.join(_TEST_DIR, "data", "ts", "QUEM.pos")


@pytest.fixture(scope="module")
def ts_quem():
    """Load QUEM.pos once per module."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found at %s" % QUEM_POS)
    return Gts.read(QUEM_POS, fmt="pos", verbose=False)


def _short_ts(ts: Gts, n: int = 80) -> Gts:
    """Return a shallow copy of ts with at most n epochs for faster tests."""
    ts2 = ts.copy()
    if ts2.data is not None and ts2.data.shape[0] > n:
        ts2.data = ts2.data[:n].copy()
    return ts2


# -------------------------------------------------------------------------
# Median filter
# -------------------------------------------------------------------------


def test_filters_median_filter_returns_gts_and_smooths(ts_quem):
    """median_filter(n) returns new Gts with same shape and reduced variance."""
    pytest.importorskip("scipy")
    from pyacs.gts.lib.filters.median import median_filter

    ts = _short_ts(ts_quem)
    out = median_filter(ts, n=5, in_place=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape

    # Median filter should not increase variance strongly; allow equality for degenerate cases.
    var_in = np.var(ts.data[:, 1])
    var_out = np.var(out.data[:, 1])
    assert var_out <= var_in * 1.01


# -------------------------------------------------------------------------
# Savitzky-Golay filter
# -------------------------------------------------------------------------


def test_filters_savitzky_golay_returns_gts(ts_quem):
    """savitzky_golay() returns new Gts with same shape."""
    pytest.importorskip("scipy")
    from pyacs.gts.lib.filters.savitzky_golay import savitzky_golay

    ts = _short_ts(ts_quem)
    out = savitzky_golay(ts, window_length=7, polyorder=2, deriv=0)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


# -------------------------------------------------------------------------
# Wiener filter
# -------------------------------------------------------------------------


def test_filters_wiener_returns_gts(ts_quem):
    """wiener() returns new Gts with same shape (default in_place=False)."""
    pytest.importorskip("scipy")
    from pyacs.gts.lib.filters.wiener import wiener

    ts = _short_ts(ts_quem)
    out = wiener(ts, my_size=7, noise=None, in_place=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


def test_filters_wiener_in_place_modifies_self(ts_quem):
    """wiener(in_place=True) modifies the input Gts and returns None."""
    pytest.importorskip("scipy")
    from pyacs.gts.lib.filters.wiener import wiener

    ts = _short_ts(ts_quem)
    before = ts.data.copy()
    out = wiener(ts, my_size=7, noise=None, in_place=True)

    assert out is None
    assert ts.data.shape == before.shape
    # Data should have changed in at least one component.
    assert not np.allclose(ts.data[:, 1], before[:, 1])


# -------------------------------------------------------------------------
# Vondrak filter
# -------------------------------------------------------------------------


def test_filters_vondrak_returns_gts(ts_quem):
    """vondrak(fc) returns new Gts with same shape."""
    pytest.importorskip("scipy")
    from pyacs.gts.lib.filters.vondrak import vondrak

    ts = _short_ts(ts_quem)
    out = vondrak(ts, fc=0.5, in_place=False, component="NEU")

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


# -------------------------------------------------------------------------
# Minimum component filter
# -------------------------------------------------------------------------


def test_filters_minimum_component_returns_gts(ts_quem):
    """minimum_component() returns new Gts and preserves masked periods."""
    from pyacs.gts.lib.filters.minimum_component import minimum_component

    ts = _short_ts(ts_quem)
    # Mask a small central window
    t_min, t_max = ts.data[0, 0], ts.data[-1, 0]
    t_mid = 0.5 * (t_min + t_max)
    mask_period = [t_mid - 0.05, t_mid + 0.05]

    out = minimum_component(ts, mask_period=mask_period, in_place=False, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape

    # Values inside masked window should be unchanged.
    idx = np.where((ts.data[:, 0] >= mask_period[0]) & (ts.data[:, 0] <= mask_period[1]))[0]
    if idx.size > 0:
        np.testing.assert_allclose(out.data[idx, 1:4], ts.data[idx, 1:4], atol=1e-9)


# -------------------------------------------------------------------------
# disp2vel (requires trendfilter: l1trend + ivel)
# -------------------------------------------------------------------------


def test_filters_disp2vel_returns_gts_when_trendfilter_available(ts_quem):
    """disp2vel() returns new Gts (velocity) when trendfilter is installed."""
    if not hasattr(Gts, "l1trend"):
        pytest.skip("trendfilter not installed")
    from pyacs.gts.lib.filters.disp2vel import disp2vel

    ts = _short_ts(ts_quem, n=60)
    out = disp2vel(ts, alpha="auto", ndays_mf=1)

    assert isinstance(out, Gts)
    assert out is not ts
    # ivel returns one fewer epoch (derivative at mid dates)
    assert out.data.shape[0] == ts.data.shape[0] - 1
    assert out.data.shape[1] == ts.data.shape[1]


# -------------------------------------------------------------------------
# edge_filter (requires trendfilter)
# -------------------------------------------------------------------------


def test_filters_edge_filter_returns_gts_when_trendfilter_available(ts_quem):
    """edge_filter() returns new Gts when trendfilter is installed."""
    if not hasattr(Gts, "edge_filter"):
        pytest.skip("edge_filter not available")
    ts = _short_ts(ts_quem, n=60)
    out = ts.edge_filter(alpha="auto", ndays_mf=1)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


# -------------------------------------------------------------------------
# ivel (instantaneous velocity)
# -------------------------------------------------------------------------


def test_filters_ivel_returns_gts_one_fewer_epoch(ts_quem):
    """ivel() returns new Gts with one fewer epoch (velocity at mid-dates)."""
    from pyacs.gts.lib.filters.ivel import ivel

    ts = _short_ts(ts_quem)
    out = ivel(ts)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape[0] == ts.data.shape[0] - 1
    assert out.data.shape[1] == ts.data.shape[1]


# -------------------------------------------------------------------------
# l1trend (requires trendfilter)
# -------------------------------------------------------------------------


def test_filters_l1trend_returns_gts_when_trendfilter_available(ts_quem):
    """l1trend() returns new Gts when trendfilter is installed."""
    if not hasattr(Gts, "l1trend"):
        pytest.skip("trendfilter not installed")
    ts = _short_ts(ts_quem, n=50)
    out = ts.l1trend(alpha=1.0, ndays_mf=1, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


# -------------------------------------------------------------------------
# l1trendi (same as l1trend when trendfilter installed)
# -------------------------------------------------------------------------


def test_filters_l1trendi_returns_gts_when_trendfilter_available(ts_quem):
    """l1trendi() returns new Gts when trendfilter is installed."""
    if not hasattr(Gts, "l1trendi"):
        pytest.skip("trendfilter not installed")
    ts = _short_ts(ts_quem, n=50)
    out = ts.l1trendi(alpha=1.0, ndays_mf=1, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


def test_filters_l1trendi_function_from_l1trend_package(ts_quem):
    """l1trendi(ts, ...) from pyacs.gts.lib.l1trend returns new Gts when trendfilter installed."""
    try:
        from pyacs.gts.lib.l1trend import l1trendi
    except (ImportError, ModuleNotFoundError):
        pytest.skip("trendfilter not installed")
    ts = _short_ts(ts_quem, n=500)
    out = l1trendi(ts, alpha='auto', ndays_mf=1, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape
    assert out.code == ts.code


# -------------------------------------------------------------------------
# smooth
# -------------------------------------------------------------------------


def test_filters_smooth_returns_gts(ts_quem):
    """smooth() returns new Gts with same shape."""
    from pyacs.gts.lib.filters.smooth import smooth

    ts = _short_ts(ts_quem)
    # window_len must be odd and < number of points
    out = smooth(ts, window_len=5, window="hanning", in_place=False, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape
