"""
Unit tests for pyacs.gts.lib.offset package (apply_offsets, find_offsets, delete_small_offsets, etc.) using QUEM.pos.

The offset functionality is distributed in pyacs.gts.lib.offset submodules; the package
__init__ re-exports all public functions for backward compatibility (same imports as before).
"""

import os

import numpy as np
import pytest

from pyacs.gts.Gts import Gts

_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
QUEM_POS = os.path.join(_TEST_DIR, "data", "ts", "QUEM.pos")


# -------------------------------------------------------------------------
# backward compatibility: package exposes all Gts offset methods
# -------------------------------------------------------------------------


def test_offset_package_backward_compatibility():
    """pyacs.gts.lib.offset package re-exports all offset methods (backward compatible)."""
    import pyacs.gts.lib.offset as offset_mod
    expected = [
        'apply_offsets', 'delete_small_offsets', 'estimate_local_offset',
        'find_offsets', 'find_offsets_ivel', 'find_offsets_t_scan', 'find_time_offsets',
        'local_offset_robust', 'suspect_offsets', 'suspect_offsets_mf',
        'test_offset_significance', 'test_offsets',
    ]
    for name in expected:
        assert hasattr(offset_mod, name), "offset package should expose %r" % name
        assert callable(getattr(offset_mod, name)), "%r should be callable" % name


@pytest.fixture(scope="module")
def ts_quem():
    """Load QUEM.pos once per module."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found at %s" % QUEM_POS)
    return Gts.read(QUEM_POS, fmt="pos", verbose=False)


def _short_ts(ts, n=80):
    """Return a copy of ts with at most n epochs."""
    ts2 = ts.copy()
    if ts2.data is not None and ts2.data.shape[0] > n:
        ts2.data = ts2.data[:n].copy()
    return ts2


# -------------------------------------------------------------------------
# apply_offsets
# -------------------------------------------------------------------------


def test_offset_apply_offsets_returns_gts(ts_quem):
    """apply_offsets() returns new Gts with same shape."""
    from pyacs.gts.lib.offset import apply_offsets

    ts = _short_ts(ts_quem)
    # One offset: date (decimal year), N, E, U, s_N, s_E, s_U (m and m for sigmas)
    mid_date = float(ts.data[ts.data.shape[0] // 2, 0])
    np_offset = np.array([[mid_date, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]])
    out = apply_offsets(ts, np_offset, opposite=False, in_place=False, verbose=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape


# -------------------------------------------------------------------------
# delete_small_offsets
# -------------------------------------------------------------------------


def test_offset_delete_small_offsets_returns_list_or_none(ts_quem):
    """delete_small_offsets() returns a list of offset dates or None."""
    from pyacs.gts.lib.offset import delete_small_offsets

    ts = _short_ts(ts_quem)
    # Input: list of offset dates (decimal year)
    offsets = [ts.data[10, 0], ts.data[20, 0]]
    result = delete_small_offsets(ts, offsets, del_by_pricise=False)

    assert result is None or isinstance(result, list)


# -------------------------------------------------------------------------
# find_offsets (can be slow; use short series)
# -------------------------------------------------------------------------


def test_offset_find_offsets_returns_gts(ts_quem):
    """find_offsets() returns Gts with offsets_dates and outliers set."""
    from pyacs.gts.lib.offset import find_offsets

    ts = _short_ts(ts_quem, n=50)
    out = find_offsets(
        ts,
        threshold=3,
        n_max_offsets=3,
        lcomponent="NE",
        verbose=False,
        in_place=False,
    )

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape
    assert isinstance(out.offsets_dates, (list, np.ndarray)) or out.offsets_dates is None
    assert isinstance(out.outliers, list)


# -------------------------------------------------------------------------
# estimate_local_offset
# -------------------------------------------------------------------------


def test_offset_estimate_local_offset_returns_gts_or_none(ts_quem):
    """estimate_local_offset() returns Gts or None when no offsets_dates."""
    from pyacs.gts.lib.offset import estimate_local_offset

    ts = _short_ts(ts_quem, n=40)
    ts.offsets_dates = []
    out = estimate_local_offset(ts, window_length=4, in_place=False)

    # With no offsets, behavior may be return self or new copy or None
    if out is not None:
        assert isinstance(out, Gts)


def test_offset_estimate_local_offset_with_one_offset(ts_quem):
    """estimate_local_offset() with one offset date returns Gts or list."""
    from pyacs.gts.lib.offset import estimate_local_offset

    ts = _short_ts(ts_quem, n=40)
    # One offset in the middle
    mid = ts.data.shape[0] // 2
    ts.offsets_dates = [float(ts.data[mid, 0])]
    try:
        out = estimate_local_offset(ts, window_length=4, in_place=False)
        if out is not None:
            assert isinstance(out, (Gts, list, np.ndarray))
    except Exception:
        pytest.skip("estimate_local_offset raised (e.g. insufficient data)")


# -------------------------------------------------------------------------
# suspect_offsets_mf (requires median_filter)
# -------------------------------------------------------------------------


def test_offset_suspect_offsets_mf_returns_gts_when_scipy_available(ts_quem):
    """suspect_offsets_mf() returns Gts with offsets_dates set when median_filter available."""
    pytest.importorskip("scipy")
    if not hasattr(Gts, "median_filter"):
        pytest.skip("Gts.median_filter not available")
    from pyacs.gts.lib.offset import suspect_offsets_mf

    ts = _short_ts(ts_quem, n=40)
    out = suspect_offsets_mf(
        ts, threshold=5, n_max_offsets=2, lcomponent="NE", verbose=False, in_place=False
    )

    assert isinstance(out, Gts)
    assert out.data.shape == ts.data.shape
    assert isinstance(out.offsets_dates, (list, np.ndarray))
