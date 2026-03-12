"""
Unit tests for pyacs.gts.lib.outliers (remove_outliers, find_outliers_*, etc.) using QUEM.pos.
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


def _short_ts(ts, n=80):
    """Return a copy of ts with at most n epochs."""
    ts2 = ts.copy()
    if ts2.data is not None and ts2.data.shape[0] > n:
        ts2.data = ts2.data[:n].copy()
    return ts2


# -------------------------------------------------------------------------
# remove_outliers
# -------------------------------------------------------------------------


def test_outliers_remove_outliers_no_outliers_returns_copy(ts_quem):
    """remove_outliers() with no outliers returns new Gts with same data."""
    from pyacs.gts.lib.outliers.remove_outliers import remove_outliers

    ts = _short_ts(ts_quem)
    ts.outliers = []
    out = remove_outliers(ts, in_place=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape == ts.data.shape
    assert out.outliers == []


def test_outliers_remove_outliers_with_outliers_removes_rows(ts_quem):
    """remove_outliers() with outliers set returns Gts with fewer rows."""
    from pyacs.gts.lib.outliers.remove_outliers import remove_outliers

    ts = _short_ts(ts_quem, n=50)
    # Mark first two indices as outliers
    ts.outliers = [0, 1]
    out = remove_outliers(ts, in_place=False)

    assert isinstance(out, Gts)
    assert out is not ts
    assert out.data.shape[0] == ts.data.shape[0] - 2
    assert out.data.shape[1] == ts.data.shape[1]
    assert out.outliers == []


# -------------------------------------------------------------------------
# find_outliers_simple
# -------------------------------------------------------------------------


def test_outliers_find_outliers_simple_returns_gts(ts_quem):
    """find_outliers_simple() returns Gts and may set outliers attribute."""
    from pyacs.gts.lib.outliers.find_outliers_simple import find_outliers_simple

    ts = _short_ts(ts_quem)
    # High threshold so we likely get no or few outliers on clean data
    out = find_outliers_simple(
        ts, threshold=100, window_length=10, in_place=False, verbose=False
    )

    assert isinstance(out, Gts)
    assert out.data.shape == ts.data.shape
    assert isinstance(out.outliers, list)


# -------------------------------------------------------------------------
# find_outliers_percentage
# -------------------------------------------------------------------------


def test_outliers_find_outliers_percentage_returns_gts(ts_quem):
    """find_outliers_percentage() returns Gts with outliers list."""
    from pyacs.gts.lib.outliers.find_outliers_percentage import find_outliers_percentage

    ts = _short_ts(ts_quem, n=60)
    out = find_outliers_percentage(
        ts, percentage=0.03, in_place=False, verbose=False, component="NEU"
    )

    assert isinstance(out, Gts)
    assert out.data.shape == ts.data.shape
    assert isinstance(out.outliers, list)


# -------------------------------------------------------------------------
# find_outliers_sliding_window
# -------------------------------------------------------------------------


def test_outliers_find_outliers_sliding_window_returns_gts(ts_quem):
    """find_outliers_sliding_window() returns Gts."""
    from pyacs.gts.lib.outliers.find_outliers_sliding_window import (
        find_outliers_sliding_window,
    )

    ts = _short_ts(ts_quem, n=50)
    out = find_outliers_sliding_window(
        ts, threshold=5, window_len=10, in_place=False, verbose=False
    )

    assert isinstance(out, Gts)
    assert out.data.shape == ts.data.shape
    assert isinstance(out.outliers, list)
