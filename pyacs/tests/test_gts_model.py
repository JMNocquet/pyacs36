"""
Unit tests for pyacs.gts.lib.model (Gts detrend, make_model, mmodel, remove_pole, frame, trajectory, etc.)
using pyacs/tests/data/ts/QUEM.pos as toy data.
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


# -----------------------------------------------------------------------------
# detrend
# -----------------------------------------------------------------------------


def test_model_detrend_returns_gts(ts_quem):
    """detrend(method='L2') returns new Gts with velocity set; original unchanged."""
    out = ts_quem.detrend(method="L2")
    assert out is not None
    assert isinstance(out, Gts)
    assert out is not ts_quem
    assert out.velocity is not None
    assert len(out.velocity) >= 3
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# detrend_annual
# -----------------------------------------------------------------------------


def test_model_detrend_annual_returns_gts(ts_quem):
    """detrend_annual(method='L2') returns Gts with velocity and annual set."""
    out = ts_quem.detrend_annual(method="L2")
    assert out is not None
    assert out.velocity is not None
    assert out.annual is not None
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# detrend_seasonal
# -----------------------------------------------------------------------------


def test_model_detrend_seasonal_returns_gts(ts_quem):
    """detrend_seasonal(method='L2') returns Gts with velocity, annual, semi_annual set."""
    out = ts_quem.detrend_seasonal(method="L2")
    assert out is not None
    assert out.velocity is not None
    assert out.annual is not None
    assert out.semi_annual is not None
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# remove_pole
# -----------------------------------------------------------------------------


def test_model_remove_pole_returns_gts(ts_quem):
    """remove_pole(pole, pole_type='euler') returns Gts (velocity removed)."""
    pole = np.array([-132.21, -18.83, 0.121])
    out = ts_quem.remove_pole(pole, pole_type="euler")
    assert out is not None
    assert isinstance(out, Gts)
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# frame
# -----------------------------------------------------------------------------


def test_model_frame_returns_gts(ts_quem):
    """frame(frame='soam') returns Gts in requested frame."""
    out = ts_quem.frame(frame="soam")
    assert out is not None
    assert isinstance(out, Gts)
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# make_model
# -----------------------------------------------------------------------------


def test_model_make_model_detrend(ts_quem):
    """make_model(option='detrend', method='L2') returns Gts with residuals and velocity."""
    out = ts_quem.make_model(option="detrend", method="L2")
    assert out is not None
    assert out.velocity is not None
    assert out.data.shape[0] == ts_quem.data.shape[0]
    assert out.data.shape[1] == ts_quem.data.shape[1]


def test_model_make_model_detrend_annual(ts_quem):
    """make_model(option='detrend_annual') returns Gts with annual term."""
    out = ts_quem.make_model(option="detrend_annual", method="L2")
    assert out is not None
    assert out.annual is not None
    assert out.velocity is not None


def test_model_make_model_detrend_seasonal(ts_quem):
    """make_model(option='detrend_seasonal') returns Gts with semi_annual."""
    out = ts_quem.make_model(option="detrend_seasonal", method="L2")
    assert out is not None
    assert out.semi_annual is not None


# -----------------------------------------------------------------------------
# mmodel
# -----------------------------------------------------------------------------


def test_model_mmodel_returns_gts(ts_quem):
    """mmodel() after make_model returns model time series (same shape as data)."""
    ts_with_vel = ts_quem.make_model(option="detrend", method="L2")
    model_ts = ts_with_vel.mmodel()
    assert model_ts is not None
    assert isinstance(model_ts, Gts)
    assert model_ts.data.shape[0] == ts_quem.data.shape[0]
    assert model_ts.data.shape[1] == ts_quem.data.shape[1]


# -----------------------------------------------------------------------------
# detrend_median
# -----------------------------------------------------------------------------


def test_model_detrend_median_returns_gts_or_none(ts_quem):
    """detrend_median() returns Gts or None (None if series too short for 1-year pairs)."""
    out = ts_quem.detrend_median(delta_day=None, verbose=False)
    if out is not None:
        assert isinstance(out, Gts)
        assert out.velocity is not None
    assert out is None or out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# detrend_seasonal_median
# -----------------------------------------------------------------------------


def test_model_detrend_seasonal_median_returns_gts(ts_quem):
    """detrend_seasonal_median(wl=11) returns Gts (unchanged if <3 years)."""
    out = ts_quem.detrend_seasonal_median(wl=11, verbose=False)
    assert out is not None
    assert isinstance(out, Gts)
    assert out.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# trajectory
# -----------------------------------------------------------------------------


def test_model_trajectory_trend_returns_tuple(ts_quem):
    """trajectory(model_type='trend') returns (results_dict, model_Gts, residual_Gts, daily_Gts)."""
    pytest.importorskip("scipy")
    result = ts_quem.trajectory(model_type="trend", verbose=False)
    assert isinstance(result, (tuple, list))
    assert len(result) == 4
    results_dict, model_gts, residual_gts, daily_gts = result
    assert model_gts is not None and isinstance(model_gts, Gts)
    assert residual_gts is not None and isinstance(residual_gts, Gts)
    assert residual_gts.data.shape[0] == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# detrend_hectorp (optional: requires hectorp)
# -----------------------------------------------------------------------------


def test_model_detrend_hectorp_skip_if_not_installed(ts_quem):
    """detrend_hectorp skips or runs; if runs, returns Gts with velocity."""
    import subprocess
    try:
        subprocess.run(
            ["estimatetrend", "--help"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=5,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired, Exception):
        pytest.skip("hectorp (estimatetrend) not installed or not runnable")
    out = ts_quem.extract_periods([2008.,2015.]).detrend_hectorp(component="NEU")
    assert out is not None
    assert isinstance(out, Gts)
    assert out.velocity is not None


# -----------------------------------------------------------------------------
# detrend_pytrf (optional: requires pytrf)
# -----------------------------------------------------------------------------


def test_model_detrend_pytrf_skip_if_not_installed(ts_quem):
    """detrend_pytrf skips if pytrf not installed; else returns Gts."""
    try:
        from pytrf.ts import ts as pytrf_ts
    except ImportError:
        pytest.skip("pytrf not installed")
    out = ts_quem.extract_periods([2008.,2015.]).detrend_pytrf()
    assert out is not None
    assert isinstance(out, Gts)
