"""
Unit tests for pyacs.gts.lib.noise (wrms, add_vel_sigma, realistic_sigma, etc.) using QUEM.pos.
"""

import os
import shutil
import subprocess

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
# wrms
# -------------------------------------------------------------------------


def test_noise_wrms_returns_ndarray_of_three(ts_quem):
    """wrms() returns array of shape (3,) for NEU."""
    from pyacs.gts.lib.noise import wrms

    ts = _short_ts(ts_quem)
    result = wrms(ts)

    assert isinstance(result, np.ndarray)
    assert result.shape == (3,)
    assert np.all(np.isfinite(result))
    assert np.all(result >= 0)


# -------------------------------------------------------------------------
# add_vel_sigma
# -------------------------------------------------------------------------


def test_noise_add_vel_sigma_requires_velocity(ts_quem):
    """add_vel_sigma() returns None when velocity is not set."""
    from pyacs.gts.lib.noise import add_vel_sigma

    ts = _short_ts(ts_quem)
    ts.velocity = None
    out = add_vel_sigma(ts, verbose=False)

    assert out is None


def test_noise_add_vel_sigma_returns_gts_when_velocity_set(ts_quem):
    """add_vel_sigma() returns new Gts with velocity sigmas when velocity set."""
    from pyacs.gts.lib.noise import add_vel_sigma

    ts = _short_ts(ts_quem)
    # Set velocity via detrend
    ts = ts.make_model(option="detrend", method="L2")
    assert ts.velocity is not None
    out = add_vel_sigma(ts, verbose=False)

    assert out is not None
    assert isinstance(out, Gts)
    assert out is not ts
    assert out.velocity is not None
    assert len(out.velocity) >= 6
    # Sigma components (indices 3,4,5) should be positive
    assert np.all(np.array(out.velocity[3:6]) > 0)


# -------------------------------------------------------------------------
# tsfit availability
# -------------------------------------------------------------------------


def test_tsfit_command_available():
    """System command tsfit is available in PATH and runnable."""
    tsfit_path = shutil.which("tsfit")
    if tsfit_path is None:
        pytest.skip("tsfit not found in PATH")
    try:
        result = subprocess.run(
            ["tsfit", "--help"],
            capture_output=True,
            timeout=5,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        pytest.fail("tsfit found but failed to run")
    # Many implementations exit 0 for --help; some may exit 1 and still print to stderr
    assert result.returncode in (0, 1) or result.stdout or result.stderr, (
        "tsfit did not run successfully (returncode=%d)" % result.returncode
    )


# -------------------------------------------------------------------------
# realistic_sigma (tsfit option may need external program; test structure only)
# -------------------------------------------------------------------------


def test_noise_realistic_sigma_returns_gts_or_self(ts_quem):
    """realistic_sigma(option='tsfit') returns Gts or None (skip if tsfit fails)."""
    from pyacs.gts.lib.noise import realistic_sigma

    ts = _short_ts(ts_quem)
    # tsfit may not be installed; accept return or skip
    try:
        out = realistic_sigma(ts, option="tsfit", in_place=False, verbose=False)
        if out is not None:
            assert isinstance(out, Gts)
    except (FileNotFoundError, OSError):
        pytest.skip("tsfit not available")
