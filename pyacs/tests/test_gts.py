"""
Unit tests for pyacs.gts.Gts using pyacs/tests/data/ts/QUEM.pos as toy data.
"""

import os
import tempfile

import numpy as np
import pytest

from pyacs.gts.Gts import Gts

# Path to toy data (relative to this test file)
_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
QUEM_POS = os.path.join(_TEST_DIR, "data", "ts", "QUEM.pos")


# -----------------------------------------------------------------------------
# Fixture
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ts_quem():
    """Load QUEM.pos once per module."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found at %s" % QUEM_POS)
    return Gts.read(QUEM_POS, fmt="pos", verbose=False)


# -----------------------------------------------------------------------------
# Read and basic attributes
# -----------------------------------------------------------------------------


def test_gts_read_pos_returns_gts():
    """Gts.read(QUEM.pos, fmt='pos') returns a Gts instance."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found")
    ts = Gts.read(QUEM_POS, fmt="pos", verbose=False)
    assert ts is not None
    assert isinstance(ts, Gts)


def test_gts_read_pos_code(ts_quem):
    """Loaded time series has code QUEM."""
    assert ts_quem.code == "QUEM"


def test_gts_read_pos_data_shape(ts_quem):
    """data array has shape (n_epochs, 10) with dec_year, N, E, U, sigmas, correlations."""
    assert ts_quem.data is not None
    assert ts_quem.data.ndim == 2
    assert ts_quem.data.shape[1] == 10
    assert ts_quem.data.shape[0] >= 1


def test_gts_read_pos_lon_lat_set(ts_quem):
    """Lon, lat, h are set from header reference position."""
    assert ts_quem.lon is not None
    assert ts_quem.lat is not None
    assert ts_quem.h is not None
    # QUEM is in Ecuador (from file header)
    assert -90 < ts_quem.lon < 0
    assert -10 < ts_quem.lat < 10


def test_gts_read_pos_t0(ts_quem):
    """t0 is first epoch decimal year."""
    assert ts_quem.t0 is not None
    assert np.isclose(ts_quem.t0, ts_quem.data[0, 0], atol=1e-9)


def test_gts_read_pos_dates_ordered(ts_quem):
    """Decimal years in data are strictly increasing."""
    dy = ts_quem.data[:, 0]
    assert np.all(np.diff(dy) > 0)


def test_gts_read_pos_reference_xyz(ts_quem):
    """X0, Y0, Z0 reference position are set."""
    assert ts_quem.X0 is not None
    assert ts_quem.Y0 is not None
    assert ts_quem.Z0 is not None
    # Rough check: QUEM reference from file header
    assert abs(ts_quem.X0 - 1272483.44) < 1.0
    assert abs(ts_quem.Y0 - (-6252975.34)) < 1.0


# -----------------------------------------------------------------------------
# Copy
# -----------------------------------------------------------------------------


def test_gts_copy(ts_quem):
    """copy() returns a new Gts with same code and data shape."""
    ts2 = ts_quem.copy()
    assert ts2 is not ts_quem
    assert ts2.code == ts_quem.code
    assert ts2.data.shape == ts_quem.data.shape
    np.testing.assert_allclose(ts2.data, ts_quem.data, atol=1e-12)


# -----------------------------------------------------------------------------
# to_pandas_df (if used by n_obs)
# -----------------------------------------------------------------------------


def test_gts_to_pandas_df(ts_quem):
    """to_pandas_df() returns a DataFrame with expected index and columns."""
    pytest.importorskip("pandas")
    df = ts_quem.to_pandas_df()
    assert df is not None
    assert len(df) == ts_quem.data.shape[0]
    # Index should be datetime-like or date
    assert len(df.index) == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# n_obs (optional, may depend on pandas)
# -----------------------------------------------------------------------------


def test_gts_n_obs_period(ts_quem):
    """n_obs with decimal-year period returns count in range."""
    pytest.importorskip("pandas")
    dy = ts_quem.data[:, 0]
    t1, t2 = float(dy.min()), float(dy.max())
    n = ts_quem.n_obs([t1, t2])
    assert n == ts_quem.data.shape[0]


# -----------------------------------------------------------------------------
# Read with wrong format / missing file
# -----------------------------------------------------------------------------


def test_gts_read_missing_file():
    """Gts.read(missing_file, fmt='pos') raises or returns None."""
    missing = os.path.join(_TEST_DIR, "data", "ts", "NONEXISTENT.pos")
    assert not os.path.isfile(missing)
    try:
        ts = Gts.read(missing, fmt="pos", verbose=False)
        assert ts is None
    except Exception:
        pass  # acceptable to raise


# -----------------------------------------------------------------------------
# pyacs.gts.lib.format (Gts format I/O)
# -----------------------------------------------------------------------------


def test_format_write_pos_read_round_trip(ts_quem):
    """write_pos then read_pos recovers same code and data (within tolerance)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ts_quem.write_pos(idir=tmpdir, verbose=False)
        out_path = os.path.join(tmpdir, "QUEM.pos")
        assert os.path.isfile(out_path)
        ts_back = Gts.read(out_path, fmt="pos", verbose=False)
        assert ts_back is not None
        assert ts_back.code == ts_quem.code
        assert ts_back.data.shape == ts_quem.data.shape
        np.testing.assert_allclose(ts_back.data, ts_quem.data, atol=1e-6, rtol=1e-5)


def test_format_force_daily_returns_gts(ts_quem):
    """force_daily(in_place=False) returns new Gts with same number of epochs."""
    ts_daily = ts_quem.force_daily(in_place=False)
    assert ts_daily is not ts_quem
    assert ts_daily.data.shape[0] == ts_quem.data.shape[0]
    assert ts_daily.code == ts_quem.code


def test_format_force_daily_dates_at_noon(ts_quem):
    """force_daily places decimal years at 12:00 (MJD.5)."""
    import pyacs.lib.astrotime as at
    ts_daily = ts_quem.force_daily(in_place=False)
    mjd = at.decyear2mjd(ts_daily.data[:, 0])
    # Should be integer + 0.5 (noon)
    frac = np.abs(mjd - np.round(mjd))
    np.testing.assert_allclose(frac, 0.5, atol=1e-6)


def test_format_force_daily_in_place(ts_quem):
    """force_daily(in_place=True) modifies self and returns self."""
    ts = ts_quem.copy()
    out = ts.force_daily(in_place=True)
    assert out is ts
    import pyacs.lib.astrotime as at
    mjd = at.decyear2mjd(ts.data[:, 0])
    frac = np.abs(mjd - np.round(mjd))
    np.testing.assert_allclose(frac, 0.5, atol=1e-6)


def test_format_to_pandas_df_with_uncertainty(ts_quem):
    """to_pandas_df(uncertainty=True) includes SN, SE, SU columns."""
    pytest.importorskip("pandas")
    df = ts_quem.to_pandas_df(uncertainty=True)
    assert "SN" in df.columns and "SE" in df.columns and "SU" in df.columns
    assert len(df) == ts_quem.data.shape[0]


def test_format_to_pandas_df_data_xyz(ts_quem):
    """to_pandas_df(data_xyz=True) has X, Y, Z columns."""
    pytest.importorskip("pandas")
    df = ts_quem.to_pandas_df(data_xyz=True)
    assert list(df.columns[:3]) == ["X", "Y", "Z"]
    assert len(df) == ts_quem.data.shape[0]


def test_format_to_pytrf(ts_quem):
    """to_pytrf() returns pytrf ts object when pytrf is installed."""
    try:
        from pytrf.ts import ts as pytrf_ts
    except ImportError:
        pytest.skip("pytrf not installed")
    r = ts_quem.to_pytrf()
    assert r is not None
    assert hasattr(r, "y") or hasattr(r, "t")
    assert len(r.t) == ts_quem.data.shape[0]


def test_format_write_cats(ts_quem):
    """write_cats(idir=...) creates cats_CODE.dat file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ts_quem.write_cats(idir=tmpdir, add_key="")
        path = os.path.join(tmpdir, "cats_QUEM.dat")
        assert os.path.isfile(path)
        data = np.loadtxt(path)
        assert data.shape[0] == ts_quem.data.shape[0]
        assert data.shape[1] == 7


def test_format_write_cats_add_key(ts_quem):
    """write_cats(add_key='test') creates cats_CODE_test.dat."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ts_quem.write_cats(idir=tmpdir, add_key="test")
        path = os.path.join(tmpdir, "cats_QUEM_test.dat")
        assert os.path.isfile(path)


def test_format_read_pos_on_empty_gts():
    """read_pos(tsfile=path) on new Gts() loads from file."""
    if not os.path.isfile(QUEM_POS):
        pytest.skip("QUEM.pos not found")
    ts = Gts()
    ts.read_pos(tsfile=QUEM_POS, verbose=False)
    assert ts.code == "QUEM"
    assert ts.data is not None
    assert ts.data.shape[1] == 10


def test_format_read_pos_xyz_true(ts_quem):
    """read_pos(tsfile=..., xyz=True) populates data_xyz."""
    assert ts_quem.data_xyz is not None
    assert ts_quem.data_xyz.shape[0] == ts_quem.data.shape[0]
    assert ts_quem.data_xyz.shape[1] >= 10


def test_format_read_cats_file_missing_raises():
    """read_cats_file with missing file raises or fails (missing file)."""
    gts = Gts(code="QUEM")
    path_missing = os.path.join(_TEST_DIR, "data", "ts", "NONEXISTENT_cats.dat")
    assert not os.path.isfile(path_missing)
    with pytest.raises(Exception):
        gts.read_cats_file(ifile=path_missing, gmt=False, verbose=False)


def test_format_read_cats_file_existing_from_pos_data(ts_quem):
    """write_cats then read_cats_file recovers same number of rows (cats format)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ts_quem.write_cats(idir=tmpdir)
        path_cats = os.path.join(tmpdir, "cats_QUEM.dat")
        gts = Gts(code="QUEM")
        gts.read_cats_file(ifile=path_cats, gmt=False, verbose=False)
        assert gts.data is not None
        assert gts.data.shape[0] == ts_quem.data.shape[0]
        assert gts.data.shape[1] >= 7

