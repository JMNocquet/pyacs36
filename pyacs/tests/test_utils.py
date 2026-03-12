"""
Unit tests for pyacs.lib.utils (string parsing, grids, recarrays, run_cmd, etc.).
"""

import os
import tempfile

import numpy as np
import pytest

import pyacs.lib.utils as U


# -----------------------------------------------------------------------------
# str2list_float
# -----------------------------------------------------------------------------


def test_str2list_float_simple():
    """str2list_float parses '[0, 2.5, 1E9]' to list of floats."""
    out = U.str2list_float("[0, 2.5, 1E9]")
    assert out == [0.0, 2.5, 1e9]


def test_str2list_float_single():
    """str2list_float parses single value '[1.5]'."""
    out = U.str2list_float("[1.5]")
    assert out == [1.5]


def test_str2list_float_no_brackets():
    """str2list_float handles string without brackets (replace removes them)."""
    out = U.str2list_float("1, 2, 3")
    assert out == [1.0, 2.0, 3.0]


def test_str2list_float_negative():
    """str2list_float parses negative and scientific notation."""
    out = U.str2list_float("[-1.5, 2e-3]")
    assert np.isclose(out[0], -1.5)
    assert np.isclose(out[1], 0.002)


# -----------------------------------------------------------------------------
# __ensure_list_of_list
# -----------------------------------------------------------------------------


def test_ensure_list_of_list_single_list():
    """Single list [a, b] becomes [[a, b]]."""
    out = U.__ensure_list_of_list([1, 2])
    assert out == [[1, 2]]


def test_ensure_list_of_list_already_list_of_list():
    """List of lists is returned unchanged."""
    inp = [[1, 2], [3, 4]]
    out = U.__ensure_list_of_list(inp)
    assert out == inp


def test_ensure_list_of_list_not_list_raises():
    """Non-list raises TypeError."""
    with pytest.raises(TypeError):
        U.__ensure_list_of_list((1, 2))


def test_ensure_list_of_list_single_nested():
    """[[1, 2]] stays [[1, 2]]."""
    out = U.__ensure_list_of_list([[1, 2]])
    assert out == [[1, 2]]


# -----------------------------------------------------------------------------
# make_grid
# -----------------------------------------------------------------------------


def test_make_grid_nx_ny():
    """make_grid with nx, ny returns (n_points, 2) array."""
    pytest.importorskip("colors")
    grid = U.make_grid(0.0, 1.0, 0.0, 1.0, nx=3, ny=2)
    assert grid.shape == (3 * 2, 2)
    assert grid[:, 0].min() == 0.0
    assert grid[:, 0].max() == 1.0
    assert grid[:, 1].min() == 0.0
    assert grid[:, 1].max() == 1.0


def test_make_grid_ny_defaults_to_nx():
    """make_grid with nx only uses ny=nx."""
    pytest.importorskip("colors")
    grid = U.make_grid(0.0, 1.0, 0.0, 1.0, nx=2)
    assert grid.shape == (4, 2)


def test_make_grid_step_x_step_y():
    """make_grid with step_x, step_y returns grid with that spacing."""
    pytest.importorskip("colors")
    grid = U.make_grid(0.0, 1.0, 0.0, 1.0, step_x=0.5, step_y=0.5)
    assert grid.shape[1] == 2
    assert 0.0 in grid[:, 0]
    assert 0.5 in grid[:, 0]
    assert 0.0 in grid[:, 1]
    assert 0.5 in grid[:, 1]


def test_make_grid_neither_nx_nor_step_raises():
    """make_grid without nx or step_x raises TypeError."""
    pytest.importorskip("colors")
    with pytest.raises(TypeError):
        U.make_grid(0.0, 1.0, 0.0, 1.0)


def test_make_grid_both_nx_and_step_raises():
    """make_grid with both nx and step_x raises TypeError."""
    pytest.importorskip("colors")
    with pytest.raises(TypeError):
        U.make_grid(0.0, 1.0, 0.0, 1.0, nx=2, step_x=0.5)


def test_make_grid_outfile():
    """make_grid with outfile writes file when format is default."""
    pytest.importorskip("colors")
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = os.path.join(tmpdir, "grid.txt")
        grid = U.make_grid(0.0, 1.0, 0.0, 1.0, nx=2, outfile=outpath)
        assert os.path.isfile(outpath)
        data = np.loadtxt(outpath)
        assert data.shape[0] == 4
        assert data.shape[1] == 8  # psvelo format


# -----------------------------------------------------------------------------
# run_cmd
# -----------------------------------------------------------------------------


def test_run_cmd_success():
    """run_cmd runs a simple command and returns CompletedProcess."""
    result = U.run_cmd("echo 1", shell=True)
    assert result.returncode == 0
    assert "1" in result.stdout


def test_run_cmd_capture_stdout():
    """run_cmd captures stdout when capture_output=True."""
    result = U.run_cmd("echo hello", shell=True, capture_output=True)
    assert "hello" in result.stdout


def test_run_cmd_fail_raises():
    """run_cmd with check=True raises on non-zero exit."""
    with pytest.raises(Exception):  # CalledProcessError
        U.run_cmd("exit 1", shell=True, check=True)


def test_run_cmd_check_false():
    """run_cmd with check=False does not raise on non-zero exit."""
    result = U.run_cmd("exit 1", shell=True, check=False)
    assert result.returncode == 1


# -----------------------------------------------------------------------------
# numpy_array_2_numpy_recarray / numpy_recarray_2_numpy_array
# -----------------------------------------------------------------------------


def test_numpy_array_2_numpy_recarray_shape_and_names():
    """numpy_array_2_numpy_recarray produces recarray with given field names."""
    A = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
    names = ["x", "y"]
    rec = U.numpy_array_2_numpy_recarray(A, names)
    assert rec.dtype.names == ("x", "y")
    assert len(rec) == 2


def test_numpy_recarray_2_numpy_array_round_trip():
    """array -> recarray -> array recovers shape and values."""
    A = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
    names = ["a", "b"]
    rec = U.numpy_array_2_numpy_recarray(A, names)
    B = U.numpy_recarray_2_numpy_array(rec)
    assert B.shape == A.shape
    # Order may be row-major; compare values
    np.testing.assert_allclose(B, A, atol=1e-10)


def test_numpy_recarray_2_numpy_array_structured():
    """numpy_recarray_2_numpy_array converts structured array to plain array (homogeneous dtype)."""
    # Function uses view(dtype[0]), so use same dtype for all fields
    arr = np.array([(1.0, 2.0), (3.0, 4.0)], dtype=[("f1", float), ("f2", float)])
    out = U.numpy_recarray_2_numpy_array(arr)
    assert out.shape == (2, 2)
    assert np.isclose(out[0, 0], 1.0)
    assert np.isclose(out[0, 1], 2.0)
    assert np.isclose(out[1, 0], 3.0)
    assert np.isclose(out[1, 1], 4.0)


# -----------------------------------------------------------------------------
# save_np_array_with_string
# -----------------------------------------------------------------------------


def test_save_np_array_with_string_writes_file():
    """save_np_array_with_string writes rows with string column."""
    A = np.array([[1, 2], [3, 4]])
    S = np.array(["a", "b"])
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = os.path.join(tmpdir, "out.dat")
        U.save_np_array_with_string(A, S, "%d %d %s", outpath)
        assert os.path.isfile(outpath)
        with open(outpath) as f:
            lines = f.readlines()
        assert lines[0].startswith("#")
        assert "1 2 a" in lines[1] or "1  2 a" in lines[1]
        assert "3 4 b" in lines[2] or "3  4 b" in lines[2]


def test_save_np_array_with_string_length_mismatch_raises():
    """save_np_array_with_string raises when len(S) != A.shape[0]."""
    A = np.array([[1, 2], [3, 4]])
    S = np.array(["a", "b", "c"])
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = os.path.join(tmpdir, "out.dat")
        with pytest.raises(TypeError):
            U.save_np_array_with_string(A, S, "%d %d %s", outpath)


def test_save_np_array_with_string_comment():
    """save_np_array_with_string writes comment line."""
    A = np.array([[1]])
    S = np.array(["x"])
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = os.path.join(tmpdir, "out.dat")
        U.save_np_array_with_string(A, S, "%d %s", outpath, comment="header")
        with open(outpath) as f:
            first = f.readline()
        assert first.strip().startswith("#")
