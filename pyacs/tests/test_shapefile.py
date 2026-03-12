"""
Unit tests for pyacs.lib.shapefile (converters to shapefile).
"""

import os
import tempfile
import warnings

import pytest

import pyacs.lib.shapefile as SHF


# -----------------------------------------------------------------------------
# Module exports
# -----------------------------------------------------------------------------


def test_shapefile_module_exports():
    """Module exposes psvelo_to_shapefile, pyeblock_fault."""
    assert hasattr(SHF, "psvelo_to_shapefile")
    assert hasattr(SHF, "pyeblock_fault")
    assert callable(SHF.psvelo_to_shapefile)
    assert callable(SHF.pyeblock_fault)


# -----------------------------------------------------------------------------
# psvelo_to_shapefile
# -----------------------------------------------------------------------------


def test_psvelo_to_shapefile_creates_shp_and_prj():
    """psvelo_to_shapefile creates .shp and .prj from a valid psvelo file."""
    pytest.importorskip("shapefile")
    # Minimal psvelo: lon lat Ve Vn S_Ve S_Vn S_Ven name (8 columns)
    content = "0.0 45.0 1.0 0.0 0.1 0.1 0.0 SITE1\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        psvelo_path = os.path.join(tmpdir, "vel.psvelo")
        shp_base = os.path.join(tmpdir, "out")
        with open(psvelo_path, "w") as f:
            f.write(content)
        SHF.psvelo_to_shapefile(psvelo_path, shp_base, verbose=False)
        assert os.path.isfile(shp_base + ".shp")
        assert os.path.isfile(shp_base + ".prj")
        with open(shp_base + ".prj") as f:
            prj = f.read()
        assert "WGS" in prj or "GEOGCS" in prj


def test_psvelo_to_shapefile_two_points():
    """psvelo_to_shapefile with two lines writes two point features."""
    pytest.importorskip("shapefile")
    content = (
        "0.0 45.0 1.0 0.0 0.1 0.1 0.0 A\n"
        "1.0 46.0 0.0 1.0 0.1 0.1 0.0 B\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        psvelo_path = os.path.join(tmpdir, "vel.psvelo")
        shp_base = os.path.join(tmpdir, "out")
        with open(psvelo_path, "w") as f:
            f.write(content)
        SHF.psvelo_to_shapefile(psvelo_path, shp_base, verbose=False)
        import shapefile
        r = shapefile.Reader(shp_base)
        assert len(r.shapes()) == 2
        assert len(r.records()) == 2
        r.close()


def test_psvelo_to_shapefile_empty_file_returns_early():
    """psvelo_to_shapefile on empty file returns without writing (no exception)."""
    pytest.importorskip("shapefile")
    with tempfile.TemporaryDirectory() as tmpdir:
        psvelo_path = os.path.join(tmpdir, "empty.psvelo")
        shp_base = os.path.join(tmpdir, "out")
        with open(psvelo_path, "w") as f:
            f.write("")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # genfromtxt empty file
            SHF.psvelo_to_shapefile(psvelo_path, shp_base, verbose=False)
        # Empty file: np_vel.size == 0, function returns; .shp may not exist
        assert not os.path.isfile(shp_base + ".shp")


def test_psvelo_to_shapefile_missing_file_raises():
    """psvelo_to_shapefile on missing file raises (ERROR with exit=True or IOError)."""
    pytest.importorskip("shapefile")
    with tempfile.TemporaryDirectory() as tmpdir:
        missing = os.path.join(tmpdir, "nonexistent.psvelo")
        shp_base = os.path.join(tmpdir, "out")
        with pytest.raises((OSError, SystemExit, Exception)):
            SHF.psvelo_to_shapefile(missing, shp_base, verbose=False)


# -----------------------------------------------------------------------------
# pyeblock_fault
# -----------------------------------------------------------------------------


def test_pyeblock_fault_creates_shp_and_prj():
    """pyeblock_fault creates .shp and .prj from a valid fault file."""
    pytest.importorskip("shapefile")
    # 8 numeric cols + 2 string cols (left_block, right_block)
    content = "0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 left right\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        fault_path = os.path.join(tmpdir, "faults.dat")
        shp_base = os.path.join(tmpdir, "faults")
        with open(fault_path, "w") as f:
            f.write(content)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PendingDeprecationWarning)  # np.matrix
            SHF.pyeblock_fault(fault_path, shp_base, verbose=False)
        assert os.path.isfile(shp_base + ".shp")
        assert os.path.isfile(shp_base + ".prj")
        assert os.path.isfile(shp_base + ".gmt")
        with open(shp_base + ".prj") as f:
            prj = f.read()
        assert "WGS" in prj or "GEOGCS" in prj


def test_pyeblock_fault_polyline_count():
    """pyeblock_fault writes one polyline per fault line."""
    pytest.importorskip("shapefile")
    content = (
        "0 0.0 0.0 1.0 1.0 10.0 20.0 90.0 block1 block2\n"
        "1 2.0 2.0 3.0 3.0 10.0 20.0 90.0 block2 block3\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        fault_path = os.path.join(tmpdir, "faults.dat")
        shp_base = os.path.join(tmpdir, "faults")
        with open(fault_path, "w") as f:
            f.write(content)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PendingDeprecationWarning)  # np.matrix
            SHF.pyeblock_fault(fault_path, shp_base, verbose=False)
        import shapefile
        r = shapefile.Reader(shp_base)
        assert len(r.shapes()) == 2
        assert r.shapeType == shapefile.POLYLINE
        r.close()


def test_pyeblock_fault_missing_file_raises():
    """pyeblock_fault on missing file raises or propagates exception."""
    pytest.importorskip("shapefile")
    with tempfile.TemporaryDirectory() as tmpdir:
        missing = os.path.join(tmpdir, "nonexistent_fault.dat")
        shp_base = os.path.join(tmpdir, "out")
        with pytest.raises(Exception):
            SHF.pyeblock_fault(missing, shp_base, verbose=False)


