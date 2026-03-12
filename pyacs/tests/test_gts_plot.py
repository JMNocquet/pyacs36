"""
Unit tests for Gts.plot().
"""

import numpy as np
import pytest

from pyacs.gts.Gts import Gts


def _minimal_gts(n=20):
    """Return a minimal Gts with synthetic N/E/U data for plotting."""
    gts = Gts(code="TEST")
    # data: decimal year, N, E, U, sN, sE, sU (meters)
    t = np.linspace(2020.0, 2022.0, n)
    data = np.column_stack([
        t,
        np.cumsum(np.random.randn(n) * 1e-3),
        np.cumsum(np.random.randn(n) * 1e-3),
        np.cumsum(np.random.randn(n) * 1e-3),
        np.full(n, 1e-4),
        np.full(n, 1e-4),
        np.full(n, 1e-4),
    ])
    gts.data = data
    return gts


# -------------------------------------------------------------------------
# plot
# -------------------------------------------------------------------------


def test_gts_plot_callable():
    """Gts has a plot method and it is callable."""
    assert hasattr(Gts, "plot")
    assert callable(Gts.plot)


def test_gts_plot_runs_with_show_false():
    """Gts().plot(show=False) runs without opening a window (no GUI)."""
    pytest.importorskip("matplotlib")
    gts = _minimal_gts()
    # show=False avoids blocking / display
    result = gts.plot(show=False)
    assert isinstance(result, Gts)
