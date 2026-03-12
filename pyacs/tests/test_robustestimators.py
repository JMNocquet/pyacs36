"""
Unit tests for pyacs.lib.robustestimators (Dikin L1, Dik_m LP solver).
"""

import numpy as np
import pytest

import pyacs.lib.robustestimators as RE


# -----------------------------------------------------------------------------
# Exceptions
# -----------------------------------------------------------------------------


def test_error_base_class():
    """Error is an Exception subclass."""
    assert issubclass(RE.Error, Exception)


def test_unbounded_function_error():
    """UnboundedFunctionError is an Exception subclass and can be raised."""
    assert issubclass(RE.UnboundedFunctionError, Exception)
    with pytest.raises(RE.UnboundedFunctionError):
        raise RE.UnboundedFunctionError("test message")


# -----------------------------------------------------------------------------
# Dik_m (linear program solver)
# -----------------------------------------------------------------------------


def test_dik_m_simple_feasible():
    """Dik_m solves min c'x s.t. Ax = b and returns solution of correct shape."""
    # min x1 + x2  s.t.  x1 + x2 = 1  (solution x1=x2=0.5 for min at 0.5, but LP is bounded)
    A = np.array([[1.0, 1.0]])
    b = np.array([[1.0]])
    c = np.array([[1.0], [1.0]])
    X = RE.Dik_m(c, A, b, er=1e-5)
    assert X.shape == (2, 1)
    # Constraint: A @ X ≈ b
    np.testing.assert_allclose((A @ X).flatten(), b.flatten(), atol=1e-4)


def test_dik_m_single_variable():
    """Dik_m with one variable: A (1,1), b scalar."""
    A = np.array([[1.0]])
    b = np.array([[2.0]])
    c = np.array([[1.0]])
    X = RE.Dik_m(c, A, b, er=1e-5)
    assert X.shape == (1, 1)
    np.testing.assert_allclose((A @ X).flatten(), b.flatten(), atol=1e-4)


# -----------------------------------------------------------------------------
# Dikin (L1 regression)
# -----------------------------------------------------------------------------


def test_dikin_trivial_system():
    """Dikin with A=[[1]], y=[1]: solution x≈1, residual e≈0."""
    A = np.array([[1.0]])
    y = np.array([1.0])
    x, e = RE.Dikin(A, y, W=None, eps=1e-5)
    assert x.shape == (1,)
    assert e.shape == (1,)
    np.testing.assert_allclose(x, [1.0], atol=1e-3)
    np.testing.assert_allclose(e, [0.0], atol=1e-3)


def test_dikin_residual_consistency():
    """Dikin: A @ x + e ≈ y (residual consistency)."""
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    y = np.array([1.0, 2.0, 3.0])
    x, e = RE.Dikin(A, y, W=None, eps=1e-5)
    assert x.shape == (2,)
    assert e.shape == (3,)
    np.testing.assert_allclose(A @ x + e, y, atol=1e-3)


def test_dikin_with_weights():
    """Dikin with diagonal weight W (1D) returns consistent A@x + e ≈ y."""
    A = np.array([[1.0, 0.0], [1.0, 1.0]])
    y = np.array([1.0, 2.0])
    W = np.array([1.0, 1.0])  # diagonal weights (one per row)
    x, e = RE.Dikin(A, y, W=W, eps=1e-5)
    assert x.shape == (2,)
    assert e.shape == (2,)
    np.testing.assert_allclose(A @ x + e, y, atol=1e-3)


def test_dikin_y_1d_accepted():
    """Dikin accepts 1D y (reshape to column internally)."""
    A = np.array([[1.0], [2.0]])
    y = np.array([1.0, 2.0])
    x, e = RE.Dikin(A, y, W=None, eps=1e-5)
    assert x.shape == (1,)
    assert e.shape == (2,)
    np.testing.assert_allclose(A.flatten() * x[0] + e, y, atol=1e-3)


def test_dikin_overdetermined():
    """Dikin on overdetermined system: solution minimizes L1 residual."""
    np.random.seed(42)
    A = np.random.randn(10, 2)
    x_true = np.array([1.0, -0.5])
    y = A @ x_true + 0.01 * np.random.randn(10)
    x, e = RE.Dikin(A, y, W=None, eps=1e-5)
    assert x.shape == (2,)
    assert e.shape == (10,)
    np.testing.assert_allclose(A @ x + e, y, atol=1e-3)
    np.testing.assert_allclose(x, x_true, atol=0.5)
