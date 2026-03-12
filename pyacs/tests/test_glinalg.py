"""
Unit tests for pyacs.lib.glinalg (linear algebra for geodesy).
"""

import numpy as np

import pyacs.lib.glinalg as GL


# -----------------------------------------------------------------------------
# corr_to_cov / cov_to_corr
# -----------------------------------------------------------------------------


def test_cov_to_corr_diagonal():
    """Covariance with zero off-diagonal gives correlation I and sqrt(diag)."""
    cov = np.diag([4.0, 9.0, 1.0])
    corr, sigma = GL.cov_to_corr(cov)
    assert np.allclose(corr, np.eye(3))
    assert np.allclose(sigma, [2.0, 3.0, 1.0])


def test_corr_to_cov_round_trip():
    """corr_to_cov then cov_to_corr recovers correlation and sigmas."""
    corr = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, -0.3], [0.2, -0.3, 1.0]])
    sigma = np.array([1.5, 2.0, 0.5])
    cov = GL.corr_to_cov(corr, sigma)
    corr_back, sigma_back = GL.cov_to_corr(cov)
    assert np.allclose(corr_back, corr, atol=1e-10)
    assert np.allclose(sigma_back, sigma, atol=1e-10)


def test_cov_to_corr_zero_division():
    """Zero variance on diagonal: corresponding correlation row/col zero."""
    cov = np.array([[1.0, 0.0], [0.0, 0.0]])
    corr, sigma = GL.cov_to_corr(cov)
    assert np.allclose(sigma, [1.0, 0.0])
    assert np.allclose(corr[1, :], 0.0)
    assert np.allclose(corr[:, 1], 0.0)


# -----------------------------------------------------------------------------
# syminv
# -----------------------------------------------------------------------------


def test_syminv_inverse():
    """syminv(M) is the inverse of M for symmetric positive definite M."""
    M = np.array([[4.0, 1.0, 0.0], [1.0, 3.0, 0.5], [0.0, 0.5, 2.0]])
    Minv = GL.syminv(M)
    assert np.allclose(M @ Minv, np.eye(3), atol=1e-10)
    assert np.allclose(Minv @ M, np.eye(3), atol=1e-10)


def test_syminv_symmetric_result():
    """syminv returns a symmetric matrix."""
    M = np.array([[5.0, 2.0], [2.0, 3.0]])
    Minv = GL.syminv(M)
    assert np.allclose(Minv, Minv.T)


# -----------------------------------------------------------------------------
# cov_to_invcov
# -----------------------------------------------------------------------------


def test_cov_to_invcov_inverse():
    """cov_to_invcov(C) gives C^{-1} for SPD C."""
    C = np.array([[2.0, 0.5], [0.5, 1.0]])
    Cinv = GL.cov_to_invcov(C)
    assert np.allclose(C @ Cinv, np.eye(2), atol=1e-10)


# -----------------------------------------------------------------------------
# ls
# -----------------------------------------------------------------------------


def test_ls_simple_system():
    """ls solves G x = d (overdetermined, unique min-norm)."""
    G = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    d = np.array([1.0, 2.0, 3.0])
    x = GL.ls(G, d)
    assert x.shape == (2,)
    # Normal equations: G'T G x = G'T d  => [2 1; 1 2] x = [4; 5] => x = [1, 2]
    assert np.allclose(x, [1.0, 2.0], atol=1e-10)


def test_ls_square_system():
    """ls on square full-rank system gives exact solution."""
    G = np.array([[1.0, 2.0], [3.0, 4.0]])
    d = np.array([5.0, 11.0])
    x = GL.ls(G, d)
    assert np.allclose(G @ x, d, atol=1e-10)


# -----------------------------------------------------------------------------
# lscov_full
# -----------------------------------------------------------------------------


def test_lscov_full_solution_residuals():
    """lscov_full returns X such that residuals d - G@X have correct weighted chi2."""
    np.random.seed(42)
    G = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
    d = np.array([1.0, 2.0, 3.0], dtype=float)
    cov = np.eye(3)
    X, COV, V, chi2 = GL.lscov_full(G, d, cov)
    assert X.shape == (2,)
    assert np.allclose(G @ X, d, atol=1e-10)
    assert np.allclose(V, d - G @ X, atol=1e-10)
    assert np.allclose(chi2, np.dot(V, V), atol=1e-10)


def test_lscov_full_posterior_covariance():
    """lscov_full posterior covariance is (G' W G)^{-1} for W = cov^{-1}."""
    G = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float)
    d = np.array([1.0, 2.0], dtype=float)
    cov = np.eye(2)
    X, COV, V, chi2 = GL.lscov_full(G, d, cov)
    # For G=I, d=[1,2], solution X=[1,2], posterior cov = I
    assert np.allclose(COV, np.eye(2), atol=1e-10)


# -----------------------------------------------------------------------------
# lsw / lsw_full
# -----------------------------------------------------------------------------


def test_lsw_weighted_solution():
    """lsw with std weights gives same solution as lscov with diag(sigma^2)."""
    G = np.array([[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float)
    d = np.array([1.0, 2.0, 1.5], dtype=float)
    std = np.array([1.0, 2.0, 1.0], dtype=float)
    cov = np.diag(std ** 2)
    x_lsw = GL.lsw(G, d, std)
    x_lscov = GL.lscov_full(G, d, cov)[0]
    assert np.allclose(x_lsw, x_lscov, atol=1e-10)


def test_lsw_full_returns_covariance_and_residuals():
    """lsw_full returns X, COV, V with V = d - G@X."""
    G = np.array([[1.0, 1.0], [1.0, 0.0]], dtype=float)
    d = np.array([3.0, 1.0], dtype=float)
    std = np.array([1.0, 1.0], dtype=float)
    X, COV, V = GL.lsw_full(G, d, std)
    assert np.allclose(V, d - G @ X, atol=1e-10)
    assert COV.shape == (2, 2)
    assert np.allclose(COV, COV.T)


# -----------------------------------------------------------------------------
# dot
# -----------------------------------------------------------------------------


def test_dot_matrix_vector():
    """dot(A, b) equals A @ b."""
    A = np.array([[1.0, 2.0], [3.0, 4.0]], order='C')
    b = np.array([1.0, -1.0])
    assert np.allclose(GL.dot(A, b), A @ b, atol=1e-10)


def test_dot_vector_vector():
    """dot(a, b) equals inner product for 1D arrays."""
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([-1.0, 0.5, 2.0])
    assert np.isclose(GL.dot(a, b), np.dot(a, b), atol=1e-10)


def test_dot_matrix_matrix():
    """dot(A, B) equals A @ B for 2D arrays."""
    A = np.array([[1.0, 2.0], [3.0, 4.0]], order='C')
    B = np.array([[0.0, 1.0], [1.0, 0.0]], order='C')
    assert np.allclose(GL.dot(A, B), A @ B, atol=1e-10)


# -----------------------------------------------------------------------------
# symmetrize
# -----------------------------------------------------------------------------


def test_symmetrize_triu():
    """symmetrize with triu builds symmetric matrix from upper triangle."""
    a = np.array([[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [0.0, 0.0, 6.0]])
    S = GL.symmetrize(a, 'triu')
    assert np.allclose(S, S.T)
    assert np.allclose(np.triu(S), np.triu(np.array([[1,2,3],[2,4,5],[3,5,6]])))


def test_symmetrize_tril():
    """symmetrize with tril builds symmetric matrix from lower triangle."""
    a = np.array([[1.0, 0.0, 0.0], [2.0, 3.0, 0.0], [4.0, 5.0, 6.0]])
    S = GL.symmetrize(a, 'tril')
    assert np.allclose(S, S.T)
    assert np.allclose(S, np.array([[1,2,4],[2,3,5],[4,5,6]]))


# -----------------------------------------------------------------------------
# sympinv
# -----------------------------------------------------------------------------


def test_sympinv_full_rank():
    """sympinv of full-rank symmetric matrix matches syminv."""
    M = np.array([[4.0, 1.0], [1.0, 3.0]])
    Mp = GL.sympinv(M)
    Minv = GL.syminv(M)
    assert np.allclose(Mp, Minv, atol=1e-10)


def test_sympinv_rank_deficient():
    """sympinv of rank-deficient matrix gives pseudo-inverse (M @ Mp @ M ≈ M)."""
    # Rank-1 matrix
    v = np.array([1.0, 2.0, 3.0])
    M = np.outer(v, v)
    Mp = GL.sympinv(M)
    assert np.allclose(M @ Mp @ M, M, atol=1e-8)


# -----------------------------------------------------------------------------
# dot_and_sum
# -----------------------------------------------------------------------------


def test_dot_and_sum():
    """dot_and_sum(LX, a) equals sum of a[i] * LX[i]."""
    LX = [np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([1.0, 1.0])]
    a = [1.0, 2.0, -1.0]
    out = GL.dot_and_sum(LX, a)
    expected = 1.0 * LX[0] + 2.0 * LX[1] - 1.0 * LX[2]
    assert np.allclose(out, expected, atol=1e-10)


# -----------------------------------------------------------------------------
# odot
# -----------------------------------------------------------------------------


def test_odot_scale_rows():
    """odot(a, G) scales block-rows of G by a."""
    # One row per block: a has 2 elements, G has 2 rows
    a = np.array([2.0, 3.0])
    G = np.array([[1.0, 1.0], [1.0, 1.0]])
    R = GL.odot(a, G)
    assert np.allclose(R, [[2.0, 2.0], [3.0, 3.0]], atol=1e-10)


def test_odot_block_scaling():
    """odot with 2 blocks of 2 rows scales each 2-row block."""
    a = np.array([1.0, -1.0])
    G = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
    R = GL.odot(a, G)
    expected = np.array([[1, 0], [1, 0], [0, -1], [0, -1]], dtype=float)
    assert np.allclose(R, expected, atol=1e-10)
