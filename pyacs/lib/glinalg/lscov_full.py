"""Solve least-squares with data covariance and return posterior covariance."""

import numpy as np
import scipy.linalg.lapack


def lscov_full(G, d, cov, verbose=False):
    """Solve least-squares with data covariance and return posterior covariance.

    Parameters
    ----------
    G : numpy.ndarray
        m x n design matrix.
    d : numpy.ndarray
        m observation vector.
    cov : numpy.ndarray
        m x m covariance matrix for d.
    verbose : bool, optional
        If True, print lapack info. Default is False.

    Returns
    -------
    X : numpy.ndarray
        Solution vector.
    COV : numpy.ndarray
        Posterior covariance of the solution.
    V : numpy.ndarray
        Residuals d - G*X.
    chi2 : float
        Weighted residual sum of squares.
    """
    L = np.linalg.cholesky(cov)
    SQRT_Wd = np.linalg.inv(L)

    GG = np.dot(SQRT_Wd, G)
    dd = np.dot(SQRT_Wd, d)

    N = np.dot(GG.T, GG)

    L, info = scipy.linalg.lapack.dpotrf(N, False, False)

    if info < 0:
        print('-- the i-th argument had an illegal value. i=%d' % info)
    if info > 0:
        print('-- the the leading minor of order i is not positive definite' +
              ' and the factorization could not be completed. i=%d' % info)
    if info == 0 and verbose:
        print('-- dpotrf ok')

    L_inv, info = scipy.linalg.lapack.dpotri(L)

    if info < 0:
        print('-- the i-th argument had an illegal value. i=%d' % info)
    if info > 0:
        print('-- if INFO = i, the (i,i) element of the factor U or L  is zero, and the inverse could not be computed. i=%d' % info)
    if info == 0 and verbose:
        print('-- dpotri ok')

    COV = np.triu(L_inv) + np.triu(L_inv, k=1).T

    X = np.dot(COV, np.dot(GG.T, dd))

    V = d - np.dot(G, X)

    vv = np.dot(SQRT_Wd, V)
    chi2 = np.dot(vv.T, vv)

    return X, COV, V, chi2
