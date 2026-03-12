"""Solve weighted least-squares and return solution and posterior covariance."""

import numpy as np
import scipy.linalg.lapack


def lsw_full(G, d, std, verbose=False):
    """Solve weighted least-squares and return solution and posterior covariance.

    Parameters
    ----------
    G : numpy.ndarray
        m x n design matrix.
    d : numpy.ndarray
        m observation vector.
    std : array_like
        Standard deviations for d (length m).
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
    """
    GG = (G.T / std).T
    dd = d / std

    N = np.dot(GG.T, GG)

    L, info = scipy.linalg.lapack.dpotrf(N, False, False)

    if verbose:
        print('-- info from scipy.linalg.lapack.dpotrf in lsw_full')
        print('-- info: %d', info)

    L_inv, info = scipy.linalg.lapack.dpotri(L)

    if verbose:
        print('-- info from scipy.linalg.lapack.dpotri in lsw_full')
        print('-- info: %d', info)

    COV = np.triu(L_inv) + np.triu(L_inv, k=1).T

    X = np.dot(COV, np.dot(GG.T, dd))

    V = d - np.dot(G, X)

    return X, COV, V
