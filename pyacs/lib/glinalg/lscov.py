"""Solve least-squares with data covariance matrix."""

import numpy as np

from .ls import ls


def lscov(G, d, cov, method='chol'):
    """Solve least-squares with data covariance matrix.

    Parameters
    ----------
    G : numpy.ndarray
        m x n design matrix.
    d : numpy.ndarray
        m observation vector.
    cov : numpy.ndarray
        m x m covariance matrix for d.
    method : str, optional
        'chol' for Cholesky. Default is 'chol'.

    Returns
    -------
    numpy.ndarray
        Solution vector.
    """
    if method == 'chol':
        L = np.linalg.cholesky(cov)
        SQRT_Wd = np.linalg.inv(L)
        GG = np.dot(SQRT_Wd, G)
        dd = np.dot(SQRT_Wd, d)
        return ls(GG, dd)
