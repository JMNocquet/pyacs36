"""Form the normal system for least-squares."""

import numpy as np


def make_normal_system(A, d, inv_Cd):
    """Form the normal system for A x = d with data covariance Cd.

    Normal system: (A.T inv(Cd) A) x = A.T inv(Cd) d.

    Parameters
    ----------
    A : numpy.ndarray
        Design matrix.
    d : array_like
        Observation vector.
    inv_Cd : numpy.ndarray
        Inverse of data covariance matrix.

    Returns
    -------
    N : numpy.ndarray
        A.T inv(Cd) A.
    Nd : numpy.ndarray
        A.T inv(Cd) d (column vector).
    """
    TMP = np.dot(A.T, inv_Cd)
    return np.dot(TMP, A), np.dot(TMP, d.reshape(-1, 1))
