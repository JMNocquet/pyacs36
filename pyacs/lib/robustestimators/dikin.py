"""L1-norm estimation using Dikin's method."""

import numpy

from .dik_m import Dik_m


def Dikin(A, y, W, eps=3.E-3):
    """L1-norm estimation using Dikin's method (linear program y = Ax + e).

    Parameters
    ----------
    A : numpy.ndarray
        Design matrix.
    y : numpy.ndarray
        Observation vector.
    W : numpy.ndarray, optional
        Diagonal weight matrix for observables. Can be None.
    eps : float, optional
        Stopping criterion for iteration. Default is 3e-3.

    Returns
    -------
    x : numpy.ndarray
        Estimated parameters.
    e : numpy.ndarray
        Residuals.

    Notes
    -----
    Reference: A. Khodabandeh and A. R. Amiri-Simkooei, "Recursive Algorithm for
    L1 Norm Estimation in Linear Models", J. Surveying Engineering ASCE (2011).
    doi:10.1061/(ASCE)SU.1943-5428.0000031. Translated from Matlab by J.-M. Nocquet (2011).
    """
    if y.ndim == 1:
        y = y.reshape(-1, 1)

    if W is not None:
        w = numpy.sqrt(W)
        A = (A.T * w).T
        y = (y.T * w).T

    (m, n) = numpy.shape(A)
    At = numpy.hstack((A, -A, numpy.eye(m), -numpy.eye(m)))
    f = numpy.vstack((numpy.zeros((2 * n, 1)), numpy.ones((2 * m, 1))))

    dlin = Dik_m(f, At, y, eps)
    x = dlin[0:n] - dlin[n:2 * n]
    e = dlin[2 * n:2 * n + m] - dlin[2 * n + m:2 * (n + m)]

    if W is not None:
        e = (e.T / w).T

    return (x.flatten(), e.flatten())
