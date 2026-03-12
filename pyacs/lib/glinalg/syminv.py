"""Invert a symmetric positive-definite matrix."""

from scipy import linalg


def syminv(M):
    """Invert a symmetric positive-definite matrix.

    Parameters
    ----------
    M : numpy.ndarray
        Symmetric matrix (positive definite).

    Returns
    -------
    numpy.ndarray
        Inverse matrix.
    """
    Q = linalg.lapack.dpotri(linalg.lapack.dpotrf(M)[0])[0]
    Q = Q + Q.T
    n = len(Q)
    Q[range(n), range(n)] = Q[range(n), range(n)] / 2
    return Q
