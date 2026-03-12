"""Return inverse of a covariance matrix."""

from scipy import linalg


def cov_to_invcov(M):
    """Return the inverse of a covariance matrix.

    Parameters
    ----------
    M : numpy.ndarray
        2D covariance matrix (symmetric positive definite).

    Returns
    -------
    numpy.ndarray
        Inverse covariance matrix.
    """
    Q = linalg.lapack.dpotri(linalg.lapack.dpotrf(M)[0])[0]
    Q = Q + Q.T
    n = len(Q)
    Q[range(n), range(n)] = Q[range(n), range(n)] / 2
    return Q
