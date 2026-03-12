"""Convert covariance matrix to correlation and standard deviations."""

import numpy as np


def cov_to_corr(Cov):
    """Convert covariance matrix to correlation and standard deviations.

    Parameters
    ----------
    Cov : numpy.ndarray
        Covariance matrix.

    Returns
    -------
    correlation : numpy.ndarray
        Correlation matrix.
    v : numpy.ndarray
        Vector of standard deviations (sqrt of diagonal).
    """
    v = np.sqrt(np.diag(Cov))
    outer_v = np.outer(v, v)
    correlation = Cov / outer_v
    correlation[Cov == 0] = 0
    return correlation, v
