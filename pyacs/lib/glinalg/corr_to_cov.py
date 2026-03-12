"""Convert correlation matrix and standard deviations to covariance matrix."""

import numpy as np


def corr_to_cov(corr, sigma_m):
    """Convert correlation matrix and standard deviations to covariance matrix.

    Parameters
    ----------
    corr : numpy.ndarray
        Correlation matrix.
    sigma_m : array_like
        Vector of standard deviations (sqrt of diagonal of covariance).

    Returns
    -------
    numpy.ndarray
        Covariance matrix.
    """
    outer_v = np.outer(sigma_m, sigma_m)
    return corr * outer_v
