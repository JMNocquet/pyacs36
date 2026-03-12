"""Solve weighted least-squares with data standard deviations."""

from .ls import ls


def lsw(G, d, std):
    """Solve weighted least-squares with data standard deviations.

    Parameters
    ----------
    G : numpy.ndarray
        m x n design matrix.
    d : numpy.ndarray
        m observation vector.
    std : array_like
        Standard deviations for d (length m).

    Returns
    -------
    numpy.ndarray
        Solution vector.

    Notes
    -----
    System is normalized so that G <- (G.T/std).T, d <- d/std, then solved by ordinary LS.
    """
    GG = (G.T / std).T
    dd = d / std
    return ls(GG, dd)
