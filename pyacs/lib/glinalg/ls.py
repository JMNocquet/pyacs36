"""Solve least-squares problem."""

import numpy.linalg


def ls(G, d, verbose=False):
    """Solve the least-squares problem min |Gx - d|**2.

    Parameters
    ----------
    G : numpy.ndarray
        m x n design matrix.
    d : numpy.ndarray
        m observation vector.
    verbose : bool, optional
        If True, print chi2, rank, singular values. Default is False.

    Returns
    -------
    numpy.ndarray
        Solution vector of length n. Note: variable name in code is m.
    float
        Chi-square (residual sum of squares). Only first return value is used in code.

    Notes
    -----
    Solved via numpy.linalg.lstsq.
    """
    (m, chi2, rank, s) = numpy.linalg.lstsq(G, d, rcond=-1)

    if verbose:
        print('-- ls info:')
        print('-- chi2: ', chi2)
        print('-- rank: ', rank)
        print('-- s   : ', s)

    return m
