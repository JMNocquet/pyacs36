"""Pseudo-inverse of a positive semi-definite symmetric matrix."""

import numpy
import sys
from scipy import linalg

from .dot import dot


def sympinv(M, verbose=False):
    """Pseudo-inverse of a positive semi-definite symmetric matrix.

    Parameters
    ----------
    M : numpy.ndarray
        Symmetric matrix (positive semi-definite).
    verbose : bool, optional
        If True, print info. Default is False.

    Returns
    -------
    numpy.ndarray
        Pseudo-inverse matrix.
    """
    (l, v, info) = linalg.lapack.dsyev(M)

    if verbose:
        print('-- info from linalg.lapack.dsyev in sympinv')
        print('-- info: %d', info)

    t = sys.float_info.epsilon * numpy.max(l)
    l[numpy.nonzero(l <= t)[0]] = numpy.inf

    return dot(v / l, v.T)
