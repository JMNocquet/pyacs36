"""Compute weighted sum of arrays."""

import numpy as np


def dot_and_sum(LX, a):
    """Compute weighted sum of arrays (matrix-by-scalar product then sum).

    Parameters
    ----------
    LX : list of array_like
        Arrays with the same shape.
    a : list of float
        Scalars (same length as LX).

    Returns
    -------
    numpy.ndarray
        Sum of a[i] * LX[i].
    """
    np_LX = np.array(LX)
    np_a = np.array(a)
    return np.dot(np_LX.T, np_a).T
