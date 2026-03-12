"""Repeat a matrix vertically n times."""

import numpy as np


def repeat_matrix_in_col(G, n):
    """Repeat a matrix vertically n times (stack n copies).

    Parameters
    ----------
    G : numpy.ndarray
        2D matrix.
    n : int
        Number of repetitions.

    Returns
    -------
    numpy.ndarray
        Matrix of shape (n * G.shape[0], G.shape[1]).
    """
    R = np.empty((n, G.shape[0], G.shape[1]))
    R[:] = G
    return R.reshape(n * G.shape[0], G.shape[1])
