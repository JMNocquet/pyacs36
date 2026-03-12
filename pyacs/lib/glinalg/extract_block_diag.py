"""Extract block diagonal from a matrix."""

import numpy as np


def extract_block_diag(a, n, k=0):
    """Extract block diagonal from a 2D array.

    Parameters
    ----------
    a : array_like
        2D array.
    n : int
        Block size.
    k : int, optional
        Diagonal offset. Default is 0.

    Returns
    -------
    numpy.ndarray
        Block diagonal elements as 3D array of shape (n_blocks, n, n).
    """
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError("Only 2-D arrays handled")
    if not (n > 0):
        raise ValueError("Must have n >= 0")

    if k > 0:
        a = a[:, n * k:]
    else:
        a = a[-n * k:]

    n_blocks = min(a.shape[0] // n, a.shape[1] // n)

    new_shape = (n_blocks, n, n)
    new_strides = (n * a.strides[0] + n * a.strides[1],
                   a.strides[0], a.strides[1])

    return np.lib.stride_tricks.as_strided(a, new_shape, new_strides)
