"""Form a symmetric matrix from the upper or lower triangle."""

import numpy as np


def symmetrize(a, type_matrix):
    """Form a symmetric matrix from the upper or lower triangle.

    Parameters
    ----------
    a : numpy.ndarray
        2D square array.
    type_matrix : str
        'triu' or 'tril'; which triangle is used to build the symmetric matrix.

    Returns
    -------
    numpy.ndarray
        Symmetric matrix.
    """
    (r, c) = a.shape
    if r != c:
        raise TypeError('!!! ERROR: matrix must be square. a has ghape: %d ', str(a.shape))

    if type_matrix not in ['triu', 'tril']:
        raise TypeError('!!! ERROR: type_matrix must be triu or tril. type_matrix=%s', type_matrix)

    if type_matrix == 'triu':
        T = np.triu(a)

    if type_matrix == 'tril':
        T = np.tril(a)

    return T + T.T - np.diag(a.diagonal())
