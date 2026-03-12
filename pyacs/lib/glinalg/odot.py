"""Scale block rows of G by vector a."""

import numpy as np


def odot(a, G):
    """Scale block rows of G by vector a (element-wise then reshape).

    G is treated as a stack of submatrices G_1,...,G_n; result is [a_1*G_1; ...; a_n*G_n].
    Implemented with numpy broadcasting. If G.shape = (n, m) and a.shape = (l,),
    then each block has n/l rows.

    Parameters
    ----------
    a : numpy.ndarray
        1D array of scalars (multipliers).
    G : numpy.ndarray
        2D array; number of rows must be divisible by len(a).

    Returns
    -------
    numpy.ndarray
        Matrix with same shape as G; each block of rows scaled by corresponding a[i].
    """
    if G.shape[0] % a.shape[0] != 0:
        raise TypeError("!!! ERROR: Input parameters are not compatible")

    nr = G.shape[0] // a.shape[0]

    GG = G.reshape(a.shape[0], nr, G.shape[1])

    aa = np.expand_dims(np.expand_dims(a, axis=0), axis=0).T

    R = GG * aa

    return R.reshape(G.shape[0], G.shape[1])
