"""Form the Tarantola-style normal system (with prior)."""

import numpy as np


def make_normal_system_tarantola(G, d, m0, inv_Cd, inv_Cm):
    """Form the Tarantola-style normal system (with prior).

    Returns N = G.T inv(Cd) G + inv(Cm), b = G.T inv(Cd) d + inv(Cm) m0.

    Parameters
    ----------
    G : numpy.ndarray
        Design matrix.
    d : array_like
        Observation vector.
    m0 : array_like or float
        Prior model.
    inv_Cd : numpy.ndarray
        Inverse data covariance.
    inv_Cm : numpy.ndarray or float
        Inverse prior covariance.

    Returns
    -------
    tuple
        (N, b) normal matrix and right-hand side.
    """
    TMP = np.dot(G.T, inv_Cd)

    if np.ndim(m0) == 0:
        return np.dot(TMP, G) + inv_Cm, np.dot(TMP, d.reshape(-1, 1)) + inv_Cm * m0
    else:
        return np.dot(TMP, G) + inv_Cm, np.dot(TMP, d.reshape(-1, 1)) + np.dot(inv_Cm, m0.reshape(-1, 1))
