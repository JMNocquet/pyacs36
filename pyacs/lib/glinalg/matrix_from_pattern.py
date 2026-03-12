"""Build a block matrix from a pattern scaled by a structure matrix."""

import numpy as np


def matrix_from_pattern(pattern, structure):
    """Build a block matrix from a pattern scaled by a structure matrix.

    Parameters
    ----------
    pattern : numpy.ndarray
        2D pattern array to be duplicated and scaled.
    structure : numpy.ndarray
        2D array of scaling factors (block structure).

    Returns
    -------
    numpy.ndarray
        Block matrix as 2D array.
    """
    npr, npc = pattern.shape
    nsr, nsc = structure.shape

    pattern_3D = np.expand_dims(np.expand_dims(pattern, axis=0), axis=0)
    structure_3D = np.expand_dims(np.expand_dims(structure, axis=0), axis=0)

    R_4D = structure_3D.T * pattern_3D

    return np.moveaxis(R_4D, [0, 1, 2, 3], [2, 0, 1, 3]).reshape(npr * nsr, npc * nsc)
