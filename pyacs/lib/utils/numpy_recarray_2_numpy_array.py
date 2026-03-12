"""Convert recarrays/structured arrays to plain numpy arrays."""


def numpy_recarray_2_numpy_array(A):
    """Convert a structured array (homogeneous dtype) to a numpy array.

    Parameters
    ----------
    A : numpy structured array
        Input array.

    Returns
    -------
    numpy.ndarray
        Plain array.
    """
    import numpy as np

    return np.array(A.view(A.dtype[0]).reshape(A.shape + (-1,)))

