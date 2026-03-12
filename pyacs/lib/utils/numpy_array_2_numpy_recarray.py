"""Convert numpy arrays to recarrays."""


def numpy_array_2_numpy_recarray(A, names):
    """Convert a numpy array to a numpy recarray.

    Parameters
    ----------
    A : numpy.ndarray
        2D array.
    names : list of str
        Field names for each column.

    Returns
    -------
    numpy.recarray
        Record array with the given names.
    """
    import numpy as np

    return np.rec.array(
        np.core.records.array(
            list(tuple(A[:, :].T)),
            dtype={'names': names, 'formats': list(map(np.dtype, A[0, :]))},
        )
    )

