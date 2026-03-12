"""Matrix/matrix, matrix/vector, or vector/vector multiplication (BLAS)."""

from scipy import linalg


def dot(A, B):
    """Matrix/matrix, matrix/vector, or vector/vector multiplication (BLAS).

    Parameters
    ----------
    A : numpy.ndarray
        Matrix or vector.
    B : numpy.ndarray
        Matrix or vector (compatible shape).

    Returns
    -------
    numpy.ndarray
        A * B.
    """
    if (A.ndim == 2) and (B.ndim == 2):
        if (A.flags['C_CONTIGUOUS'] and B.flags['C_CONTIGUOUS']):
            C = linalg.blas.dgemm(1.0, A, B)
        elif A.flags['C_CONTIGUOUS']:
            C = linalg.blas.dgemm(1.0, A, B.T, trans_b=True)
        elif B.flags['C_CONTIGUOUS']:
            C = linalg.blas.dgemm(1.0, A.T, B, trans_a=True)
        else:
            C = linalg.blas.dgemm(1.0, A.T, B.T, trans_a=True, trans_b=True)

    elif A.ndim == 2:
        if A.flags['C_CONTIGUOUS']:
            C = linalg.blas.dgemv(1.0, A, B)
        else:
            C = linalg.blas.dgemv(1.0, A.T, B, trans=True)

    else:
        C = linalg.blas.ddot(A, B)

    return C
