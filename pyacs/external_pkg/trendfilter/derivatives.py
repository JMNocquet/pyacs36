"""
Sparse derivative matrices
"""

import numpy as np
from scipy.sparse import spdiags, coo_matrix, dia_matrix
from itertools import chain
import cvxpy


def first_derivative_matrix(n):
    """
    A sparse matrix representing the first derivative operator
    :param n: a number
    :return: a sparse matrix that applies the derivative operator
             to a numpy array or list to yield a numpy array
    """
    e = np.ones((1, n))
    return spdiags(np.vstack((-1*e, e)), range(2), n-1, n)


def first_derivative_matrix_circular(n):
    """
    A matrix representing the first derivative operator
    Circular so it wraps around. Not sparse matrix.
    :param n: a number
    :return: a numpy array that applies the derivative operator
             to a numpy array or list to yield a numpy array
    """
    return np.eye(n) - np.roll(np.eye(n), 1, axis=0)


def second_derivative_matrix(n):
    """
    A sparse matrix representing the second derivative operator
    :param n: a number
    :return: a sparse matrix that applies the second derivative operator
             to a numpy array or list to yield a numpy array
    """
    e = np.ones((1, n))
    return spdiags(np.vstack((e, -2*e, e)), range(3), n-2, n)


def first_derivative(x):
    """
    Derivative as a functional operator
    :param x: array with (N, n_time_periods) dimensions
    :return: first derivative iterable
             backward difference so will have 1 fewer elements
             derivative not defined for first element
    """
    return x[:, 1:] - x[:, 0:-1]


def second_derivative_matrix_nes(x, a_min=0.0, a_max=None, scale_free=False):
    """
    Get the second derivative matrix for non-equally spaced points
    :param : x numpy array of x-values
    :param : a_min, float. Allows for modification where the x-values
             can't get any closer than this. (Not technically the sec derv)
    :param : a_max, float. Same but max.
    :return: A matrix D such that if x.size == (n,1), D * x is the second derivative of x
    assumes points are sorted
    """
    n = len(x)
    m = n - 2

    values = []
    for i in range(1, n-1):
        # These are all positive if sorted
        a0 = float(x[i+1] - x[i])
        a2 = float(x[i] - x[i-1])

        assert (a0 >= 0) and (a2 >= 0), "Points do not appear to be sorted"

        # And neither cab be zero
        assert (a0 > 0) and (a2 > 0), "Second derivative doesn't exist for repeated points"

        # Now allow for a min and max on the differences
        # of the x separations
        if a_max is not None:
            a0 = min(a0, a_max)
            a2 = min(a2, a_max)

        a0 = max(a0, a_min)
        a2 = max(a2, a_min)
        a1 = a0 + a2

        if scale_free:
            # Just subtract derivatives
            # don't divide again by the length scale
            scf = a1/2.0
        else:
            scf = 1.0

        vals = [2.0*scf/(a1*a2), -2.0*scf/(a0*a2), 2.0*scf/(a0*a1)]
        values.extend(vals)

    i = list(chain(*[[_] * 3 for _ in range(m)]))
    j = list(chain(*[[_, _ + 1, _ + 2] for _ in range(m)]))

    d2 = coo_matrix((values, (i, j)), shape=(m, n))

    return dia_matrix(d2)


def first_derv_nes(x, y):
    n = len(x)
    ep = 1e-9
    idx = 1.0 / (x[1:] - x[0:-1] + ep)
    matrix = first_derivative_matrix(n)
    return idx * (matrix @ y)


def first_derv_nes_cvxpy(x, y):
    n = len(x)
    ep = 1e-9
    idx = 1.0 / (x[1:] - x[0:-1] + ep)
    matrix = first_derivative_matrix(n)
    return cvxpy.multiply(idx, matrix @ y)


def cumulative_matrix(n):
    """
    Matrix of nxn dimension that when multiplied with a vector
        results in the cumulative sum
    :param n: dimension
    :return: cumulative sum matrix
    """
    return np.tril(np.ones((n, n)))
