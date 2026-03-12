"""Solve standard-form linear program by Dikin's method."""

import numpy


def Dik_m(c, A, b, er):
    """Solve standard-form linear program min c'x s.t. Ax = b by Dikin's method.

    Parameters
    ----------
    c : numpy.ndarray
        Objective coefficients (n-vector).
    A : numpy.ndarray
        Constraint matrix (m x n).
    b : numpy.ndarray
        Constraint right-hand side (m-vector).
    er : float
        Stopping criterion (max discrepancy between iterations).

    Returns
    -------
    X : numpy.ndarray
        Solution vector.
    z : float
        Optimal objective value.
    D : numpy.ndarray
        Centering diagonal matrix D.
    """
    from .errors import UnboundedFunctionError

    (_m, n) = numpy.shape(A)
    c = numpy.vstack((c, 1000))

    A = numpy.hstack((A, b - numpy.dot(A, numpy.ones((n, 1)))))
    (_m, n) = numpy.shape(A)
    xo = numpy.ones((n, 1))
    D = numpy.diagflat(xo)
    A1 = numpy.dot(A, D)
    c1 = numpy.dot(D, c)
    P1 = numpy.eye(n)
    NN = numpy.linalg.inv(numpy.dot(A1, A1.T))
    P2 = -numpy.dot(numpy.dot(A1.T, NN), A1)
    P = P1 + P2
    d = -numpy.dot(P, c1)
    t = -numpy.min(d)

    if t <= 0:
        raise UnboundedFunctionError('!!! Error: the objective function(z) is unbounded!!')
    else:
        x1 = numpy.ones((n, 1)) + (.9 / t) * d
        x = numpy.dot(D, x1)
        xx = 100 * numpy.ones((n, 1))
    i = 0
    numpy.set_printoptions(precision=4, linewidth=150)

    while numpy.max(numpy.abs(x - xx)) > er:
        i = i + 1
        xx = x
        D = numpy.diagflat(x)
        A1 = numpy.dot(A, D)
        c1 = numpy.dot(D, c)
        NN = numpy.linalg.inv(numpy.dot(A1, A1.T))

        P = numpy.eye(n) - numpy.dot(numpy.dot(A1.T, NN), A1)
        d = -numpy.dot(P, c1)
        t = -numpy.min(d)
        if t <= 0:
            raise UnboundedFunctionError('!!! Error: the objective function(z) is unbounded!!')

        x1 = numpy.ones((n, 1)) + (.9 / t) * d
        x = numpy.dot(D, x1)
    X = x[0:n - 1]
    D = numpy.diagflat(X)

    return X
