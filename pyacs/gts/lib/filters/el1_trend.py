"""
Taken from https://github.com/wudingcheng/TSAnalyzer/blob/master/TSAnalyzer/algorithms/l1extensive.py
"""

import cvxpy as cvx
import numpy as np
from scipy.sparse import eye, csr_matrix, hstack, linalg
from scipy.signal import argrelextrema


###############################################################################
def el1_trend(self, lam, rho, periods=None, in_place=False, return_offset=False, return_periodic=False, verbose=True, component='NEU'):
###############################################################################
    """
    extensive l1 trend filtering

    :param lam: weight of regularization of filtered data
    :param rho: weight of regularization of offsets
    :param period: tuple, periods to be estimated (1.,0.5) will estimate annual and semi-annual terms
    :param in_place: if True then replace the current time series
    :param return_offset: if True also return offset time serie
    :param return_periodic: if True also return periodic time serie
    :param verbose: boolean, verbose mode
    :param component: string. Default 'NEU'
    :return: the filtered time series

    """
    # import
    from pyacs.gts.lib.message import message

    ### copy
    new_gts_f = self.copy(data_xyz=None)
    new_gts_o = self.copy(data_xyz=None)
    new_gts_s = self.copy(data_xyz=None)

    new_gts_f.data[:,4:8] = 1.E-3
    new_gts_o.data[:, 4:8] = 1.E-3
    new_gts_s.data[:, 4:8] = 1.E-3

### filter
    if 'N' in component:
        if verbose:
            message('Computing extensive L1 trend filter for component North')
        new_gts_f.data[:, 1],new_gts_o.data[:, 1],new_gts_s.data[:, 1] = l1filter(self.data[:, 0], self.data[:, 1], lam=lam,rho=rho, periods=periods, solver=cvx.MOSEK, verbose=verbose)

    if 'E' in component:
        if verbose:
            message('Computing extensive L1 trend filter for component East')
        new_gts_f.data[:, 2],new_gts_o.data[:, 2],new_gts_s.data[:, 2] = l1filter(self.data[:, 0], self.data[:, 2], lam=lam,rho=rho, periods=periods, solver=cvx.MOSEK, verbose=verbose)

    if 'U' in component:
        if verbose:
            message('Computing extensive L1 trend filter for component Up')
        new_gts_f.data[:, 3],new_gts_o.data[:, 3],new_gts_s.data[:, 3] = l1filter(self.data[:, 0], self.data[:, 3], lam=lam,rho=rho, periods=periods, solver=cvx.MOSEK, verbose=verbose)

    ### return
    if in_place:
        self = new_gts_f

    if return_offset and return_periodic:
        return new_gts_f,new_gts_o, new_gts_s

    if return_offset and not return_periodic:
        return new_gts_f, new_gts_o

    if not return_offset and return_periodic:
        return new_gts_f, new_gts_s

    return new_gts_f

def gen_d2(n):
    """
    Generate the 2nd difference matrix.
    :param n: int, length of time series
    :return: csr_matrix, sparse matrix
    """
    I2 = eye(n - 2)
    O2 = csr_matrix((n - 2, 1))
    return hstack((I2, O2, O2)) + hstack((O2, -2 * I2, O2)) + hstack((O2, O2, I2))


def gen_d1(n):
    """
    Generate the 1st difference matrix
    :param n: int, length of time series
    :return: csr_matrix, sparse matrix
    """
    I1 = eye(n - 1)
    O1 = csr_matrix((n - 1, 1))
    return hstack((I1, O1)) - hstack((O1, I1))


def get_max_lam(y):
    """
    Calculate the max lambda value for given time series y
    :param y: np.array, time series given
    :return: float, max lambda value
    """
    D = gen_d2(len(y))
    ddt = D.dot(D.T)
    dy = D.dot(y)
    return np.linalg.norm(linalg.spsolve(ddt, dy), np.inf)


def get_max_rho(y):
    """
    Calculate the max rho value for given time series y,
    :param y: np.array, time series given
    :return: float, max rho value
    """
    D = gen_d1(len(y))
    ddt = D.dot(D.T)
    dy = D.dot(y)
    return np.linalg.norm(linalg.spsolve(ddt, dy), np.inf)


def l1filter(t, y,
             lam=1200,
             rho=80,
             periods=(365.25, 182.625),
             solver=cvx.MOSEK,
             verbose=False):
    """
    Do l1 regularize for given time series.
    :param t: np.array, time
    :param y: np.array, time series value
    :param lam: lambda value
    :param rho: rho value
    :param periods: list, periods, same unit as t
    :param solver: cvx.solver
    :param verbose: bool, show verbose or not
    :return: x, w, s, if periods is not None, else return x, w
    """
    t = np.asarray(t, dtype=np.float64)
    t = t - t[0]
    y = np.asarray(y, dtype=np.float64)

    assert y.shape == t.shape

    n = len(t)
    D = gen_d2(n)

    x = cvx.Variable(n)
    w = cvx.Variable(n)
    errs = y - x - w
    seasonal = None
    if periods:
        tpi_t = 2 * np.pi * t
        for period in periods:
            a = cvx.Variable()
            b = cvx.Variable()
            temp = a * np.sin(tpi_t / period) + b * np.cos(tpi_t / period)
            if seasonal is None:
                seasonal = temp
            else:
                seasonal += temp
        errs = errs - seasonal
    obj = cvx.Minimize(0.5 * cvx.sum_squares(errs) +
#                       lam * cvx.sum_squares(D * x) +
#                       lam * cvx.norm(D * x, 2) +
                       lam * cvx.norm(D * x, 1) +
                       rho * cvx.tv(w))
    prob = cvx.Problem(obj)
    prob.solve(solver=solver, verbose=verbose)
    if periods:
        return np.array(x.value), np.array(w.value), np.array(seasonal.value)
    else:
        return np.array(x.value), np.array(w.value), None

    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(t)
    x = cvx.Variable(n)
    w = cvx.Variable(n)
    dx = cvx.mul_elemwise(1.0 / np.diff(t), cvx.diff(x))
    x_term = cvx.tv(dx)
    dw = cvx.mul_elemwise(1.0 / np.diff(t), cvx.diff(w))
    w_term = cvx.norm(dw, 1)
    errs = y - x - w
    seasonal = None
    if periods:
        tpi_t = 2 * np.pi * t
        for period in periods:
            a = cvx.Variable()
            b = cvx.Variable()
            temp = a * np.sin(tpi_t / period) + b * np.cos(tpi_t / period)
            if seasonal is None:
                seasonal = temp
            else:
                seasonal += temp
        errs = errs - seasonal

    obj = cvx.Minimize(0.5 * cvx.sum_squares(errs) + lam * x_term + rho * w_term)
    prob = cvx.Problem(obj)
    print('solving')
    prob.solve(solver=solver, verbose=verbose)
    print(x)
    if periods:
        return np.array(x.value), np.array(w.value), np.array(seasonal.value)
    else:
        return np.array(x.value), np.array(w.value), None




