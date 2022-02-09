#-------------------------------------------------------------------------------
# Module   : mathutils
# Purpose  : Mathematical routines
# Author   : P. Rebischung
# Created  : 22-Aug-2014
#
# Changes  :
#
# Routines : - dot           : Matrix/matrix, matrix/vector and vector/vector multiplication
#            - syminv        : Invert a symmetric matrix
#            - sympinv       : Pseudo-invert a positive semi-definite symmetric matrix
#            - rms           : Compute RMS of a vector
#            - trend         : Fit linear model to a time series
#            - detrend       : Remove linear trend from a time series
#            - cart2geo      : Transform cartesian coordinates into geographical coordinates
#            - geo2cart      : Transform geographical coordinates into cartesian coordinates
#            - xyz2enh       : Compute rotation matrices from geocentric to topocentric coordinates
#            - vondrak       : Pass a Vondrak filter through a time series
#            - lombscargle   : Compute Lomb-Scargle periodogram of a time series
#            - kolmo2d       : 2D Kolmogorov-Smirnov test wrt uniform distribution over [0,1] x [0,1]
#            - q2c           : Compute correlation matrix from covariance matrix
#            - chebychev     : Compute Chebychev polynomials
#            - spline_smooth : Cubic spline smoothing of data
#-------------------------------------------------------------------------------

# LIBRARIES
#-------------------------------------------------------------------------------
import sys
from math import *
import numpy
from scipy import linalg
from scipy import interpolate as interp

from .constants import *


#-------------------------------------------------------------------------------
# Routine : dot
# Purpose : Matrix/matrix, matrix/vector and vector/vector multiplication
# Author  : P. Rebischung
# Created : 02-Aug-2013
#
# Changes :
#
# Input   : - A : Matrix or vector
#           - B : Matrix or vector
# Output  : - C : A*B
#-------------------------------------------------------------------------------
def dot(A, B):

  # Matrix/matrix multiplication
  if ((A.ndim == 2) and (B.ndim == 2)):
    if (A.flags['C_CONTIGUOUS'] and B.flags['C_CONTIGUOUS']):
      C = linalg.blas.dgemm(1.0, A, B)
    elif (A.flags['C_CONTIGUOUS']):
      C = linalg.blas.dgemm(1.0, A, B.T, trans_b=True)
    elif (B.flags['C_CONTIGUOUS']):
      C = linalg.blas.dgemm(1.0, A.T, B, trans_a=True)
    else:
      C = linalg.blas.dgemm(1.0, A.T, B.T, trans_a=True, trans_b=True)

  # Matrix/vector multiplication
  elif (A.ndim == 2):
    if (A.flags['C_CONTIGUOUS']):
      C = linalg.blas.dgemv(1.0, A, B)
    else:
      C = linalg.blas.dgemv(1.0, A.T, B, trans=True)

  # Vector/vector multiplication
  else:
    C = linalg.blas.ddot(A, B)

  return C


#-------------------------------------------------------------------------------
# Routine : syminv
# Purpose : Invert a symmetric matrix
# Author  : P. Rebischung
# Created : 02-Aug-2013
#
# Changes :
#
# Input   : - M : Matrix
# Output  : - Q : Inverse matrix
#-------------------------------------------------------------------------------
def syminv(M):

  Q = linalg.lapack.dpotri(linalg.lapack.dpotrf(M)[0])[0]
  Q = Q + Q.T
  n = len(Q)
  Q[list(range(n)), list(range(n))] = Q[list(range(n)), list(range(n))] / 2
  
  return Q


#-------------------------------------------------------------------------------
# Routine : sympinv
# Purpose : Pseudo-invert a positive semi-definite symmetric matrix
# Author  : P. Rebischung
# Created : 02-Aug-2013
#
# Changes :
#
# Input   : - M : Matrix
# Output  : - Q : Pseudo-inverse
#-------------------------------------------------------------------------------
def sympinv(M):

  (l, v, info) = linalg.lapack.dsyev(M)
  t = sys.float_info.epsilon * numpy.max(l)
  l[numpy.nonzero(l <= t)[0]] = numpy.inf

  return dot(v / l, v.T)


#-------------------------------------------------------------------------------
# Routine : rms
# Purpose : Compute RMS of a vector
# Author  : P. Rebischung
# Created : 07-Jul-2011
#
# Changes :
#
# Input   : - x : Data vector
# Output  : - s : RMS
#-------------------------------------------------------------------------------
def rms(x):

  return sqrt(dot(x, x) / len(x))


#-------------------------------------------------------------------------------
# Routine : trend
# Purpose : Fit linear model to a time series
# Author  : P. Rebischung
# Created : 17-Nov-2011
#
# Changes :
#
# Input   : - t : Dates
#           - x : Time series
# Output  : - a : Coefficients of the linear fit
#-------------------------------------------------------------------------------
def trend(t, x):

  n = len(t)
  St = sum(t)
  Sx = sum(x)
  Stx = sum(t * x)
  Stt = sum(t ** 2)
  d = n * Stt - St ** 2

  a = numpy.zeros(2)
  a[0] = (Stt * Sx - St * Stx) / d
  a[1] = (n * Stx - St * Sx) / d

  return a


#-------------------------------------------------------------------------------
# Routine : detrend
# Purpose : Remove linear trend from a time series
# Author  : P. Rebischung
# Created : 17-Nov-2011
#
# Changes :
#
# Input   : - t : Dates
#           - x : Time series
# Output  : - d : Detrended time series
#-------------------------------------------------------------------------------
def detrend(t, x):

  a = trend(t, x)
  d = x - a[0] - a[1] * t

  return d


#-------------------------------------------------------------------------------
# Routine : cart2geo
# Purpose : Transform cartesian coordinates into geographical coordinates
# Author  : P. Rebischung
# Created : 25-May-2011
#
# Changes :
#
# Input   : - X   : [X, Y, Z] (cartesian coordinates in m)
# Output  : - phi : Latitudes  (rad)
#           - lam : Longitudes (rad)
#           - h   : Ellipsoidal heights (m)
#-------------------------------------------------------------------------------
def cart2geo(X):

  reshape = False
  if (X.shape == (3,)):
    reshape = True
    X.resize(1, 3)

  p = numpy.sqrt(X[..., 0] ** 2 + X[..., 1] ** 2)
  r = numpy.sqrt(X[..., 0] ** 2 + X[..., 1] ** 2 + X[..., 2] ** 2)
  u = numpy.arctan2(X[..., 2] / p * (1 - fe + ee ** 2 * ae / r), 1)

  lam = 2 * numpy.arctan2(X[..., 1], X[..., 0] + p)
  phi = numpy.arctan2(X[..., 2] * (1 - fe) + ee ** 2 * ae * numpy.sin(u) ** 3, (1 - fe) * (p - ee ** 2 * ae * numpy.cos(u) ** 3))
  h = p * numpy.cos(phi) + X[..., 2] * numpy.sin(phi) - ae * numpy.sqrt(1 - ee ** 2 * numpy.sin(phi) ** 2)
  
  # lam = numpy.arctan2(X[...,1], X[...,0])
  # phi = numpy.arctan2(X[...,2], numpy.sqrt(X[...,0]**2 + X[...,1]**2))
  # h = numpy.sqrt(X[...,0]**2 + X[...,1]**2 + X[...,2]**2) - ae

  if (reshape):
    X.resize(3)
    phi = phi[0]
    lam = lam[0]
    h = h[0]

  return (phi, lam, h)


#-------------------------------------------------------------------------------
# Routine : geo2cart
# Purpose : Transform geographical coordinates into cartesian coordinates
# Author  : P. Rebischung
# Created : 25-May-2011
#
# Changes :
#
# Input   : - phi : Latitudes  (rad)
#           - lam : Longitudes (rad)
#           - h   : Ellipsoidal heights (m)
# Output  : - X   : [X, Y, Z] (cartesian coordinates in m)
#-------------------------------------------------------------------------------
def geo2cart(phi, lam, h):

  N = ae / numpy.sqrt(1 - (ee * numpy.sin(phi)) ** 2)

  if (type(phi) == type(float())):
    X = numpy.zeros(3,)
  else:
    X = numpy.zeros(phi.shape + (3,))

  X[..., 0] = (N + h) * numpy.cos(phi) * numpy.cos(lam)
  X[..., 1] = (N + h) * numpy.cos(phi) * numpy.sin(lam)
  X[..., 2] = (N * (1 - ee ** 2) + h) * numpy.sin(phi)

  # X[...,0] = ae * numpy.cos(phi) * numpy.cos(lam)
  # X[...,1] = ae * numpy.cos(phi) * numpy.sin(lam)
  # X[...,2] = ae * numpy.sin(phi)

  return X


#-------------------------------------------------------------------------------
# Routine : xyz2enh
# Purpose : Compute rotation matrices from geocentric to topocentric coordinates
# Author  : P. Rebischung
# Created : 25-May-2011
#
# Changes :
#
# Input   : - X : [X, Y, Z] (cartesian coordinates in m)
# Output  : - R : Rotation matrices from geocentric to topocentric coordinates
#-------------------------------------------------------------------------------
def xyz2enh(X):

  reshape = False
  if (X.shape == (3,)):
    reshape = True
    X.resize(1, 3)

  (phi, lam, h) = cart2geo(X)

  R = numpy.zeros(phi.shape + (3, 3))
  R[..., 0, 0] = -numpy.sin(lam)
  R[..., 0, 1] = numpy.cos(lam)
  R[..., 1, 0] = -numpy.sin(phi) * numpy.cos(lam)
  R[..., 1, 1] = -numpy.sin(phi) * numpy.sin(lam)
  R[..., 1, 2] = numpy.cos(phi)
  R[..., 2, 0] = numpy.cos(phi) * numpy.cos(lam)
  R[..., 2, 1] = numpy.cos(phi) * numpy.sin(lam)
  R[..., 2, 2] = numpy.sin(phi)

  if (reshape):
    X.resize(3)
    R.resize((3, 3))

  return R


#-------------------------------------------------------------------------------
# Routine : enh2uvh
# Purpose : Compute rotation matrices from topocentric frames of two stations to their UVH frame
# Author  : P. Rebischung
# Created : 04-Oct-2016
#
# Changes :
#
# Input   : - X1 : Cartesian coordinates of first station (m)
#           - X2 : Cartesian coordinates of second station (m)
# Output  : - R1 : Rotation matrix from ENH frame to UVH frame of first station
#           - R2 : Rotation matrix from ENH frame to UVH frame of second station
#-------------------------------------------------------------------------------
def enh2uvh(X1, X2):

  (p1, l1, h1) = cart2geo(X1)
  (p2, l2, h2) = cart2geo(X2)
    
  az = atan2(sin(l2 - l1), cos(p1) * tan(p2) - sin(p1) * cos(l2 - l1))
  s = sin(az)
  c = cos(az)
  
  R1 = numpy.zeros((3, 3))
  R1[0, 0] = s
  R1[0, 1] = c
  R1[1, 0] = -c
  R1[1, 1] = s
  R1[2, 2] = 1
  
  az = atan2(sin(l1 - l2), cos(p2) * tan(p1) - sin(p2) * cos(l1 - l2))
  s = sin(az + pi)
  c = cos(az + pi)
  
  R2 = numpy.zeros((3, 3))
  R2[0, 0] = s
  R2[0, 1] = c
  R2[1, 0] = -c
  R2[1, 1] = s
  R2[2, 2] = 1
  
  return (R1, R2)


#-------------------------------------------------------------------------------
# Routine : vondrak
# Purpose : Pass a Vondrak filter through a time series
# Author  : P. Rebischung
# Created : 17-Dec-2011
#
# Changes :
#
# Input   : - t               : Dates
#           - x               : Time series
#           - fc              : Cutoff frequency
#           - return_partials : True if the partial derivatives dd/dy should be
#                               returned. Default is False.
# Output  : - xs              : Filtered time series
#-------------------------------------------------------------------------------
def vondrak(t, x, fc, return_partials=False):

  import numpy
  import scipy.linalg as linalg
  
  eps = (7.223147119819503 * fc) ** 6 / (len(t) - 3)
  num = 6 * numpy.sqrt(t[2:-1] - t[1:-2])
  den = numpy.sqrt(t[-1] - t[0])

  a = numpy.hstack((0, 0, 0, num / den / ((t[0:-3] - t[1:-2]) * (t[0:-3] - t[2:-1]) * (t[0:-3] - t[3:])), 0, 0, 0))
  b = numpy.hstack((0, 0, 0, num / den / ((t[1:-2] - t[0:-3]) * (t[1:-2] - t[2:-1]) * (t[1:-2] - t[3:])), 0, 0, 0))
  c = numpy.hstack((0, 0, 0, num / den / ((t[2:-1] - t[0:-3]) * (t[2:-1] - t[1:-2]) * (t[2:-1] - t[3:])), 0, 0, 0))
  d = numpy.hstack((0, 0, 0, num / den / ((t[3:] - t[0:-3]) * (t[3:] - t[1:-2]) * (t[3:] - t[2:-1])), 0, 0, 0))

  d0 = eps + a[3:] ** 2 + b[2:-1] ** 2 + c[1:-2] ** 2 + d[0:-3] ** 2
  d1 = a[3:-1] * b[3:-1] + b[2:-2] * c[2:-2] + c[1:-3] * d[1:-3]
  d2 = a[3:-2] * c[3:-2] + b[2:-3] * d[2:-3]
  d3 = a[3:-3] * d[3:-3]

  A = numpy.diag(d0) + numpy.diag(d1, 1) + numpy.diag(d1, -1) + numpy.diag(d2, 2) + numpy.diag(d2, -2) + numpy.diag(d3, 3) + numpy.diag(d3, -3)

  if (return_partials):
    A = eps * syminv(A)
    xs = dot(A, x)
    return (xs, A)
  else:
    xs = linalg.solve(A, eps * x)
    return xs


#-------------------------------------------------------------------------------
# Routine : lombscargle
# Purpose : Compute Lomb-Scargle periodogram of a time series
# Author  : P. Rebischung
# Created : 17-Dec-2011
#
# Changes :
#
# Input   : - t         : Dates
#           - x         : Time series
#           - sf        : Oversampling factor. Default is 4.
#           - f         : Frequencies. Default is None (automatically set).
#           - normalize : True for a normalized periodogram. Default is False.
# Output  : - f         : Frequencies
#           - p         : Powers
#-------------------------------------------------------------------------------
def lombscargle(t, x, sf=4, f=None, normalize=False):

  # Length and time span of the series
  n = len(x)
  T = t[-1] - t[0]
  
  # Normalize time series if needed
  if (normalize):
    x = x - numpy.mean(x)
    x = x / numpy.std(x)

  # Build list of frequencies
  if (type(f).__name__ == 'NoneType'):
    f0 = (n - 1) / (n * T * sf)
    fc = (n - 1) / (2 * T)
    f = numpy.arange(f0, fc + f0, f0)

  # Initialize periodogram
  p = numpy.zeros(len(f))

  # Loop over frequencies
  for i in range(len(f)):
    w = 2 * pi * f[i]
    
    # Compute some sums
    xc = numpy.sum(x * numpy.cos(w * t))
    xs = numpy.sum(x * numpy.sin(w * t))
    cc = numpy.sum(numpy.cos(w * t) ** 2)
    ss = numpy.sum(numpy.sin(w * t) ** 2)
    cs = numpy.sum(numpy.cos(w * t) * numpy.sin(w * t))

    # Offset
    tau = atan2(2 * cs, cc - ss) / (2 * w)
    ctau = cos(w * tau)
    stau = sin(w * tau)

    # Power at frequency i
    p[i] = (ctau * xc + stau * xs) ** 2 / (ctau ** 2 * cc + 2 * ctau * stau * cs + stau ** 2 * ss) + (ctau * xs - stau * xc) ** 2 / (ctau ** 2 * ss - 2 * ctau * stau * cs + stau ** 2 * cc)
    
  return (f, p)


#-------------------------------------------------------------------------------
# Routine : kolmo2d
# Purpose : 2D Kolmogorov-Smirnov test wrt uniform distribution over [0,1] x [0,1]
# Author  : P. Rebischung
# Created : 08-Feb-2012
#
# Changes :
#
# Input   : - x, y : Coordinates of data points
# Output  : - D    : Maximal absolute difference between theoretical and empirical
#                    cumulative distributions
#-------------------------------------------------------------------------------
def kolmo2d(x, y):

  # Initializations
  D = 0
  n = len(x)

  # Loop over data points
  for i in range(n):
    
    # Count number of other data points in each quadrant
    n1a = numpy.sum(numpy.logical_and(x < x[i], y < y[i]))
    n2a = numpy.sum(numpy.logical_and(x < x[i], y > y[i]))
    n3a = numpy.sum(numpy.logical_and(x > x[i], y < y[i]))
    n4a = numpy.sum(numpy.logical_and(x > x[i], y > y[i]))
    n1b = numpy.sum(numpy.logical_and(x <= x[i], y <= y[i]))
    n2b = numpy.sum(numpy.logical_and(x <= x[i], y >= y[i]))
    n3b = numpy.sum(numpy.logical_and(x >= x[i], y <= y[i]))
    n4b = numpy.sum(numpy.logical_and(x >= x[i], y >= y[i]))

    # Theoretical probabilities for each quadrant
    p1 = x[i] * y[i]
    p2 = x[i] * (1 - y[i])
    p3 = (1 - x[i]) * y[i]
    p4 = (1 - x[i]) * (1 - y[i])

    # Compare empirical and theoretical distributions
    d1a = abs(p1 - float(n1a) / float(n))
    d1b = abs(p1 - float(n1b) / float(n))
    d2a = abs(p2 - float(n2a) / float(n))
    d2b = abs(p2 - float(n2b) / float(n))
    d3a = abs(p3 - float(n3a) / float(n))
    d3b = abs(p3 - float(n3b) / float(n))
    d4a = abs(p4 - float(n4a) / float(n))
    d4b = abs(p4 - float(n4b) / float(n))

    # Update max absolute difference if neccesary
    d = max([d1a, d1b, d2a, d2b, d3a, d3b, d4a, d4b])
    if (d > D):
      D = d

  return D


#-------------------------------------------------------------------------------
# Routine : Q2C
# Purpose : Compute correlation matrix from covariance matrix
# Author  : P. Rebischung
# Created : 21-Feb-2012
#
# Changes :
#
# Input   : - Q : Covariance matrix
# Output  : - C : Correlation matrix
#-------------------------------------------------------------------------------
def Q2C(Q):

  f = 1. / numpy.sqrt(numpy.diag(Q))

  return f * (Q * f).T


#-------------------------------------------------------------------------------
# Routine : chebychev
# Purpose : Compute Chebychev polynomials
# Author  : P. Rebischung
# Created : 14-Mar-2012
#
# Changes :
#
# Input   : - x : Array of abscissas in [-1,1]
#           - n : Number of polynomials to compute
# Output  : - y : Values of polynomials in x
#-------------------------------------------------------------------------------
def chebychev(x, n):

  # Initialization
  y = numpy.zeros((len(x), n))

  # Degree 0 and 1 polynomials
  y[:, 0] = 1
  if (n > 1):
    y[:, 1] = x

  # Degree >=2 polynomials
  for i in range(2, n):
    y[:, i] = 2 * x * y[:, i - 1] - y[:, i - 2]

  return y


#-------------------------------------------------------------------------------
# Routine : spline_smooth
# Purpose : Cubic spline smoothing of data
# Author  : P. Rebischung
# Created : 26-Sep-2013
#
# Changes :
#
# Input   : - x               : Abscissas
#           - y               : Ordinates
#           - s               : Sigmas
#           - l               : Smoothing parameter
#           - return_partials : True if the partial derivatives dd/dy should be
#                               returned. Default is False.
# Output  : - d               : Smoothed ordinates
#-------------------------------------------------------------------------------
def spline_smooth(x, y, s, l, return_partials=False):

  # Useful things
  n = len(x) - 1
  mu = 2 * (1 - l) / (3 * l)
  h = x[1:] - x[:-1]
  p = 2 * (h[:-1] + h[1:])
  r = 3 / h
  f = -(r[:-1] + r[1:])

  # Some matrices
  S = numpy.diag(s)
  R = numpy.diag(p) + numpy.diag(h[1:-1], 1) + numpy.diag(h[1:-1], -1)
  Q = numpy.zeros((n + 1, n - 1))
  Q[0:n - 1, :] = Q[0:n - 1, :] + numpy.diag(r[:-1])
  Q[1:n + 0, :] = Q[1:n + 0, :] + numpy.diag(f)
  Q[2:n + 1, :] = Q[2:n + 1, :] + numpy.diag(r[1:])

  # Compute the partial derivatives dd/dy
  D = syminv(mu * dot(Q.T, dot(S, Q)) + R)
  D = numpy.eye(n + 1) - mu * dot(S, dot(Q, dot(D, Q.T)))

  # Smoothed data
  d = dot(D, y)

  if (return_partials):
    return (d, D)
  else:
    return d
