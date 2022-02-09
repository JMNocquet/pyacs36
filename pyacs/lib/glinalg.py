"""
Linear algebra for Geodesy problems
"""

###############################################################################
# LEAST-SQUARES
###############################################################################
"""
Routines to solve linear systems using least-squares 
"""

###############################################################################
def ls(G,d, verbose=False):
###############################################################################
    """
    Solve the least-squares (LS) problem m so that  |Gx-d|**2 is minimum.
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param verbose: verbose mode
    
    :return: x,chi2: m (1D numpy array of dim m), chi2 (chi-square)
    
    :note: solved through numpy.linalg.lstsq
    """

    # numpy linalg lstsq
    
    import numpy.linalg

    (m,chi2,rank,s)=numpy.linalg.lstsq(G,d,rcond=-1)

    if verbose:
        print('-- ls info:')
        print('-- chi2: ', chi2 )
        print('-- rank: ', rank )
        print('-- s   : ', s )

    return(m)

###############################################################################
def lsw(G,d,std):
###############################################################################
    """
    Solve the least-squares (LS) with data uncertainties provided as a vector
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param std: standard deviation vector for d

    :note: the system is modified to be solved by ordinary LS by the change G<- (G.T/std).T and d<- d/std
    """

    # full solution
    
    GG=(G.T/std).T
    dd=d/std

    return ls(GG,dd)

###############################################################################
def lscov(G,d,cov,method='chol'):
###############################################################################
    """
    Solve the least-squares (LS) problem with data covariance

    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param cov: covariance matrix for d

    """
    
    import numpy as np

    if method=='chol':
    
        # cholesky decomposition of cov
    
        L=np.linalg.cholesky(cov)
        SQRT_Wd=np.linalg.inv(L)
      
        # linear system
        
        GG=np.dot(SQRT_Wd,G)
        dd=np.dot(SQRT_Wd,d)
        
        return ls(GG,dd)
    
###############################################################################
def lsw_full(G,d,std,verbose=False):
###############################################################################
    """
    Solve the least-squares (LS) with data uncertainties provided as a vector and returns the posterior covariance
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param std: standard deviation vector for d

    """

    import numpy as np
    import scipy.linalg.lapack

    # normalize system

    GG=(G.T/std).T
    dd=d/std

    N=np.dot(GG.T,GG)
    
    # solve

    L , info = scipy.linalg.lapack.dpotrf(N, False, False)

    if verbose:
        print('-- info from scipy.linalg.lapack.dpotrf in lsw_full')
        print('-- info: %d' , info)
    
    L_inv , info = scipy.linalg.lapack.dpotri(L)

    if verbose:
        print('-- info from scipy.linalg.lapack.dpotri in lsw_full')
        print('-- info: %d' , info)
    
    
    COV = np.triu(L_inv) + np.triu(L_inv, k=1).T
    
    X   = np.dot( COV , np.dot( GG.T , dd ) )

    # residuals
    V = d - np.dot(G,X)

    return X,COV,V

###############################################################################
def lscov_full(G,d,cov,verbose=False):
###############################################################################
    """
    Solve the least-squares (LS) with data covariance and returns the posterior covariance
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param cov: covariance matrix for d

    """

    import numpy as np
    import scipy.linalg.lapack

    # normalize system

    # cholesky decomposition of cov

    L=np.linalg.cholesky(cov)
    SQRT_Wd=np.linalg.inv(L)
  
    # linear system
    
    GG=np.dot(SQRT_Wd,G)
    dd=np.dot(SQRT_Wd,d)

    N=np.dot(GG.T,GG)
    

    # solve

    L , info = scipy.linalg.lapack.dpotrf(N, False, False)

    if info < 0:
        print('-- the i-th argument had an illegal value. i=%d' % info)
    if info > 0:
        print('-- the the leading minor of order i is not positive definite'+\
              ' and the factorization could not be completed. i=%d' % info)
    if info == 0:
        if verbose:
            print('-- dpotrf ok')

    
    L_inv , info = scipy.linalg.lapack.dpotri(L)

    if info < 0:
        print('-- the i-th argument had an illegal value. i=%d' % info)
    if info > 0:
        print('-- if INFO = i, the (i,i) element of the factor U or L  is zero, and the inverse could not be computed. i=%d' % info)
    if info == 0:
        if verbose:
            print('-- dpotri ok')

    
    COV = np.triu(L_inv) + np.triu(L_inv, k=1).T
    
    X   = np.dot( COV , np.dot( GG.T , dd ) )

    # residuals
    V = d - np.dot(G,X)

    # chi2
    
    vv = np.dot(SQRT_Wd,V)
    
    chi2 = np.dot(vv.T,vv)

    return X,COV,V, chi2


###############################################################################
# COVARIANCE
###############################################################################

"""
Routines to manipulate covariance and correlations matrices 
"""


###############################################################################
def cov_to_invcov(M):
###############################################################################
    """
    Inverse a covariance matrix
    
    :param cov: 2D numpy array covariance matrix
    
    :return: 2D numpy array inverse covariance matrix
    """
    
    from scipy import linalg
    Q = linalg.lapack.dpotri(linalg.lapack.dpotrf(M)[0])[0]
    Q = Q + Q.T
    n = len(Q)
    Q[range(n),range(n)] = Q[range(n),range(n)] / 2
    
    return Q


###############################################################################
def corr_to_cov(corr,sigma_m):
###############################################################################
    """
    Correlation to covariance matrix

    :param corr: correlation matrix
    :param sigma_m: vector of standard deviation = sqrt(diag(Cov))
    
    :return covariance matrix
    """

    import numpy as np
    outer_v = np.outer(sigma_m, sigma_m)
    return corr * outer_v


###############################################################################
def cov_to_corr(Cov):
###############################################################################
    """
    Covariance to correlation transformation
    
    :param cov: covariance matrix
    
    :return corr,sigma_m: correlation matrix and standard deviation matrix
    """
    import numpy as np

    v = np.sqrt(np.diag(Cov))
    outer_v = np.outer(v, v)
    correlation = Cov / outer_v
    correlation[ Cov == 0] = 0
    return correlation , v


###############################################################################
def symmetrize(a,type_matrix):
###############################################################################
    """
    Extract the upper or lower triangle and make a symmetric matrix
    
    :parameter a: numpy 2D array, must be a square matrix
    :parameter type: 'triu' or 'tril', triangle used to form the matrix
    """
    
    import numpy as np
    
    (r,c) = a.shape
    if r != c:
        raise TypeError('!!! ERROR: matrix must be square. a has ghape: %d ' , str(a.shape) )
        
    if type_matrix not in ['triu','tril']:
        raise TypeError('!!! ERROR: type_matrix must be triu or tril. type_matrix=%s' , type_matrix )
    
    if type_matrix == 'triu':
        T = np.triu(a)
        
    if type_matrix == 'tril':
        T = np.tril(a)
    
    return T + T.T - np.diag(a.diagonal())



###############################################################################
# VARIOUS
###############################################################################

"""
Various usefull routines to build linear systems 
"""

###############################################################################
def dot_and_sum(LX,a):
###############################################################################
    """
    does a matrix by scalar product and sum it
    :param LX : list of arrays with the same dimension
    :param a  : list of scalars
    """
    import numpy as np
    
    np_LX = np.array(LX)
    np_a  = np.array(a)
    
    return( np.dot(np_LX.T,np_a).T )

###############################################################################
def repeat_matrix_in_col(G,n):
###############################################################################
    """
    
    Repeat a matrix along a column
    
    R = np.empty( ( n,G.shape[0], G.shape[1] ) )
    R[:] = G
    return( R.reshape( n*G.shape[0], G.shape[1])  )    
           
    """
    import numpy as np
    
    R = np.empty( ( n,G.shape[0], G.shape[1] ) )
    R[:] = G
    
    return( R.reshape( n*G.shape[0], G.shape[1])  )    

###############################################################################
def odot(a,G):
###############################################################################
    """
    
    Customize vector by matrix product.
    
    Let's assume we have a matrix made of equal dimension sub-matrices $G_1,\cdots, G_n$
    \begin{equation}
    \left[
    \begin{array}{c}
    G_1\\
    G_2\\
    \vdots \\
    G_n
    \end{array}
    \right]
    \end{equation}
    and a vector of scalars
    \begin{equation}
    \left[
    \begin{array}{c}
    a_1\\
    a_2\\
    \vdots \\
    a_n
    \end{array}
    \right]
    \end{equation}
    and we want
    \begin{equation}
    \left[
    \begin{array}{c}
    a_1 G_1\\
    a_2 G_2\\
    \vdots \\
    a_n G_n
    \end{array}
    \right]
    \end{equation}
    We do this with numpy broadcasting

    :param a: 1D numpy array of multipliers (scalars)
    :param G: 2-D numpy arrays of submatrix G_i.
    
    :return: result of the multiplication
    :note: if G.shape = (n,m), and a.shape=(l), then the shape of the G_i is (n/l,m) 
    
    """
    
    # import
    
    import numpy as np

    # check dimension
    
    if G.shape[0] % a.shape[0] != 0:
        raise TypeError("!!! ERROR: Input parameters are not compatible")

    # number of row for the individual matrices to be multiplied

    nr = G.shape[0] // a.shape[0]

    # adding a dimension to G
    
    GG = G.reshape( a.shape[0] , nr, G.shape[1] )
    
    # adding a dimension to a
    
    aa = np.expand_dims(np.expand_dims(a,axis=0), axis=0).T

    R = GG * aa
    
    # reshape
    
    return( R.reshape( G.shape[0], G.shape[1] ) )

# ###############################################################################
# def dot_array_scalar_in_array(X,A):
# ###############################################################################
#     """
#     Assuming X.shape = (m,k,k) and A.shape = (m,n,n)
#     
#     Returns a 2D numpy array of the form:
# 
#     [[ (sum X[i]*A[i,0,0]), (sum X[i]*A[i,1,0]), ...., (sum X[i]*A[i,n,0])],
#      [ (sum X[i]*A[i,1,0]), (sum X[i]*A[i,1,1]), ...., (sum X[i]*A[i,n,1]),
#      ..
#      [ (sum X[i]*A[i,n,0]), (sum X[i]*A[i,n,1]), ...., (sum X[i]*A[i,n,n])]
#      
#     :param X : 2D numpy array
#     :param a  : 2D numpy array
#     """
#     
#     import numpy as np
#     
#     # TEST
#     
#     X = np.array([[[1,2],[2,1]],[[3,4],[4,3]]]) 
#     A = np.arange(18).reshape(2,3,3)
#     
#     
#     SOL = np.zeros((A.shape[1]*X.shape[1],A.shape[1]*X.shape[1]))
# 
#     for i in np.arange( A.shape[0]+1 ):
#         for j in np.arange( A.shape[0]+1 ):
# 
#             XX = np.zeros((X.shape[0],X.shape[0]))
#             for k in np.arange(X.shape[0]):
#                 for l in np.arange(A.shape[0]):
#                     XX += X[k] * A[l,i,j]
# 
#             SOL[i*X.shape[0]:(i+1)*X.shape[0],j*X.shape[0]:(j+1)*X.shape[0] ] = XX
#     
#     print np.sum( X.reshape(1,2,X.shape[0],X.shape[0]) * A.reshape(A.shape[1]*A.shape[2],A.shape[0],1,1),axis=1)
#    
#    return( None )

###############################################################################
# PAUL REBISCHUNG ROUTINES
###############################################################################

"""
Various usefull routines provided by Paul Rebischung (IGN/LAREG Paul.Rebischung@ign.fr) 
"""
###############################################################################
def dot(A, B):
###############################################################################
    """
    Matrix/matrix, matrix/vector and vector/vector multiplication
    :author  : P. Rebischung
    :date:Created : 02-Aug-2013
    :changes :
    :param A: Matrix or vector
    :param B : Matrix or vector
    :return:  A*B
    """

    from scipy import linalg
    
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

###############################################################################
def syminv(M):
###############################################################################
    """
    Invert a symmetric matrix
    :author  : P. Rebischung
    :created : 02-Aug-2013
    :param M : Matrix
    :return:Inverse matrix
    """
    
    from scipy import linalg

    
    Q = linalg.lapack.dpotri(linalg.lapack.dpotrf(M)[0])[0]
    Q = Q + Q.T
    n = len(Q)
    Q[range(n),range(n)] = Q[range(n),range(n)] / 2
    
    return Q

###############################################################################
def sympinv(M, verbose=False):
###############################################################################
    """
    Pseudo-invert a positive semi-definite symmetric matrix
    :author: P. Rebischung
    :created : 02-Aug-2013
    :param M : Matrix
    :return: Pseudo-inverse
    """
    
    import numpy , sys 
    from scipy import linalg
    
    
    (l, v, info) = linalg.lapack.dsyev(M)
    
    if verbose:
        print('-- info from linalg.lapack.dsyev in sympinv')
        print('-- info: %d' , info)

    t = sys.float_info.epsilon * numpy.max(l)
    l[numpy.nonzero(l <= t)[0]] = numpy.inf
    
    return dot(v/l, v.T)

###############################################################################
def make_normal_system( A , d , inv_Cd ):
###############################################################################
    """
    Given the linear system A x = d with Cd the covariance matrix of d
    the associated normal system is A.T Cd-1 A x = A.T Cd-1 d.
    returns N= A.T Cd-1 A, Nd=A.T Cd-1 d
    """
    
    # import
    
    import numpy as np
    
    TMP = np.dot( A.T , inv_Cd )
    
    return( np.dot( TMP , A ) , np.dot( TMP , d.reshape(-1,1) ) )


###############################################################################
def make_normal_system_tarantola( G , d , m0, inv_Cd, inv_Cm ):
###############################################################################
    """
    returns Gt inv_Cd G + inv_Cm , Gt inv_Cd d + inv_Cm m0
    """
    
    # import
    
    import numpy as np
    
    TMP = np.dot( G.T , inv_Cd )

    if np.ndim(m0) == 0:
        return( np.dot( TMP , G ) + inv_Cm  , np.dot( TMP , d.reshape(-1,1) ) + inv_Cm * m0  )
    else:
        return( np.dot( TMP , G ) + inv_Cm  , np.dot( TMP , d.reshape(-1,1) ) + np.dot(inv_Cm , m0.reshape(-1,1)  ) )


###############################################################################
def matrix_from_pattern( pattern , structure ):
###############################################################################
    """
    creates a matrix made of a pattern according to a structure
    
    :param pattern : pattern 2D numpy array to be duplicated
    :param structure: 2D numpy array giving the structure as multiplicating factors
    :return: the block matrix as a 2D numpy array
    
    Examples
    --------
    >>> pattern = np.arange(6).reshape(2,3)
    >>> structure = np.arange(6).reshape(3,2)

    >>> print(pattern)
    >>> print('--')
    >>> print(structure)
    >>> print('--')
    >>> print( matrix_from_pattern( pattern , structure ) )

        [[0 1 2]
        [3 4 5]]
        --
        [[0 1]
         [2 3]
         [4 5]]
        --
        [[ 0  0  0  0  1  2]
         [ 0  0  0  3  4  5]
         [ 0  2  4  0  3  6]
         [ 6  8 10  9 12 15]
         [ 0  4  8  0  5 10]
         [12 16 20 15 20 25]]
    
    
    
    """
    
    import numpy as np

    npr , npc = pattern.shape
    nsr , nsc = structure.shape
    
    pattern_3D = np.expand_dims( np.expand_dims( pattern , axis=0 ) , axis=0)
    structure_3D = np.expand_dims( np.expand_dims( structure , axis=0 ), axis = 0 )

    R_4D =  structure_3D.T * pattern_3D
    
    return np.moveaxis(R_4D ,[0,1,2,3],[2, 0, 1, 3]).reshape( npr*nsr , npc*nsc )


###############################################################################
def extract_block_diag(a, n, k=0):
###############################################################################
    import numpy as np
    
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError("Only 2-D arrays handled")
    if not (n > 0):
        raise ValueError("Must have n >= 0")

    if k > 0:
        a = a[:,n*k:] 
    else:
        a = a[-n*k:]

    n_blocks = min(a.shape[0]//n, a.shape[1]//n)

    new_shape = (n_blocks, n, n)
    new_strides = (n*a.strides[0] + n*a.strides[1],
                   a.strides[0], a.strides[1])

    return np.lib.stride_tricks.as_strided(a, new_shape, new_strides)
