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
    import pyacs.glinalg

    if method=='chol':
    
        # cholesky decomposition of cov
    
        L=np.linalg.cholesky(cov)
        SQRT_Wd=np.linalg.inv(L)
      
        # linear system
        
        GG=np.dot(SQRT_Wd,G)
        dd=np.dot(SQRT_Wd,d)
        
        return pyacs.glinalg.ls(GG,dd)
