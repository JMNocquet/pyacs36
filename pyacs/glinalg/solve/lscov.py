###############################################################################
def lscov(G,d,cov,method='chol'):
###############################################################################
    """
    Solve the least-squares problem with data covariance matrix.

    Parameters
    ----------
    G : ndarray
        m x n model matrix (2D).
    d : ndarray
        m observation vector (1D).
    cov : ndarray
        m x m covariance matrix for d.
    method : str, optional
        'chol' for Cholesky (default).

    Returns
    -------
    ndarray
        Solution (same as ls on whitened system).
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
