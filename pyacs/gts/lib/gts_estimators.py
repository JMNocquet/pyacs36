
import numpy as np

###################################################################
## LEAST-SQUARES - old version written by Trong Tran
###################################################################

def least_square(A,L,P=None):
    """
    Least-squares estimation for AX + L = 0 with optional weight matrix P.

    Parameters
    ----------
    A : ndarray
        Design matrix.
    L : ndarray
        Observation vector.
    P : ndarray, optional
        Weight matrix for L (default identity).

    Returns
    -------
    tuple
        X (unknowns), s_X (parameter std), V (residuals), s_V (residual std),
        std (sigma_0), S (model A*X).

    Raises
    ------
    ValueError
        If N is singular or has negative diagonal in Q.
    """
    if isinstance(P,np.ndarray):
        if len(P.shape)==2:
            ATP = np.dot(np.transpose(A),P)
            N = np.dot(ATP,A)
            M = np.dot(ATP,L)
        else:
            N = np.dot(A.T * P , A )
            M = np.dot (A.T * P , L)
    else:
        N = np.dot(np.transpose(A),A)
        M = np.dot(np.transpose(A),L)

    if np.linalg.det(N) == 0: 
        X, V, std, s_X, s_V, S = None, None, None, None, None, None
        raise ValueError("!!! Null determinant from np.linalg.det" )
    else:
        Q = np.linalg.inv(N)
        if len(np.argwhere(np.diag(Q)<0))!=0: 
            X, V, std, s_X, s_V, S = None, None, None, None, None, None
            raise ValueError("!!! Singular matrix from np.linalg.inv")
        else:               
            # unknowns
            X = np.dot(Q,M)
            # residuals
            V = L - np.dot(A,X)
            # model
            S = np.dot(A,X)

            #reduced
            if isinstance(P,np.ndarray): 
                if len(P.shape)==2:
                    chi2 = np.dot(np.dot(np.transpose(V),P),V)
                else:
                    chi2 = np.dot(V.T * P,V)
                    
            else: chi2 = np.dot(np.transpose(V),V)

            #standard deviation
            if  len(V) != len(X):              
                std = np.sqrt(chi2/(len(V)-len(X)))
            else:
                std = 0.0
            #variance of unknowns
            s_X = std*np.sqrt(np.diag(Q))
            
            #variance of residuals: s_V = std*sqrt(C-A(Q)AT)
            if isinstance(P,np.ndarray): 
                if len(P.shape)==2:
                    s_V = std*np.sqrt(1./np.diag(P)-np.diag(np.dot(np.dot(A,Q), np.transpose(A))))
                else:
                    s_V = std*np.sqrt(1./P)
                    
            else: s_V = std*np.ones(len(L))
    return X, s_X, V, s_V, std, S
