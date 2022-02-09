"""
Useful operations on vectors (1-D numpy array)
"""

###################################################################
def np_array(X):
###################################################################

    """
    Converts to 1-D numpy array
    
    :param X: list or tuple to be converted
    
    :return : 1D numpy array
    
    """
    
    import numpy as np
    
    return np.array(X) 

###################################################################
def norm(X):
###################################################################

    """
    Returns the norm of a vector

    :param X: list, tuple or 1D numpy array
    
    :return : norm (float)


    """

    import numpy as np
    
    return np.sqrt(np.sum(np_array(X)**2))


###################################################################
def normalize(X):
###################################################################

    """
    Normalize a vector
    
    :param X: list, tuple or 1D numpy array
    
    :return : 1D numpy array
    
    """

    return np_array(X) / norm(X)

###################################################################
def scal_prod(X,Y):
###################################################################

    """
    Returns the scalar products
    
    :param X,Y: list, tuple or 1D numpy array
    
    :return : scalar product as float
    
    """
   
    import numpy as np
    
    return np.sum( np_array(X) * np_array(Y) )

###################################################################
def vector_product(A,B):
###################################################################
    """
    Vector product
    
                          A[1]*B[2] - B[1]*A[2]
    vector_product(A,B)=  A[2]*B[0] - B[2]*A[0]
                          A[0]*B[1] - B[0]*A[1]
                          
    :param A,B: 1D numpy arrays
    
    :return AxB: 1D numpy array
    """
    
    import numpy as np
    
    C=np.zeros([3])
    C[0]=A[1]*B[2]-B[1]*A[2]
    C[1]=A[2]*B[0]-B[2]*A[0]
    C[2]=A[0]*B[1]-B[0]*A[1]
    
    return(C)



    