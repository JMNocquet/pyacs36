"""
RobustEstimators.py includes a number of robust estimators
to solve linear problems.
"""

class Error(Exception):
    """Base class for exceptions in module robustestimators.py"""
    pass

class UnboundedFunctionError(Exception):
    """Exception raised for unbounded objective function."""
    pass

######################################################################################################
def Dikin(A,y,W,eps=3.E-3):
######################################################################################################
    """
    L1-norm estimation of parameters x using the Dikin's method
    using Linear optimization in linear model  y=Ax+e
    
    :param A: Model matrix
    :param y: observed values (observation vector)
    :param W: Weight matrix of observables (of DIAGONAL type)
    :param eps: small value for iteration (default: eps = 1e-6)
    
    :return : x,e: vector of estimated parameters and residuals
    
    :note: translated from Matlab code kindly provided by Amir Khodabandeh june.2009
    reference:Recursive Algorithm for L1 Norm Estimation in Linear Models,
    A. Khodabandeh and A. R. Amiri-Simkooei, JOURNAL OF SURVEYING ENGINEERING ASCE / FEBRUARY 2011 / 1
    doi:10.1061/ASCESU.1943-5428.0000031
    Translated to python from Matlab original by J.-M. Nocquet July 2011
    """
    
    import numpy

    if y.ndim == 1:
        y=y.reshape(-1,1)


#%%%%%%%Transformation%%%%%%%%%%%%%%
    
    if W is not None:
        w=numpy.sqrt(W)
        
        A=(A.T*w).T
        y=(y.T*w).T
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    (m,n)=numpy.shape(A)
    
    At=numpy.hstack((A,-A,numpy.eye(m),-numpy.eye(m)))
    f=numpy.vstack((numpy.zeros((2*n,1)),numpy.ones((2*m,1))))
    
#%%%%%% Call to main function %%%%%%%   
    dlin=Dik_m(f,At,y,eps)
    x=dlin[0:n]-dlin[n:2*n]
    e=dlin[2*n:2*n+m]-dlin[2*n+m:2*(n+m)]

#%%%%%%% transform back to the original residuals
    if W is not None:
        e=(e.T/ w).T

    # return
    
    return(x.flatten(),e.flatten())

######################################################################################################
def Dik_m(c,A,b,er):
######################################################################################################
    """
    subject:solve the standard form of linear programming by affine/Dikin's
    method("an interior point method") 
    minimize z=c'*x; subject to Ax=b;
    input:(c):coefficients of objective function(z) as n-vector(hint:a column
    vector)
    (A):matrix of constraint set with size m*n
    (b):m-vector of constraint set 
    (er):maximum discrepancy between two iteration.("stopping criterion")
    output:(X):unknown variables
    (z):optimal value of objective function
    (D):Centering transformer "D"(a diagonal matrix)
    Amir khodabandeh Oct.2008
    """
    import numpy
    
    (_m,n)=numpy.shape(A)
    c=numpy.vstack((c,1000))
    
    A=numpy.hstack((A,b-numpy.dot(A,numpy.ones((n,1)))))
    (_m,n)=numpy.shape(A);
    xo=numpy.ones((n,1))
    D=numpy.diagflat(xo)
    A1=numpy.dot(A,D)
    c1=numpy.dot(D,c)
    P1=numpy.eye(n)
#    NN=numpy.dot(A1,A1.T).I
    
    NN=numpy.linalg.inv(numpy.dot(A1,A1.T))
    P2=-numpy.dot(numpy.dot(A1.T,NN),A1)
    P=P1+P2
    d=-numpy.dot(P,c1)
    t=-numpy.min(d)

    if t<=0:
            raise UnboundedFunctionError('!!! Error: the objective function(z) is unbounded!!')
    else:
        x1=numpy.ones((n,1))+(.9/t)*d
        x=numpy.dot(D,x1)
        xx=100*numpy.ones((n,1))
    i=0
    numpy.set_printoptions(precision=4,linewidth=150   )

    while numpy.max(numpy.abs(x-xx))>er:
        i=i+1
        xx=x
        D=numpy.diagflat(x)
        A1=numpy.dot(A,D)
        c1=numpy.dot(D,c)
        NN=numpy.linalg.inv(numpy.dot(A1,A1.T))

        P=numpy.eye(n)-numpy.dot(numpy.dot(A1.T,NN),A1)
        d=-numpy.dot(P,c1)
        t=-numpy.min(d)
        if t<=0:
            raise UnboundedFunctionError('!!! Error: the objective function(z) is unbounded!!')
        
        x1=numpy.ones((n,1))+(.9/t)*d
        x=numpy.dot(D,x1)
#    z=numpy.dot(c.T,x)
    X=x[0:n-1]
    D=numpy.diagflat(X)
    
    return(X)

# ######################################################################################################
# def ransac():
# ######################################################################################################
#     """
#     Robust linear model estimation using sklearn/RANSAC
#     http://scikit-learn.org/stable/auto_examples/linear_model/plot_ransac.html
#     """
#     import numpy as np
#     from sklearn import linear_model, datasets
