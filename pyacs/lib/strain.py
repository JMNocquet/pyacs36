"""
Strain rate calculation library
"""

###############################################################################
def vgrad(X,Y,VE,VN,SVE=None,SVN=None,CVEN=None,CODE=None, method='WLS',verbose=False):
###############################################################################
    """
    Linear estimates of a horizontal velocity gradient
    
    :param X,Y: 1D numpy array of longitudes and latitudes
    :param VE,VN: 1D numpy array of east and north velocity. Unit is assummed to be mm/yr
    :param SVE,SVN,CVEN: 1D numpy array of east and north velocity standard deviation (mm/yr) and correlation.
    :param CODE: 1D numpy of string with site codes
    :param method: estimator in 'L1','WLS' (default: weighted least-squares 'WLS')
    
    :return: DV, the velocity gradient tensor as a 2x2 2D-numpy array and a 2-columns (East and North) residuals 2D numpy array 
     
    """
    
    import numpy as np
    
    ###########################################################################
    def __vgrad_obs_eq__(l0,p0,l,p,ve,vn,sve,svn,cven):
    ###########################################################################
        """
        observation equation for the horizontal velocity gradient
        we assume that observation are in mm/yr and result will be in nstrain / yr
        """

        Rt = 6371.E3 # Earth's radius in metres

        Rt3nstrain = Rt / 1.E6

        import pyacs.lib.glinalg
        
        Ai = np.zeros((2,6))
        Bi = np.zeros(2)


        Corri = np.zeros((2,2))

        # velocity at the barycenter
                
        Ai[0,0] = 1.
        Ai[1,1] = 1.
        
        # ve
        
        Ai[0,2] = np.radians(l-l0) * Rt3nstrain * np.cos( np.radians( (p+p0)/2. ) )
        Ai[0,3] = np.radians(p-p0) * Rt3nstrain
        
        # vn
        
        Ai[1,4] = np.radians(l-l0) * Rt3nstrain * np.cos( np.radians( (p+p0)/2. ) )
        Ai[1,5] = np.radians(p-p0) * Rt3nstrain
        
        # Observation vector
        
        Bi[0] = ve
        Bi[1] = vn
        
        # Covariance matrix
        
        Corri[0,0] = 1.
        Corri[1,1] = 1.
        Corri[1,0] = cven
        Corri[0,1] = cven
        
        Cvi = pyacs.lib.glinalg.corr_to_cov( Corri,np.array([sve,svn]) )
                                
        return(Ai,Bi,Cvi)
    
    ###########################################################################
    def __barycenter__(X,Y):
    ###########################################################################
        """
        the barycenter of a set of point given by geographical coordinates
        """
        
        import pyacs.lib.coordinates
        
        XYZ = np.array( list(map(pyacs.lib.coordinates.geo2xyz,np.radians(X),np.radians(Y),np.zeros(X.size)) ))

        [mx,my,mz] = np.mean( XYZ, axis=0 )

        [mlong,mlat, _ ] = pyacs.lib.coordinates.xyz2geo(mx, my, mz, unit='dec_deg')

        return(mlong,mlat)
    
    # MAIN
    
    # barycenter of network
    
    (l0,p0) = __barycenter__(X,Y)
    
    if verbose:
        print("-- Barycenter at (%10.5lf , %10.5lf )" % (l0,p0))


    # build the linear system

    A=np.zeros((2*X.size,6))
    B=np.zeros((2*X.size))
    CV=np.zeros((2*X.size,2*X.size))

    for i in np.arange(X.size):
        (Ai,Bi,Cvi) = __vgrad_obs_eq__(l0,p0,X[i],Y[i],VE[i],VN[i],SVE[i],SVN[i],CVEN[i])
        
        A[2*i:2*i+2,:] = Ai
        B[2*i:2*i+2] = Bi
        CV[2*i:2*i+2,2*i:2*i+2] = Cvi 
    
    # Solve
    
    import pyacs.lib.glinalg
    
    SOL,COV,RESIDUALS, chi2 = pyacs.lib.glinalg.lscov_full(A,B,CV)
    
    
    return(l0,p0,SOL,COV,RESIDUALS, chi2)

        
