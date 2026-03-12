"""Helmert (similarity) transformation estimation and matrix construction."""

from pyacs.sol.errors import PyacsSol_HelmertError
import inspect


###############################################################################
def helmert_matrix(x,y,z,tran=True,rot=True,scale=True,equilibrate=True):
###############################################################################
    """Generate Helmert transformation matrix for a single point.

    Builds the design matrix for translation/rotation/scale transformation.

    Parameters
    ----------
    x : float
        Point geocentric X coordinate (m).
    y : float
        Point geocentric Y coordinate (m).
    z : float
        Point geocentric Z coordinate (m).
    tran : bool, optional
        Include translation. Default is True.
    rot : bool, optional
        Include rotation. Default is True.
    scale : bool, optional
        Include scale. Default is True.
    equilibrate : bool, optional
        If True (default), use km for position and mas for scale for better conditioning.

    Returns
    -------
    numpy.matrix
        Transformation design matrix (2D).
    """    
    
    import numpy as np
    
    # equilibrate
    if equilibrate:
        c=np.pi/180./3600. # factor used to equilibrate Helmert's equation
        m=1.E-6
    else:
        c=1.
        m=1.
    
    # initialization
    A = np.zeros([3,7], float)
    
    if tran:
        A[0,0] = 1.
        A[1,1] = 1.
        A[2,2] = 1.
    
    if rot:
        A[0,5] = z*c
        A[0,6] = -y*c
        A[1,4] = -z*c
        A[1,6] = x*c
        A[2,4] = y*c
        A[2,5] = -x*c
        
    if scale:
        A[0,3] = x*m
        A[1,3] = y*m
        A[2,3] = z*m

    return np.matrix(A)

###############################################################################
def observation_equation_helmert_local(xf,yf,zf,xr,yr,zr,tran=True,rot=True,scale=True):
###############################################################################
    """Generate linear observation equations for Helmert transformation (E,N,U).

    Parameters
    ----------
    xf, yf, zf : float
        Geocentric coordinates in the initial ('free') frame (m).
    xr, yr, zr : float
        Geocentric coordinates in the final (reference) frame (m).
    tran : bool, optional
        Include translation. Default is True.
    rot : bool, optional
        Include rotation. Default is True.
    scale : bool, optional
        Include scale. Default is True.

    Returns
    -------
    A : numpy.matrix
        Design matrix.
    B : numpy.matrix
        Observation vector.

    Notes
    -----
    VCV (variance-covariance) is not accounted for yet.
    """

    import numpy as np
    import pyacs.lib.coordinates
    
    B = np.zeros(3, float)
    B[0] = xr-xf
    B[1] = yr-yf
    B[2] = zr-zf
    B = np.transpose(np.matrix(B))
    
    (lam,phi,_h)=pyacs.lib.coordinates.xyz2geo(xr, yr, zr)
    R = pyacs.lib.coordinates.mat_rot_general_to_local(lam,phi)
    
    B = np.dot(R,B)
    A = helmert_matrix(xr,yr,zr,tran=True,rot=True,scale=True)

    A = np.dot(R,A)   
    
    return A, B

###############################################################################
def estimate_helmert(H_Gpoint_ref,H_Gpoint_free,
            vcv_xyz=None,psd=None,
            tran=True,rot=True,scale=True,
            lexclude=[],
            method='Dikin_LS', threshold=5.,
            verbose=True):
###############################################################################
    
    """Estimate Helmert (or similarity) parameters between two coordinate sets.

    Default is 7-parameter transformation (translation, rotation, scale).

    Parameters
    ----------
    H_Gpoint_ref : dict
        Dictionary of Gpoint instances used as reference.
    H_Gpoint_free : dict
        Dictionary of Gpoint instances used as free (to be transformed).
    vcv_xyz : array-like, optional
        Variance-covariance of coordinates. Default is None.
    psd : object, optional
        Post-seismic deformation object from sinex library.
    tran : bool, optional
        Include translation. Default is True.
    rot : bool, optional
        Include rotation. Default is True.
    scale : bool, optional
        Include scale. Default is True.
    lexclude : list, optional
        Site codes to exclude. Default is [].
    method : str, optional
        'LS' (least-squares) or 'Dikin_LS' (L1 then L2). Default is 'Dikin_LS'.
    threshold : float, optional
        Outlier rejection threshold. Default is 5.
    verbose : bool, optional
        If True, print progress. Default is True.

    Returns
    -------
    object
        Result of the estimation (transformation parameters and related info).
    """

    from icecream import ic

    import numpy as np
    import pyacs.lib.glinalg
    import pyacs.lib.robustestimators as RobustEstimators


    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    VERBOSE("Helmert parameters estimation using method: %s" % method)

    n_points= len(list(H_Gpoint_ref.keys()))

    A=np.zeros((3*n_points,7))
    B=np.zeros((3*n_points,1))

    i=0
    for code,soln in list(H_Gpoint_ref.keys()):

        Xref=H_Gpoint_ref[code,soln]

        Xref.assign_index(i)
        (xr,yr,zr)=Xref.posxyz()
        epoch_ref=Xref.epoch
        
        Xfree= H_Gpoint_free[code,soln]
        (xf,yf,zf)=Xfree.posxyz()
        #(sxf, syf, szf) = Xfree.staxyz()
        epoch_free=Xfree.epoch
        
        (vxr,vyr,vzr)=Xref.velxyz()
        delta_t=epoch_free-epoch_ref
        x=xr+delta_t*vxr
        y=yr+delta_t*vyr
        z=zr+delta_t*vzr

        if psd is not None:
            from pyacs.sinex import snxutils
            
            import pyacs.lib.astrotime
            import pyacs.lib.coordinates

            epoch = pyacs.lib.astrotime.decyear2epoch(Xfree.epoch)
            # Compute ENH post-seismic deformations
            (denh, senh) = snxutils.compute_psd(psd, code, epoch)
            if np.any(denh):
                VERBOSE("psd correction (mm, ENU): %s %10.2lf %10.2lf %10.2lf " % (code,denh[0] * 1000., denh[1] * 1000., denh[2] * 1000.)  )
            # Compute XYZ post-seismic deformations
            (l,p,he)=pyacs.lib.coordinates.xyz2geo(xf, yf, zf)
            R = pyacs.lib.coordinates.mat_rot_local_to_general(l, p)
            dxyz = np.dot(R, denh.reshape(3,-1)).flatten()
            # Add post-seismic deformations
            x=x+dxyz[0]
            y=y+dxyz[1]
            z=z+dxyz[2]


        (Ai,Bi)=observation_equation_helmert_local(xf,yf,zf,x,y,z)

        A[3*i:3*i+3,:]=Ai
        B[3*i:3*i+3,:]=Bi

        i=i+1
        
    if method == 'Dikin_LS':
    
    ### Robust estimators of Dikin + LS
        VERBOSE("Estimating Helmert parameters - L1 norm (Dikins)")
        
        
        X,V=RobustEstimators.Dikin(A,B,W=None)

        X=X.reshape(-1,1)
        V=V.reshape(-1,1)

    else:
        VERBOSE("Estimating Helmert parameters - L2 norm")

        X=pyacs.lib.glinalg.ls(A,B)
        V = B-np.dot(A,X)

        X=X.reshape(-1,1)
        V=V.reshape(-1,1)


    [median_east,median_north,median_up]=np.median(np.sqrt((V.reshape(-1,3)*1000.0)**2),axis=0)

### reject the outliers et re-adjust by LS

    VERBOSE("Testing outliers; Threshold for outliers detection : %.1lf " % threshold)

    np_outliers = find_outliers_Helmert(V,threshold=threshold)

    l_np_outliers=np_outliers.tolist()

    VERBOSE("# of outliers detected: %d" % len(l_np_outliers))

    # remove outliers for LS Helmert and flag

    FLAG_VV=V.reshape(-1,3)*0.0

    if len(l_np_outliers)!=0:


        loutliers=[]

        for [i,j] in l_np_outliers:

            FLAG_VV[i,j]=1
            loutliers.append(3*i+j)

        A_rejected=np.delete(A,loutliers,axis=0)
        B_rejected=np.delete(B,loutliers,axis=0)

        # raise PyacsSol_HelmertError if problem
        if A_rejected.shape[0] < 7:
            msg = 'Too many outliers. Cannot calculate Helmert Transformation'
            WARNING( msg )
            raise PyacsSol_HelmertError(inspect.stack()[0][3],__name__,msg)

        VERBOSE("Doing LS estimation without outliers")

        LS_X=pyacs.lib.glinalg.ls(A_rejected,B_rejected)
        LS_V = B-np.dot(A,X)

    #
    # else:
    # ### LS estimation, no outlier detection
    #
    #     #VERBOSE("Doing LS estimation without outlier detection")
    #     VERBOSE("Using L2 estimates with outliers rejection")
    #
    #     X=pyacs.lib.glinalg.ls(A,B)
    #     V = B-np.dot(A,X)
    #     l_np_outliers=[]
    #     FLAG_VV = V.reshape(-1, 3) * 0.0

### Creating residuals, accounting for rejected
    
    LS_V = B - np.dot(A,X)
      
    VV=LS_V.reshape(-1,3)
    
    S_VV=np.array(VV*1000.0,dtype=str)

    
    # adding code,soln
    
    S_VV=np.insert(S_VV,0,np.array(list(H_Gpoint_ref.keys()))[:,1],axis=1)
    S_VV=np.insert(S_VV,0,np.array(list(H_Gpoint_ref.keys()))[:,0],axis=1)

    # Adding flag
    
    S_VV=np.hstack((S_VV,FLAG_VV))
    
    # Statistics
    
    wrms=np.sqrt( np.sum ( (VV*1000.0)**2 / (FLAG_VV*1.E12 +1.),axis=0) / np.sum ( 1. / (FLAG_VV*1.E12 +1.),axis=0) )

    H_stat={}
    
    H_stat['rejection_threshold']=threshold
    H_stat['n_site']=len(list(H_Gpoint_ref.keys()))
    H_stat['n_obs_raw']=H_stat['n_site']*3
    H_stat['n_obs_rejected']=len(l_np_outliers)
    H_stat['percent_rejected']=float(len(l_np_outliers))/float(H_stat['n_obs_raw'])*100.0

    if method == 'Dikin_LS':
        H_stat['median_L1']=[median_east,median_north,median_up]
    else:
        H_stat['median_L1']=[0.0,0.0,0.0]
    
    H_stat['wrms']=wrms


    return (X,S_VV,H_stat)



###############################################################################
def find_outliers_Helmert(R,threshold=5):
###############################################################################
    """
    Returns index outliers points
    """
    from icecream import ic

    import numpy as np
    np.set_printoptions(precision=2,suppress=True)
    
    V=R.reshape(-1,3)*1000.0

    median=np.median(np.sqrt(V**2),axis=0)

    X=V/(threshold*median)
    
    np_outliers=np.argwhere(X**2>1)


    # If only East or North have been flagged then flag both
    
    for [i,j] in np_outliers:
        if j==0:
            np_outliers=np.append(np_outliers,np.array([[i,1]]),axis=0)
        if j==1:
            np_outliers=np.append(np_outliers,np.array([[i,0]]),axis=0)
    
    # accepted from Numpy >= 1.13
    #np_outliers=np.unique( np_outliers, axis=0 )

    if np.size(np_outliers) !=  0:

        np_outliers=np_outliers[np.argsort(np_outliers[:, 0])]
        #np_outliers = np.vstack({tuple(row) for row in np_outliers })
        #ic(np_outliers)
        np_outliers = np.vstack([tuple(row) for row in np_outliers ])
        #ic(np_outliers)

    return np_outliers


# ###############################################################################
# def LS_adjustment(A,L,P=None):
# ###############################################################################
#     if P!=None:
#         ATP = np.dot(np.transpose(A),P)
#         N = np.dot(ATP,A)
#         M = np.dot(ATP,L)
#     else:
#         N = np.dot(np.transpose(A),A)
#         M = np.dot(np.transpose(A),L)
#         
# 
#     Q=np.linalg.inv(N)
#     X = np.dot(Q,M)
#     V = L-np.dot(A,X)
#     
#     rows,cols=np.shape(A)
#     if P!=None:
#         sigma=np.sqrt(np.dot(np.dot(np.transpose(V),P),V)/(rows-cols))[0,0] #variance_unit_weight
#     else:
#         if rows != cols:
#             sigma=np.sqrt(np.dot(np.transpose(V),V)/(rows-cols)) #variance_unit_weight
#         else:
#             sigma=np.sqrt(np.dot(np.transpose(V),V)) #variance_unit_weight
#             
#     SX = sigma*np.sqrt(np.diag(Q))
#   
#     return X,SX,V,sigma

