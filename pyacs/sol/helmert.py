"""
calculates Helmert transformation 
"""

from pyacs.sol.errors import PyacsSol_HelmertError
import inspect


###############################################################################
def helmert_matrix(x,y,z,tran=True,rot=True,scale=True,equilibrate=True):
###############################################################################
    """
    Generates a Helmert transformation matrix for a single point
    Can also create translation/rotation/scale transformation matrix
    
    :param x,y,z: point geocentric coordinates (m)
    :param tran,rot,scale: boolean telling that translation/rotation/scale will be included
    :param equilibrate: (default True) changes the unit of the matrix so that x are in km and scale in mas for better conditionning of the matrix
    
    :returns H: the transformation matrix as 2D numpy matrix
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
    """
    Generate the linear transformation equation system for coordinate of a sites in two reference system.
    The equation system is generated for E,N,U components
    
    :param xf,yf,zf: geocentric coordinates in the initial reference (also called 'free' frame) 
    :param xr,yr,zr: geocentric coordinates in the final reference (the reference frame) 
    :param tran,rot,scale: boolean telling that translation/rotation/scale will be included

    :return A,B: the design matrix and observation vector
    
    :note: VCV not accounted yet 
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
    
    """
    Estimates Helmert or transformation parameters between two sets of coordinates
    Default is a 7-parameters transformation (translation, rotation & scale)
    
    
    :param H_Gpoint_ref : dictionary of Geodetic Points (Gpoint instance) used as reference
    :param H_Gpoint_free: dictionary of Geodetic Points (Gpoint instance) used as free
    :param psd: psd object from sinex library for post-seismic deformation
    :param tran,rot,scale: boolean telling if the transformation will use translation/rotation/scale 
    :param method: 'LS' (Least-Squares) or 'Dikin_LS' (a sequence of L1 and L2 norm) estimation
    :param threshold: threshold value for outlier rejection (default 5.)
    :param verbose: boolean
    
    """
    
    import numpy as np
    import pyacs.lib.glinalg
    import pyacs.lib.robustestimators as RobustEstimators 

    if verbose:print("-- Helmert parameters estimation using method: ",method)

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
            if np.any(denh) and verbose:
                print('-- psd correction: ',code,denh * 1000. )
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
        if verbose:
            print("-- Doing Dikin estimation")
        
        
        X,V=RobustEstimators.Dikin(A,B,W=None)

        X=X.reshape(-1,1)
        V=V.reshape(-1,1)

        [median_east,median_north,median_up]=np.median(np.sqrt((V.reshape(-1,3)*1000.0)**2),axis=0)
    
    ### reject the outliers et re-adjust by LS
    
        if verbose:
            print("-- Testing outliers; Threshold for outliers detection: ",threshold)
        
        np_outliers = find_outliers_Helmert(V,threshold=threshold)

        l_np_outliers=np_outliers.tolist()

        if verbose:
            print("-- # of outliers detected: ",len(l_np_outliers))

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
                raise PyacsSol_HelmertError(inspect.stack()[0][3],__name__,msg)

            if verbose:
                print("-- Doing LS estimation without outliers")
            
            LS_X=pyacs.lib.glinalg.ls(A_rejected,B_rejected)
            LS_V = B-np.dot(A,X)

                
    else:
    ### LS estimation, no outlier detection

        if verbose:print("-- Doing LS estimation without outlier detection")

        LS_X=pyacs.lib.glinalg.ls(A,B)
        V = B-np.dot(A,X)
        l_np_outliers=[]

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
        np_outliers = np.vstack({tuple(row) for row in np_outliers })
    
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

