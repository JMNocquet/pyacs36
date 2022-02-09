"""

Euler poles manipulation

"""

###############################################################################
def rot2euler(Wx,Wy,Wz):
###############################################################################
    """Converts a rotation rate vector Wx Wy Wz in radians/yr into 
    an Euler pole (long. lat. deg/Myr) into geographical (spherical) coordinates and angular velocity
 
    :param Wx,Wy,Wz: rotation rate vector in geocentric cartesian coordinates with units of radians/yr
    :returns longitude,latitude,omega: longitude and latitude of the Euler pole in radians and the
    angular velocity in decimal degrees per Myr.
    
    :note: longitude and latitude are relative to the sphere, not the ellipsoid. This is because 
    Euler pole and rigid rotations only have sense on a sphere. 
    
    """  

    import numpy as np
    
    W = np.sqrt(Wx**2+Wy**2+Wz**2)
    W_lat = 90.0-np.arccos(Wz/W)*180.0/np.pi

    if (Wx>0):
        W_long = np.arctan(Wy/Wx)*180.0/np.pi
    elif (Wx<0):
        if (Wy>0): 
            W_long = np.arctan(Wy/Wx)*180.0/np.pi+180.0
        else:
            W_long = np.arctan(Wy/Wx)*180.0/np.pi-180.0
    else:
        # Wx = 0
        if (Wy>0): 
            W_long = 90.
        else:
            W_long = -90.
        
    
    W_omega=W/np.pi*180*1.0E6

    return(W_long,W_lat,W_omega)

###############################################################################
def euler2rot(lon, lat, omega):
###############################################################################
    """Converts Euler poles (long., lat., deg/Myr) into cartesian geocentric 
       rotation rate vector Wx Wy Wz in radians/yr
       
    :parameters longitude,latitude,omega: longitude and latitude of the Euler pole in decimal degrees and the
    angular velocity in decimal degrees per Myr.
    :returns Wx Wy Wz: in radians/yr 
   
    :note: longitude and latitude are relative to the sphere, not the ellipsoid.
    """

    import numpy as np

    wx=np.cos(np.radians(lat))*np.cos(np.radians(lon))*np.radians(omega)*1.0E-6
    wy=np.cos(np.radians(lat))*np.sin(np.radians(lon))*np.radians(omega)*1.0E-6
    wz=np.sin(np.radians(lat))*np.radians(omega)*1.0E-6

    return (wx,wy,wz)

###############################################################################
def euler_uncertainty( w, vcv ):
###############################################################################
    """
    Calculates Euler pole parameters uncertainty
    
    :param w: rotation vector in XYZ coordinates as numpy 1D array
    :param vcv: covariance of w as numpy 2D array
    
    :return: vcv_euler as numpy 2D array
    """

    import pyacs.lib.coordinates as Coordinates
    import numpy as np

    nw = np.sqrt( w[0]**2 + w[1]**2 + w[2]**2 )
        
    (llambda,phi,omega) = rot2euler( w[0],w[1],w[2] )

    # Propagates vcv into the local frame to get the Euler pole uncertainties
    
    R=Coordinates.mat_rot_general_to_local(np.radians(llambda),np.radians(phi))

    VCV_POLE_ENU=np.dot(np.dot(R, vcv),R.T)

    vcv11=VCV_POLE_ENU[0,0]
    vcv12=VCV_POLE_ENU[0,1]
    vcv22=VCV_POLE_ENU[1,1]
    vcv33=VCV_POLE_ENU[2,2]

    # Error ellipse parameters
    
    max_sigma=0.5*(vcv11+vcv22+np.sqrt((vcv11-vcv22)**2+4*vcv12**2))
    max_sigma=np.degrees(np.arctan(np.sqrt(max_sigma)/ nw ))
    min_sigma=0.5*(vcv11+vcv22-np.sqrt((vcv11-vcv22)**2+4*vcv12**2))
    min_sigma=np.degrees(np.arctan(np.sqrt(min_sigma)/ nw ))
    azimuth=np.degrees(2*np.arctan(2*vcv12/(vcv11-vcv22)))
    s_omega=np.degrees(np.sqrt(vcv33))*1.E06

    return( max_sigma, min_sigma, azimuth, s_omega )


###############################################################################
def vel_from_euler(lonp,latp,lon_euler,lat_euler,omega_euler):
###############################################################################
    """
    Return the horizontal velocity predicted at lonp,latp from an Euler pole
    
    :param lonp,latp: longitude,latitude in decimal degrees where velocity will be predicted
    :param lon_euler,lat_euler_omega_euler: longitude and latitude of the Euler pole in decimal degrees and the
    angular velocity in decimal degrees per Myr.
    
    :return: ve,vn in mm/yr
    """

    import numpy as np
    import pyacs.lib.coordinates
    
    # convert to rotation rate vector
    (wx,wy,wz) = euler2rot(lon_euler, lat_euler, omega_euler)
    
    # W
    W=np.zeros((3,1))
    W[0,0]=wx
    W[1,0]=wy
    W[2,0]=wz


    # get spherical coordinates

    he = 0.
    (x,y,z)=pyacs.lib.coordinates.geo2xyz( lonp , latp , he , unit='dec_deg')
    (l,p,_)=pyacs.lib.coordinates.xyz2geospheric(x,y,z,unit='radians')
    R=pyacs.lib.coordinates.mat_rot_general_to_local(l,p)

    # Observation equation in local frame
            
    Ai=np.zeros([3,3],float)
    
    Ai[0,1]= z
    Ai[0,2]=-y
    Ai[1,0]=-z
    Ai[1,2]= x
    Ai[2,0]= y
    Ai[2,1]=-x
            
    RAi=np.dot(R,Ai)
        
    # Predicted velocity
            
    V = np.dot(RAi,W) * 1.E3 # to get mm/yr
    
    return( V[0,0] , V[1,0])
    
###############################################################################
def pole_matrix( coor ):
###############################################################################
    """
    Calculates the matrix relating the horizontal velocity to a rotation rate vector.
    Given a 2D-numpy array of n positions [lonp , latp] in decimal degrees the return matrix is W so that
    np.dot( W , w ) gives a 2D-numpy array of [ve1,vn1,ve2,vn2,....] expressed in m/yr for a rotation rate vector in rad/yr 
    
    :param coor: 2D numpy array of [lon, lat] in decimal degrees
    
    :return: the pole matrix as a 2D-numpy array
    """
    
    # import 
    import numpy as np
    import pyacs.lib.coordinates

    # ensures dimension
    
    coor = coor.reshape(-1,2)

    # pole_matrix

    pole_matrix = None

    # loop if points
    
    for i in np.arange( coor.shape[0] ):
        
        # get spherical coordinates
    
        he = 0.
        (x,y,z) = pyacs.lib.coordinates.geo2xyz( coor[i,0] , coor[i,1] , he , unit='dec_deg')
        (l,p,_) = pyacs.lib.coordinates.xyz2geospheric(x,y,z,unit='radians')
        R       = pyacs.lib.coordinates.mat_rot_general_to_local(l,p)
    
        # Observation equation in local frame
                
        Ai=np.zeros([3,3],float)
        
        Ai[0,1]= z
        Ai[0,2]=-y
        Ai[1,0]=-z
        Ai[1,2]= x
        Ai[2,0]= y
        Ai[2,1]=-x
                
        RAi=np.dot(R,Ai)[:2,:]
        
        # Predicted velocity
                
        if pole_matrix is None:
            pole_matrix = RAi
        else:
            pole_matrix = np.vstack(( pole_matrix , RAi ))
        

    return pole_matrix

###############################################################################
def pole_matrix_fault( coor , strike , order=None ):
###############################################################################
    """
    Calculates the matrix relating the along strike and normal slip components of a fault to a rotation rate vector.
    Given a 2D-numpy array of n positions [lonp , latp] in decimal degrees and strike counter-clockwise from north the return matrix is W so that
    np.dot( W , w ) gives a 2D-numpy array of [ss1,ns1,ss2,ns2,....] expressed in m/yr for a rotation rate vector in rad/yr 
    
    :param coor  : 2D numpy array of [lon, lat] in decimal degrees
    :param strike: 1D numpy array of [strike] in decimal degrees
    
    :return: the pole matrix as a 2D-numpy array
    :note: np.dot( W , w ).reshape(-1,2) gives the along strike,and normal components in two columns
    :note: the value are given for the hanging-wall block (right-polygon)  
    """
    
    # import 
    import numpy as np
    import pyacs.lib.coordinates

    # ensures dimension
    
    coor = coor.reshape(-1,2)

    # strike in radians
    
    #strike = strike.reshape(-1,1) 
    rstrike = np.radians( strike )

    if order is not None:
        pmf_SS = np.zeros( ( strike.shape[0] , 3  ))
        pmf_DS = np.zeros( ( strike.shape[0] , 3  ))

    # pole_matrix

    pole_matrix = None

    # loop if points
    
    for i in np.arange( coor.shape[0] ):
        
        # get spherical coordinates
    
        he = 0.
        (x,y,z)=pyacs.lib.coordinates.geo2xyz( coor[i,0] , coor[i,1] , he , unit='dec_deg')
        (l,p,_)=pyacs.lib.coordinates.xyz2geospheric(x,y,z,unit='radians')
        R=pyacs.lib.coordinates.mat_rot_general_to_local(l,p)
    
        # Observation equation in local frame
                
        Ai=np.zeros([3,3],float)
        
        Ai[0,1]= z
        Ai[0,2]=-y
        Ai[1,0]=-z
        Ai[1,2]= x
        Ai[2,0]= y
        Ai[2,1]=-x
                
        RAi=np.dot(R,Ai)[:2,:]
        
        # rotates into the fault coordinates

        RF = np.array( [ [  np.sin( rstrike[i] )  ,  np.cos( rstrike[i] ) ] , \
                         [ -np.cos( rstrike[i] )  ,  np.sin( rstrike[i] ) ] ] )

        # RW
        RW = np.dot( RF , RAi  )
        
        # Predicted velocity

        if order is not None:
            pmf_SS[i] = RW[0,:]
            pmf_DS[i] = RW[1,:]

                
        if pole_matrix is None:
            pole_matrix = RW
        else:
            pole_matrix = np.vstack(( pole_matrix , RW ))
    
    if order is not None:
        return np.vstack( (pmf_SS,pmf_DS) )
    else:
        return pole_matrix


