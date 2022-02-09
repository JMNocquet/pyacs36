"""
Coordinates library.

-- Local/Geocentric frame conversion

-- Geodetic/Geocentric frame conversion

-- Spherical/Geocentric conversion

-- Geodetic/Flat Earth conversion


"""

from pyacs.lib.errors import *

###################################################################
def mat_rot_general_to_local(lam,phi,unit='radians'):
###################################################################
    
    """ 
    Generates the conversion matrix R from general geocentric cartesian coordinates (XYZ) to local cartesian coordinates (ENU):
    
    :param lam,phi: longitude
    :param unit: 'radians' or 'dec_deg' (default is 'radians')
    
    :returns: R as a 2D numpy array
    
    
    :note: R = [[-sin(lambda)            cos(lambda)         0       ],
    
    [-sin(phi)*cos(lambda) , -sin(phi)*sin(lamda) , cos(phi)],
    
    [ cos(phi)*cos(lambda) , cos(phi)*sin(lamda) , sin(phi)]]
    """

    #print("unit:",unit)


    if unit not in ['radians','dec_deg']:
        raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)

    import numpy as np


    if unit == 'dec_deg':
        lam = np.radians(lam)
        phi = np.radians(phi)
    
    R = np.zeros([3,3], float)
    R[0,0] = -np.sin(lam)
    R[0,1] = np.cos(lam)

    R[1,0] = -np.sin(phi)*np.cos(lam)
    R[1,1] = -np.sin(phi)*np.sin(lam)
    R[1,2] = np.cos(phi)
    
    R[2,0] = np.cos(phi)*np.cos(lam)
    R[2,1] = np.cos(phi)*np.sin(lam)
    R[2,2] = np.sin(phi)

    return R


###################################################################
def mat_rot_local_to_general(lam,phi,unit='radians'):
###################################################################
    """
    Generates the conversion matrix R from local cartesian coordinates (ENU) to general geocentric cartesian coordinates (XYZ):

    :param lam,phi: longitude in radians
    :param unit: 'radians' or 'dec_deg' (default is 'radians')

    :returns: R as a 2D numpy array
    
    :note: Since R is orthogonal, it is the inverse and also the transpose of the conversion matrix from general geocentric cartesian coordinates (XYZ) to local cartesian coordinates (ENU)
    """
    if unit not in ['radians','dec_deg']:
        raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)
    
    return(mat_rot_general_to_local(lam,phi,unit=unit).T)


###################################################################
def denu_at_x0y0z0_to_xyz(de,dn,du,x0,y0,z0):
###################################################################
    """
    Converts a local vector [dn,de,du] at [x0,y0,z0] in geocentric cartesian coordinates
    into X,Y,Z geocentric cartesian coordinates
    
    :param de,dn,du: east,north, up components of the local vector in meters
    :param x0,x0,z0: reference point for the local vector in geocentric cartesian coordinates
    
    :returns: x,y,z in geocentric cartesian coordinates
    
    """
    
    import numpy as np
    
    (lam,phi,_)=xyz2geo(x0,y0,z0)
    R=mat_rot_local_to_general(lam,phi)
    
    DL=np.array([de,dn,du])
    X0=np.array([x0,y0,z0])
    XYZ=np.dot(R,DL)+X0
    return(XYZ[0],XYZ[1],XYZ[2])


###################################################################
def xyz2geo(x,y,z,A=6378137.,E2 = 0.006694380022903, unit='radians'):
###################################################################
    """
    Converts geocentric cartesian coordinates (XYZ) to geodetic coordinates (lambda,phi,h)
    
    :param x,y,z: XYZ coordinates in meters
    :param unit: 'radians' or 'dec_deg' (default is 'radians')
    
    :returns: long,lat,he: longitude and latitude in radians, he in meters above the ellipsoid

    :note: Default ellipsoid is GRS80 used for WGS84 with
    
    A  = 6378137.         # semi major axis = equatorial radius
    
    E2 = 0.00669438003    # eccentricity and then
    
    F = 1.0 - sqrt(1-E2)  # flattening

    :ref: Bowring, 1985, The accuracy of geodetic latitude and height equations, Survey Review, 28, 202-206.

    """    
    if unit not in ['radians','dec_deg']:
        raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)

    import numpy as np
    
    F=1.0 - np.sqrt(1-E2)
    
    TP = np.sqrt (x**2 + y**2)
    R = np.sqrt (TP**2 + z**2)
 
    TMU = np.arctan2( (z/TP)*((1. - F) + (E2*A)/R) ,1)
    RLBDA = np.arctan2(y,x)

    S3 = (np.sin(TMU))**3
    C3 = (np.cos(TMU))**3
    T1 = z*(1 - F) + E2*A*S3
    T2 = (1 - F)*( TP - E2*A*C3 )


    RPHI = np.arctan2(T1,T2)
    RHE = TP*(np.cos(RPHI)) + z*(np.sin(RPHI))
    RHE = RHE - A*( np.sqrt(1-E2*(np.sin(RPHI)**2)) )

    if unit == 'dec_deg':
        RLBDA = np.degrees(RLBDA)
        RPHI = np.degrees(RPHI)

    return(RLBDA,RPHI,RHE)



###################################################################
def geo2xyz(llambda,phi,he,unit='radians',A=6378137.,E2 = 0.006694380022903):
###################################################################
    """
    Converts geodetic coordinates (long,lat,he) to geocentric cartesian coordinates (XYZ)

    :param llambda,phi: longitude, latitude
    :param he: ellipsoidal height in meters
    :param unit: 'radians' or 'dec_deg' for llambda and phi (default is 'radians')

    :returns: x,y,z in meters
    
    :note: Default ellipsoid is GRS80 used for WGS84 with
    
    A  = 6378137.         # semi major axis = equatorial radius
    
    E2 = 0.00669438003    # eccentricity and then
    
    F = 1.0 - sqrt(1-E2)  # flattening

    """

    if unit not in ['radians','dec_deg']:
        raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)

    import numpy as np
      
    if unit=='dec_deg':
        llambda=np.radians(llambda)
        phi=np.radians(phi)

    a=A
    e2=E2

    wn=wnorm(phi, A=a, E2=e2)
   
    xx=(wn+he)*np.cos(phi)*np.cos(llambda)
    yy=(wn+he)*np.cos(phi)*np.sin(llambda)
    zz=(wn*(1.0-e2)+he)*np.sin(phi)
   
    return(xx,yy,zz)     


###################################################################
def wnorm(phi,A=6378137.,E2 = 0.006694380022903):
###################################################################
    """
    Calculates the geodetic radius normal to the ellipsoid
    """

    import numpy as np

    a=A
    e2=E2
  
    wd=np.sqrt(1.-e2*np.sin(phi)**2) 
    result=a/wd
    
    return (result)  

###################################################################
def xyz2geospheric(x,y,z,unit='radians'):
###################################################################
    """
    Converts geocentric cartesian coordinates (XYZ) to geo-spherical coordinates (longitude,latitude,radius).
    
    :param x,y,z: geocentric cartesian coordinates in meters
    :param unit: 'radians' or 'dec_deg' (default is 'radians')
    
    :returns: longitude,latitude (in radians), radius from the Earth's center in meters
    
    :note: Be aware that the obtained coordinates are not what is usually taken as spherical coordinates, which uses co-latitude
    """

    if unit not in ['radians','dec_deg']:
        raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)
    
    import numpy as np
    
    R=np.sqrt(x**2+y**2+z**2)
    RLBDA = np.arctan2(y,x)
    
    Req=np.sqrt(x**2+y**2)
    if Req !=0.:
        RPHI = np.arctan(z/Req)
    else:
        RPHI=np.pi/2.*z/np.sqrt(z**2)

    if unit == 'dec_deg':
        RLBDA = np.degrees(RLBDA)
        RPHI = np.degrees(RPHI)
    
    return(RLBDA,RPHI,R)

###################################################################
def xyz_spherical_distance(x1,y1,z1,x2,y2,z2,Rt=6.371E6):
###################################################################
    """
    Returns the spherical distance between two points of XYZ coordinates x1,y1,z1 and x2,y2,z2
     
    :param x1,y1,z1,x2,y2,z2: in meters
    :param Rt: mean Earth's radius in meters (default 6 371 000.0 m)
 
    :returns: distance in meters
     
    """
 
    import pyacs.lib.vectors
    import numpy as np
     
    sp = pyacs.lib.vectors.scal_prod(pyacs.lib.vectors.normalize([x1,y1,z1]), pyacs.lib.vectors.normalize([x2,y2,z2]))
    sp = np.around(sp,decimals=13)
    return np.fabs( np.arccos( sp ) * Rt )

###################################################################
def geo_spherical_distance(lon1,lat1,h1,lon2,lat2,h2,unit='radians',Rt=6.371E6):
###################################################################
    """
    Returns the spherical distance between two points of geodetic coordinates lon1,lat1,lon2,lat2
     
    :param lon1,lat1,h1,lon2,lat2,h2: longitude, latitude, radius
    :param unit: 'radians' or 'dec_deg' for lon & lat (default is 'radians')
    :param Rt: mean Earth's radius in meters (default 6 371 000.0 m)
 
    :returns: distance in meters
     
    """
    if unit not in ['radians','dec_deg']:
       raise OptionError("unit option must be in [radians,dec_deg];unit=",unit)

    (x1,y1,z1)=geo2xyz(lon1, lat1, 0.0, unit=unit)
    (x2,y2,z2)=geo2xyz(lon2, lat2, 0.0, unit=unit)
 
 
    return xyz_spherical_distance(x1,y1,z1,x2,y2,z2)


###############################################################################
def azimuth_to_en(azimuth):
###############################################################################
    """
    converts an azimuth (clockwise from north) to east and north unit vector

    :param azimuth: azimuth in decimal degrees
    :return east,north: unit vector  
    :note: calculation is made using the sphere
    """

    import numpy as np

    r_az=np.radians(azimuth)
    return(np.sin(r_az),np.cos(r_az))


###################################################################
def geo2flat_earth(longitude,latitude):
##################################################################
    """
    Converts geographical coordinates to flat earth approximation
    uses pyproj and web Mercator projection.
    
    :param longitude,latitude: geographical (ellipsoid) coordinates in decimal degrees. Also works with 1D numpy arrays
    
    :returns x,y: in km
    """
    # import
    from pyproj import Transformer
    import numpy as np

    longitude=np.array(longitude)
    latitude = np.array(latitude)

    TRAN_4326_TO_3857 = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    (x,y) = TRAN_4326_TO_3857.transform( longitude , latitude )
    return x * 1.E-3, y * 1.E-3


###################################################################
def flat_earth2geo(x, y):
##################################################################
    """
    Converts web Mercator coordinates (units km) to geographical coordinates
    uses pyproj and web Mercator projection.

    :param longitude,latitude: geographical (ellipsoid) coordinates in decimal degrees. Also works with 1D numpy arrays

    :returns x,y: in km
    """
    # import
    import numpy as np
    from pyproj import Transformer

    x = np.array(x)
    y = np.array(y)

    TRAN_3857_TO_4326 = Transformer.from_crs("EPSG:3857","EPSG:4326", always_xy=True )
    return TRAN_3857_TO_4326.transform(x * 1.E3 , y *1.E3,  )


###################################################################
def spherical_baseline_length_rate(slon, slat, sve, svn, elon, elat, eve, evn, sigma=None, verbose=False):
###################################################################
    """
    calculates the baseline (great circle) length rate of change between two points
    
    :param slon,slat: coordinates of profile start point (decimal degrees) 
    :param sve,svn: east and north components of start point (mm/yr) 
    :param elon,elat: coordinates of profile end   point (decimal degrees)
    :param eve,evn: east and north components of end point (mm/yr) 
    :param sigma: list of velocities uncertainties: [sig_sve, sig_svn, corr_sven, sig_eve, sig_evn, corr_even] 
    
    :return         :  length, rate, 1D strain rate
    """
    # Rt in meters
    
    Rt = 6371.0E3
    
    # 
    #print('slon, slat, sve, svn, elon, elat, eve, evn',slon, slat, sve, svn, elon, elat, eve, evn)
    #print('sigma ', sigma )
    
    # import
    
    import numpy as np
    from pyacs.lib import vectors
    from pyacs.lib.gmtpoint import GMT_Point
    from pyacs.lib import glinalg
    from pyacs.lib import euler

    (xi, yi, zi) = geo2xyz(np.radians(slon), np.radians(slat), 0.0 , unit='radians')
    Xi = np.array([xi, yi, zi])
    #print('Xi' , Xi)

    (xs, ys, zs) = geo2xyz(np.radians(elon), np.radians(elat), 0.0 , unit='radians')
    Xs = np.array([xs, ys, zs])
    #print('Xs' , Xs)

    POLE = vectors.vector_product(Xi, Xs)
    #print( euler.rot2euler(POLE[0], POLE[1], POLE[2]) )
    POLE = vectors.normalize(POLE)

    # along / transverse great circle component start point
    R = mat_rot_general_to_local(np.radians(slon), np.radians(slat), unit='radians')
    OM = np.array([xi, yi, zi])
    OM = vectors.normalize(OM)

    unit_parallele_xyz = vectors.vector_product(POLE, OM)
    unit_parallele_enu = np.dot(R, unit_parallele_xyz)
    #print('unit_parallele_enu' , unit_parallele_enu*10. )

    sv_parallele = sve * unit_parallele_enu[0] + svn * unit_parallele_enu[1]
    sv_perpendicular = sve * unit_parallele_enu[1] - svn * unit_parallele_enu[0]

    # covariance
    if sigma is not None:
        sig = np.array([ sigma[0] , sigma[1] ])
        corr = np.eye(2)
        corr[0,1] = sigma[2]
        corr[1,0] = sigma[2]
        cov = glinalg.corr_to_cov(corr, sig )
        # rotation
        angle_rad = np.arctan2( unit_parallele_enu[1] , unit_parallele_enu[2] )
        #print('angle_rad in deg ' , np.degrees( angle_rad ))
        rot = np.zeros((2,2))
        rot[0,0] = np.cos( angle_rad )
        rot[0,0] = np.cos( angle_rad )
        rot[1,1] = np.cos( angle_rad )
        rot[0,1] = -np.sin( angle_rad )
        rot[1,0] = np.sin( angle_rad )
        # apply rotation
        vcv_en = np.dot( rot , np.dot( cov , rot.T) )
        
        s_sv_parallel = np.sqrt( vcv_en[0,0] )

    # along / transverse great circle component end point
    R = mat_rot_general_to_local(np.radians(elon), np.radians(elat), unit='radians')
    OM = np.array([xs, ys, zs])
    OM = vectors.normalize(OM)

    unit_parallele_xyz = vectors.vector_product(POLE, OM)
    unit_parallele_enu = np.dot(R, unit_parallele_xyz)
    #print('unit_parallele_enu' , unit_parallele_enu*10. )
    ev_parallele = eve * unit_parallele_enu[0] + evn * unit_parallele_enu[1]
    ev_perpendicular = eve * unit_parallele_enu[1] - evn * unit_parallele_enu[0]

    # covariance
    if sigma is not None:
        sig = np.array([ sigma[3] , sigma[4] ])
        corr = np.eye(2)
        corr[0,1] = sigma[5]
        corr[1,0] = sigma[5]
        cov = glinalg.corr_to_cov(corr, sig )
        # rotation
        angle_rad = np.arctan2( unit_parallele_enu[1] , unit_parallele_enu[2] )
        #print('angle_rad in deg ' , np.degrees( angle_rad ))
        rot = np.zeros((2,2))
        rot[0,0] = np.cos( angle_rad )
        rot[0,0] = np.cos( angle_rad )
        rot[1,1] = np.cos( angle_rad )
        rot[0,1] = -np.sin( angle_rad )
        rot[1,0] = np.sin( angle_rad )
        # apply rotation
        vcv_en = np.dot( rot , np.dot( cov , rot.T) )
        
        s_ev_parallel = np.sqrt( vcv_en[0,0] )
        
        
    # great circle arc length
    MS = GMT_Point(lon=slon, lat=slat)
    ME = GMT_Point(lon=elon, lat=elat)
    length_arc = MS.spherical_distance(ME)

    # lengthening in mm/yr
    length_rate = ev_parallele - sv_parallele

    # 1D strain rate in 1/yr
    strain_rate = length_rate * 1E-3 / length_arc
    
    # case covariance
    if sigma is None:
        return  length_arc , length_rate , strain_rate
    else:
        sig_v = np.sqrt( s_sv_parallel**2 + s_ev_parallel**2 )
        sig_strain_rate = sig_v * 1.E-3 / length_arc
        return  length_arc , length_rate , strain_rate,  sig_v, sig_strain_rate







