
###############################################################################
def v_to_rake(ve,vn,strike,dip,style='inverse'):
###############################################################################

    """
    for a given relative horizontal velocity (ve,vn) between two blocks separated by a fault of a given strike & dip
    , returns the expected rake for the eq.
    
    :param ve,vn: east and north components of motion (any unit)
    :param strike,dip: in decimal degrees
    :param style: string among {'inverse', 'normal', 'leftlateral','rightlateral'}
    :returns rake: in decimal degrees 
        
    :note: Because only providing ve,vn is ambiguous (we don't know whether if ve,vn is hanging wall motion wrt the footwall or the opposite) \
    an additional information (style) in the form of one of {'inverse', 'normal', 'leftlateral','rightlateral'} must be provided.
    """
    
    import numpy as np

    if dip == 90.:
        import sys
        print("!!! ERROR: dip is 90 degrees. Cannot calculate rake")
        sys.exit()
        
        

    r_dip=np.radians(dip)
    ralpha=np.radians(90.-strike)
    v=np.array([ve,vn])
    v_norm=v/np.sqrt(v[0]**2+v[1]**2)
    
    R=np.array([[np.cos(ralpha),-np.sin(ralpha)],[np.sin(ralpha),np.cos(ralpha)]])
    
    VRAKE=np.dot(R.T,v_norm)
    VRAKE[1]=VRAKE[1]/np.cos(r_dip)
    
    rake_1=np.degrees(np.arctan2(VRAKE[1],VRAKE[0]))
    
    rake_2=rake_1+180.0
    if rake_2>180:rake_2=rake_2-360.
    
    rakes=np.array([rake_1,rake_2])
    
    if style is not None:
        
        (rake_negative,rake_positive)=(np.min(rakes),np.max(rakes))
        if np.abs(rake_1)<=90:
            rake_leftlateral=rake_1;rake_rightlateral=rake_2 
        else:
            rake_leftlateral=rake_2;rake_rightlateral=rake_1 
        
        if style=='inverse':rake=rake_positive
        if style=='normal':rake=rake_negative
        if style=='leftlateral':rake=rake_leftlateral
        if style=='rightlateral':rake=rake_rightlateral
    
    else:
        rake = rake_1

    return (rake)

###############################################################################
def v_to_n_ss(ve,vn,strike):
###############################################################################
    """
    
    for a given relative horizontal velocity (ve,vn) between two blocks separated by a fault of a given strike
    , returns normal and strike-slip component of motion.
    The convention is that the point is in the left-domain with respect to the fault.
    Shortening & right-lateral are positive, extension and left-lateral negative
    
    :param ve,vn: east and north components of motion (any unit)
    :param strike: in decimal degrees
    :returns normal,strike-slip: same units as ve,vn 
    
    """

    import numpy as np
    
    azimuth_radian=np.arctan2(ve,vn)
    radian_strike=np.radians(strike)
    
    angle=azimuth_radian-radian_strike
    
    v_i=np.sqrt(ve**2+vn**2)
    ss=np.cos(angle)*v_i
    normal=np.sin(angle)*v_i

    return (normal,ss)

###############################################################################
def strike_dip_rake_to_dir(strike,dip,rake):
###############################################################################
    """
    for a given (strike,dip,rake)  returns the azimuth of the horizontal motion
    :param strike,dip,rake: in decimal degrees
    :returns: direction: slip direction (modulo 180 degrees)
    """


    import numpy as np
    r_rake=np.radians(rake)
    r_dip=np.radians(dip)
    dir_slip=strike-np.degrees(np.arctan(np.tan(r_rake)*np.cos(r_dip)))

    if dir_slip > 180.:dir_slip=dir_slip-180.0
    if dir_slip > 180.:dir_slip=dir_slip-180.0
    if dir_slip<0:dir_slip=dir_slip+180.0

    return(dir_slip)

###################################################################
def rake_from_euler(longitude,latitude,strike, dip, euler):
###################################################################
    """
    predicts rake for a given fault according to an Euler pole, position and strike of the fault

    :param longitude,latitude: in decimal degrees
    :param strike: fault strike from north in decimal degrees
    :param dip: fault dip from north in decimal degrees
    :param euler: Euler pole as a string '/long/lat/w/style' (style among 'inverse', 'normal', 'leftlateral','rightlateral') 
    
    :return rake: in decimal degrees
    :note: style needs to be provided to ensure the correct sense of slip.
    """

    import pyacs.lib.gmtpoint
    import numpy as np
    import pyacs.lib.faultslip

    M=pyacs.lib.gmtpoint.GMT_Point(code='XXXX',lon=longitude,lat=latitude,he=0., Ve=0.,Vn=0.,Vu=0.,\
                  SVe=0,SVn=0,SVu=0,SVen=0,\
                  Cv_xyz=None, Cv_enu=None, index=None)

    _tmp,elon,elat,ew,motion_type=euler.split('/')

    W=np.array( list(map(float, [elon,elat,ew])))

    N=M.pole(W=W,type_euler='euler',option='predict')
    
    (ve,vn)=(N.Ve,N.Vn)
    
    rake = v_to_rake(ve,vn,strike,dip,motion_type)
    
    return(rake)

###################################################################
def slip_rake_2_ds_ss(slip,rake):
###################################################################
    """
    Converts a slip vector provided as slip and rake to dip slip and strike-slip components

    :param slip: slip
    :param rake: rake in degrees

    :return ( slip*np.sin( np.radians(rake) ) ,  slip*np.cos( np.radians(rake) ))

    :note:
    ds will be positive for rake in [0,180] that is reverse motion
    ss will be positive for left-lateral motion
    
    """
    
    import numpy as np

    return ( slip*np.sin( np.radians(rake) ) ,  slip*np.cos( np.radians(rake) ))

###################################################################
def ss_ns_2_ve_vn(ss,ns,strike):
###################################################################
    """
    Converts strike-slip and normal slip components of a fault to ve, vn

    :param ss: vector of strike-slip
    :param ns: vector of normal-slip
    :param strike: vector of fault strike, counter-clockwise from north in degrees

    :return: ve, vn, same unit as ss and ds

    :note:
    ss is assume to be positive for left-lateral motion
    ds will be positive for reverse motion
    """
    
    import numpy as np

    # strike in radians
    rstrike = np.radians( strike )
    # rcos and rsin
    rcos = np.cos( rstrike )
    rsin = np.sin( rstrike )

    return ss*rsin - ns*rcos , ss*rcos + ns*rsin 


###############################################################################
def geo_to_strike(ilon,ilat,elon,elat):
###############################################################################

    """
    for a given fault segment starting at ilon,lat and ending at elon, elat
    , returns the strike.
    
    :param ilon,ilat: geographical coordinates of fault segment start point in decimal degrees
    :param elon,elat: geographical coordinates of fault segment end point in decimal degrees
    :returns strike: in decimal degrees clockwise from north
    
    :note: strike here is taken as the initial bearing. Be cautious with long segments
    """
    
    import numpy as np
    
    # radian conversion
    [lon1,lat1,lon2,lat2] = np.radians( [ ilon,ilat,elon,elat ] ) 
    
    
    strike_rad_from_east_ccwise = np.arctan2(np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1) , \
                        np.sin(lon2-lon1)*np.cos(lat2)) 

    strike_deg_from_north_cwise = 90. - np.degrees( strike_rad_from_east_ccwise )
    
    if strike_deg_from_north_cwise > 180.:
        strike_deg_from_north_cwise = strike_deg_from_north_cwise - 360.

    return( strike_deg_from_north_cwise )
