

###################################################################
def spherical_baseline_length_rate(slon, slat, sve, svn, elon, elat, eve, evn, sigma=None, verbose=False):
###################################################################
    """Baseline (great circle) length rate between two points.

    Parameters
    ----------
    slon : float
        Longitude of profile start point, decimal degrees.
    slat : float
        Latitude of profile start point, decimal degrees.
    sve : float
        East velocity at start point, mm/yr.
    svn : float
        North velocity at start point, mm/yr.
    elon : float
        Longitude of profile end point, decimal degrees.
    elat : float
        Latitude of profile end point, decimal degrees.
    eve : float
        East velocity at end point, mm/yr.
    evn : float
        North velocity at end point, mm/yr.
    sigma : array_like, optional
        Velocity uncertainties [sig_sve, sig_svn, corr_sven, sig_eve, sig_evn, corr_even].
    verbose : bool, optional
        If True, enable verbose output.

    Returns
    -------
    length_arc : float
        Great circle arc length in meters.
    length_rate : float
        Lengthening rate in mm/yr.
    strain_rate : float
        1D strain rate in 1/yr.
    If sigma is provided, also returns sig_v, sig_strain_rate.
    """
    # Rt in meters

    import numpy as np
    from pyacs.lib.gmtpoint import GMT_Point
    from pyacs.lib import glinalg
    from pyacs.lib import euler
    from .mat_rot_general_to_local import mat_rot_general_to_local
    from .geo2xyz import geo2xyz
    
    Rt = 6371.0E3
    
    # 
    #print('slon, slat, sve, svn, elon, elat, eve, evn',slon, slat, sve, svn, elon, elat, eve, evn)
    #print('sigma ', sigma )
    
    (xi, yi, zi) = geo2xyz(np.radians(slon), np.radians(slat), 0.0, unit="radians")
    Xi = np.array([xi, yi, zi])
    #print('Xi' , Xi)

    (xs, ys, zs) = geo2xyz(np.radians(elon), np.radians(elat), 0.0, unit="radians")
    Xs = np.array([xs, ys, zs])
    #print('Xs' , Xs)

    POLE = np.cross(Xi, Xs)
    #print( euler.rot2euler(POLE[0], POLE[1], POLE[2]) )
    POLE = POLE / np.linalg.norm(POLE)

    # along / transverse great circle component start point
    R = mat_rot_general_to_local(np.radians(slon), np.radians(slat), unit="radians")
    OM = np.array([xi, yi, zi])
    OM = OM / np.linalg.norm(OM)

    unit_parallele_xyz = np.cross(POLE, OM)
    unit_parallele_enu = np.dot(R, unit_parallele_xyz)
    #print('unit_parallele_enu' , unit_parallele_enu*10. )

    sv_parallele = sve * unit_parallele_enu[0] + svn * unit_parallele_enu[1]
    sv_perpendicular = sve * unit_parallele_enu[1] - svn * unit_parallele_enu[0]

    # covariance
    if sigma is not None:
        sig = np.array([sigma[0], sigma[1]])
        corr = np.eye(2)
        corr[0, 1] = sigma[2]
        corr[1, 0] = sigma[2]
        cov = glinalg.corr_to_cov(corr, sig)
        # rotation
        angle_rad = np.arctan2(unit_parallele_enu[1], unit_parallele_enu[2])
        #print('angle_rad in deg ' , np.degrees( angle_rad ))
        rot = np.zeros((2, 2))
        rot[0, 0] = np.cos(angle_rad)
        rot[0, 0] = np.cos(angle_rad)
        rot[1, 1] = np.cos(angle_rad)
        rot[0, 1] = -np.sin(angle_rad)
        rot[1, 0] = np.sin(angle_rad)
        # apply rotation
        vcv_en = np.dot(rot, np.dot(cov, rot.T))
        
        s_sv_parallel = np.sqrt(vcv_en[0, 0])

    # along / transverse great circle component end point
    R = mat_rot_general_to_local(np.radians(elon), np.radians(elat), unit="radians")
    OM = np.array([xs, ys, zs])
    OM = OM / np.linalg.norm(OM)

    unit_parallele_xyz = np.cross(POLE, OM)
    unit_parallele_enu = np.dot(R, unit_parallele_xyz)
    #print('unit_parallele_enu' , unit_parallele_enu*10. )
    ev_parallele = eve * unit_parallele_enu[0] + evn * unit_parallele_enu[1]
    ev_perpendicular = eve * unit_parallele_enu[1] - evn * unit_parallele_enu[0]

    # covariance
    if sigma is not None:
        sig = np.array([sigma[3], sigma[4]])
        corr = np.eye(2)
        corr[0, 1] = sigma[5]
        corr[1, 0] = sigma[5]
        cov = glinalg.corr_to_cov(corr, sig)
        # rotation
        angle_rad = np.arctan2(unit_parallele_enu[1], unit_parallele_enu[2])
        #print('angle_rad in deg ' , np.degrees( angle_rad ))
        rot = np.zeros((2, 2))
        rot[0, 0] = np.cos(angle_rad)
        rot[0, 0] = np.cos(angle_rad)
        rot[1, 1] = np.cos(angle_rad)
        rot[0, 1] = -np.sin(angle_rad)
        rot[1, 0] = np.sin(angle_rad)
        # apply rotation
        vcv_en = np.dot(rot, np.dot(cov, rot.T))
        
        s_ev_parallel = np.sqrt(vcv_en[0, 0])
        
        
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
        sig_v = np.sqrt(s_sv_parallel ** 2 + s_ev_parallel ** 2)
        sig_strain_rate = sig_v * 1.0e-3 / length_arc
        return  length_arc , length_rate , strain_rate,  sig_v, sig_strain_rate


