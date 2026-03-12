"""Profile projection methods for Velocity_Field."""

import numpy as np
from pyacs.lib import coordinates as Coordinates
from pyacs.lib.gmtpoint import GMT_Point


def proj_profile(self, slon, slat, elon, elat, d, save=None, verbose=False):
    """Project velocity components along a great-circle profile.

    Parameters
    ----------
    slon : float
        Longitude of profile start (decimal degrees).
    slat : float
        Latitude of profile start (decimal degrees).
    elon : float
        Longitude of profile end (decimal degrees).
    elat : float
        Latitude of profile end (decimal degrees).
    d : float
        Maximum distance (km) from profile for a point to be included.
    save : str, optional
        If provided, write results to this file.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    tuple
        np_code, np_distance_along_profile, np_distance_to_profile,
        np_Ve, np_Vn, np_SVe, np_SVn, np_v_parallele, np_v_perpendicular,
        np_sigma_v_parallele, np_sigma_v_perpendicular, np_lazimuth.
    """
    lcode = []
    ldistance_along_profile = []
    ldistance_to_profile = []
    lVe = []
    lVn = []
    lSVe = []
    lSVn = []
    lv_parallele = []
    lv_perpendicular = []
    lsigma_v_parallele = []
    lsigma_v_perpendicular = []
    lazimuth = []

    Rt = 6371.0E3

    (xi, yi, zi) = Coordinates.geo2xyz(np.radians(slon), np.radians(slat), 0.0)
    Xi = np.array([xi, yi, zi])
    MS = GMT_Point(lon=slon, lat=slat)

    (xs, ys, zs) = Coordinates.geo2xyz(np.radians(elon), np.radians(elat), 0.0)
    Xs = np.array([xs, ys, zs])
    ME = GMT_Point(lon=elon, lat=elat)

    length_arc = MS.spherical_distance(ME) / 1000.0

    POLE = np.cross(Xi, Xs)
    POLE = POLE / np.linalg.norm(POLE)

    H_distance = {}

    Xi_unit = Xi / np.linalg.norm(Xi)

    for M in self.sites:
        if verbose:
            print('-- projecting ', M.code)
        (x, y, z) = Coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), 0.0)
        R = Coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat), unit='radians')
        OM = np.array([x, y, z])
        OM = OM / np.linalg.norm(OM)
        unit_parallele_xyz = np.cross(POLE, OM)
        unit_parallele_enu = np.dot(R, unit_parallele_xyz)
        v_parallele = M.Ve * unit_parallele_enu[0] + M.Vn * unit_parallele_enu[1]
        v_perpendicular = M.Ve * unit_parallele_enu[1] - M.Vn * unit_parallele_enu[0]
        azimuth = np.arctan2(unit_parallele_enu[0], unit_parallele_enu[1])
        x_m = np.dot(OM, Xi_unit)
        Xi_cross = np.cross(Xi, -POLE)
        Xi_cross = Xi_cross / np.linalg.norm(Xi_cross)
        y_m = np.dot(OM, Xi_cross)
        longitud_m = np.arctan2(y_m, x_m)
        distance_along_profile = longitud_m * Rt
        z_m = np.dot(OM, POLE)
        latitud_m = np.arcsin(z_m)
        distance_to_profile = latitud_m * Rt
        local_R = np.array([[np.cos(azimuth + np.pi / 2.), np.sin(azimuth + np.pi / 2.)],
                            [-np.sin(azimuth + np.pi / 2.), np.cos(azimuth + np.pi / 2.)]])
        cov = M.SVe * M.SVn * M.SVen
        vcv_en = np.array([[M.SVe ** 2, cov], [cov, M.SVn ** 2]])
        vcv_profile = np.dot(np.dot(local_R, vcv_en), local_R.T)
        sigma_v_parallele = np.sqrt(vcv_profile[0, 0])
        sigma_v_perpendicular = np.sqrt(vcv_profile[1, 1])

        if (np.fabs(distance_to_profile / 1000.0) < d) and (distance_along_profile / 1000.0 >= 0.0) and (distance_along_profile / 1000.0 <= length_arc):
            lcode.append(M.code)
            ldistance_along_profile.append(distance_along_profile / 1000.0)
            ldistance_to_profile.append(distance_to_profile / 1000.0)
            lVe.append(M.Ve)
            lVn.append(M.Vn)
            lSVe.append(M.SVe)
            lSVn.append(M.SVn)
            lv_parallele.append(v_parallele)
            lv_perpendicular.append(v_perpendicular)
            lsigma_v_parallele.append(sigma_v_parallele)
            lsigma_v_perpendicular.append(sigma_v_perpendicular)
            lazimuth.append(np.degrees(azimuth))
            H_distance[distance_along_profile] = (
                "%4s           %10.3lf               %10.3lf        %10.2lf %10.2lf %10.2lf %10.2lf    %10.2lf          %10.2lf              %10.2lf            %10.2lf          %10.1lf\n" %
                (M.code, distance_along_profile / 1000.0, distance_to_profile / 1000.0, M.Ve, M.Vn, M.SVe, M.SVn, v_parallele, v_perpendicular, sigma_v_parallele, sigma_v_perpendicular, np.degrees(azimuth)))

    if save is not None:
        print("-- writing results to ", save)
        fs = open(save, 'w')
        fs.write("#site  profile_distance_to_A (km)  distance_to_profile (km)       Ve         Vn        SVe        SVn    V_parallele     V_perpendicular            SV_parallele            SV_perpendicular        azimuth \n")
        fs.write("#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
        for distance in sorted(H_distance.keys()):
            fs.write("%s" % H_distance[distance])
        fs.close()

    np_distance_along_profile = np.array(ldistance_along_profile)
    lindex = np.argsort(np_distance_along_profile)
    np_code = np.array(lcode, dtype=str)[lindex]
    np_distance_along_profile = np_distance_along_profile[lindex]
    np_distance_to_profile = np.array(ldistance_to_profile)[lindex]
    np_Ve = np.array(lVe)[lindex]
    np_Vn = np.array(lVn)[lindex]
    np_SVe = np.array(lSVe)[lindex]
    np_SVn = np.array(lSVn)[lindex]
    np_v_parallele = np.array(lv_parallele)[lindex]
    np_v_perpendicular = np.array(lv_perpendicular)[lindex]
    np_sigma_v_parallele = np.array(lsigma_v_parallele)[lindex]
    np_sigma_v_perpendicular = np.array(lsigma_v_perpendicular)[lindex]
    np_lazimuth = np.array(lazimuth)[lindex]

    return (np_code, np_distance_along_profile, np_distance_to_profile,
            np_Ve, np_Vn, np_SVe, np_SVn,
            np_v_parallele, np_v_perpendicular,
            np_sigma_v_parallele, np_sigma_v_perpendicular, np_lazimuth)
