import numpy as np


###################################################################
def vel_xyz_to_enu(vel_xyz, approx_xyz=None, approx_lon_lat=None):
###################################################################
    """
    Converts a velocity expressed in geocentric cartesian coordinates (X, Y, Z)
    into local North, East, Up components (dN, dE, dU) at a reference point.

    :param vel_xyz: velocity in geocentric cartesian coordinates (meters)
    :param approx_xyz: reference point for the local frame in geocentric cartesian coordinates
    :param approx_lon_lat: reference point for the local frame in geographical coordinates (decimal degree)

    :returns: Venu
    """

    from .xyz2geo import xyz2geo
    from .mat_rot_general_to_local import mat_rot_general_to_local
    from .geo2xyz import geo2xyz

    # check if approx_xyz or approx_lon_lat is provided
    if approx_xyz is None and approx_lon_lat is None:
        raise ValueError("Either approx_xyz or approx_lon_lat must be provided")
    if approx_xyz is not None and approx_lon_lat is not None:
        raise ValueError("Only one of approx_xyz or approx_lon_lat must be provided")
    if approx_xyz is None:
        approx_xyz = geo2xyz(approx_lon_lat[0], approx_lon_lat[1], 0)
    else:
        lon, lat, _ = approx_lon_lat

    # rotation matrix from global to local (returns components ordered as E, N, U)
    lam, phi, _ = xyz2geo(approx_xyz[0], approx_xyz[1], approx_xyz[2])
    rot = mat_rot_general_to_local(lam, phi)

    dxyz = np.array([vel_xyz[0], vel_xyz[1], vel_xyz[2]])
    enu = np.dot(rot, dxyz)

    # reorder to return in N, E, U order
    return enu


