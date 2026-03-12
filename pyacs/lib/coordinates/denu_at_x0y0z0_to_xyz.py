import numpy as np


###################################################################
def denu_at_x0y0z0_to_xyz(de, dn, du, x0, y0, z0):
###################################################################
    """Convert local vector (dE, dN, dU) at a reference point to geocentric XYZ.

    Parameters
    ----------
    de : float or array_like
        East component of local vector, meters.
    dn : float or array_like
        North component of local vector, meters.
    du : float or array_like
        Up component of local vector, meters.
    x0 : float
        X of reference point in geocentric cartesian, meters.
    y0 : float
        Y of reference point in geocentric cartesian, meters.
    z0 : float
        Z of reference point in geocentric cartesian, meters.

    Returns
    -------
    x, y, z : float or ndarray
        Resulting point in geocentric cartesian coordinates, meters.
    """
    from .xyz2geo import xyz2geo
    from .mat_rot_local_to_general import mat_rot_local_to_general
    
    (lam, phi, _) = xyz2geo(x0, y0, z0)
    R = mat_rot_local_to_general(lam, phi)
    
    DL = np.array([de, dn, du])
    X0 = np.array([x0, y0, z0])
    XYZ = np.dot(R, DL) + X0
    return (XYZ[0], XYZ[1], XYZ[2])


