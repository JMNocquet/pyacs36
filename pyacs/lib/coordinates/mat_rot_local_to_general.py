

###################################################################
def mat_rot_local_to_general(lam, phi, unit="radians"):
###################################################################
    """Rotation matrix from local (ENU) to geocentric (XYZ).

    Parameters
    ----------
    lam : float
        Longitude.
    phi : float
        Latitude.
    unit : {'radians', 'dec_deg'}, optional
        Units for lam and phi. Default is 'radians'.

    Returns
    -------
    R : ndarray, shape (3, 3)
        Rotation matrix. R is orthogonal; inverse equals transpose (same as
        mat_rot_general_to_local transposed).
    """
    from pyacs.lib.errors import OptionError
    from .mat_rot_general_to_local import mat_rot_general_to_local

    if unit not in ["radians", "dec_deg"]:
        raise OptionError("unit option must be in [radians,dec_deg];unit=", unit)
    
    return mat_rot_general_to_local(lam, phi, unit=unit).T


