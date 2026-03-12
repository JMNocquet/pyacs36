from pyacs.lib.errors import OptionError
import numpy as np


###################################################################
def mat_rot_general_to_local(lam, phi, unit="radians"):
###################################################################

    """Rotation matrix from geocentric (XYZ) to local (ENU).

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
        Rotation matrix. R = [[-sin(lam), cos(lam), 0],
        [-sin(phi)*cos(lam), -sin(phi)*sin(lam), cos(phi)],
        [cos(phi)*cos(lam), cos(phi)*sin(lam), sin(phi)]].
    """

    if unit not in ["radians", "dec_deg"]:
        raise OptionError("unit option must be in [radians,dec_deg];unit=", unit)

    if unit == "dec_deg":
        lam = np.radians(lam)
        phi = np.radians(phi)
    
    R = np.zeros([3, 3], float)
    R[0, 0] = -np.sin(lam)
    R[0, 1] = np.cos(lam)

    R[1, 0] = -np.sin(phi) * np.cos(lam)
    R[1, 1] = -np.sin(phi) * np.sin(lam)
    R[1, 2] = np.cos(phi)
    
    R[2, 0] = np.cos(phi) * np.cos(lam)
    R[2, 1] = np.cos(phi) * np.sin(lam)
    R[2, 2] = np.sin(phi)

    return R


