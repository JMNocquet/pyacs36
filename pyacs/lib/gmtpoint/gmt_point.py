"""GMT_Point class: position and optional velocity (GMT psvelo style)."""


class GMT_Point:
    """Simple point with position and optional velocity (GMT psvelo style).

    Attributes
    ----------
    lon : float
        Longitude in decimal degrees.
    lat : float
        Latitude in decimal degrees.
    code : str, optional
        4-letter site code.
    he : float, optional
        Height above ellipsoid (m). Default is 0.
    Ve, Vn : float, optional
        East and north velocity components (mm/yr).
    Vu : float, optional
        Up component (mm/yr). Default is 0.
    SVe, SVn, SVu : float, optional
        Standard deviations on velocity components. Default is 0.
    SVen : float, optional
        Correlation or covariance term for Ve/Vn. Default is 0.
    """

    def __init__(self, code=None, lon=None, lat=None, he=0, Ve=None, Vn=None, Vu=0,
                 SVe=0, SVn=0, SVu=0, SVen=0,
                 Cv_xyz=None, Cv_enu=None, index=None):
        try:
            self.code = code
            self.Ve = Ve
            self.Vn = Vn
            self.Vu = Vu
            self.SVe = SVe
            self.SVn = SVn
            self.SVen = SVen
            self.SVu = SVu
            self.lon = lon
            self.lat = lat
            self.he = he
            self.Cv_xyz = Cv_xyz
            self.Cv_enu = Cv_enu
            self.index = index
        except Exception:
            raise ValueError('!!! Error. Could not init GMT_Point.')


# Attach methods from submodules
from . import copy as _copy_mod
from . import get_info as _get_info_mod
from . import magaz as _magaz_mod
from . import assign_index as _assign_index_mod
from . import get_index as _get_index_mod
from . import spherical_distance as _spherical_distance_mod
from . import pole as _pole_mod
from . import add_to_gmt_psvelo as _add_to_gmt_psvelo_mod
from . import rotate_vel as _rotate_vel_mod
from . import midpoint as _midpoint_mod

GMT_Point.copy = _copy_mod.copy
GMT_Point.get_info = _get_info_mod.get_info
GMT_Point.magaz = _magaz_mod.magaz
GMT_Point.assign_index = _assign_index_mod.assign_index
GMT_Point.get_index = _get_index_mod.get_index
GMT_Point.spherical_distance = _spherical_distance_mod.spherical_distance
GMT_Point.pole = _pole_mod.pole
GMT_Point.add_to_gmt_psvelo = _add_to_gmt_psvelo_mod.add_to_gmt_psvelo
GMT_Point.rotate_vel = _rotate_vel_mod.rotate_vel
GMT_Point.midpoint = _midpoint_mod.midpoint
