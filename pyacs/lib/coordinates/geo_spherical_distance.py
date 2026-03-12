from pyacs.lib.errors import OptionError


###################################################################
def geo_spherical_distance(lon1, lat1, h1, lon2, lat2, h2, unit="radians", Rt=6.371E6):
###################################################################
    """Spherical distance between two points given in geodetic coordinates.

    Parameters
    ----------
    lon1 : float or array_like
        Longitude of first point.
    lat1 : float or array_like
        Latitude of first point.
    h1 : float or array_like
        Height of first point (used with geo2xyz; radius implied).
    lon2 : float or array_like
        Longitude of second point.
    lat2 : float or array_like
        Latitude of second point.
    h2 : float or array_like
        Height of second point.
    unit : {'radians', 'dec_deg'}, optional
        Units for longitude and latitude. Default is 'radians'.
    Rt : float, optional
        Mean Earth radius in meters. Default 6.371e6.

    Returns
    -------
    float or ndarray
        Distance in meters.
    """

    from .geo2xyz import geo2xyz
    from .xyz_spherical_distance import xyz_spherical_distance

    if unit not in ["radians", "dec_deg"]:
       raise OptionError("unit option must be in [radians,dec_deg];unit=", unit)

    (x1, y1, z1) = geo2xyz(lon1, lat1, 0.0, unit=unit)
    (x2, y2, z2) = geo2xyz(lon2, lat2, 0.0, unit=unit)
 
 
    return xyz_spherical_distance(x1, y1, z1, x2, y2, z2)


