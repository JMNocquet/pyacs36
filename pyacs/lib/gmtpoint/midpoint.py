"""Midpoint between two GMT_Points."""

import pyacs.lib.coordinates


def midpoint(self, point, code='Mm'):
    """Return the midpoint between this point and another.

    Parameters
    ----------
    point : GMT_Point
        Other point.
    code : str, optional
        Code for the new point. Default is 'Mm'.

    Returns
    -------
    GMT_Point
        Midpoint (lon, lat, he set; velocity not set).
    """
    (Xi, Yi, Zi) = pyacs.lib.coordinates.geo2xyz(self.lon, self.lat, 0.0, 'dec_deg')
    (Xe, Ye, Ze) = pyacs.lib.coordinates.geo2xyz(point.lon, point.lat, 0.0, 'dec_deg')

    X = (Xi + Xe) / 2.
    Y = (Yi + Ye) / 2.
    Z = (Zi + Ze) / 2.

    (lonMm, latMm, heMm) = pyacs.lib.coordinates.xyz2geo(X, Y, Z, unit='dec_deg')

    Mm = self.__class__(code=code)
    Mm.lon = lonMm
    Mm.lat = latMm
    Mm.he = heMm

    return Mm
