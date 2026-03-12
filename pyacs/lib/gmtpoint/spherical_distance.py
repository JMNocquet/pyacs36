"""Spherical distance for GMT_Point."""

import pyacs.lib.coordinates
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def spherical_distance(self, M):
    """Return spherical distance between this point and another (meters).

    Parameters
    ----------
    M : GMT_Point
        Other point.

    Returns
    -------
    float
        Distance along the sphere in meters.
    """
    DEBUG(("%.5lf %.5lf %.5lf %.5lf ") % (self.lon, self.lat, M.lon, M.lat))

    return pyacs.lib.coordinates.geo_spherical_distance(
        self.lon, self.lat, 0.0, M.lon, M.lat, 0.0, unit='dec_deg')
