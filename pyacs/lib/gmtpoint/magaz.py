"""Magnitude and azimuth of velocity for GMT_Point."""

import math


def magaz(self):
    """Return magnitude and azimuth of velocity for the current GMT_Point.

    Returns
    -------
    mag : float
        Velocity magnitude in mm/yr.
    az : float
        Azimuth in decimal degrees, clockwise from North.

    Notes
    -----
    Uncertainty on magnitude/azimuth is not implemented.
    """
    mag = math.sqrt(self.Ve**2 + self.Vn**2)
    az = math.degrees(math.atan2(self.Ve, self.Vn))
    return (mag, az)
