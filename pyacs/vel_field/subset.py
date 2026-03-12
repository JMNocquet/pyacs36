"""Subset and transform methods for Velocity_Field."""

import copy
import numpy as np


def subset(self, lonly=None, lexclude=None):
    """Return a new Velocity_Field from a subset of sites.

    Parameters
    ----------
    lonly : list, optional
        If provided, only these site codes are included.
    lexclude : list, optional
        If provided, these site codes are excluded (ignored if lonly is set).

    Returns
    -------
    Velocity_Field
        New velocity field with the subset of sites.
    """
    lGpoint = []
    for M in self.sites:
        if lonly is not None:
            if M.code in lonly:
                lGpoint.append(M)
        else:
            if lexclude is not None and M.code not in lexclude:
                lGpoint.append(M)
    return self.__class__(lgmt_points=lGpoint)


def radial(self, center):
    """Return a velocity field with radial and tangential components about a center.

    Parameters
    ----------
    center : array_like
        [longitude, latitude] of center in decimal degrees.

    Returns
    -------
    Velocity_Field
        New velocity field with Ve = radial, Vn = tangential.
    """
    lnew_point = []
    for code in self.lcode():
        gmt_point = self.site(code)
        radial_unit_vector = np.array([gmt_point.lon - center[0], gmt_point.lat - center[1]])
        radial_unit_vector = radial_unit_vector / np.sqrt(np.sum(radial_unit_vector ** 2))
        tangential_unit_vector = np.array([radial_unit_vector[1], -radial_unit_vector[0]])
        new_gmt_point = copy.copy(gmt_point)
        new_gmt_point.Ve = np.sum(np.array([gmt_point.Ve, gmt_point.Vn]) * radial_unit_vector)
        new_gmt_point.Vn = np.sum(np.array([gmt_point.Ve, gmt_point.Vn]) * tangential_unit_vector)
        lnew_point.append(new_gmt_point)
    return self.__class__(lgmt_points=lnew_point)
