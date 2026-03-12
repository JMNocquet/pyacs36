"""Rotate velocity components for GMT_Point."""

import numpy as np


def rotate_vel(self, angle, unit='radians'):
    """Rotate velocity components by an angle (clockwise).

    Parameters
    ----------
    angle : float
        Rotation angle.
    unit : str, optional
        'radians' or 'degrees'. Default is 'radians'.

    Returns
    -------
    GMT_Point
        New point with rotated Ve, Vn.
    """
    if unit == 'degrees':
        angle = np.radians(angle)

    N = self.copy()

    def __rotate__(V, rad_angle):
        cos = np.cos(rad_angle)
        sin = np.sin(rad_angle)
        R = np.array([[cos, -sin], [sin, cos]])
        return np.dot(R, V)

    V = np.array([self.Ve, self.Vn])
    rotated = __rotate__(V, angle)
    (N.Ve, N.Vn) = (rotated[0], rotated[1])

    return N
