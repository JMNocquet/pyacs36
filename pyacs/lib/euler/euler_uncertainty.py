"""Euler pole parameter uncertainties from rotation vector covariance."""

import numpy as np
import pyacs.lib.coordinates as Coordinates

from .rot2euler import rot2euler


def euler_uncertainty(w, vcv):
    """Calculate Euler pole parameter uncertainties from rotation vector covariance.

    Propagates the covariance of the rotation vector (radians/yr) into
    uncertainties of the Euler pole (semi-major, semi-minor, azimuth, omega).

    Parameters
    ----------
    w : array_like, shape (3,)
        Rotation vector in XYZ (geocentric) coordinates, radians/yr.
    vcv : array_like, shape (3, 3)
        Covariance matrix of `w`.

    Returns
    -------
    max_sigma : float
        Semi-major axis of error ellipse, degrees.
    min_sigma : float
        Semi-minor axis of error ellipse, degrees.
    azimuth : float
        Azimuth of the error ellipse, degrees.
    s_omega : float
        Uncertainty of angular velocity, deg/Myr.
    """
    nw = np.sqrt(w[0]**2 + w[1]**2 + w[2]**2)

    (llambda, phi, omega) = rot2euler(w[0], w[1], w[2])

    R = Coordinates.mat_rot_general_to_local(np.radians(llambda), np.radians(phi))

    VCV_POLE_ENU = np.dot(np.dot(R, vcv), R.T)

    vcv11 = VCV_POLE_ENU[0, 0]
    vcv12 = VCV_POLE_ENU[0, 1]
    vcv22 = VCV_POLE_ENU[1, 1]
    vcv33 = VCV_POLE_ENU[2, 2]

    max_sigma = 0.5 * (vcv11 + vcv22 + np.sqrt((vcv11 - vcv22)**2 + 4 * vcv12**2))
    max_sigma = np.degrees(np.arctan(np.sqrt(max_sigma) / nw))
    min_sigma = 0.5 * (vcv11 + vcv22 - np.sqrt((vcv11 - vcv22)**2 + 4 * vcv12**2))
    min_sigma = np.degrees(np.arctan(np.sqrt(min_sigma) / nw))
    azimuth = np.degrees(2 * np.arctan(2 * vcv12 / (vcv11 - vcv22)))
    s_omega = np.degrees(np.sqrt(vcv33)) * 1.E06

    return (max_sigma, min_sigma, azimuth, s_omega)
