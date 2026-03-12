"""Convert radians to milliarcseconds."""


def rad2mas(rad):
    """Convert radians to milliarcseconds.

    Parameters
    ----------
    rad : float or array_like
        Angle(s) in radians.

    Returns
    -------
    float or numpy.ndarray
        Angle(s) in milliarcseconds.
    """
    import numpy as np

    return (rad / 2. / np.pi * 360. * 60. * 60. * 1000.)

