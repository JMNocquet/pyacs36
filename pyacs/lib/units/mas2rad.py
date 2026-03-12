"""Convert milliarcseconds to radians."""


def mas2rad(mas):
    """Convert milliarcseconds to radians.

    Parameters
    ----------
    mas : float or array_like
        Angle(s) in milliarcseconds.

    Returns
    -------
    float or numpy.ndarray
        Angle(s) in radians.

    Examples
    --------
    >>> import pyacs
    >>> pyacs.mas2rad(1.)
    4.84813681109536e-09
    >>> pyacs.mas2rad(np.array([1.0, 2.0]))
    array([  4.84813681e-09,   9.69627362e-09])
    """
    import numpy as np

    if not isinstance(mas, np.ndarray):
        np_mas = np.array(mas)
    else:
        np_mas = mas

    return (np_mas * 2. * np.pi / 360. / 60. / 60. / 1000.)

