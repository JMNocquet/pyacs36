"""Convert moment magnitude to seismic moment."""


def magnitude_to_moment(magnitude):
    """Convert moment magnitude Mw to seismic moment in N.m.

    Parameters
    ----------
    magnitude : float or array_like
        Moment magnitude Mw.

    Returns
    -------
    float or numpy.ndarray
        Seismic moment in N.m (10^(1.5*Mw + 9.05)).
    """
    import numpy as np

    moment = np.power(10, (1.5 * magnitude + 9.05))
    return moment

