"""Convert seismic moment to moment magnitude."""


def moment_to_magnitude(moment):
    """Convert seismic moment (N.m) to moment magnitude Mw.

    Parameters
    ----------
    moment : float or array_like
        Seismic moment in N.m.

    Returns
    -------
    float or numpy.ndarray
        Moment magnitude Mw (2/3 * (log10(M0) - 9.05)).
    """
    import numpy as np

    magnitude = 2. / 3. * (np.log10(moment) - 9.05)
    return magnitude

