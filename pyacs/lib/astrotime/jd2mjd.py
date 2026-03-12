"""
Convert Julian day to modified Julian day.
"""

def jd2mjd(jd):
    """
    Convert Julian day to modified Julian day.

    Parameters
    ----------
    jd : float
        Julian day.

    Returns
    -------
    float
        jd - 2400000.5.
    """

    return jd-2400000.5
