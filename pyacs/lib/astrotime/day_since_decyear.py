"""
Calculate the number of days since a reference decimal year.
"""

def day_since_decyear(decyear,ref_decyear):
    """
    Calculate the number of days since a reference decimal year.

    Parameters
    ----------
    decyear : float or array-like
        Decimal year(s).
    ref_decyear : float
        Reference decimal year.

    Returns
    -------
    float or ndarray
        Days elapsed. Negative means before reference date. Useful for
        postseismic analysis.

    Examples
    --------
    >>> import pyacs
    >>> ref_decyear=pyacs.cal2decyear(16,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 33))
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(17,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 32)),ref_decyear)
    >>> 0.9999884259232203
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(17,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 33)),ref_decyear)
    >>> 0.9999999999854481
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(16,4,2017,ut=pyacs.hmsmicros2ut(23, 58, 33)),ref_decyear)
    >>> 364.9999999999345
    >>> ld=pyacs.mjd2decyear(np.arange(pyacs.cal2mjd(13,4,2016),pyacs.cal2mjd(19,4,2016)))
    >>> pyacs.day_since_decyear(ld,ref_decyear)
    >>> array([-3.499, -2.499, -1.499, -0.499,  0.501,  1.501])
    """

    from .decyear2mjd import decyear2mjd

    mjd=decyear2mjd(decyear)
    
    ref_mjd=decyear2mjd(ref_decyear)

    return mjd-ref_mjd
