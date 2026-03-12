"""
Convert fractional day to seconds since 00:00:00.0.
"""

def ut2uts(ut):
    """
    Convert fractional day to seconds since 00:00:00.0.

    Parameters
    ----------
    ut : float
        Fractional day in [0., 1.[.

    Returns
    -------
    float
        ut * 60.*60.*24. (seconds in day).
    """
    
    nsec_day=60.*60.*24.
    return(ut*nsec_day)
