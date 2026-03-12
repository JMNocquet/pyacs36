"""
Convert seconds since 00:00:00.0 to fractional day.
"""

def uts2ut(uts):
    """
    Convert seconds since 00:00:00.0 to fractional day.

    Parameters
    ----------
    uts : float or array-like
        Seconds in [0., 86400.[.

    Returns
    -------
    float or ndarray
        uts / (60.*60.*24.) (fractional day).
    """
    
    from ._common import __utsOK
    import numpy as np

    if isinstance(uts, list):
        uts=np.array(uts)
    
    if isinstance(uts, np.ndarray):
        ut=np.array(list(map(uts2ut,uts)))

    else:
        
        __utsOK(uts)
        nsec_day=60.*60.*24.
        ut=uts/nsec_day
    
    return ut
