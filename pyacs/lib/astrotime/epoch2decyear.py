"""
Convert SINEX-format epoch string 'YY:DOY:SOD' (e.g. '17:100:43185') to decimal year.
"""

def epoch2decyear(epoch):
    """
    Convert SINEX-format epoch string 'YY:DOY:SOD' (e.g. '17:100:43185') to decimal year.

    Parameters
    ----------
    epoch : str or array-like
        Epoch string(s).

    Returns
    -------
    float or ndarray
        Decimal year.
    """

    from .yr2year import yr2year
    from .dayno2decyear import dayno2decyear
    from .uts2ut import uts2ut
    import numpy as np

    if isinstance(epoch, list):
        epoch=np.array(epoch)

    if isinstance(epoch, np.ndarray):
        decyear=np.array(list(map(epoch2decyear,epoch))).T

    else:
    
        lepoch = epoch.split(':')
        yr = int(lepoch[0])
        doy = int(lepoch[1])
        sod = int(lepoch[2])
    
        year=yr2year(yr)
    
        decyear=dayno2decyear(doy,year,ut=uts2ut(sod))

    return(decyear)
