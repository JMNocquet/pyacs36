"""
Convert day of year and year to a Python datetime instance.
"""

def dayno2datetime(dayno,year,ut=0.5):
    """
    Convert day of year and year to a Python datetime instance.

    Parameters
    ----------
    dayno : int or array-like
        Day of year.
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction. Default 0.5 (12:00).

    Returns
    -------
    datetime or ndarray of datetime
        Datetime instance(s). At 12:00 by default; specify ut for other times.
    """

    from datetime import datetime
    from ._common import __daynoOK, __utOK
    from .dayno2cal import dayno2cal
    from .uts2hmsmicros import uts2hmsmicros
    from .ut2uts import ut2uts
    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)
    
    if isinstance(dayno, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=dayno*0.0+ut
        return_datetime=np.array(list(map(dayno2datetime,dayno,year,ut)))

    else:
    
        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
        __utOK(ut) 

        (day, month)=dayno2cal(dayno,year)
        (h,mn,s,microsecond)=uts2hmsmicros(ut2uts(ut))
        return_datetime=datetime(year,month,day, h, mn, s, int(microsecond))
    
    return return_datetime
