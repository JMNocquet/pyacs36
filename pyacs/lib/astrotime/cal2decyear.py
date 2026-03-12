"""
Convert a calendar date to decimal year.
"""

def cal2decyear(day,month,year,ut=0.5):
    """
    Convert a calendar date to decimal year.

    Parameters
    ----------
    day : int or array-like
        Day of month.
    month : int or array-like
        Month (1-12).
    year : int or array-like
        Year.
    ut : float, optional
        Day fraction in [0., 1.[. Default 0.5 (middle of day).

    Returns
    -------
    float or ndarray
        Decimal year.

    Notes
    -----
    When ut is not provided, the middle of the day is used.
    """

    from ._common import __monthOK, __dayOK
    from .cal2dayno import cal2dayno
    from .dayno2decyear import dayno2decyear
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)
   
    if isinstance(day, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=day*0.0+ut
        decyear=np.array(list(map(cal2decyear,day, month, year, ut)))

    else:
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year) 

        dayno=cal2dayno(day,month,year)
        decyear=dayno2decyear(dayno,year,ut=ut)
        
    return decyear
