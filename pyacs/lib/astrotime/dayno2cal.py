"""
Return the day and month corresponding to day of year.
"""

def dayno2cal(dayno, year):
    """
    Return the day and month corresponding to day of year.

    Parameters
    ----------
    dayno : int or array-like
        Day of year (1-366).
    year : int or array-like
        Year.

    Returns
    -------
    day : int or ndarray
        Day of month.
    month : int or ndarray
        Month (1-12).
    """

    from ._common import __daynoOK, days
    from .leap_year import leap_year
    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)

 
    if isinstance(dayno, np.ndarray):
        [day,month]=np.array(list(map(dayno2cal,dayno,year))).T

    else:

        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
    
        if (leap_year(year)):days[1] = 29
        else:days[1] = 28
    
        month = 0
        end = days[month]
        while dayno > end:
            month=month+1
            end= end + days[month]
      
        end = end - days[month]
        day = dayno - end
        month=month+1

    return( np.array( day, dtype=int)+0 ,\
            np.array( month, dtype=int)+0)
