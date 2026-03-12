"""
Return True if year is a leap year.
"""

def leap_year(year):
    """
    Return True if year is a leap year.

    Parameters
    ----------
    year : int or array-like
        Year (YYYY).

    Returns
    -------
    bool or ndarray
        True if leap year, False otherwise.
    """

    import numpy as np

    if isinstance(year, np.ndarray):
        OK=np.array(list(map(leap_year,year))).T

    else:
    
        year=int(year)
        OK=(( not (year%4))and(year%100))or(not(year%400))
        
    return OK
