"""
Convert modified Julian day to calendar date (universal time).
"""

def mjd2cal(mjd):
    """
    Convert modified Julian day to calendar date (universal time).

    Based on the slalib routine sla_djcl.

    Parameters
    ----------
    mjd : float or array-like
        Modified Julian day (JD - 2400000.5).

    Returns
    -------
    day : int or ndarray
        Day of month.
    month : int or ndarray
        Month (1-12).
    year : int or ndarray
        Year.
    ut : float or ndarray
        Day fraction in [0., 1.[.
    """

    from ._common import __utOK
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)

    if isinstance(mjd, np.ndarray):
        [day,month,year,ut]=np.array(list(map(mjd2cal,mjd))).T

    else:
    
    
        ut = mjd-int(mjd)

        # check arguments and raise an Error if not OK
        __utOK(ut) 
    
        mmjd=int(mjd)
    
        jd = mmjd + 2400001
    
        # Do some rather cryptic calculations
        # For Python3.6 a/b returns a float even if a & b are integer
        # for integer operation, / must be changed to //
        

        temp1 = 4*(jd+((6*(((4*jd-17918) // 146097))) // 4+1) // 2-37)
        temp2 = 10*(((temp1-237)%1461) // 4)+5

#        temp1 = 4 * ( jd + int( ( int( (6 * ( ( int( (4*jd-17918) // 146097 ) ) ) ) // 4)  + 1 )/2 ) - 37 )
#        temp2 = 10*(((temp1-237)%1461)/4)+5
    
        year = temp1 // 1461-4712
        month =((temp2 // 306+2)%12)+1
        day = (temp2%306) // 10+1

    return day, month, year, ut
