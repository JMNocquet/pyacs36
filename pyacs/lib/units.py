"""
Various routines dealing with units conversion.
"""


###################################################################
def mas2rad(mas):
###################################################################
    """Converts milliarcseconds to radians

    :param mas: angle in millarcseconds
    
    :returns: rad: in radians
    
    :Example:
    >>> import pyacs
    >>> pyacs.mas2rad(1.)
    >>> 4.84813681109536e-09
    >>> pyacs.mas2rad(np.array([1.0,2.0]))
    >>> array([  4.84813681e-09,   9.69627362e-09])
    
    """
    import numpy as np
    
    if not isinstance(mas, np.ndarray):
        np_mas = np.array( mas )
    else:
        np_mas = mas
    
    return( np_mas*2.*np.pi/360./60./60./1000.)

###################################################################
def rad2mas(rad):
###################################################################
    """
    Converts radians to milliarcseconds
    
    :param rad: angle in radians
    
    :returns mas: angle in millarcseconds: 
    
    """

    import numpy as np

    return(rad/2./np.pi*360.*60.*60.*1000.)

###################################################################
def radians2deg_mn_sec(rad,angle_range='-180-180'):
###################################################################
    """
    Converts an angle from radians to deg, minutes, seconds 

    :param rad: angle in radians 
    :param angle_range: either '-180-180' or '0-360'

    :returns: deg, minutes, seconds 
    :rtype: int,int,float
    
    """

    import numpy as np

    
    # to degrees
    decimal_degrees=np.degrees(rad)
    # within 0-360
    decimal_degrees=np.remainder(decimal_degrees,360.)
    # -180 - 180 if asked
    if angle_range=='-180-180':
        if decimal_degrees>180.:
            decimal_degrees = decimal_degrees - 360.
    
    angle_sign=np.sign(decimal_degrees)

    # we work only with positive angles
    decimal_degrees=decimal_degrees * angle_sign
    
            
    deg=int(decimal_degrees)
    
    deg_frac=decimal_degrees-float(deg)
    
    minutes=int(deg_frac*60.)
    
    minutes_frac=deg_frac-float(minutes)/60.
    
    seconds=minutes_frac*3600.
    
    return(angle_sign*deg,minutes,seconds)

###################################################################
def moment_to_magnitude(moment):
###################################################################
    """
    Converts a moment in N.m to Mw
    
    :param moment: in N.m
    :return magnitude: 2./3.*(math.log10(M0)-9.05)
    """
    
    import numpy as np
    magnitude=2./3.*(np.log10(moment)-9.05)
    return( magnitude )

###################################################################
def magnitude_to_moment(magnitude):
###################################################################
    """
    Converts a Mw magnitude to moment in N.m
    
    :param magnitude: Mw
    :return moment: in N.m, 10^(1.5*magnitude+9.5)
    """
    
    import numpy as np
    moment = np.power(10,(1.5*magnitude+9.05))
    return( moment )

