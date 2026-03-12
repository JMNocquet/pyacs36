"""Convert radians to degrees/minutes/seconds."""


def radians2deg_mn_sec(rad, angle_range='-180-180'):
    """Convert an angle from radians to degrees, minutes, seconds.

    Parameters
    ----------
    rad : float
        Angle in radians.
    angle_range : str, optional
        '0-360' or '-180-180'. Default is '-180-180'.

    Returns
    -------
    deg : int
        Degrees (signed if angle_range is '-180-180').
    minutes : int
        Minutes.
    seconds : float
        Seconds.
    """
    import numpy as np

    # to degrees
    decimal_degrees = np.degrees(rad)
    # within 0-360
    decimal_degrees = np.remainder(decimal_degrees, 360.)
    # -180 - 180 if asked
    if angle_range == '-180-180':
        if decimal_degrees > 180.:
            decimal_degrees = decimal_degrees - 360.

    angle_sign = np.sign(decimal_degrees)

    # we work only with positive angles
    decimal_degrees = decimal_degrees * angle_sign

    deg = int(decimal_degrees)

    deg_frac = decimal_degrees - float(deg)

    minutes = int(deg_frac * 60.)

    minutes_frac = deg_frac - float(minutes) / 60.

    seconds = minutes_frac * 3600.

    return (angle_sign * deg, minutes, seconds)

