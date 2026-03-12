"""
Round a Python datetime instance to the nearest second.
"""

def datetime_round_second(my_datetime):
    """
    Round a Python datetime instance to the nearest second.

    Parameters
    ----------
    my_datetime : datetime
        Datetime instance to round.

    Returns
    -------
    datetime
        New rounded datetime (microsecond set to 0).
    """

    import datetime
    

    if my_datetime.microsecond <= 500000:
        my_datetime = my_datetime.replace( microsecond = 0 )
    else:
        my_datetime = my_datetime.replace( microsecond = 0 )
        deltatime = datetime.timedelta(seconds=1)
        my_datetime = my_datetime + deltatime
    
    return(my_datetime)
