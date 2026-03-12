"""
Guess a date from input (string or float) and return decimal year.
"""

def guess_date(date):
    """
    Guess a date from input (string or float) and return decimal year.

    Parameters
    ----------
    date : str or float
        Date as string (e.g. '2016/04/16', '2016 107', '57500.0') or float
        (MJD or decimal year).

    Returns
    -------
    float
        Decimal year.

    Examples
    --------
    Calendar date formats (year/month/day or day/month/year):
    >>> import pyacs
    >>> pyacs.guess_date('2016/04/16')
    2016.2909836065573
    >>> pyacs.guess_date('2016-04-16')
    2016.2909836065573

    Day of year format (year doy or doy year):
    >>> pyacs.guess_date('2016 107')
    2016.2909836065573
    >>> pyacs.guess_date('107 2016')
    2016.2909836065573

    Decimal year format:
    >>> pyacs.guess_date('2016.29')
    2016.29
    >>> pyacs.guess_date(2016.29)
    2016.29

    MJD (Modified Julian Date) format:
    >>> pyacs.guess_date('57500.0')
    2016.2909836065573
    >>> pyacs.guess_date(57500.0)
    2016.2909836065573
    """
    
    from .mjd2decyear import mjd2decyear
    from .dayno2decyear import dayno2decyear
    from .cal2decyear import cal2decyear
    
    decyear='UNKNOWN'
    if isinstance(date,float):
        # it is then mjd or decimal year
        if date > 3000.0:
            # certainly a date is mjd type
            decyear=mjd2decyear(date)
        else:
            # certainly a date is decyear type
            decyear=date

    if isinstance(date,str):
        # cal,dayno or decyear passed as a string
        import re

        ldate=re.findall(r"[\w'.]+", date)
        
        if len(ldate)==1: # again decyear or mjd
            if float(date) > 3000.0:decyear=mjd2decyear(float(date))
            else:decyear=float(date)

        if len(ldate)==2: # dayno
            arg1=int(ldate[0])
            arg2=int(ldate[1])
            if arg1> 1980:
                year=arg1
                doy=arg2
            else:
                year=arg2
                doy=arg1
            decyear=dayno2decyear(doy,year,ut=0.5)
        if len(ldate)==3: # cal
            arg1=int(ldate[0])
            arg2=int(ldate[1])
            arg3=int(ldate[2])
            if arg1>1980:
                year=arg1
                month=arg2
                mday=arg3
            if arg3>1980:
                year=arg3
                month=arg2
                mday=arg1
            decyear=cal2decyear(mday, month, year, ut=0.5)
    
    if decyear == 'UNKNOWN':
        # we were not able to guess
        raise ValueError( ('!!! guess_date was unable to guess the format of %s ' % ( str(date) ) ) )
    
    return(decyear)
