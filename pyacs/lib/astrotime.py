
"""

AstroTime.py is a traduction from perl module Astro::Time distributed by CPAN 
and written by Chris Phillips (Chris.Phillips@csiro.au). Conversion to python made by Jean-Mathieu Nocquet (Geoazur - CNRS-IRD-OCA-Univ. Nice, nocquet@geoazur.unice.fr)

AstroTime contains a set of Python routines for time based
conversions, such as conversion between calendar dates and Modified
Julian day and conversion from UT to local sidereal time. Included are
routines for conversion between numerical and string representation of
angles. All string functions removed.

Conversion from and to datetime objects have also been added for convenience.

:note: This module is intended dates manipulation, but leap seconds for instance are not handled. \
       See the GPSTime module for GPS, UTC times conversions. 
:note: python 3.6 : all / operator changed to // to force integer operations 

:warning: !!! all conversions not specifying ut or uts use 12h00mn00s as a default. ut is the decimal fraction of day
       and uts is the number of seconds since 00h00m00.0s

"""



days = [31,28,31,30,31,30,31,31,30,31,30,31]

###############################################################################
# PRIVATE ARGUMENT CHECK FUNCTIONS
###############################################################################

# Is the dayno valid?
def __daynoOK(dayno,year):
    
    if (dayno<1 or dayno>366 or (dayno>365 and not leap_year(year))):
        raise ValueError( ('!!! doy out of range: (doy/year)=(%d,%d) ' % ( dayno,year ) ) )

# Is the month valid?
def __monthOK(month):
    import numpy as np 
    if (not isinstance(month,int)) and (not isinstance(month,np.int64)):
        raise TypeError( ('!!! bad type for month. Must be integer month= ' , month ) )
        

    if (month > 12 or month < 1):
        raise ValueError( ('!!! month out of range. Must be in the range 1-12 month=%d ' % month ) )
        

# IS the day of month OK? (assumes month IS ok - should be checked first)
def __dayOK (day,month,year):
    month=month-1             # For array indexing
    if (leap_year(year)):days[1] = 29
    else:days[1] = 28

    if (day < 1 or day > days[month]):
        raise ValueError( ('!!! day of month out of range: (day/month/year)=(%d,%d,%d) ' % ( day,month,year ) ) )

# Is the day fraction OK?
def __utOK (ut):
    if (ut < 0.0 or ut >= 1.0):
        raise ValueError( ('!!! fraction day ut out of range [0;1[: ut=%.9lf  ' % ( ut ) ) )

# Is the number of seconds  OK?
def __utsOK (uts):
    if (uts < 0.0 or uts >= 24.0*60.0*60):
        raise ValueError( ('!!! numer of seconds out of range [0;24.0*60.0*60[: uts=%.9lf  ' % ( uts ) ) )


###############################################################################
# CALENDAR DATES CONVERSIONS
###############################################################################

###############################################################################
def cal2dayno (day, month, year):
###############################################################################
    """
    Returns the day of year (doy).
    
    :param day,month,year:
    :type: int,int,int
    :returns: doy of year as int

    :Example:
    >>> import pyacs
    >>> pyacs.cal2dayno(1,1,2000)
    >>> 1
    >>> pyacs.cal2dayno(29,2,2000)
    >>> 60
    >>> pyacs.cal2dayno(29.,2,2000)
    >>> 60.
    >>> pyacs.cal2dayno(29,2,2001)
    >>> !!!  29  out of range
    >>> !!! bad day
    """
    
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        dayno=np.array(list(map(cal2dayno,day, month, year)))

    else:
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year) 
    
        month=month-1 # For array indexing
    
        if (leap_year(year)):days[1] = 29
        else:days[1] = 28
    
    
        dayno = day
        mon=0
        while mon<month:
            dayno = dayno + days[mon]
            mon=mon+1

    return( np.array(dayno, dtype=int)+0 )

###############################################################################
def cal2mjd(day, month, year, ut=0.5):
###############################################################################
    """  
    Converts a calendar date (universal time) into modified Julian day number.
    mjd     Modified Julian day (JD-2400000.5)

    :param day, month, year:
    :param ut : day fraction ([0.,1.[)
    :returns: mjd (modified Julian day (julian day -2400000.5))
    :rtype: float
    :note: Be aware that when ut is not provided, then the middle of the day is used. See example.

    :Example:
    >>> import pyacs
    >>> pyacs.cal2mjd(29,2,2000)
    >>> 51603.5
    >>> pyacs.cal2mjd(29,2,2000,ut=0.0)
    >>> 51603.0
    """

    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=day*0.0+ut
        mjd=np.array(list(map(cal2mjd,day, month, year, ut)))

    else:
        
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year)
        __utOK(ut)
    
        if (month <= 2):
            m = int(month+9)
            y = int(year-1)
        else:
            m = int(month-3)
            y = int(year)
    
        c = int(y/100)
        y = y-c*100
        x1 = int(146097.0*c/4.0)
        x2 = int(1461.0*y/4.0)
        x3 = int((153.0*m+2.0)/5.0)

        mjd=x1+x2+x3+day-678882+ut

    return(mjd)

###############################################################################
def cal2decyear(day,month,year,ut=0.5):
###############################################################################
    """
    converts a calendar date to decimal year

    :param day, month, year:
    :param ut : day fraction ([0.,1.[)
    :returns: decimal year
    :rtype: float
    :note: Be aware that when ut is not provided, then the middle of the day is used.
    """

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

###############################################################################
def cal2datetime(day,month,year,ut=0.5):
###############################################################################
    """
    converts a calendar date to a datime python instance

    :param day, month, year:
    :param ut : day fraction ([0.,1.[), optional

    :note:!!! the returned datetime is at 12:00, that is ut=0.5 by default
    """
    
    import numpy as np

    if isinstance(day, list):
        day=np.array(day)
        month=np.array(month)
        year=np.array(year)

    
    if isinstance(day, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=day*0.0+ut
        return_datetime=np.array(list(map(cal2datetime,day,month,year,ut)))

    else:
    
        # check arguments and raise an Error if not OK
        __monthOK(month)
        __dayOK(day, month, year) 

        dayno=cal2dayno(day,month,year)
        return_datetime=dayno2datetime(dayno,year,ut=ut)
        
    return return_datetime


###############################################################################
# DAY OF YEAR CONVERSIONS
###############################################################################

###############################################################################
def dayno2cal(dayno, year):
###############################################################################

    """        
    Returns the day and month corresponding to dayno of year.
    
    :param dayno, year:
    :returns: day, month
    
    """

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

###############################################################################
def dayno2mjd(dayno, year, ut=0.5):
###############################################################################
    """
    converts a dayno and year to modified Julian day

    :param dayno, year:
    :param ut: day fraction, optional
    :returns: modified Julian day (julian day - 2400000.5)
    :note: Be aware that when ut is not provided, then the middle of the day is used. See example.
    :Example:
    >>> import pyacs
    >>> pyacs.dayno2mjd(60,2000)
    >>> 51603.5
    >>> pyacs.dayno2mjd(60.,2000,ut=0)
    >>> 51603.0
    """

    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)
   
    if isinstance(dayno, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=dayno*0.0+ut
        mjd=np.array(list(map(dayno2mjd,dayno, year, ut)))

    else:

        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
        __utOK(ut)
    
        (day, month) = dayno2cal(dayno, year)
        mjd=cal2mjd(day, month, year, ut)
    
    return mjd 

###############################################################################
def dayno2decyear(dayno,year,ut=0.5):
###############################################################################
    """
    converts julian day to decimal year

    :param dayno, year:
    :param ut: day fraction, optional

    :returns: decimal year
    :note: !!! accounts for leap years with 366 days
    :note: !!! unless ut is specified, the returned decimal years is at 12:00, that is ut=0.5 by default

    :Example:
    >>> import pyacs
    >>> pyacs.dayno2decyear(60,1999)
    >>> 1999.16301369863
    >>> pyacs.dayno2decyear(60,2000)
    >>> 2000.1625683060108
    >>> pyacs.dayno2decyear(60,2001)
    >>> 2001.16301369863
    >>> pyacs.dayno2decyear(60,2001,ut=0.)
    >>> 2001.1616438356164
    
    """
    
    import numpy as np

    if isinstance(dayno, list):
        dayno=np.array(dayno)
        year=np.array(year)

    
    if isinstance(dayno, np.ndarray):
        if not isinstance(ut, np.ndarray):ut=dayno*0.0+ut
        decyear=np.array(list(map(dayno2decyear,dayno, year, ut)))

    else:
    
        # check arguments and raise an Error if not OK
        __daynoOK(dayno, year) 
        __utOK(ut) 
    
        nday=cal2dayno(31,12,year)
        dayno=dayno-1
        
        dayno=dayno+ut

        decyear = float(year)+float(dayno)/float(nday)
    
    return decyear

###############################################################################
def dayno2datetime(dayno,year,ut=0.5):
###############################################################################
    """
    julian day to datetime python instance

    :param dayno, year:
    :param ut: day fraction, optional

    :returns: datetime

    :note:!!! the returned datetime is at 12:00, that is uts=24.0*60.0*60/2.0 by default

    :note:!!! specify uts for other time of the day

    """

    from datetime import datetime

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


###############################################################################
# MODIFIED JULIAN DAY CONVERSIONS
###############################################################################

###############################################################################
def mjd2cal(mjd):
###############################################################################

    """
    Converts a modified Julian day number into calendar date (universal time). (based on the slalib routine sla_djcl).
    
    :param mjd: modified Julian day (JD-2400000.5)
    :returns: day,month,year,ut
    :rtype: int,int,int,float
    :note: ut is the the day fraction in [0., 1.[.
    """

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


###############################################################################
def mjd2dayno(mjd):
###############################################################################
    """  
    converts a modified Julian day into year and dayno (universal time).
    
    :param mjd: modified julian day
    :returns: dayno,year,ut
    :rtype: int,int,float

    :Example:
    >>> import pyacs
    >>> pyacs.mjd2dayno(51603.0)
    >>> (60, 2000, 0.0)
    >>> pyacs.mjd2dayno(51603.5)
    >>> (60, 2000, 0.5)
    """
    
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
   
    if isinstance(mjd, np.ndarray):
        [dayno,year,ut]=np.array(list(map(mjd2dayno,mjd))).T
    
    else:
        
        (day, month, year, ut) = mjd2cal(mjd)
        dayno=cal2dayno(day,month,year)
    return  dayno, year, ut

###############################################################################
def mjd2decyear(mjd):
###############################################################################
    """
    converts a modified julian day to decimal year conversion

    :param mjd: modified julian day
    :returns: decimal year

    """
    
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
 
    if isinstance(mjd, np.ndarray):
        decyear=np.array(list(map(mjd2decyear,mjd)))
    
    else:
        (dayno,year,ut)=mjd2dayno(mjd)
        decyear = dayno2decyear(dayno,year,ut=ut)
    
    return decyear

###############################################################################
def mjd2datetime(mjd):
###############################################################################
    """

    converts Modified Julian Day date to datime python instance

    :param mjd: modified julian day
    :returns: datetime instance

    """

    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
  
    if isinstance(mjd, np.ndarray):
        return_datetime=np.array(list(map(mjd2datetime,mjd)))
                                 
    else:
    
        (dayno,year,ut)=mjd2dayno(mjd)
        return_datetime=dayno2datetime(dayno,year,ut=ut)
    
    return return_datetime

###############################################################################
def mjd2gpsweek(mjd):
###############################################################################

    """
    converts a modified julian day to gps week and gps day of week
    
    :param mjd: modified julian day
    :returns: gps_week,day_of_week
    
    """
    
    import numpy as np

    if isinstance(mjd, list):
        mjd=np.array(mjd)
    
    if isinstance(mjd, np.ndarray):
        [gweek,dow]=np.array(list(map(mjd2gpsweek,mjd))).T
    
    else:
    
        mjd010580 = 44243
        mjd_d = ( int(mjd) - mjd010580 ) - 1
        gweek = mjd_d // 7
        dow = mjd_d % 7

    return(gweek,dow)


###############################################################################
# DECIMAL YEAR CONVERSION
###############################################################################

###############################################################################
def decyear2cal(decyear):
###############################################################################
    """
    decimal year to calendar date conversion

    :param decyear: decimal year
    :returns: day of month,month,ut
    :note: the input decimal year is assumed to account for leap years, that is the day of year is the \
    decimal of decyear * number of days in the year.

    """

    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)

    if isinstance(decyear, np.ndarray):
        [mday,month,ut]=np.array(list(map(decyear2cal,decyear))).T
    
    else:

        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)
        (mday,month)=dayno2cal(noday,year)

    # this trick ensures that a single decyear provided as float
    # returns int,int,float for monthday, month and ut
    # while a list or np.array of decyear will return
    # int 1D np array for monthday and month and float 1D np array for ut
    
    return( np.array(mday,dtype=int)+0,
            np.array(month,dtype=int)+0,
            ut)

###############################################################################
def decyear2dayno(decyear):
###############################################################################
    """
    converts decimal year to the day of year
    
    :param decyear: decimal year
    :returns: day of year,ut
    :note: the input decimal year is assumed to account for leap years, that is the day of year is the \
    decimal of decyear * number of days in the year.

    """

    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)
    
    if isinstance(decyear, np.ndarray):
        [noday,ut]=np.array(list(map(decyear2dayno,decyear))).T

    else:
    
        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)

    return( np.array(noday, dtype=int)+0,ut)

###############################################################################
def decyear2mjd(decyear):
###############################################################################
    """
    decimal year to modified julian day conversion
    
    :param decyear: decimal year
    :returns: modified julian day
    """

    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array(decyear)
    
    if isinstance(decyear, np.ndarray):
        mjd=np.array(list(map(decyear2mjd,decyear)))

    else:

        year=int(decyear)
        frac_year=decyear-year
        nday=cal2dayno(31,12,year)
        doy=frac_year*float(nday)
        noday=int(doy)+1
        ut=doy-int(doy)
        mjd=dayno2mjd(noday,year,ut)
        
    return mjd

###############################################################################
def decyear2datetime(decyear):
###############################################################################
    """
    converts a decimal date to a datime python instance
    
    :param decyear: decimal year
    :returns: datetime instance

    """

    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array( decyear )

    if isinstance(decyear, np.ndarray):
        return_datetime=np.array(list(map(decyear2datetime,decyear)))

    else:

        year=int(decyear)
        (dayno,ut)=decyear2dayno(decyear)
        return_datetime=dayno2datetime(dayno,year,ut=ut)

    return return_datetime

###############################################################################
def decyear2epoch(decyear):
###############################################################################
    """
    converts decimal year to formatted epoch %02d:%03d:%05d
    
    :param decyear: decimal year
    :returns : 'yy:doy:sec'
    
    :note: the input decimal year is assumed to account for leap years, that is the day of year is the \
    decimal of decyear * number of days in the year.
    
    :Example:
    >>> import pyacs
    >>> pyacs.cal2decyear(1,3,2001)
    >>> 2001.16301369863
    >>> pyacs.dec
    >>> pyacs.decyear2epoch(2001.16301369863)
    >>> '01:060:43200'
    >>> pyacs.decyear2epoch(2000.16301369863)
    >>> '00:060:57284'
    >>> pyacs.decyear2epoch(pyacs.cal2decyear(1,3,2000))
    >>> '00:061:43200'
    
    """
    
    import numpy as np

    if isinstance(decyear, list):
        decyear=np.array( decyear )

    if isinstance(decyear, np.ndarray):
        epoch=np.array(list(map(decyear2epoch,decyear)))

    else:
    
        (doy,ut)=decyear2dayno(decyear)
        seconds=ut*(60*60*24)
        yr=year2yr(decyear)
        epoch=("%02d:%03d:%05d" %(yr,doy,round(seconds)))

    return(epoch)



###############################################################################
# DATETIME CONVERSIONS
###############################################################################

###############################################################################
def datetime_from_calarray(calarray, hour=12, minute=0, sec=0,ms=0):
###############################################################################
    """
    Build a datetime array from an array of calendar dates
    
    :param calarray: 1D or 2D numpy array. Columns order are year,month,mday,hour,minute,second,microsecond
    
    :return: numpy 1D array of datime instances
    :note: year,month,mday are mandatory. If not provided then hour=12, minute=0, sec=0,ms=0 (middle of the day) 
    """

    #import
    import numpy as np
    from datetime import datetime

    # convert to numpy array if a list of list was provided
    calarray = np.array( calarray )
    

    # case array is 1D
    if calarray.ndim == 1:
        calarray = calarray.reshape(-1,calarray.size)
        
    # prepare a formatted array
    
    fmt_array = np.zeros( (calarray.shape[0] , 7) , dtype = int )
    
    # default value for hour
    
    fmt_array[:,3] = 12
    
    # fill the year, month, mday columns
    
    fmt_array[:,:3] = calarray[:,:3].astype(int)
    
    # if calarray is not fully integer, then we assume that microseconds are not provided, but seconds are decimal
    if  calarray.dtype is not np.dtype('int'):
        # if there are 7 columns, then microseconds are provided and dtype was accidentally of float type
        if calarray.shape[1] == 7:
            fmt_array = calarray.astype(int)
        # else, the decimal part of seconds are converted to microseconds
        else:
            seconds = calarray[:,5].astype(int)
            microseconds = ( np.modf(calarray[:,5])[0] * 1.E6).astype(int) 
            
            fmt_array[:,3] = calarray[:,3].astype(int)
            fmt_array[:,4] = calarray[:,4].astype(int)
            fmt_array[:,5] = seconds
            fmt_array[:,6] = microseconds
            
    return( np.array(list(map(datetime,fmt_array[:,0],fmt_array[:,1],fmt_array[:,2],fmt_array[:,3],fmt_array[:,4],fmt_array[:,5],fmt_array[:,6]))) )
        

###############################################################################
def datetime2cal(datetime):
###############################################################################
    """
    converts from python datetime to calendar date
    
    :param datetime: datetime instance
    :returns: year,month,monthday,ut
    
    """
    
    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)

    if isinstance(datetime, np.ndarray):
        [year,month,mday,ut]=np.array(list(map(datetime2cal,datetime))).T

    else:
    
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        uts=hmsmicros2uts(datetime.hour,datetime.minute,datetime.second,datetime.microsecond)
        ut=uts2ut(uts)

    return( np.array(year, dtype=int)+0,\
            np.array(month, dtype=int)+0,\
            np.array(mday, dtype=int)+0,\
            ut )


###############################################################################
def datetime2dayno(datetime):
###############################################################################
    """
    Converts python datetime to dayno,year,ut
    
    :param datetime: pythob datetime instance
    :returns: year,dayno,ut
    
    """

    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)
    
    if isinstance(datetime, np.ndarray):
        [year,dayno,ut]=np.array(list(map(datetime2dayno,datetime))).T
    
    else:
    
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        uts=hmsmicros2uts(datetime.hour,datetime.minute,datetime.second,datetime.microsecond)
        ut=uts2ut(uts)
        dayno=cal2dayno(mday, month, year)

    return( np.array(year , dtype=int)+0,\
            np.array(dayno, dtype=int)+0,\
            ut)

###############################################################################
def datetime2mjd(datetime):
###############################################################################
    """
    converts a python datetime instance to a modified julian day
    
    :param datetime: python datetime instance
    :returns: modified julian day
    
    """

    import numpy as np

    if isinstance(datetime, list):
        datetime=np.array(datetime)
    
    if isinstance(datetime, np.ndarray):
        mjd=np.array(list(map(datetime2mjd,datetime)))
    
    else:

        
        year=datetime.year
        month=datetime.month
        mday=datetime.day
        h=datetime.hour
        mn=datetime.minute
        s=datetime.second
        microsecond=datetime.microsecond
        uts=hmsmicros2uts(h,mn,s=s,microsecond=microsecond)
        ut=uts2ut(uts)
        mjd=cal2mjd(mday, month, year, ut=ut)
        
    return mjd


###############################################################################
def datetime2decyear(datetime):
###############################################################################
    """
    converts a python datetime instance to decimal year
    
    :param datetime: python datetime instance or array of it
    :returns: decimal year
    
    """

    return mjd2decyear( datetime2mjd(datetime) )


###############################################################################
def datetime_round_second(my_datetime):
###############################################################################
    """
    rounds a python datetime instance to the nearest second
    
    :param datetime: python datetime instance or array of it
    :returns: new rounded datetime
    
    """

    import datetime
    

    if my_datetime.microsecond <= 500000:
        my_datetime = my_datetime.replace( microsecond = 0 )
    else:
        my_datetime = my_datetime.replace( microsecond = 0 )
        deltatime = datetime.timedelta(seconds=1)
        my_datetime = my_datetime + deltatime
    
    return(my_datetime)
 

###############################################################################
# GPS WEEK CONVERSIONS
###############################################################################

###############################################################################
def gpsweek2mjd(gpsweek,dow):
###############################################################################

    """
    converts a GPS week and dow to modified julian day
    
    :param gpsweek,dow: gps week, day of week (0 is sunday, 6 saturday)
    :returns: modified julian day
    
    """

    import numpy as np

    if isinstance(gpsweek, list):
        gpsweek=np.array(gpsweek)
        dow=np.array(dow)


    if isinstance(gpsweek, np.ndarray):
        mjd=np.array(list(map(gpsweek2mjd,gpsweek,dow))).T
    
    else:
    
        mjd010580 = 44243
        mjd = ( gpsweek * 7 ) + dow + mjd010580 + 1

    return( mjd )


###############################################################################
# EPOCH CONVERSIONS
###############################################################################

###############################################################################
def epoch2decyear(epoch):
###############################################################################
    """
    Converts an epoch string used in the SINEX format 'YY:DOY:SOD' e.g.'17:100:43185' into a decimal year
    
    :param epoch: epoch formatted string or list of epochs
    
    :returns: decimal year or array of decimal years
    
    """

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

###############################################################################
# JULIAN DAY, LEAP, DAY FRACTIONS, HOUR/MIN/SEC/MICROSECONDS CONVERSIONS
###############################################################################

###############################################################################
def jd2mjd(jd):
###############################################################################
    """
    Converts a Julian day to a Modified Julian day
    
    :param jd: Julian day
    
    :returns: jd-2400000.5
    
    """

    return jd-2400000.5

###############################################################################
def mjd2jd(mjd):
###############################################################################
    """

    converts a Modified Julian day to Julian day
    
    :param mjd:     Modified Julian day
    
    :returns: mjd + 2400000.5
    
    """
    
    return mjd+2400000.5


###############################################################################
def yr2year(yr):
###############################################################################
    """
    converts two digits year to four digits year
    
    :param yr: 2-digits year (e.g. 01)
    :returns: 4-digit years (2001). This is yr+2000 if yr<80 yr+1900 otherwise.
    
    :note: This is an heritage from a bad practice in the geodesy community assuming that origin of the universe started with GPS. See example.

    :Example:
    >>> import pyacs
    >>> pyacs.yr2year(96)
    >>> 1996
    >>> pyacs.yr2year(11)
    >>> 2011
    """

    import numpy as np

    if isinstance(yr, list):
        yr=np.array(yr)

    if isinstance(yr, np.ndarray):
        year=np.array(list(map(yr2year,yr)))

    else:
    
        if yr<80:year=yr+2000
        else:year=yr+1900

    return(year)

###############################################################################
def year2yr(year):
###############################################################################
    """
    converts 4-digits year to 2-digits year

    :param yr: 4-digits year (e.g. 2001)
    :returns: 2-digit years (01). This is yr-2000 if yr<80 yr-1900 otherwise.
    
    :note: This is an heritage from a bad practice in the geodesy community assuming that origin of the universe started with GPS in 1980. See example.

    :Example:
    >>> import pyacs
    >>> pyacs.year2yr(1996)
    >>> 96
    >>> pyacs.year2yr(2011)
    >>> 11

    """

    import numpy as np

    if isinstance(year, list):
        year=np.array(year)

    if isinstance(year, np.ndarray):
        yr=np.array(list(map(year2yr,year)))

    else:
    
        
        if year > 1999:yr=year-2000
        else:yr=year-1900

    return(yr)

###############################################################################
def leap_year(year):
###############################################################################
    
    """
    Returns true if year is a leap year.
    
    :param year: year in YYYY
    
    :returns: True if year is leap, False otherwise 
    """

    import numpy as np

    if isinstance(year, np.ndarray):
        OK=np.array(list(map(leap_year,year))).T

    else:
    
        year=int(year)
        OK=(( not (year%4))and(year%100))or(not(year%400))
        
    return OK

###############################################################################
def uts2hmsmicros(uts):
###############################################################################
    """
    converts decimal seconds on a day (0-86400) to hours, min, seconde, microseconds
    
    :param uts: decimal seconds since 00:00:00.0. Must be in the [0.,86400.[ range.
    :returns: hours,minutes,seconds,microseconds
    :rtype: int,int,int,float
    :Example:
    >>> import pyacs
    >>> pyacs.uts2hmsmicros(86399.999999)
    >>> (23, 59, 59, 999999.0000069374)
    >>> pyacs.uts2hmsmicros(86400)
    >>> !!!  86400 uts out of range [0- 86400.0 [
    
    """

    import numpy as np
    
    if isinstance(uts, list):
        uts=np.array(uts)

    if isinstance(uts, np.ndarray):
        [h,mn,s,microsecond]=np.array(list(map(uts2hmsmicros,uts))).T

    else:

        __utsOK(uts)
        
        nsec_hour= 60 * 60
        nsec_min= 60
        
        int_sec=int(uts)
        
        h = int_sec // nsec_hour
        r = int_sec % nsec_hour

        mn = r // nsec_min
        s =  r % nsec_min

        microsecond= (uts - float(int_sec)) * 1.E6

        return( np.array(h, dtype=int)+0,\
                np.array(mn, dtype=int)+0,\
                np.array(s,dtype=int)+0,\
                microsecond)


###############################################################################
def hmsmicros2uts(h,mn,s=0,microsecond=0):
###############################################################################
    """
    converts hour minute seconds to uts
    
    :param hours,minutes,seconds,microseconds
    :returns: uts
    """

    nsec_hour=60 * 60
    nsec_min=60

    return((h*nsec_hour+mn*nsec_min+s)+(microsecond)/1.E6)    

###############################################################################
def hmsmicros2ut(h,mn,s=0,microsecond=0):
###############################################################################
    """
    converts hour minute seconds to fractional day
    
    :param hours,minutes,seconds,microseconds
    :returns: ut
    """

    uts=hmsmicros2uts(h,mn,s=s,microsecond=microsecond)

    return(uts2ut(uts))    


###############################################################################
def ut2uts(ut):
###############################################################################

    """
    converts fractional day to seconds
    
    :param ut: fractional day in [0.,1.[
    :returns: ut * 60.*60.*24.
    """
    
    nsec_day=60.*60.*24.
    return(ut*nsec_day)

###############################################################################
def uts2ut(uts):
###############################################################################
    """
    converts uts in seconds to fractional day
    
    :param uts: seconds in [0.,86400.[
    :returns: uts / (60.*60.*24.)
    """
    
    import numpy as np

    if isinstance(uts, list):
        uts=np.array(uts)
    
    if isinstance(uts, np.ndarray):
        ut=np.array(list(map(uts2ut,uts)))

    else:
        
        __utsOK(uts)
        nsec_day=60.*60.*24.
        ut=uts/nsec_day
    
    return ut


###############################################################################
def day_since_decyear(decyear,ref_decyear):
###############################################################################
    """
    Calculates the number of days since a reference date.
    
    :param decyear: decimal year
    :returns: days elapsed since decyear (float)
    
    :note: negative number means before the reference date
    :note: a useful function for postseismic analysis.
    
    :Example:
    >>> import pyacs
    >>> ref_decyear=pyacs.cal2decyear(16,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 33)) # Ecuador, Mw 7.8 EQ
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(17,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 32)),ref_decyear)
    >>> 0.9999884259232203
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(17,4,2016,ut=pyacs.hmsmicros2ut(23, 58, 33)),ref_decyear)
    >>> 0.9999999999854481
    >>> pyacs.day_since_decyear(pyacs.cal2decyear(16,4,2017,ut=pyacs.hmsmicros2ut(23, 58, 33)),ref_decyear)
    >>> 364.9999999999345
    >>> ld=pyacs.mjd2decyear(np.arange(pyacs.cal2mjd(13,4,2016),pyacs.cal2mjd(19,4,2016)))
    >>> ld
    >>> array([ 2016.283,  2016.286,  2016.288,  2016.291,  2016.294,  2016.296])
    >>> pyacs.day_since_decyear(ld,ref_decyear)
    >>> array([-3.499, -2.499, -1.499, -0.499,  0.501,  1.501])

    """

    mjd=decyear2mjd(decyear)
    
    ref_mjd=decyear2mjd(ref_decyear)

    return mjd-ref_mjd

###############################################################################
def guess_date(date):
###############################################################################
    """
    Tries to guess a date from input parameter
    returns a decimal year

    :param date: date as string or float
    :returns: decimal year

    :Example:
    >>> import pyacs
    >>> pyacs.guess_date('2016/04/16')
    >>> 2016.2909836065573


    """
    
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

##############################################################################
def decyear2seconds( np_decyear , rounding='day' ):
##############################################################################
    """
    converts decimal year to a seconds since 1980/1/1/0/0/0
    
    :param rounding: controls rounding. 'day' means at 12:00:00, 
    'hour' at 00:00, 'minute' at the current minute with 00 seconds, 'second' at the integer of the current second.
    :return: 1D integer (np.int64) numpy array 
    """

    # import  
    import numpy as np
    from datetime import datetime, timedelta
    
    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)
    
    np_datetime = decyear2datetime( np_decyear )

    if not( isinstance(np_datetime, np.ndarray) ):
        np_datetime = np.array([ np_datetime ])

    # handle round
    
    for i in np.arange( np_datetime.shape[0] ):
        
        if rounding == 'day':
            np_datetime[i] = np_datetime[i].replace(hour= 12, minute = 0, second = 0, microsecond = 0)
        if rounding == 'hour':
            if np_datetime[i].minute >= 30:
                np_datetime[i] = np_datetime[i] + timedelta( minutes=30 ) 
            np_datetime[i] = np_datetime[i].replace( minute = 0, second = 0, microsecond = 0)
        if rounding == 'minute':
            if np_datetime[i].second >= 30:
                np_datetime[i] = np_datetime[i] + timedelta( seconds=30 ) 
            np_datetime[i] = np_datetime[i].replace( second = 0, microsecond = 0)
        if rounding == 'second':
            np_datetime[i] = np_datetime[i].replace( microsecond = 0)

    
    # convert to timedelta and then to seconds

    np_dates_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)

    # return
    
    return( np_dates_s ) 

##############################################################################
def datetime2seconds( np_datetime ):
##############################################################################
    """
    converts dattime to seconds since 1980/1/1/0/0/0
    :param datetime:
    :return: 1D numpy array of seconds
    """


    # import
    import numpy as np
    from datetime import datetime, timedelta

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)


    if not (isinstance(np_datetime, np.ndarray)):
        np_datetime = np.array([np_datetime])


    # convert to timedelta and then to seconds

    np_dates_s = np.array(list(map(int, [x.total_seconds() for x in (np_datetime - ref_date_time)])), dtype=np.int64)

    # return

    return (np_dates_s)


##############################################################################
def seconds2datetime( seconds ):
##############################################################################
    """
    converts seconds since 1980/1/1/0/0/0 to a datetime object
    
    :return:a 1D numpy array of datetime object
    """
    # import  
    import numpy as np
    from datetime import datetime , timedelta
    
    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)

    if isinstance( seconds , np.ndarray):
        np_datetime =   np.array(list(map(lambda x: ref_date_time + timedelta(seconds=x) ,  seconds.astype(int).tolist() ))) 
    else:
        np_datetime =   ref_date_time + timedelta(seconds=int(seconds)) 
    return np_datetime
    

