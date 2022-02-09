#-------------------------------------------------------------------------------
# Module   : date
# Purpose  : Handle dates in various formats
# Author   : P. Rebischung
# Created  : 22-May-2011
#
# Changes  :
#
# Routines : - __init__      : Initialize a date object
#            - __str__       : Convert a date object to a string
#            - add_s         : Add (or substract) n seconds to a date object
#            - add_h         : Add (or substract) n hours to a date object
#            - add_d         : Add (or substract) n days to a date object
#            - from_mjd      : Create a date object from a modified Julian date
#            - from_ymdhms   : Create a date object from year, month, day, hour, minute and second
#            - from_weekdow  : Create a date object from GPS week and day of week
#            - from_snxepoch : Create a date object from a date in SINEX format
#            - from_ydec     : Create a date object from a decimal year
#            - from_wdec     : Create a date object from a decimal GPS week
#            - snxepoch      : Write a date object in SINEX date format
#            - ydec          : Compute decimal year of a date object
#            - isoepoch      : Write a date object in ISO date format
#            - wdec          : Compute decimal GPS week of a date object
#-------------------------------------------------------------------------------

# LIBRARIES
#-------------------------------------------------------------------------------
from time import *
import datetime
import calendar

from .constants import *

# CONSTANTS
#-------------------------------------------------------------------------------
gps0 = calendar.timegm(datetime.datetime(1980, 1, 6, 0, 0, 0).timetuple())
j2000 = calendar.timegm(datetime.datetime(2000, 1, 1, 12, 0, 0).timetuple())


# DATE CLASS
#-------------------------------------------------------------------------------
class date:

#-------------------------------------------------------------------------------
# Routine : __init__
# Purpose : Initialize a date object
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : tsys : System time (seconds since 1970-01-01 00:00:00). Default is "now".
# Output  : t    : date object
#-------------------------------------------------------------------------------
    def __init__(self, *args):
    
        t=self
        # If no argument, t = now!
        if (len(args) == 0):
            t.tsys = time()
    
        # If one argument, it should be the system time
        elif (len(args) == 1):
            t.tsys = args[0]
    
        # Intermediate struct_time object
        st = gmtime(t.tsys)
    
        # Define "standard" attributes of the date object
        t.yyyy = strftime('%Y', st)
        t.yy = strftime('%y', st)
        t.mm = strftime('%m', st)
        t.dd = strftime('%d', st)
        t.hour = strftime('%H', st)
        t.min = strftime('%M', st)
        t.sec = strftime('%S', st)
        t.doy = strftime('%j', st)
        t.dow = strftime('%w', st)
    
        # GPS week
        t.week = int((t.tsys - gps0) / 86400 / 7)
        t.week = '{0:0>4}'.format(t.week)
    
        # Week of year
        t.wk = '{0:0>2}'.format(int(float(t.doy) / 7) + 1)
    
        # MJD
        t.mjd = mjd2000 + float(t.tsys - j2000) / 86400

#-------------------------------------------------------------------------------
# Routine : __str__
# Purpose : Convert a date object to a string (for print)
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : 
# Output  : String in the format 'YYYY-MM-DD HH:MM:SS'
#-------------------------------------------------------------------------------
    def __str__(self):
        t=self
        return t.yyyy + '-' + t.mm + '-' + t.dd + ' ' + t.hour + ':' + t.min + ':' + t.sec

#-------------------------------------------------------------------------------
# Routine : add_s
# Purpose : Add (or substract) n seconds to a date object
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : n : Number of seconds to add (or to substract if n < 0)
# Output  :
#-------------------------------------------------------------------------------
    def add_s(self, n):
        t=self
        t.__init__(t.tsys + n)

#-------------------------------------------------------------------------------
# Routine : add_h
# Purpose : Add (or substract) n hours to a date object
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : n : Number of hours to add (or to substract if n < 0)
# Output  :
#-------------------------------------------------------------------------------
    def add_h(self, n):
        t=self
        t.__init__(t.tsys + 60 * n)

#-------------------------------------------------------------------------------
# Routine : add_d
# Purpose : Add (or substract) n days to a date object
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : n : Number of days to add (or to substract if n < 0)
# Output  :
#-------------------------------------------------------------------------------
    def add_d(self, n):
        t=self
        t.__init__(t.tsys + 86400 * n)

#-------------------------------------------------------------------------------
# Routine : from_mjd
# Purpose : Create a date object from a modified Julian date
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : d : Modified Julian date
# Output  : t : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_mjd(cls, d):
        tsys = j2000 + (d - mjd2000) * 86400
        return date(tsys)

#-------------------------------------------------------------------------------
# Routine : from_ymdhms
# Purpose : Create a date object from year, month, day, hour, minute and second
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : - y    : Year
#           - m    : Month
#           - d    : Day
#           - hour : Hour.   Default is 0.
#           - min  : Minute. Default is 0.
#           - sec  : Second. Default is 0.
# Output  : - t    : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_ymdhms(cls, y, m, d, hour=0, min=0, sec=0):
        tsys = calendar.timegm(datetime.datetime(y, m, d, hour, min, sec).timetuple())
        return date(tsys)

#-------------------------------------------------------------------------------
# Routine : from_weekdow
# Purpose : Create a date object from GPS week and day of week
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : - w : GPS week
#           - d : Day in week
# Output  : - t : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_weekdow(cls, w, d):
        tsys = gps0 + (7 * w + d) * 86400
        return date(tsys)

#-------------------------------------------------------------------------------
# Routine : from_snxepoch
# Purpose : Create a date object from a date in SINEX format
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : s : Date in SINEX format (e.g. '08:265:73400')
# Output  : t : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_snxepoch(cls, s):
        tsys = calendar.timegm(strptime(s[0:6], '%y:%j'))
        tsys += int(s[7:12])
        return date(tsys)

#-------------------------------------------------------------------------------
# Routine : from_ydec
# Purpose : Create a date object from a decimal year
# Author  : P. Rebischung
# Created : 20-Sep-2013
#
# Changes :
#
# Input   : y : Decimal year
# Output  : t : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_ydec(cls, y):
        yint = int(y)
        yfra = y - yint
        t = date.from_ymdhms(yint, 1, 1, 0, 0)
        if (((y % 4 == 0) and (y % 100 != 0)) or (y % 400 == 0)):
          t.add_d(yfra * 366)
        else:
          t.add_d(yfra * 365)
        return t
  
#-------------------------------------------------------------------------------
# Routine : from_wdec
# Purpose : Create a date object from a decimal GPS week
# Author  : P. Rebischung
# Created : 22-Oct-2014
#
# Changes :
#
# Input   : w : Decimal GPS week
# Output  : t : date object
#-------------------------------------------------------------------------------
    @classmethod
    def from_wdec(cls, w):
        wint = int(w)
        wfra = w - wint
        t = date.from_weekdow(wint, wfra * 7)
        return t

#-------------------------------------------------------------------------------
# Routine : snxepoch
# Purpose : Write a date object in SINEX date format
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   :
# Output  : s : Date in SINEX format (e.g. '08:265:73400')
#-------------------------------------------------------------------------------
    def snxepoch(self):
        t=self
        sec = 3600 * int(t.hour) + 60 * int(t.min) + int(t.sec)
        return t.yy + ':' + t.doy + ':' + '{0:0>5}'.format(sec)

#-------------------------------------------------------------------------------
# Routine : ydec
# Purpose : Compute decimal year of a date object
# Author  : P. Rebischung
# Created : 09-Feb-2012
#
# Changes :
#
# Input   :
# Output  : y : Decimal year
#-------------------------------------------------------------------------------
    def ydec(self):
            t=self
            y = float(t.yyyy)
            if (((y % 4 == 0) and (y % 100 != 0)) or (y % 400 == 0)):
                return y + (float(t.doy) - 1 + (float(t.hour) + (float(t.min) + float(t.sec) / 60.) / 60.) / 24.) / 366.
            else:
                return y + (float(t.doy) - 1 + (float(t.hour) + (float(t.min) + float(t.sec) / 60.) / 60.) / 24.) / 365.

#-------------------------------------------------------------------------------
# Routine : isoepoch
# Purpose : Write a date object in ISO date format
# Author  : P. Rebischung
# Created : 27-Aug-2012
#
# Changes :
#
# Input   :
# Output  : s : Date in ISO format (e.g. '08-12-31T23:59:59')
#-------------------------------------------------------------------------------
    def isoepoch(self):
        t=self
        return t.yyyy + '-' + t.mm + '-' + t.dd + 'T' + t.hour + ':' + t.min + ':' + t.sec

#-------------------------------------------------------------------------------
# Routine : wdec
# Purpose : Compute decimal GPS week of a date object
# Author  : P. Rebischung
# Created : 27-Aug-2012
#
# Changes :
#
# Input   :
# Output  : w : Decimal GPS week
#-------------------------------------------------------------------------------
    def wdec(self):
        t=self
        return float(t.week) + (float(t.dow) + (float(t.hour) + (float(t.min) + float(t.sec) / 60.) / 60.) / 24.) / 7.
