"""
Note: Leap seconds still needs to be checked an improve
A Python implementation of GPS related time conversions.

Copyright 2002 by Bud P. Bruegger, Sistema, Italy
mailto:bud@sistema.it
http://www.sistema.it

Modifications for GPS seconds by Duncan Brown

PyUTCFromGpsSeconds added by Ben Johnson

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc., 59
Temple Place, Suite 330, Boston, MA  02111-1307  USA

GPS Time Utility functions

This file contains a Python implementation of GPS related time conversions.

The two main functions convert between UTC and GPS time (GPS-week, time of
week in seconds, GPS-day, time of day in seconds).  The other functions are
convenience wrappers around these base functions.  

A good reference for GPS time issues is:
http://www.oc.nps.navy.mil/~jclynch/timsys.html

Note that python time types are represented in seconds since (a platform
dependent Python) Epoch.  This makes implementation quite straight forward
as compared to some algorigthms found in the literature and on the web.  
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date: 2006/02/16 04:36:09 $'
__version__ = '$Revision: 1.6 $'[11:-2]

import time, math

secsInWeek = 604800
secsInDay = 86400
gpsEpoch = (1980, 1, 6, 0, 0, 0)  # (year, month, day, hh, mm, ss)

def dayOfWeek(year, month, day):
    """Return GPS day of week: 0=Sun, 1=Mon, ..., 6=Sat.

    Parameters
    ----------
    year, month, day : int
        Date components.

    Returns
    -------
    int
        Day of week (0-6).
    """
    hr = 12  #make sure you fall into right day, middle is save
    t = time.mktime((year, month, day, hr, 0, 0.0, 0, 0, -1))
    pyDow = time.localtime(t)[6]
    gpsDow = (pyDow + 1) % 7
    return gpsDow

def gpsWeek(year, month, day):
    """Return full GPS week number for given UTC date.

    Parameters
    ----------
    year, month, day : int
        UTC date.

    Returns
    -------
    int
        GPS week number.
    """
    hr = 12  #make sure you fall into right day, middle is save
    #print "toto",year
    #return year
    #print year, month, day, hr, 0, 0.0
    return gpsFromUTC(year, month, day, hr, 0, 0.0)[0]


def julianDay(year, month, day):
    """Return day of year (1-366) for given date.

    Parameters
    ----------
    year, month, day : int
        Date components.

    Returns
    -------
    int
        Day since Jan 1 of year.
    """
    hr = 12  #make sure you fall into right day, middle is save
    t = time.mktime((year, month, day, hr, 0, 0.0, 0, 0, -1))
    julDay = time.localtime(t)[7]
    return julDay

def mkUTC(year, month, day, hour, minute, sec):
    """Convert UTC date/time to seconds since epoch (like mktime but for UTC).

    Parameters
    ----------
    year, month, day, hour, minute, sec : int or float
        UTC components.

    Returns
    -------
    float
        Seconds since epoch.
    """
    spec = [year, month, day, hour, minute, sec] + [0, 0, 0]
    utc = time.mktime(spec) - time.timezone
    return utc

def ymdhmsFromPyUTC(pyUTC):
    """Return (year, month, day, hour, minute, sec) from Python UTC time.

    Parameters
    ----------
    pyUTC : float
        Seconds since epoch (e.g. from time.time()).

    Returns
    -------
    tuple
        (year, month, day, hour, minute, sec).
    """
    ymdhmsXXX = time.gmtime(pyUTC)
    return ymdhmsXXX[:-3]

def wtFromUTCpy(pyUTC, leapSecs=14):
    """Convert Python UTC time to GPS week and time of week.

    Parameters
    ----------
    pyUTC : float
        Seconds since epoch (UTC).
    leapSecs : int, optional
        Leap seconds offset. Default is 14.

    Returns
    -------
    tuple
        (gpsWeek, secsOfWeek).
    """
    ymdhms = ymdhmsFromPyUTC(pyUTC)
    wSowDSoD = gpsFromUTC(*ymdhms + (leapSecs,))
    return wSowDSoD[0:2]

def gpsFromUTC(year, month, day, hour, minute, ssec, leapSecs=30):
    """Convert UTC to GPS week, seconds of week, GPS day, seconds of day.

    Parameters
    ----------
    year, month, day, hour, minute : int
        UTC date/time components.
    ssec : float
        Seconds (and fraction).
    leapSecs : int, optional
        Leap seconds (GPS - UTC). Default is 30; update as leap seconds change.

    Returns
    -------
    tuple
        (gpsWeek, secsOfWeek, gpsDay, secsOfDay).

    Notes
    -----
    GPS time is seconds since 1980-01-06 00:00:00. Week starts Saturday midnight.
    See http://www.oc.nps.navy.mil/~jclynch/timsys.html. Python uses integer
    seconds; fractional seconds in ssec are preserved in the return value.
    """
    secFract = ssec % 1
    sec=int(math.floor(ssec))
    epochTuple = gpsEpoch + (-1, -1, 0)
    t0 = time.mktime(epochTuple)
    
    t = time.mktime((year, month, day, hour, minute, sec, -1, -1, 0))
    
    # Note: time.mktime strictly works in localtime and to yield UTC, it should be
    #       corrected with time.timezone
    #       However, since we use the difference, this correction is unnecessary.
    # Warning:  trouble if daylight savings flag is set to -1 or 1 !!!
    t = t + leapSecs   
    tdiff = t - t0
    gpsSOW = (tdiff % secsInWeek)  + secFract
    gpsWeek = int(math.floor(tdiff/secsInWeek)) 
    gpsDay = int(math.floor(gpsSOW/secsInDay))
    gpsSOD = (gpsSOW % secsInDay) 
    return (gpsWeek, gpsSOW, gpsDay, gpsSOD)


def UTCFromGps(gpsWeek, SOW, leapSecs=14):
    """Convert GPS week and seconds of week to UTC.

    Parameters
    ----------
    gpsWeek : int
        Full GPS week number (not modulo 1024).
    SOW : float
        Seconds of week.
    leapSecs : int, optional
        Leap seconds (GPS - UTC). Default is 14.

    Returns
    -------
    tuple
        (year, month, day, hour, minute, sec).
    """
    secFract = SOW % 1
    epochTuple = gpsEpoch + (-1, -1, 0) 
    t0 = time.mktime(epochTuple) - time.timezone  #mktime is localtime, correct for UTC
    tdiff = (gpsWeek * secsInWeek) + SOW - leapSecs
    t = t0 + tdiff
    (year, month, day, hh, mm, ss, dayOfWeek, julianDay, _daylightsaving) = time.gmtime(t)
    #use gmtime since localtime does not allow to switch off daylighsavings correction!!!
    return (year, month, day, hh, mm, ss + secFract)

def GpsSecondsFromPyUTC(pyUTC, _leapSecs=14):
    """Convert Python epoch time to GPS seconds (since GPS epoch).

    Parameters
    ----------
    pyUTC : float
        Seconds since Python epoch (e.g. time.time()).
    _leapSecs : int, optional
        Unused; kept for API compatibility.

    Returns
    -------
    int
        GPS seconds since 1980-01-06 00:00:00.
    """
    t = gpsFromUTC(*ymdhmsFromPyUTC( pyUTC ))
    return int(t[0] * 60 * 60 * 24 * 7 + t[1])

#def PyUTCFromGpsSeconds(gpsseconds):
#    """converts gps seconds to the
#    python epoch. That is, the time
#    that would be returned from time.time()
#    at gpsseconds.
#    """
#    pyUTC
    
#===== Tests  =========================================

def testTimeStuff():
    print("-"*20)
    print()
    print("The GPS Epoch when everything began (1980, 1, 6, 0, 0, 0, leapSecs=0)")
    (w, sow, d, sod) = gpsFromUTC(1980, 1, 6, 0, 0, 0, leapSecs=0)
    print("**** week: %s, sow: %s, day: %s, sod: %s" % (w, sow, d, sod))
    print("     and hopefully back:")
    print("**** %s, %s, %s, %s, %s, %s\n" % UTCFromGps(w, sow, leapSecs=0))

    print("The time of first Rollover of GPS week (1999, 8, 21, 23, 59, 47)")
    (w, sow, d, sod) = gpsFromUTC(1999, 8, 21, 23, 59, 47)
    print("**** week: %s, sow: %s, day: %s, sod: %s" % (w, sow, d, sod))
    print("     and hopefully back:")
    print("**** %s, %s, %s, %s, %s, %s\n" % UTCFromGps(w, sow, leapSecs=14))

    print("Today is GPS week 1186, day 3, seems to run ok (2002, 10, 2, 12, 6, 13.56)")
    (w, sow, d, sod) = gpsFromUTC(2002, 10, 2, 12, 6, 13.56)
    print("**** week: %s, sow: %s, day: %s, sod: %s" % (w, sow, d, sod))
    print("     and hopefully back:")
    print("**** %s, %s, %s, %s, %s, %s\n" % UTCFromGps(w, sow))

def testJulD():
    print('2002, 10, 11 -> 284  ==??== ', julianDay(2002, 10, 11))

def testGpsWeek():
    print('2002, 10, 11 -> 1187  ==??== ', gpsWeek(2002, 10, 11))

def testDayOfWeek():
    print('2002, 10, 12 -> 6  ==??== ', dayOfWeek(2002, 10, 12))
    print('2002, 10, 6  -> 0  ==??== ', dayOfWeek(2002, 10, 6))

def testPyUtilties():
    ymdhms = (2002, 10, 12, 8, 34, 12.3)
    print("testing for: ", ymdhms)
    pyUtc = mkUTC(*ymdhms)
    back =  ymdhmsFromPyUTC(pyUtc)
    print("yields     : ", back)
#*********************** !!!!!!!!    
    #assert(ymdhms == back)
    #! TODO: this works only with int seconds!!! fix!!!
    (w, t) = wtFromUTCpy(pyUtc)
    print("week and time: ", (w,t))


#===== Main =========================================
if __name__ == "__main__":
    pass
    testTimeStuff()
    testGpsWeek()
    testJulD()
    testDayOfWeek()
    testPyUtilties()
