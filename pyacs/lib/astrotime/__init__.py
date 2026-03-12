"""
AstroTime: time conversions (calendar, MJD, decimal year, datetime).

Traduction from perl module Astro::Time (CPAN, Chris Phillips).
Conversion to Python by Jean-Mathieu Nocquet (Geoazur).
Routines for calendar dates, Modified Julian day, UT/local sidereal time,
and numerical/string angles. String functions removed. Datetime support added.

Notes
-----
This module is for date manipulation; leap seconds are not handled.
See the GPSTime module for GPS/UTC conversions.
In Python 3.6+, / was changed to // for integer operations.
All conversions not specifying ut or uts use 12h00mn00s as default; ut is the
decimal fraction of day, uts is seconds since 00h00m00.0s.
"""

# Import all public functions to maintain backward compatibility
from .cal2dayno import cal2dayno
from .cal2mjd import cal2mjd
from .cal2decyear import cal2decyear
from .cal2datetime import cal2datetime
from .dayno2cal import dayno2cal
from .dayno2mjd import dayno2mjd
from .dayno2decyear import dayno2decyear
from .dayno2datetime import dayno2datetime
from .mjd2cal import mjd2cal
from .mjd2dayno import mjd2dayno
from .mjd2decyear import mjd2decyear
from .mjd2datetime import mjd2datetime
from .mjd2gpsweek import mjd2gpsweek
from .gpsweek2mjd import gpsweek2mjd
from .decyear2cal import decyear2cal
from .decyear2dayno import decyear2dayno
from .decyear2mjd import decyear2mjd
from .decyear2datetime import decyear2datetime
from .decyear2epoch import decyear2epoch
from .datetime_from_calarray import datetime_from_calarray
from .datetime2cal import datetime2cal
from .datetime2dayno import datetime2dayno
from .datetime2mjd import datetime2mjd
from .datetime2decyear import datetime2decyear
from .datetime_round_second import datetime_round_second
from .epoch2decyear import epoch2decyear
from .jd2mjd import jd2mjd
from .mjd2jd import mjd2jd
from .yr2year import yr2year
from .year2yr import year2yr
from .leap_year import leap_year
from .uts2hmsmicros import uts2hmsmicros
from .hmsmicros2uts import hmsmicros2uts
from .hmsmicros2ut import hmsmicros2ut
from .ut2uts import ut2uts
from .uts2ut import uts2ut
from .day_since_decyear import day_since_decyear
from .guess_date import guess_date
from .decyear2seconds import decyear2seconds
from .datetime2seconds import datetime2seconds
from .seconds2datetime import seconds2datetime
from .jpldate_2_datetime import jpldate_2_datetime
from .jpldate_2_decyear import jpldate_2_decyear

__all__ = [
    'cal2dayno',
    'cal2mjd',
    'cal2decyear',
    'cal2datetime',
    'dayno2cal',
    'dayno2mjd',
    'dayno2decyear',
    'dayno2datetime',
    'mjd2cal',
    'mjd2dayno',
    'mjd2decyear',
    'mjd2datetime',
    'mjd2gpsweek',
    'gpsweek2mjd',
    'decyear2cal',
    'decyear2dayno',
    'decyear2mjd',
    'decyear2datetime',
    'decyear2epoch',
    'datetime_from_calarray',
    'datetime2cal',
    'datetime2dayno',
    'datetime2mjd',
    'datetime2decyear',
    'datetime_round_second',
    'epoch2decyear',
    'jd2mjd',
    'mjd2jd',
    'yr2year',
    'year2yr',
    'leap_year',
    'uts2hmsmicros',
    'hmsmicros2uts',
    'hmsmicros2ut',
    'ut2uts',
    'uts2ut',
    'day_since_decyear',
    'guess_date',
    'decyear2seconds',
    'datetime2seconds',
    'seconds2datetime',
    'jpldate_2_datetime',
    'jpldate_2_decyear',
]
