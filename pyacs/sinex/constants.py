#-------------------------------------------------------------------------------
# Module   : constants
# Purpose  : Define general constants
# Author   : P. Rebischung
# Created  : 04-Oct-2013
#
# Changes  :
#
# Routines : 
#-------------------------------------------------------------------------------



# LIBRARIES
#-------------------------------------------------------------------------------
from math import *
import numpy



# CONSTANTS
#-------------------------------------------------------------------------------

# Agency
agency = 'IGN'

# Conversion factor from as, mas and uas to rad
as2rad  = pi / 180 / 3600
mas2rad = as2rad / 1000
uas2rad = mas2rad / 1000

# Conversion factor from s, ms and us to rad
s2rad  = pi / 12 / 3600
ms2rad = s2rad / 1000
us2rad = ms2rad / 1000

# MJD at epoch J2000
mjd2000 = 51544.5

# Offset between GPS time and TT
tt_gps = (19. + 32.184) / 86400.

# GPS-UTC leap seconds
gps_utc  = numpy.arange(-9., 19.)
mjd_leap = numpy.array([41317., 41499., 41683., 42048., 42413., 42778., 43144.,
                        43509., 43874., 44239., 44786., 45151., 45516., 46247.,
                        47161., 47892., 48257., 48804., 49169., 49534., 50083.,
                        50630., 51179., 53736., 54832., 56109., 57204., 57754.])

# GRS80 parameters
ae = 6378137.
fe = 0.00335281068118
ee = sqrt(2*fe - fe**2)
be = ae * sqrt(1 - ee**2)

# Gravitational constant
G = 6.67428e-11

# Earth gravitational constant
GM = 3.986004418e14

# Mean equatorial gravity
ge = 9.7803278

# Mean Earth angular velocity
Oe = 7.292115e-5

# Astronomical unit (m)
AU = 149597870691.0

# Rate of advance of ERA
dera_dt = 1.00273781191135448



# DUMMY CLASS
#-------------------------------------------------------------------------------
class record:
    pass
