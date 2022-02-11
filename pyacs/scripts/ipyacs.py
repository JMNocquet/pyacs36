#!/usr/bin/env python
import pyacs.message.message as MESSAGE
import logging
logging.getLogger("my_logger").setLevel(logging.WARNING)

MESSAGE("Welcome to pyacs interactive environment")
import sys, os
MESSAGE("Importing pyacs core module")
import pyacs
MESSAGE("pyacs veriosn: %s" % pyacs.__version__)
MESSAGE("Importing pyacs.gts module")
from pyacs.gts.Sgts import Sgts
from pyacs.gts.Gts import Gts
MESSAGE("Importing class Velocity_Field from pyacs.lib.vel_field module as vf")
from pyacs.lib.vel_field import Velocity_Field as vf
MESSAGE("Importing numpy as np")
import numpy as np
MESSAGE("Importing matplotlib.pyplot as plt")
import matplotlib.pyplot as plt
MESSAGE("Importing pyacs.lib.astrotime as at")
import pyacs.lib.astrotime as at
MESSAGE("Importing pyacs.lib.coordinates as coo")
import pyacs.lib.coordinates as coo
MESSAGE('Trying to read time series files')
ts=Sgts()


