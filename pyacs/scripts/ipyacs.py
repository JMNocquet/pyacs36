#!/usr/bin/env python

version='0.02'
print("-- Welcome to pyacs interactive environment -- version ",version)
import sys, os
print("- Importing pyacs core module")
import pyacs
print("- Importing pyacs.gts module")
from pyacs.gts.Sgts import Sgts
from pyacs.gts.Gts import Gts
print("- Importing class Velocity_Field from pyacs.lib.vel_field module as vf")
from pyacs.lib.vel_field import Velocity_Field as vf
print("- Importing numpy as np")
import numpy as np
print("- Importing matplotlib.pyplot as plt")
import matplotlib.pyplot as plt
print("- Importing pyacs.lib.astrotime as at")
import pyacs.lib.astrotime as at
print("- Importing pyacs.lib.coordinates as coo")
import pyacs.lib.coordinates as coo
print('- Trying to read time series files')
ts=Sgts()

