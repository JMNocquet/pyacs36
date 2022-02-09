# TESTS FOR gts subpackage

from pyacs.gts.Sgts import Sgts as Sgts
import numpy as np
from pyacs.gts.lib.outliers import find_outliers_sliding_window

dir_test='data/ts'

dir_output = 'output'

# Reading ts files
ts = Sgts(dir_test)

# select site

code = ts.lcode()[0]

# test plot basic
ts.__dict__[code].plot(save=dir_output+'/'+code+'_01.png',verbose=True, title='ex_01: raw time series')

# test outliers percentage
ts.__dict__[code]\
.find_outliers_percentage()\
.plot(save=dir_output+'/'+code+'_02.png',verbose=True, title='ex_02: outliers (percentage method)')

# test outliers smoothing time windows
ts.__dict__[code]\
.find_outliers_sliding_window()\
.plot(save=dir_output+'/'+code+'_03.png',verbose=True, title='ex_03: outliers (sliding window method)')

# test outliers smoothing time windows
ts.__dict__[code]\
.find_out\
.plot(save=dir_output+'/'+code+'_03.png',verbose=True, title='ex_03: outliers (sliding window method)')
