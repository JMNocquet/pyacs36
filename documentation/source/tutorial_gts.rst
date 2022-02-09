Quick start on GPS time series analysis
=======================================

.. toctree::
   :maxdepth: 1

Entering ipyacs for interactive analysis
----------------------------------------
::

	ipyacs

should return:

::

	Python 3.7.1 (default, Dec 14 2018, 13:28:58) 
	Type 'copyright', 'credits' or 'license' for more information
	IPython 7.9.0 -- An enhanced Interactive Python. Type '?' for help.
	-- Welcome to pyacs interactive environment -- version  0.02
	- Importing pyacs core module
	- Importing pyacs.gts module
	- Importing class Velocity_Field from pyacs.lib.vel_field module as vf
	- Importing numpy as np
	- Importing matplotlib.pyplot as plt
	- Importing pyacs.lib.astrotime as at
	- Importing pyacs.lib.coordinates as coo
	- Trying to read time series files
	-- Reading directory:  .
	-- No PYACS pck file found
	-- No pride pos files found
	-- No pride_files found
	-- No mb_files found
	-- No tdp_files found
	-- No kenv file found
	-- No cats file found
	-- No Gamit/Globk pos file found
	-- No pyacs t_xyz file found
	-- No Gamit/Globk track NEU file found
	-- read  0  time series in directory  .


pyacs has tried to read everything it could. The resulting loaded time series are stored 
in a Sgts (Super Geodetic Time Series) instance called *ts*. Here ts is empty (yet).


Loading a time series from `UNR <http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap_MAG.html>`_
-------------------------------------------------------------------------------------------------------

::
	
	In [1]: ts.append(Gts().get_unr('ALBH')) 


This command appends a geodetic time series to the ts. Gts() creates an empty time series, which is then fed with the data downloaded from the UNR. 

Visualizing time series
-----------------------

To visualize an individual time series: ::

	In [2]: ts.ALBH.plot()


Individual Gts (Geodetic Time Series) are store as attribute of *ts* and are accessed through ts.XXXX. plot is a method applying to Gts instances.


Detrending time series
______________________

Detrending is simply achieved applying the detrend() method to the Gts instance ts.QUEM: ::

	In [3]: detrended_ALBH=ts.ALBH.detrend()
	In [4]: detrended_ALBH.plot()

Using pyacs'help
----------------

For any function, help is available from the command line of the interactive environment: ::

	In [5]: help(ts.ALBH.plot)

Press *q* to exit from the *help* environment.
Selecting a specific period *[2008.0, 2010.0]* and highlighting two periods [[2008.1,2008.7],[2009.6,2009.8]] can be simply obtained: ::

	In [6]: ts.ALBH.plot(date=[2008.0, 2010.0],lperiod=[[2008.1,2008.7],[2009.6,2009.8]])
	
Chaining methods
----------------

In general, a method applied on a Gts instance will return a new Gts so that various methods can be successively applied: ::

	ts.ALBH.find_outliers_percentage(percentage=0.005).plot().remove_outliers().plot()
	
The line above does:

* select the 0.5% largest residuals of the detrended time series; returns a new Gts
* plot the returned Gts; returns the Gts
* remove the flagged outliers; returns a new time series with the outliers removed
* plot the new time series

