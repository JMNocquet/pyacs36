How to install pyacs?
=====================

Requirements
------------

pyacs is written in python. Since 2018, pyacs has moved to python 3.6.

I recommend to install first a Scientific Python Package like Anaconda or Enthought, which
comes with all required libraries installed.

Getting the software
--------------------

Check the latest release `here <https://github.com/JMNocquet/pyacs36/tree/master/dist/>`_


Installation
------------

Installation should be straightforward if you have Anaconda installed. For a tar bz2 ball pyacs.X.tar.bz2, just type 

::

	pip install pyacs.X.tar.bz2

Anaconda pip program should install all required packages first and then install pyacs components.

You can uninstall pyacs any time by typing

::

	pip uninstall pyacs

Setting your environment
------------------------

* First, check that env python properly launches python
* Check that env ipython properly launches ipython
* verify that Numpy, Scipy and Matplotlib can properly be imported

In your *.bashrc* or *.bash_profile*, or *.cshrc* depending on your environment , add (example for a bash environment) ::

	# PYACS ALIAS
	alias ipyacs='ipython `which ipyacs.py` -i'

Once done, source you configuration file, or open a new terminal window. Check pyacs is working ::

	ipyacs
	
should return something like 

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
	In [1]:                                                                                                                                          


try to import pyacs ::

	import pyacs
	

Additional software
-------------------

Some additional software are useful to work with PYACS and visualize results

* GMT

* Qgis

Issues
------
Sgts.show_map(tile=True) requires cartopy 0.18.0 which does not import properly. Use:

::

	CFLAGS='-stdlib=libc++' pip install cartopy --upgrade

