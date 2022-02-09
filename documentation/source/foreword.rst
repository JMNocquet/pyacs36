
.. toctree::
   :maxdepth: 2


What is pyacs?
==============

Pyacs is a set of scripts and python modules for analyzing and modeling geodetic data.

The philosophy of pyacs
-----------------------

Although pyacs includes a few ready-to-use wrappers, pyacs is rather intended to provide 
a set of tools for analysts that can be used to develop a specific application. pyacs is 
constantly evolving, as new capabilities are added
depending on the need of users. pyacs comes as it is and no warranty.

Pyacs is not intending of replacing more complex software. Rather, pyacs enables quick
and robust solutions to be derived and analyzed. Rapid solutions are achieved through neglecting the
full variance associated with solution. Robustness is ensure through extensive use of the
L1 norm.


The components of pyacs
-----------------------

pyacs includes several components:

* **lib:** library of core functions for dates, coordinates, least-squares

* **ipyacs.py:** pyacs interactive environment through Ipython

* **pyacs_make_time_series.py:** master script to generate time series from free solutions

* **gts:** package for GNSS time series analysis

* **vel_field:** GPS horizontal velocity field analysis
