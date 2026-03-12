.. toctree::
   :maxdepth: 2


What is PYACS?
==============

**PYACS** is a collection of Python modules and scripts for analyzing and
modeling geodetic data, with a focus on tectonic applications.

It provides tools for handling GNSS solutions, geodetic time series, and
velocity fields.

PYACS is designed primarily for researchers working on crustal deformation
and tectonics using geodetic observations.

Python ≥ 3.8 (current version tested with Python 3.12)

Documentation
-------------

Online documentation:

https://jmnocquet.github.io/pyacs_docs/pyacs

Source code:

https://github.com/JMNocquet/pyacs36


The philosophy of PYACS
=======================

Although PYACS includes a few ready-to-use scripts, it is primarily intended
as a collection of tools that analysts can use to develop their own
applications and workflows.

PYACS is continuously evolving as new capabilities are added depending on
the needs of users.

The software is provided **as is**, without warranty.

PYACS is not intended to replace more comprehensive geodetic processing
software. Instead, it aims to enable rapid and robust analysis of geodetic
data.

Rapid solutions are obtained by focusing on practical approximations rather
than fully propagating the complete variance of the solutions. Robustness is
achieved through extensive use of **L1-norm estimators**, which are less
sensitive to outliers.


Main components of PYACS
========================

PYACS includes several main components.

Core libraries
--------------

A set of Python modules handling:

- coordinate transformations
- time representations
- least-squares estimators
- format conversions


GNSS time series generation
---------------------------

The script ``pyacs_make_time_series.py`` converts free GNSS solutions
(typically in SINEX format) into time series expressed in a consistent
reference frame.


Geodetic Time Series module (Gts)
---------------------------------

The ``gts`` package provides a framework for the analysis of geodetic time
series, including:

- loading time series
- visualization
- filtering
- offset estimation
- trend estimation


Velocity field analysis
-----------------------

The ``vel_field`` module provides tools for analyzing horizontal GNSS
velocity fields, including:

- Euler pole estimation
- strain-rate analysis


Interactive environment
-----------------------

``ipyacs.py`` provides an interactive environment for working with PYACS
through IPython.


Related software
================

PYACS also serves as a core module for several related software tools:

- **pygeca** – processing of large GNSS networks using GAMIT,
  including execution in HPC environments

- **pyeblock** – elastic block modelling

- **pyeq / pyaks** – time-dependent slip inversion on faults


Typical workflow
================

A typical PYACS workflow for GNSS time-series analysis may involve the
following steps.

1. Generate station time series from free GNSS solutions in SINEX format::

       pyacs_make_time_series.py

2. Load time series and work interactively, for example in a Jupyter
   notebook environment::

       from pyacs.gts.Sgts import Sgts

       ts = Sgts("timeseries_directory")
       ts.CODE.plot()
       ts.CODE.add_offsets_dates([2016.29, 2018.315]).plot()

       dts = ts.gts('detrend')

3. Export velocity fields and perform Euler pole estimation or strain-rate
   analysis.