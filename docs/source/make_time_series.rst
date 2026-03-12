.. toctree::
   :maxdepth: 2


pyacs_make_time_series.py
=========================

Generates time series from a list of loosely constrained (free) daily solutions.

::

	pyacs_make_time_series.py [-h] -lsinex LSINEX_NAME -experiment EXPT
                                 [-ref_sinex REF_SINEX_NAME]
                                 [--ref_apr REF_APR_NAME]
                                 [--eq_rename EQ_RENAME]
                                 [--discontinuity DISCONTINUITY]
                                 [--codomes CODOMES] [--psd PSD]
                                 [--conf_file OPTIONS_FILE]
                                 [--dikin_outlier DIKIN_OUTLIER]
                                 [--rep MIN_REPEATABILITY] [--method METHOD]
                                 [--verbose] [--debug] [--replace]
                                 [--uncertainty] [--glred] [--no_clean]

What pyacs_make_time_series.py does
-----------------------------------

pyacs_make_time_series.py performs Helmert transformation of daily (or any position solution) solution onto a reference solution and output time series and statictics to evaluate the solution.



