# VERSION HISTORY
## 0.66.20 on 04/03/2025
- pyproject.toml added
- installation documentation revised
- README.md updated
- script make_github.sh added to copy only relevant files for distribution.
- foreword.rst / index.rst changed to reflect README.md
- trendfilter (slightly modified for compatibility) now included in pyacs/external_pkg with 
- added hectorp a package required in setup.py
- change in pyeblock_fault.py to maintain compatibility with pyshp 3.x
## 0.66.19 on 02/03/2025
- vel_field moved to pyacs level. methods are now distributed in individual files
- pyacs.lib.euler, glinalg, gmtpoint, icoshaedron, robustestimators, shafile, timeperiod, units, utils modules converted to packages
- remove lib.pygamit_module. Was saved in pygeca.
- work in progress towards documentation build
- pyacs36/make_pyacs_doc_html_sphinx.sh script now stable to build documentation. Documentation still need update.
- clean setup.py
- remove all previous implementations of l1trend using prox_tv or plwf
- pyacs.gts.lib offset methods have been refactored. This needs further testing.
## 0.66.18 on 26/02/2025
- All Docstring were reformatted to follow Numpy-style
- remove pyacs.lib.astrotime.py and pyacs.lib.coordinates.py files since they were converted to package since pyacs.00.66.16 
## 0.66.17 on 23/02/2025
- make read_pos more robust to read UGA & Geodesy Plotter (https://www.poleterresolide.fr/geodesy-plotter/#/) pos format
- implemented Gts.detrend_hectorp
- implemented Gts.detrend_pytrf
- remove check_apr, that was mistakelly in Sgts_methods. This is a pygeca script.
## 0.66.16 on 19/12/2025
- split pyacs.lib.coordinates.py into a package
## 0.66.15 on 11/12/2025
- split pyacs.lib.astrotime.py into a package
## 0.66.14 on 10/12/2025
- added Sgts.sel_radius_eq
## 0.66.13 on 28/11/2025
- change in Sgts.read_ts to use windows compatible path
## 0.66.12 on 14/11/2025
- added Gts.read_sol for SGC sol format
## 0.66.11 on 14/11/2025
- some improvement to make installation easier: contextily (used in Sgts.show_map) removed from mandatory library
## 0.66.1 on 30/09/2025 & 06/11/2025
- a few part in pyacs_make_time_series.py have been parallelized
- update pandas.read_csv for future compatibility
## 0.66.00 on 10/09/2025
- add Sgts.get_dates
- added sum_l1_trend: print a summary of a previsouly l1trend filtered time series
- get_unr changed to get IGS20
- added ncpu option in pyacs_make_time_series.py
## 0.65.99 on 19/08/2025
- gts.lib.l1trend.optimal_l1trend_workflow introduced in 0.65.98 is stabilized - still needs extensive testing
- pyacs.message has been refactored so that it also handles full logging in files.
## 0.65.98 on 13/08/2025
- l1trend has been refactored in a separate directory gts.lib.l1trend for better maintainability.
## 0.65.97 on 11/08/2025
- added functions l1trend_to_breakpoints and clean_l1trend in pyacs.gts.lib.filter.refine_l1trend. Those functions avoid too many close breakpoints without altering results.
## 0.65.96 on 20250610-20250626
- vel field has been updated to avoid np.asmatrix which is deprecated.
- solved remaing bug with date option in plot
- change in Gts.neu2xyz. Now uses lon, lat, h if X0, Y0, Z0 are missing. Allows to work with CATS format that does not include explicit position of the site.
- added save_file option in Sgts.stat_site. save_dir option kept for compatibility.
- Gts.read_cats has now Gts.data with 10 columns
## 0.65.95 on 20250513
- added ability to read new format for Pride kinematic processing 
- change in Sgts.stat_site - suppress numpy warning when no velocity have estimated
- Gts().write_pos and Gts().write_cats now use pos and cats defaults directory for output. 
## 0.65.94 on 20250425
- Improve robustness of l1trend
- Improve robustness in Sgts.gts_mp. In particular, better handles engine crash.
- Added title option in Gts.show_map
- There was an issue arising from a change in pickle/numpy interaction for python 3.12, making Sgts.code being np.str_ instread of pure str. str is now directly enforced when reading a Super time series .pck.
- change in substract_daily_ts. Code are now kept 4-characters for compatibility with integrity tests in Sgts.gts_mp
## 0.65.93 on 20241119-20250114
- Added Sgts.correct_offsets_from_file to correct offsets from a file
- Added Gts.find_offsets_ivel to estimate offsets with l1trend/ivel
- Added Gts.get_coseismic_l1trend to estimate coseismic offsets with l1trend
- Added Sgts.show_ivel_map_gmt to plot daily  velocity field with GMT
- Correct minor bug in pyacs.lib.shapefile when writing velocity field (name was not properly written)
- Added jpl dates conversions in astrotime
- Added Gts.get_values_at_date
## 0.65.92 on 20240906
- Improve verbose message for Gts.l1trend helping to evaluate the algorithm
- Gts.find_large_uncertainty was in_place. Now returns a new Gts
## 0.65.91 on 20240605
- change in Sgts.show_map to handle change in geopandas. Now download and use GMT coastlines shapefiles
## 0.65.90 on 20240515
- added Sgts.make_distance_matrix_from_sgts and Sgts.nearest
- change in remove_outliers to properly handle .data_xyz
## 0.65.89 on 20240416
- gts.refine_l1trend: test 4-segments optimization
## 0.65.88 on 20240415
- pyacs on python 3.11: compatibility tests + add pyinterp dependency
- sel_from_grid: bug correction when grid crosses the 180/-180 longitude like Alsaka & Kermadec
## 0.65.87 on 20240411
- Bug correction: .data was not created when --pck_only option alone. corrected.
## 0.65.87 on 20240321
- Gts.l1trend: added refined and option to select the component to be processed
## 0.65.86 on 20230927
- Sgts.to_tspck added: time series tensor format
## 0.65.85 on 20230607
- Sgts.delnone added
## 0.65.84 on 20230602
- new implementation of Gts.l1trend. Now allows AICc, Cp and MIX. Optimization alpha algorithm improved
## 0.65.83 on 20230531
- added citerion='AICc' in Gts.l1trend
## 0.65.82 on 20230505
- improved verbose
- added new functionality in Sgts.get_unr allowing spatial selection 
- remove Gts.l1_trend from pyacs distribution
- added test & harmonized import throughout the code; import pyacs.glinalg.solve
- made writing time series more efficient in pyacs_make_time_series.py
## 0.65.81 on 20230504
- read_igs_discontinuity: print line of error when IGS discontinuity file has format problem
- bug correction reading an apr in sinex.py when velocity are not 0.
- correct a rare bug using max_ivel option in Sgts.compute_common_mode
## 0.65.80 on 20230417
- added Gts.ivel
- compute_common_mode_l1trend robust to np.nan
## 0.65.79 on 20230412
- correction line 176-177 in pyacs.lib.shapefile to handle psvelo file with only 1 record
## 0.65.78 on 20230407
- added Gts.edge_filter which is similar to l1trend
## 0.65.77 on 20230405
- added Sgts.gts_mp. requires ipyparallel
## 0.65.76 on 20230404
- added Sgts.compute_common_mode_l1trend
## 0.65.75 on 20230331
- added Gts.read_series to read GipsyX format
## 0.65.74 on 20230328
- Gts.get_unr_loading added: gets loading prediction from UNR
## 0.65.73 on 20230306
- Sgts.save_velocity refactoring
## 0.65.72 on 20230131
- Gts.apply_offset: sets .data_xyz to None to enforce consistency between .data and .data_xyz
- bug in Gts.frame when argument passed w
## 0.65.71 on 20221117
- Sgts.sel_from_grid: selection over a grid like slab2 model

## 0.65.70 on 20220805
- Sgts.plot_component: multiple time series plot

## 0.65.69 on 202200617
- Gts.l1trend: l1 trend filter with optimal BIC criterion

## 0.65.68 on 202200608
- Gts.sub: default is now linclude=None to avoid misinterpretation of []

## 0.65.67 on 20220318 
- added exclude option in pyacs_make_time_series.py conf_file
- new verbose format loading time series Sgts.read_ts
- added Gts.extrapolate
- added Gts.n_obs
## 0.65.66 on 20220209 
- painful cartopy dependency removed thanks to change in Gts.show_map
- added Sgts.to_kml and Sgts.plot_data_sum methods (small bug corrected 20220222)
- set data_xyz = None in Gts.model methods
- add CHANGELOG.md
- Start to use logging to select the verbose level. Work not finished.
- added pyacs.verbose and pyacs.debug
- added art in setup.py
## 0.65.65 on 20211130 
corrected small bug. Sgts.sel_radius when range=[0,X] did not include provided center site
This release has been distributed to various users.
## 0.65.64 on 20211105 
- change to gts.to_pandas_df to allow xyz
- corrected date bug in gts.get_unr (not enough decimal in UNR txyz2 format for decimal year)
- accept a pck ts file as reference
- error through logging module started to be implemented
## 0.65.63 on 20211105 
added Sgts_methods.get_unr
## 0.65.62 on 20211102 
added Sgts_methods.apply_coseismic
## 0.65.61 on 20211013 
added Sgts_methods.to_obs_tensor ; relax version requirements in setup.py
## 0.65.6 on 20211013  
added Jarrin's Euler pole in Sgts.frame & nasty bug in get_coseismic which was not using the eq_date
## 0.65.5 on 20210905  
new setup.py for docker
## 0.65.4 on 20210505  
put obs_tensor2sgts and sgts2obs_tensor taken from pyeq directly in pyacs
## 0.65.3 on 20210303  
refactoring of gts.lib.__init__.py and gts.Gts. added gts.insert_ts and update --replace in pyacs_make_time_series.py
## 0.65.2 on 20210302  
refactoring of gts.lib.model
## 0.65.1 on 20210215  
added pck_only option in pacs_make_time_series.py
## 0.65.0 on 20210215  
added cvxpy as required in setup.py. Remove the png files.
## 0.64.9 on 20210113  
add min_yaxis option in ts.plot
## 0.64.8 on 20201220  
faultslip exploded into individual subroutines. geo2flat changed to use pyproj and web Mercator proj
## 0.64.7 on 20201204  
added outliers package to gts.lib and a new outliers.find_l1trend function
## 0.64.6 on 20201115  
extensive l1 trend filtering added and bug in copy in trajectory (data_xyz must be ignored)
## 0.64.5 on 20201115  
slight change in print Euler pole estimation results
## 0.64.3 on 20201029  
added l1_trend filter
## 0.64.2 on 20201001  
change in Sgts.show_map
## 0.64.1 on 20201001  
correct a bug in pyacs_make_time_series.py when there is a rename in the conf_file and uncertainties are asked
## 0.64.0 on 20200910  
write pck now use DEFAULT_PROTOCOL instead of HIGHEST for backward compatibility
## 0.63.9 on 20200728  
working version, plwf & prox-tv removed from requirement in setup.py
## 0.63.9 on 20200508  
working version
## 0.63.8 on 20200325  
refactoring version glinalg
## 0.63.7 on 20200212  
add common mode testing status
## 0.63.6 on 20200207  
add make_normal_system in pyacs.lib.glinalg
## 0.63.5 on 20200104  
Sgts reorganized
## 0.63.4 on 20200104  
handling of .data and .data_xyz is hopefully more consistent
## 0.63.3 on 20191217  
- some refactoring. gts.read constructor. Should be compatible with jupyter notebook plotting_time_series.ipynb
- turns display off by default
## 0.63.2 on 20191218  
- gts.plot improvement.
- turns display off by default
## 0.63.1 on 20191010  
add seconds conversion in pyacs.lib.astrotime
## 0.63.0 on 20190819  
- add option ref_only_all in pyacs.sol.read_conf, pyacs.sol_sinex, pyacs_make_time_series.py
- fix non optimized tsxyz approach
- pck format for Sgts added ; also added in paycs_make_time_series.py
- pyacs.gts.plot add info option , superimposed can have several gts 
## 0.62.9 on 20190814  
refactoring of gts.lib.primitive. Fix bug handling .data_xyz
## 0.62.8 on 20190805  
- added pyacs.lib.astrotime & coordinates import in pyacs.py. 
- Added pyacs.lib.utils. 
- Bug corrected in pyacs.lib.shapefile when wrting polygons (due to syntax change in pyshp).
## 0.62.7 on 20190607  
added gts.primitive.substract_ts_daily
## 0.62.6 on 20190604  
added spherical_baseline_length_rate in pyacs.lib.coordinates
## 0.62.5 on 20190520  
added scripts/pyacs_gvel_strain.py
## 0.62.4 on 20190514  
Gts.lib.format.py re-organized in Gts.lib.format
## 0.62.3 on 20190503  
reorganize filter for gts. Added total variation filters and piecewise linear fit.
## 0.62.2 on 20190411  
in gts.extract_periods if data_xyz present then data is rebuilt and X0,Y0,Z0,lon,lat,h attributes are updated. Bug in data_xyz extraction corrected.
## 0.62.1 on 20190409  
added pyacs_gvel_pole.py script
## 0.62.0 on 20190409  
small potential bugs in astrotime corrected. Enforcing output type to int for mday, month, noday
## 0.62.0 on 20190401  
pygeca separated from pyacs
## 0.61.6 on 20190312  
script for combining psvelo velocity field added. Experimental.
## 0.61.5 on 20190124  
minor bug in gts.remove_outliers when max index of outliers exceeds the length of the new time series 
## 0.61.5 on 20190124  
minor bug in paycs_make_time_series.py where the correlation coefficients were not properly filled by rotating the XYZ cov matrix 
## 0.61.5 on 20190123  
small change in Sgts.plot to handle the save option
## 0.61.5 on 20190122  
change in Sgts.medvel. Now return the Sgts instance with velocity populated
## 0.61.5 on 20181227  
small bug in gts.plot which was preventing to see the first point with date_unit='cal' option ; also dates are plot at their true time rather at the day at 0000
## 0.61.5 on 20181220  
working version for pyeq_kinematic_inversion.py python 3.7
## 0.61.4 on 20181217  
first Exception handling for gts. Still in test
## 0.61.3 on 20181211  
scaling_laws and magnitude re-introduced, pole estimates, work on median & medvel in Sgts
## 0.61.3 on 20181120  
working version for Green's function cleaning
## 0.61.2 on 20181112  
new calc_pole method in vel_field
## 0.61.2 on 20181112  
bug in Sgts.stat_site where ve and vn had been swaped
## 0.61.2 on 20181112  
bug in gts.exclude_periods corrected
## 0.61.2 on 20181005  
non linear trajectory model added to gts
## 0.61.1.20181004     
Pure python Okada's formulas implementation ; bug corrected in option periods in gts.detrend_median() 
## 0.61.1.20180927     
Major change in pyacs_make_time_series.py that now allows coordinates uncertainties to be calculated through projection
## 0.60.1.20180902     
first python 3.6 working bundle with pyacs_make_time_series and the core library ok. 
## 0.60 MAJOR CHANGE TO PYTHON 3 
- for version > 0.60, PYACS is being progressively migrated to python 3.6 & 3.7, only few features 0.51/python2.7 available for now 
## 0.51                
working code with pygeca
## 0.50                
added observation reweighting options. InSAR capability tested. 
## 0.49                
added insar capability in pyeq_static_inversion.py (quick made from Amazonia, after World Cup final, requires more carefully result printing) 
## 0.48                
add pyeq_model_to_disp.py
## 0.47                
- adding obs-model displacement psvelo gmt files in pyeq_kinematics.py
- occasional bug in apply_offset corrected (> dev3)
- bug in coordinates.xyz_spherical_distance (arcos missing) 
## 0.46                
- adding pyeq_scaling_laws.py
- adding pyacs/lib/GMT.py for backwards compatibility (> dev2)
## 0.45               
- corrected sign error for up Green function using Meade tde
- detrend_median added to Gts.model
- bug gts.file -> gts.ifile (for > dev2) 
- automatic shapefiles for cumulative slip and plots of time series in pyeq_kinematic.py/print_log (dev >3)
- datetime_from_calarray added to pyacs.lib.astrotime
- bug corrected in uts2hmsmicros for seconds calculation
- read_tdp added to gts.lib.format
- added itermax=10 option in Gts.lazy
- added force_day to gts to force daily solution at 1200
- change to gts.format.write_pos to force round at the second
- prototype of new detect_steps methods
## 0.42                
- new lazy command with offset detection using a median filter implemented
- bug correction in plot when use of date and date_unit='days'
## 0.41                
lazy_pyacs seems to be operational again. Patch in pyeq_parametrize_curve_surface_triangles.py and coordinates for GMT5 compatibility
## 0.40                
geometry and green generation scripts included, as well as the static inversion -not fully operational yet.
## 0.39                
new gts method in Sgts
## 0.38                
this version has been distributed to Vergnolle/DeChabalier/Tissandier/
## 0.35                
-cslip time dependent modeling operational again
- pyacs_qgis_model2polygon.py script added to distribution
