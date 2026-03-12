
Time series library
===================

PYACS includes two main classes:

* Gts objects are individual time series

	- Gts allow to read, write, plot, detrend, filter individual time series
	- :doc:`gts_class`: Gts structure description
	- :doc:`gts_methods`: list of Gts methods available
	- :doc:`time_series`: tutorial and jupyter notebooks to get started

* Sgts objects are a list of Gts objects

	- :doc:`sgts_class`: Gts structure description


Expert users can also use the tensor format :doc:`./pyacs.gts.tensor_ts`.


Gts class API
-------------

Full reference: :class:`pyacs.gts.Gts.Gts`

The **Gts** class holds NEU (or XYZ) time series with metadata (code, lon/lat/h, velocity, offsets, etc.). Below is a summary of available methods by category.


Summary of available methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read / write and format
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.read`
     - Read time series from file (format auto-detected or set by ``fmt``).
   * - :meth:`~pyacs.gts.Gts.Gts.read_pos`
     - Read GAMIT/GLOBK pos format.
   * - :meth:`~pyacs.gts.Gts.Gts.read_pride`
     - Read PRIDE kinematics format.
   * - :meth:`~pyacs.gts.Gts.Gts.read_cats_file`
     - Read CATS format.
   * - :meth:`~pyacs.gts.Gts.Gts.read_series`
     - Read series format.
   * - :meth:`~pyacs.gts.Gts.Gts.write_pos`
     - Write to pos format.
   * - :meth:`~pyacs.gts.Gts.Gts.write_cats`
     - Write to CATS format.
   * - :meth:`~pyacs.gts.Gts.Gts.to_pandas_df`
     - Convert to pandas DataFrame.
   * - :meth:`~pyacs.gts.Gts.Gts.to_pytrf`
     - Convert to PyTRF time series object.
   * - :meth:`~pyacs.gts.Gts.Gts.get_unr`
     - Fetch time series from UNR.
   * - :meth:`~pyacs.gts.Gts.Gts.force_daily`
     - Force daily sampling.


Primitives (extract, merge, transform)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.copy`
     - Return a copy of the time series.
   * - :meth:`~pyacs.gts.Gts.Gts.extract_periods`
     - Extract one or more time periods.
   * - :meth:`~pyacs.gts.Gts.Gts.extract_dates`
     - Extract data at given dates.
   * - :meth:`~pyacs.gts.Gts.Gts.exclude_periods`
     - Exclude given periods from the series.
   * - :meth:`~pyacs.gts.Gts.Gts.substract_ts`
     - Subtract another Gts (with date matching).
   * - :meth:`~pyacs.gts.Gts.Gts.add_obs`
     - Add NEU observation at a date.
   * - :meth:`~pyacs.gts.Gts.Gts.xyz2neu`
     - Convert XYZ data to NEU.
   * - :meth:`~pyacs.gts.Gts.Gts.neu2xyz`
     - Convert NEU data to XYZ.
   * - :meth:`~pyacs.gts.Gts.Gts.reorder`
     - Reorder by date.
   * - :meth:`~pyacs.gts.Gts.Gts.decimate`
     - Decimate to a given time step.
   * - :meth:`~pyacs.gts.Gts.Gts.displacement`
     - Compute displacement between periods.
   * - :meth:`~pyacs.gts.Gts.Gts.rotate`
     - Rotate NEU components by an angle.
   * - :meth:`~pyacs.gts.Gts.Gts.insert_gts_data`
     - Insert data from another Gts.
   * - :meth:`~pyacs.gts.Gts.Gts.interpolate`
     - Interpolate to regular dates.
   * - :meth:`~pyacs.gts.Gts.Gts.set_zero_at_date`
     - Set displacement to zero at a reference date.
   * - :meth:`~pyacs.gts.Gts.Gts.split`
     - Split at given dates.
   * - :meth:`~pyacs.gts.Gts.Gts.get_coseismic`
     - Get coseismic jump at a date.


Model and detrending
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.detrend`
     - Remove trend (velocity + optional offsets).
   * - :meth:`~pyacs.gts.Gts.Gts.detrend_annual`
     - Remove annual signal.
   * - :meth:`~pyacs.gts.Gts.Gts.detrend_seasonal`
     - Remove seasonal signal.
   * - :meth:`~pyacs.gts.Gts.Gts.detrend_median`
     - Detrend using median filter.
   * - :meth:`~pyacs.gts.Gts.Gts.trajectory`
     - Compute trajectory model.
   * - :meth:`~pyacs.gts.Gts.Gts.frame`
     - Transform to another reference frame (Euler).
   * - :meth:`~pyacs.gts.Gts.Gts.make_model`
     - Build parametric model.
   * - :meth:`~pyacs.gts.Gts.Gts.remove_pole`
     - Remove pole from velocity field.
   * - :meth:`~pyacs.gts.Gts.Gts.detrend_pytrf`
     - Detrend using PyTRF (realistic noise and sigmas).


Offsets
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.find_offsets`
     - Detect and estimate offsets.
   * - :meth:`~pyacs.gts.Gts.Gts.suspect_offsets`
     - Propose candidate offset dates.
   * - :meth:`~pyacs.gts.Gts.Gts.apply_offsets`
     - Apply offset corrections.
   * - :meth:`~pyacs.gts.Gts.Gts.test_offset_significance`
     - Test significance of an offset at a date.
   * - :meth:`~pyacs.gts.Gts.Gts.find_offsets_t_scan`
     - Find offsets using t-scan step detection.
   * - :meth:`~pyacs.gts.Gts.Gts.find_offsets_ivel`
     - Find offsets using instantaneous velocity.


Outliers
^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.remove_outliers`
     - Remove flagged outliers.
   * - :meth:`~pyacs.gts.Gts.Gts.find_outliers_simple`
     - Detect outliers with threshold and window.
   * - :meth:`~pyacs.gts.Gts.Gts.find_outliers_vondrak`
     - Detect outliers using Vondrak filter residuals.
   * - :meth:`~pyacs.gts.Gts.Gts.find_outliers_sliding_window`
     - Detect outliers in sliding windows.
   * - :meth:`~pyacs.gts.Gts.Gts.find_outlier_around_date`
     - Find outlier around a given date.


Noise and uncertainty
^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.realistic_sigma`
     - Estimate realistic uncertainties (e.g. tsfit).
   * - :meth:`~pyacs.gts.Gts.Gts.wrms`
     - Weighted RMS.
   * - :meth:`~pyacs.gts.Gts.Gts.sigma_vel_tsfit`
     - Velocity sigma from time series fit.
   * - :meth:`~pyacs.gts.Gts.Gts.add_vel_sigma`
     - Velocity sigma from white + flicker noise (Williams 2003).


Filters
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.vondrak`
     - Vondrak filter.
   * - :meth:`~pyacs.gts.Gts.Gts.wiener`
     - Wiener filter.
   * - :meth:`~pyacs.gts.Gts.Gts.spline`
     - Spline smoothing.
   * - :meth:`~pyacs.gts.Gts.Gts.median_filter`
     - Median filter.
   * - :meth:`~pyacs.gts.Gts.Gts.edge_filter`
     - Edge-preserving filter.
   * - :meth:`~pyacs.gts.Gts.Gts.disp2vel`
     - Displacement to velocity (filtering).
   * - :meth:`~pyacs.gts.Gts.Gts.ivel`
     - Instantaneous velocity.
   * - :meth:`~pyacs.gts.Gts.Gts.smooth`
     - Smooth time series.


L1-trend (when Trendfilter is installed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.l1trendi`
     - L1 trend filtering.
   * - :meth:`~pyacs.gts.Gts.Gts.refine_l1trend`
     - Refine breakpoint dates from L1 trend.
   * - :meth:`~pyacs.gts.Gts.Gts.clean_l1trend`
     - Clean L1 trend (remove spurious breakpoints).
   * - :meth:`~pyacs.gts.Gts.Gts.simplify_l1trend`
     - Simplify L1 trend (merge segments).
   * - :meth:`~pyacs.gts.Gts.Gts.l1trend_to_breakpoints`
     - Export L1 trend breakpoints.
   * - :meth:`~pyacs.gts.Gts.Gts.l1trend_optimal_workflow`
     - Full L1-trend workflow.


Metadata and I/O
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Gts.Gts.info`
     - Print or return summary info.
   * - :meth:`~pyacs.gts.Gts.Gts.save_velocity`
     - Save velocity to file.
   * - :meth:`~pyacs.gts.Gts.Gts.save_offsets`
     - Save offsets to file.
   * - :meth:`~pyacs.gts.Gts.Gts.plot`
     - Plot the time series.


Full API
~~~~~~~~

.. automodule:: pyacs.gts.Gts
   :no-members:
   :show-inheritance:

Each method in the summary tables above links to its own page under ``gts_methods/``
(e.g. :meth:`~pyacs.gts.Gts.Gts.read` → ``gts_methods/read.html``).


Sgts class API
--------------

Full reference: :class:`pyacs.gts.Sgts.Sgts`

The **Sgts** class holds a collection of Gts (geodetic time series) and provides methods to load, select, and operate on them collectively. Below is a summary of available methods by category.


Summary of available methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read and load
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.read_ts`
     - Load time series from a directory (or list of files).
   * - :meth:`~pyacs.gts.Sgts.Sgts.read_gmt`
     - Read GMT-style time series files.
   * - :meth:`~pyacs.gts.Sgts.Sgts.read_soln`
     - Read from a solution directory.
   * - :meth:`~pyacs.gts.Sgts.Sgts.read_gts_conf`
     - Read time series using a configuration file.


Selection and filtering
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.sel_period`
     - Keep time series overlapping a given period.
   * - :meth:`~pyacs.gts.Sgts.Sgts.sel_radius`
     - Select sites within a radius of a center point.
   * - :meth:`~pyacs.gts.Sgts.Sgts.sel_rectangle`
     - Select sites inside a rectangular bounds.
   * - :meth:`~pyacs.gts.Sgts.Sgts.sel_from_grid`
     - Select sites from a grid with depth range.
   * - :meth:`~pyacs.gts.Sgts.Sgts.sel_radius_eq`
     - Select sites within radius of an earthquake.
   * - :meth:`~pyacs.gts.Sgts.Sgts.sub`
     - Subset by including or excluding site codes.
   * - :meth:`~pyacs.gts.Sgts.Sgts.lGts`
     - Return list of Gts (optionally filtered).
   * - :meth:`~pyacs.gts.Sgts.Sgts.n`
     - Number of time series (optionally filtered).
   * - :meth:`~pyacs.gts.Sgts.Sgts.has_ts`
     - Check if a site code is in the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.same_site`
     - Keep only one time series per site.
   * - :meth:`~pyacs.gts.Sgts.Sgts.nearest`
     - Return the n nearest sites to a given site.


Collection and modification
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.append`
     - Append a Gts to the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.copy`
     - Return a copy of the Sgts.
   * - :meth:`~pyacs.gts.Sgts.Sgts.delts`
     - Remove time series by site code.
   * - :meth:`~pyacs.gts.Sgts.Sgts.delnone`
     - Remove entries with no data.
   * - :meth:`~pyacs.gts.Sgts.Sgts.add_offsets_dates`
     - Add offset dates to time series from a file.
   * - :meth:`~pyacs.gts.Sgts.Sgts.correct_offsets_from_file`
     - Apply offsets from a file to matching time series.
   * - :meth:`~pyacs.gts.Sgts.Sgts.remove_observations`
     - Remove specified observations from all time series.


Frame and apply Gts methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.frame`
     - Transform all time series to another reference frame (Euler).
   * - :meth:`~pyacs.gts.Sgts.Sgts.gts`
     - Call a Gts method on each time series in the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.gts_mp`
     - Call a Gts method on each time series in parallel.


Common mode and coseismic
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.common_mode`
     - Remove common mode signal from the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.compute_common_mode_l1trend`
     - Compute common mode using L1-trend approach.
   * - :meth:`~pyacs.gts.Sgts.Sgts.apply_coseismic`
     - Apply coseismic corrections from GMT files.


Export and write
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.save_velocity`
     - Save velocities to file.
   * - :meth:`~pyacs.gts.Sgts.Sgts.write_pck`
     - Write to PCK format.
   * - :meth:`~pyacs.gts.Sgts.Sgts.to_obs_tensor`
     - Convert to observation tensor.
   * - :meth:`~pyacs.gts.Sgts.Sgts.to_displacement`
     - Export displacement time series.
   * - :meth:`~pyacs.gts.Sgts.Sgts.to_kml`
     - Export to KML.
   * - :meth:`~pyacs.gts.Sgts.Sgts.to_tsnpz`
     - Save to compressed NumPy archive.
   * - :meth:`~pyacs.gts.Sgts.Sgts.to_tspck`
     - Save to TSPCK format.


Info, dates and statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.info`
     - Print or return summary info for the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.dates`
     - Return dates (e.g. decimal year) for the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.get_dates`
     - Get dates in requested format.
   * - :meth:`~pyacs.gts.Sgts.Sgts.lcode`
     - List of site codes (optionally filtered).
   * - :meth:`~pyacs.gts.Sgts.Sgts.medvel`
     - Median velocity estimation for the collection.
   * - :meth:`~pyacs.gts.Sgts.Sgts.stat_site`
     - Site statistics (e.g. velocity, RMS).
   * - :meth:`~pyacs.gts.Sgts.Sgts.make_distance_matrix_from_sgts`
     - Build distance matrix between sites.


Plot and map
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.show_map`
     - Display a map of the sites.
   * - :meth:`~pyacs.gts.Sgts.Sgts.plot_component`
     - Plot a component for all or selected time series.
   * - :meth:`~pyacs.gts.Sgts.Sgts.plot_data_sum`
     - Plot stacked or summary of time series data.
   * - :meth:`~pyacs.gts.Sgts.Sgts.show_ivel_map_gmt`
     - Show instantaneous velocity map (GMT).


Data
^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.gts.Sgts.Sgts.get_unr`
     - Fetch time series from UNR for sites in the collection.


Sgts full API
~~~~~~~~~~~~~

.. automodule:: pyacs.gts.Sgts
   :no-members:
   :show-inheritance:

.. autoclass:: pyacs.gts.Sgts.Sgts
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:


.. toctree::
   :hidden:
   :glob:

   gts_methods/*
