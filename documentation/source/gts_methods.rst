Gts methods
=============

- :ref:`Gts format methods`
- :ref:`Gts primitive methods`
- :ref:`Gts model methods`
- :ref:`Gts filter methods`
- :ref:`Gts outliers methods`
- :ref:`Gts plot methods`


Gts format methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: read_pos,write_pos,read_cats_file,write_cats,force_daily,get_unr,read_kenv,read_mb_file,write_mb_file,read_pride,read_pride_pos,read_tdp,to_pandas_df,read_track_NEU

Gts primitive methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: cdata, copy, differentiate, extract_periods, exclude_periods, extract_dates, substract_ts, substract_ts_daily, add_offsets_dates, remove_velocity , set_zero_at_date, decimate, displacement, add_obs, add_obs_xyz, xyz2neu, neu2xyz, reorder, extract_ndates_before_date, extract_ndates_after_date, extract_ndates_around_date, get_coseismic, correct_duplicated_dates, rotate, insert_gts_data, insert_ts, find_large_uncertainty,split_gap, interpolate                



Gts model methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: detrend,detrend_seasonal,detrend_annual,detrend_median,detrend_seasonal_median,frame,make_model,mmodel,remove_pole,add_vel_sigma, trajectory
	
Gts filter methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: el1_trend, l1_trend, median, minimum_component, piecewise_linear, savitzky_golay, smooth, spline, total_varation, vondrak, wiener
	

Gts outliers methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: remove_outliers, find_l1trend, find_outliers_around_daten find_outliers_percenrage, find_outlier_simple, find_outliers_sliding_window, find_outliers_vondrak

Gts plot methods
*****************************
.. autoclass:: pyacs.gts.Gts.Gts
	:members: plot