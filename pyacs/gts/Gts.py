
#from pyacs.gts.lib.offset import __fmt_date
"""
    The individual Geodetic Time Series Class

    The Gts implemented in PYACS has the following attributes:

    Mandatory attributes:
        - data:        a 2D numpy array with 10 columns: dec.year, N, E, U, S_N, S_E, S_U, S_NE, S_NU, S_EU
        - code:            station 4-letters code
    Coordinates attributes
        - lon,lat,h:       approximate longitude, latitude (geodetic, deg.dec) and ellipsoidal height (m)
        - X0,Y0,Z0         XYZ reference position in the Geocentric Frame. N,E,U are considered with respect to X0,Y0,Z0
        - t0               reference date in decimal year for X0,Y0,Z0
    Not persisting attributes
        - data_xyz:        a 2D numpy array with 10 columns: dec.year, X, Y, Z, SX, SY, SZ, S_XY, S_XZ, S_YZ
            because many Gts methods are applied on NEU components, .data_xyz is often set to None.
            it can however be rebuilt using the neu2xyz method
    Attributes populated after some analysis
        - outliers:        list of index of outliers in a time series (all components)
        - offsets_values:  a 2D numpy array with 7 columns: dec.year N, E, U, S_N, S_E, S_U
        - offsets_dates:   a list of dates for offsets
        - velocity:        a 1D numpy array with 6 columns: vel_N, vel_E, vel_U, S_vel_N, S_vel_E, S_vel_U
        - annual:          a 1D numpy array with 6 columns: Amplitude_N, Phase_N, Amplitude_E, Phase_E, Amplitude_U, Phase_U
        - semi_annual:     a 1D numpy array with 6 columns: Amplitude_N, Phase_N, Amplitude_E, Phase_E, Amplitude_U, Phase_U
    Metadata attributes
        - ifile:          original input file of the time series
        - log:            log of operations
        - metadata:       any information the analyst would like to be recorded

    Units Conventions
        - dates are in decimal year
        - coordinates are in meters
        - phases are in radians
"""


###################################################################
def get_index_from_dates(dates,data,tol=0.25):
###################################################################
    """
    returns the list of index in data corresponding to given dates within tolerance
    
    :param dates: list of dates in decimal year
    :param data: a 2D numpy array with decimal dates in the first column
    :param tol: date tolerance to decide that two dates are equal. (default 0.25 day)
    
    :return : list of index 
     
    """

    # check argument dates
    
    if isinstance( dates , list ):
        if dates == [] : 
            return( [] )
    
    # check argument data
    
    if data is None:
        return( [] )
    
    # import
    import numpy as np
    
    # sort dates - very important
     
    dates=sorted(dates)
 
    # tolerance in decimal year
    tolerance=tol/365.25
     
    if len(data.shape)>1:data=data[:,0]
     
    index=0
    lindex=[]

    for i in range(len(dates)):

        while data[index]<dates[i]+tolerance:
        
            if np.sqrt((data[index]-dates[i])**2)<tolerance:
                lindex.append(index)
                break
            index=index+1
            
            if index > data.size - 1:
                break
    
    return(lindex)
 
###################################################################
## Init Gts
###################################################################

class Gts:

    def __init__ (self,code=None,\
                  lat=None,lon=None,h=None,\
                  X0=None,Y0=None,Z0=None, t0=None,\
                  data=None,\
                  data_xyz=None,data_corr_neu=None,data_corr_xyz=None,\
                  offsets_dates=[],offsets_values=None,\
                  outliers=[],\
                  annual=None,semi_annual=None,\
                  velocity=None,\
                  ifile=None, log=None, metadata=None):
        """
        The Geodetic time series class
        """
        
        self.code=code

        self.lat = lat
        self.lon = lon
        self.h   = h

        self.X0=X0
        self.Y0=Y0
        self.Z0=Z0
        self.t0=t0


        self.data = data
        
        self.data_xyz = data_xyz
        self.data_corr_neu = data_corr_neu
        self.data_corr_xyz = data_corr_xyz
        
        
        self.outliers = outliers
        self.offsets_dates = offsets_dates
        self.offsets_values = offsets_values
        self.annual = annual
        self.semi_annual = semi_annual
        self.velocity = velocity
        
        self.ifile=ifile
        self.log=log
        self.metadata=metadata

# Old attributes from Trong in case of post-seismic
#        self.earthquakes = earthquakes
#        self.post_seismics = post_seismics
#        self.rate_change = rate_change
#        self.H_conf = H_conf


    @classmethod        
###################################################################
    def read(cls, tsfile, fmt=None, verbose=False):
###################################################################
        
        """
        Reads a time series file
        """

        # import
        import inspect

        # init
        ts = Gts()
        from pyacs.gts.lib.errors import GtsReadFileError

        # PRIDE POS
        if fmt is None or fmt == 'pride_pos':
            try:
                ts.read_pride_pos(tsfile=tsfile , verbose=verbose )
            except:
                return None
                
        # PRIDE KINEMATICS
        if fmt is None or fmt == 'pride':
            try:
                ts.read_pride(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None

        # MB_FILE
        if fmt is None or fmt == 'mb_file':
            try:
                ts.read_mb_file(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None

        # TDP
        if fmt is None or fmt == 'tdp':
            try:
                ts.read_tdp(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None
        
        # KENV
        if fmt is None or fmt == 'kenv':
            try:
                ts.read_kenv(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None
        
        # CATS
        if fmt is None or fmt == 'cats':
            try:
                ts.read_cats(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None
        
        
        # POS FORMAT
        if fmt is None or fmt == 'pos':
            try:
                ts = ts.read_pos(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None

        # TRACK FORMAT
        if fmt is None or fmt == 'track':
            try:
                ts.read_track_neu(tsfile=tsfile , verbose=verbose )
            except:
                raise GtsReadFileError( inspect.stack()[0][3],__name__,tsfile )
                return None

        return( ts )

        
        

###################################################################
##  METHODS IMPORT
###################################################################

# DATES

import pyacs.gts
import pyacs.gts.lib


import pyacs.gts.lib.gts_dates
Gts.np_yyyy_mm_dd_hh_mm_ss_2_decyear=pyacs.gts.lib.gts_dates.np_yyyy_mm_dd_hh_mm_ss_2_decyear
Gts.np_yyyy_mm_dd_hh_mm_ss_2_datetime=pyacs.gts.lib.gts_dates.np_yyyy_mm_dd_hh_mm_ss_2_datetime
Gts.np_datetime_2_eq_time=pyacs.gts.lib.gts_dates.np_datetime_2_eq_time
Gts.decyear2days=pyacs.gts.lib.gts_dates.decyear2days

# FORMATS
import pyacs.gts.lib.format.cats
import pyacs.gts.lib.format.force_daily
import pyacs.gts.lib.format.kenv
import pyacs.gts.lib.format.mb_file
import pyacs.gts.lib.format.pos
import pyacs.gts.lib.format.pride
import pyacs.gts.lib.format.tdp
import pyacs.gts.lib.format.track
import pyacs.gts.lib.format.get_unr
import pyacs.gts.lib.format.to_pandas_df



Gts.read_pos         = pyacs.gts.lib.format.pos.read_pos
Gts.read_pride       = pyacs.gts.lib.format.pride.read_pride
Gts.read_pride_pos    = pyacs.gts.lib.format.pride.read_pride_pos
Gts.read_kenv        = pyacs.gts.lib.format.kenv.read_kenv
Gts.read_mb_file     = pyacs.gts.lib.format.mb_file.read_mb_file
Gts.read_tdp         = pyacs.gts.lib.format.tdp.read_tdp
Gts.read_cats_file   = pyacs.gts.lib.format.cats.read_cats_file
Gts.read_track_NEU   = pyacs.gts.lib.format.track.read_track_NEU

Gts.write_pos        = pyacs.gts.lib.format.pos.write_pos
Gts.write_mb_file    = pyacs.gts.lib.format.mb_file.write_mb_file
Gts.write_cats       = pyacs.gts.lib.format.cats.write_cats

Gts.force_daily      = pyacs.gts.lib.format.force_daily.force_daily

Gts.get_unr          = pyacs.gts.lib.format.get_unr.get_unr
Gts.to_pandas_df     = pyacs.gts.lib.format.to_pandas_df.to_pandas_df

# PRIMITIVE
import  pyacs.gts.lib.primitive.cdata
import  pyacs.gts.lib.primitive.copy
import  pyacs.gts.lib.primitive.differentiate
import  pyacs.gts.lib.primitive.extract_periods
import  pyacs.gts.lib.primitive.exclude_periods
import  pyacs.gts.lib.primitive.extract_dates
import  pyacs.gts.lib.primitive.substract_ts
import  pyacs.gts.lib.primitive.substract_ts_daily
import  pyacs.gts.lib.primitive.add_offsets_dates
import  pyacs.gts.lib.primitive.remove_velocity
import  pyacs.gts.lib.primitive.set_zero_at_date
import  pyacs.gts.lib.primitive.decimate
import  pyacs.gts.lib.primitive.displacement
import  pyacs.gts.lib.primitive.add_obs
import  pyacs.gts.lib.primitive.add_obs_xyz
import  pyacs.gts.lib.primitive.xyz2neu
import  pyacs.gts.lib.primitive.neu2xyz
import  pyacs.gts.lib.primitive.reorder
import  pyacs.gts.lib.primitive.extract_ndates_before_date
import  pyacs.gts.lib.primitive.extract_ndates_after_date
import  pyacs.gts.lib.primitive.extract_ndates_around_date
import  pyacs.gts.lib.primitive.get_coseismic
import  pyacs.gts.lib.primitive.correct_duplicated_dates
import  pyacs.gts.lib.primitive.rotate
import  pyacs.gts.lib.primitive.insert_gts_data
import  pyacs.gts.lib.primitive.find_large_uncertainty
import  pyacs.gts.lib.primitive.split_gap
import  pyacs.gts.lib.primitive.interpolate
import  pyacs.gts.lib.primitive.insert_ts


Gts.cdata                      = pyacs.gts.lib.primitive.cdata.cdata
Gts.copy                       = pyacs.gts.lib.primitive.copy.copy
Gts.differentiate              = pyacs.gts.lib.primitive.differentiate.differentiate
Gts.extract_periods            = pyacs.gts.lib.primitive.extract_periods.extract_periods
Gts.exclude_periods            = pyacs.gts.lib.primitive.exclude_periods.exclude_periods
Gts.extract_dates              = pyacs.gts.lib.primitive.extract_dates.extract_dates
Gts.substract_ts               = pyacs.gts.lib.primitive.substract_ts.substract_ts
Gts.substract_ts_daily         = pyacs.gts.lib.primitive.substract_ts_daily.substract_ts_daily
Gts.add_offsets_dates          = pyacs.gts.lib.primitive.add_offsets_dates.add_offsets_dates
Gts.remove_velocity            = pyacs.gts.lib.primitive.remove_velocity.remove_velocity
Gts.set_zero_at_date           = pyacs.gts.lib.primitive.set_zero_at_date.set_zero_at_date
Gts.decimate                   = pyacs.gts.lib.primitive.decimate.decimate
Gts.displacement               = pyacs.gts.lib.primitive.displacement.displacement
Gts.add_obs                    = pyacs.gts.lib.primitive.add_obs.add_obs
Gts.add_obs_xyz                = pyacs.gts.lib.primitive.add_obs_xyz.add_obs_xyz
Gts.xyz2neu                    = pyacs.gts.lib.primitive.xyz2neu.xyz2neu
Gts.neu2xyz                    = pyacs.gts.lib.primitive.neu2xyz.neu2xyz
Gts.reorder                    = pyacs.gts.lib.primitive.reorder.reorder
Gts.extract_ndates_before_date = pyacs.gts.lib.primitive.extract_ndates_before_date.extract_ndates_before_date
Gts.extract_ndates_after_date  = pyacs.gts.lib.primitive.extract_ndates_after_date.extract_ndates_after_date
Gts.extract_ndates_around_date = pyacs.gts.lib.primitive.extract_ndates_around_date.extract_ndates_around_date
Gts.get_coseismic              = pyacs.gts.lib.primitive.get_coseismic.get_coseismic
Gts.correct_duplicated_dates   = pyacs.gts.lib.primitive.correct_duplicated_dates.correct_duplicated_dates
Gts.rotate                     = pyacs.gts.lib.primitive.rotate.rotate
Gts.insert_gts_data            = pyacs.gts.lib.primitive.insert_gts_data.insert_gts_data
Gts.find_large_uncertainty     = pyacs.gts.lib.primitive.find_large_uncertainty.find_large_uncertainty
Gts.split_gap                  = pyacs.gts.lib.primitive.split_gap.split_gap
Gts.interpolate                = pyacs.gts.lib.primitive.interpolate.interpolate
Gts.insert_ts                  = pyacs.gts.lib.primitive.insert_ts.insert_ts

# METADATA

import pyacs.gts.lib.metadata

Gts.read_lon_lat            = pyacs.gts.lib.metadata.read_lon_lat
Gts.save_velocity           = pyacs.gts.lib.metadata.save_velocity
Gts.save_offsets            = pyacs.gts.lib.metadata.save_offsets
Gts.read_eq_rename          = pyacs.gts.lib.metadata.read_eq_rename
Gts.save_eq_rename          = pyacs.gts.lib.metadata.save_eq_rename
Gts.save_apr                = pyacs.gts.lib.metadata.save_apr
Gts.make_dynamic_apr        = pyacs.gts.lib.metadata.make_dynamic_apr
Gts.read_offset_dates       = pyacs.gts.lib.metadata.read_offset_dates
Gts.info                    = pyacs.gts.lib.metadata.info

# PLOT
import pyacs.gts.lib.plot.plot
## plot
Gts.plot = pyacs.gts.lib.plot.plot.plot

# MODEL

import pyacs.gts.lib.model.add_vel_sigma
import pyacs.gts.lib.model.detrend
import pyacs.gts.lib.model.detrend_annual
import pyacs.gts.lib.model.detrend_seasonal
import pyacs.gts.lib.model.remove_pole
import pyacs.gts.lib.model.frame
import pyacs.gts.lib.model.make_model
import pyacs.gts.lib.model.mmodel
import pyacs.gts.lib.model.detrend_median
import pyacs.gts.lib.model.detrend_seasonal_median
import pyacs.gts.lib.model.trajectory

Gts.add_vel_sigma            = pyacs.gts.lib.model.add_vel_sigma.add_vel_sigma
Gts.detrend                  = pyacs.gts.lib.model.detrend.detrend
Gts.detrend_annual           = pyacs.gts.lib.model.detrend_annual.detrend_annual
Gts.detrend_seasonal         = pyacs.gts.lib.model.detrend_seasonal.detrend_seasonal
Gts.remove_pole              = pyacs.gts.lib.model.remove_pole.remove_pole
Gts.frame                    = pyacs.gts.lib.model.frame.frame
Gts.make_model               = pyacs.gts.lib.model.make_model.make_model
Gts.mmodel                   = pyacs.gts.lib.model.mmodel.mmodel
Gts.detrend_median           = pyacs.gts.lib.model.detrend_median.detrend_median
Gts.detrend_seasonal_median  = pyacs.gts.lib.model.detrend_seasonal_median.detrend_seasonal_median
Gts.trajectory               = pyacs.gts.lib.model.trajectory.trajectory

# OFFSET

import pyacs.gts.lib.offset

Gts.suspect_offsets            = pyacs.gts.lib.offset.suspect_offsets
Gts.suspect_offsets_mf         = pyacs.gts.lib.offset.suspect_offsets_mf
Gts.test_offset_significance   = pyacs.gts.lib.offset.test_offset_significance
Gts.find_time_offsets          = pyacs.gts.lib.offset.find_time_offsets
Gts.delete_small_offsets       = pyacs.gts.lib.offset.delete_small_offsets
#Gts.test_offsets_significance=pyacs.gts.lib.offset.test_offsets_significance
Gts.test_offsets               = pyacs.gts.lib.offset.test_offsets
Gts.estimate_local_offset      = pyacs.gts.lib.offset.estimate_local_offset
Gts.apply_offsets              = pyacs.gts.lib.offset.apply_offsets
Gts.find_offsets               = pyacs.gts.lib.offset.find_offsets
Gts.local_offset_robust        = pyacs.gts.lib.offset.local_offset_robust
Gts.find_offsets_t_scan        = pyacs.gts.lib.offset.find_offsets_t_scan

# OUTLIERS

import pyacs.gts.lib.outliers
import pyacs.gts.lib.outliers.find_l1trend
import pyacs.gts.lib.outliers.find_outliers_percentage
import pyacs.gts.lib.outliers.find_outliers_sliding_window
import pyacs.gts.lib.outliers.find_outliers_vondrak
import pyacs.gts.lib.outliers.remove_outliers
import pyacs.gts.lib.outliers.find_outliers_around_date
import pyacs.gts.lib.outliers.find_outliers_simple



Gts.remove_outliers                                   = pyacs.gts.lib.outliers.remove_outliers.remove_outliers
Gts.find_outliers_percentage                          = pyacs.gts.lib.outliers.find_outliers_percentage.find_outliers_percentage
Gts.find_outliers_simple                              = pyacs.gts.lib.outliers.find_outliers_simple.find_outliers_simple
#Gts.find_outliers_and_offsets_through_differentiation = pyacs.gts.lib.outliers_old.find_outliers_and_offsets_through_differentiation
#Gts.find_outliers_by_RMS_ts                           = pyacs.gts.lib.outliers_old.find_outliers_by_RMS_ts
#Gts.find_outliers_by_residuals                        = pyacs.gts.lib.outliers_old.find_outliers_by_residuals
Gts.find_outliers_sliding_window                      = pyacs.gts.lib.outliers.find_outliers_sliding_window.find_outliers_sliding_window
Gts.find_outlier_around_date                          = pyacs.gts.lib.outliers.find_outliers_around_date.find_outliers_around_date
Gts.find_outliers_vondrak                             = pyacs.gts.lib.outliers.find_outliers_vondrak.find_outliers_vondrak
Gts.find_outliers_l1trend                             = pyacs.gts.lib.outliers.find_l1trend.find_l1trend

# NOISE

import pyacs.gts.lib.noise

Gts.realistic_sigma   = pyacs.gts.lib.noise.realistic_sigma
Gts.wrms              = pyacs.gts.lib.noise.wrms
Gts.sigma_vel_tsfit   = pyacs.gts.lib.noise.sigma_vel_tsfit
Gts.sigma_cats        = pyacs.gts.lib.noise.sigma_cats

# FILTER

import pyacs.gts.lib.filters.median
import pyacs.gts.lib.filters.minimum_component
import pyacs.gts.lib.filters.savitzky_golay
import pyacs.gts.lib.filters.total_variation
import pyacs.gts.lib.filters.vondrak
import pyacs.gts.lib.filters.wiener
import pyacs.gts.lib.filters.piecewise_linear
import pyacs.gts.lib.step_detect_edge_filter
import pyacs.gts.lib.filters.l1_trend
#import pyacs.gts.lib.filters.el1_trend
import pyacs.gts.lib.filters.spline
import pyacs.gts.lib.filters.smooth

Gts.spline            = pyacs.gts.lib.filters.spline.spline
Gts.smooth            = pyacs.gts.lib.filters.smooth.smooth
Gts.median_filter     = pyacs.gts.lib.filters.median.median_filter
Gts.minimum_component = pyacs.gts.lib.filters.minimum_component.minimum_component
Gts.savitzky_golay    = pyacs.gts.lib.filters.savitzky_golay.savitzky_golay
Gts.edge_filter       = pyacs.gts.lib.filters.total_variation.edge
Gts.tv_l2_filter      = pyacs.gts.lib.filters.total_variation.edge_l2
Gts.vondrak           = pyacs.gts.lib.filters.vondrak.vondrak
Gts.wiener            = pyacs.gts.lib.filters.wiener.wiener
Gts.pwlf              = pyacs.gts.lib.filters.piecewise_linear.pwlf
Gts.l1_trend          = pyacs.gts.lib.filters.l1_trend.l1_trend
#Gts.el1_trend          = pyacs.gts.lib.filters.el1_trend.el1_trend
Gts.find_offsets_edge_filter = pyacs.gts.lib.step_detect_edge_filter.find_offsets_edge_filter

