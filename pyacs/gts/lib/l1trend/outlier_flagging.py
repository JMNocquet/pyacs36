"""
Outlier flagging using L1-trend analysis.

This module provides functions to identify and flag outliers in GPS time series
based on their deviation from L1-trend filtered representations.
"""

import numpy as np
import pyacs.message.verbose_message as VERBOSE


def flag_outliers_using_l1trend(self, l1trend, threshold=5):
    """
    Flag outliers in the time series based on a user provided l1trend filter representation of the time series.
    
    This function identifies outliers by comparing the original time series to its L1-trend
    representation. Outliers are flagged based on their deviation from the trend using
    the Median Absolute Deviation (MAD) method.
    
    Parameters
    ----------
    self : Gts object
        The time series to flag outliers in.
    l1trend : Gts object
        The l1trend representation of the time series.
    threshold : float, optional
        The threshold for flagging outliers. Default is 5 for 5 times the median 
        absolute deviation (MAD) from the median.
        
    Returns
    -------
    Gts
        The original time series object with outliers flagged in self.outliers
        
    Notes
    -----
    The method uses the Median Absolute Deviation (MAD) which is more robust to outliers
    than standard deviation. Outliers are identified as points that deviate more than
    threshold * MAD from the median of the absolute differences between the original
    time series and the L1-trend representation.
    """
    VERBOSE("Flagging outliers using L1-trend with threshold: %.1f" % threshold)
    
    # Calculate the difference between original time series and L1-trend
    diff_ts = self.substract_ts_daily(l1trend)
    
    # Get absolute differences for N, E, U components
    data_fabs = np.fabs(diff_ts.data[:, 1:4])  # Columns 1, 2, 3 for N, E, U
    
    # Calculate median of absolute differences for each component
    median = np.median(data_fabs, axis=0)
    
    # Calculate Median Absolute Deviation (MAD) for each component
    mad = np.median(np.fabs(data_fabs - median), axis=0)
    
    # Calculate threshold values for each component
    threshold_values = threshold * mad
    
    VERBOSE("MAD values (N, E, U): [%.3f, %.3f, %.3f] mm" % tuple(mad * 1000))
    VERBOSE("Threshold values (N, E, U): [%.3f, %.3f, %.3f] mm" % tuple(threshold_values * 1000))
    
    # Find outliers: points where any component exceeds its threshold
    # Check each component separately and combine results
    outliers_by_component = []
    for i in range(3):  # N, E, U components
        component_outliers = np.where(data_fabs[:, i] > threshold_values[i])[0]
        outliers_by_component.append(component_outliers)
    
    # Combine outliers from all components
    outliers = np.unique(np.concatenate(outliers_by_component))
    
    # Do not flag outliers that are close to identified offsets (5 days)
    if len(outliers) > 0 and len(self.offsets_dates) > 0:
        outliers_dates = self.data[outliers, 0]
        
        outliers_close_to_offsets = np.array([], dtype=int)
        for offset_date in self.offsets_dates:
            # Find outliers that are within 5 days of this offset
            close_indices = np.where(np.abs(outliers_dates - offset_date) < 5/365.25)[0]
            outliers_close_to_offsets = np.concatenate([outliers_close_to_offsets, close_indices])
        
        # Remove duplicates from outliers_close_to_offsets
        outliers_close_to_offsets = np.unique(outliers_close_to_offsets)
        
        # Remove outliers close to offsets
        outliers = np.setdiff1d(outliers, outliers[outliers_close_to_offsets])



    # create new Gts
    new_gts = self.copy()
    new_gts.outliers = np.sort(outliers).tolist()

    
    VERBOSE("Flagged %d outliers out of %d total points (%.1f%%)" % 
           (len(new_gts.outliers), len(self.data), 100.0 * len(new_gts.outliers) / len(self.data)))
    
    return new_gts
