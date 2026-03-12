"""
L1-trend summary statistics.

This module provides functions to analyze and summarize L1-trend filtered time series,
including breakpoint statistics, velocity analysis, and segment duration information.
"""

import numpy as np
import pyacs.message.verbose_message as VERBOSE
import pyacs.lib.astrotime as at


def sum_l1_trend(self):
    """
    Summarize L1-trend filtered time series with comprehensive statistics.
    
    This function analyzes a previously L1-trend filtered time series and provides
    detailed statistics including:
    - Time series period
    - Number of breakpoints per component
    - Velocity statistics (mean, median, maximum)
    - Median duration between successive breakpoints
    
    Parameters
    ----------
    self : Gts object
        The L1-trend filtered time series to analyze.
        
    Returns
    -------
    Gts
        The original time series object (unchanged)
        
    Notes
    -----
    This method is designed to work on time series that have already been processed
    with L1-trend filtering. It extracts breakpoints and computes velocity statistics
    for each component (N, E, U).
    """
    VERBOSE("Analyzing L1-trend filtered time series")
    
    # Get time series period
    start_date = self.data[0, 0]
    end_date = self.data[-1, 0]
    period_years = end_date - start_date
    
    VERBOSE("Time series period: %.3f to %.3f (%.2f years)" % 
           (start_date, end_date, period_years))
    
    # Component names
    components = ['N', 'E', 'U']
    
    # Analyze each component
    for i, component in enumerate(components):
        VERBOSE("\n--- Component %s ---" % component)
        
        # Get component data (column i+1: N=1, E=2, U=3)
        component_data = self.data[:, i+1]
        
        # Find breakpoints by detecting significant changes in slope
        breakpoints = _find_breakpoints(self.data[:, 0], component_data)
        
        VERBOSE("Number of breakpoints: %d" % len(breakpoints))
        
        if len(breakpoints) > 0:
            # Calculate velocities between breakpoints
            velocities = _calculate_velocities(self.data[:, 0], component_data, breakpoints)
            
            # Calculate durations between breakpoints
            durations = _calculate_durations(self.data[:, 0], breakpoints)
            
            # Statistics
            mean_vel = np.mean(velocities)
            median_vel = np.median(velocities)
            max_vel = np.max(np.abs(velocities))
            median_duration = np.median(durations)
            
            VERBOSE("Velocity statistics:")
            VERBOSE("  Mean velocity: %.3f mm/year" % (mean_vel * 1000))
            VERBOSE("  Median velocity: %.3f mm/year" % (median_vel * 1000))
            VERBOSE("  Maximum |velocity|: %.3f mm/year" % (max_vel * 1000))
            VERBOSE("  Median duration between breakpoints: %.1f days" % (median_duration * 365.25))
            
            # Additional statistics
            VERBOSE("Additional statistics:")
            VERBOSE("  Number of segments: %d" % (len(breakpoints) + 1))
            VERBOSE("  Velocity range: [%.3f, %.3f] mm/year" % 
                   (np.min(velocities) * 1000, np.max(velocities) * 1000))
            VERBOSE("  Duration range: [%.1f, %.1f] days" % 
                   (np.min(durations) * 365.25, np.max(durations) * 365.25))
        else:
            VERBOSE("No breakpoints detected - linear trend")
            # Calculate overall velocity
            overall_vel = (component_data[-1] - component_data[0]) / (end_date - start_date)
            VERBOSE("Overall velocity: %.3f mm/year" % (overall_vel * 1000))
    
    VERBOSE("\nL1-trend analysis completed")
    return self


def _find_breakpoints(t, y, threshold=0.1):
    """
    Find breakpoints in a time series by detecting significant changes in slope.
    
    Parameters
    ----------
    t : array_like
        Time array (decimal years)
    y : array_like
        Data array
    threshold : float
        Threshold for detecting significant slope changes
        
    Returns
    -------
    array
        Array of breakpoint indices
    """
    if len(t) < 3:
        return np.array([])
    
    # Calculate first differences (velocities)
    dt = np.diff(t)
    dy = np.diff(y)
    velocities = dy / dt
    
    # Calculate second differences (acceleration)
    d2v = np.diff(velocities)
    
    # Find significant changes in velocity
    # Use a threshold based on the standard deviation of second differences
    if len(d2v) > 0:
        std_d2v = np.std(d2v)
        significant_changes = np.abs(d2v) > threshold * std_d2v
        
        # Get indices of significant changes (add 1 because of diff)
        breakpoints = np.where(significant_changes)[0] + 1
        
        # Filter out breakpoints too close to each other (minimum 5 points apart)
        if len(breakpoints) > 1:
            min_gap = 5
            filtered_breakpoints = [breakpoints[0]]
            for bp in breakpoints[1:]:
                if bp - filtered_breakpoints[-1] >= min_gap:
                    filtered_breakpoints.append(bp)
            breakpoints = np.array(filtered_breakpoints)
    else:
        breakpoints = np.array([])
    
    return breakpoints


def _calculate_velocities(t, y, breakpoints):
    """
    Calculate velocities between breakpoints.
    
    Parameters
    ----------
    t : array_like
        Time array (decimal years)
    y : array_like
        Data array
    breakpoints : array_like
        Array of breakpoint indices
        
    Returns
    -------
    array
        Array of velocities between breakpoints
    """
    velocities = []
    
    # Add start and end points
    all_points = np.concatenate([[0], breakpoints, [len(t)-1]])
    
    for i in range(len(all_points) - 1):
        start_idx = all_points[i]
        end_idx = all_points[i + 1]
        
        if end_idx > start_idx:
            dt = t[end_idx] - t[start_idx]
            dy = y[end_idx] - y[start_idx]
            velocity = dy / dt if dt > 0 else 0
            velocities.append(velocity)
    
    return np.array(velocities)


def _calculate_durations(t, breakpoints):
    """
    Calculate durations between breakpoints.
    
    Parameters
    ----------
    t : array_like
        Time array (decimal years)
    breakpoints : array_like
        Array of breakpoint indices
        
    Returns
    -------
    array
        Array of durations between breakpoints (in years)
    """
    durations = []
    
    # Add start and end points
    all_points = np.concatenate([[0], breakpoints, [len(t)-1]])
    
    for i in range(len(all_points) - 1):
        start_idx = all_points[i]
        end_idx = all_points[i + 1]
        
        if end_idx > start_idx:
            duration = t[end_idx] - t[start_idx]
            durations.append(duration)
    
    return np.array(durations)
