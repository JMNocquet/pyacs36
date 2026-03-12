"""
Simplification functions for L1-trend filtered time series.
"""

import numpy as np
#from pyacs.gts.Gts import Gts
from scipy.stats import f

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def simplify_l1trend(self, tolerance=.5, components='ENU'):
    """
    Remove unnecessary breakpoints from an L1-trend filtered time series.
    
    This function iteratively removes breakpoints and tests if the simplified model
    still fits the original time series within a specified tolerance.
    
    Parameters
    ----------
    tolerance : float
        Maximum allowed difference (in mm) between original and simplified model.
        Default is 0.5 mm.
    components : str
        Components to process. Default is 'ENU'.
        
    Returns
    -------
    pyacs.gts.Gts.Gts
        Simplified Gts object with unnecessary breakpoints removed.
    """
    VERBOSE("Simplifying l1trend for %s with tolerance %.1f mm and components %s" % (self.code, tolerance, components))
    
    # Convert tolerance from mm to m for consistency with PYACS Gts format 
    tolerance_m = tolerance / 1000.0
    
    # Get breakpoints from the current L1-trend
    bp = self.l1trend_to_breakpoints(tol=0.1)
    
    # Count initial breakpoints for each component
    initial_bp_counts = {}
    for comp in components:
        if comp in ['N', 'E', 'U']:
            initial_bp_counts[comp] = len(bp[comp][0]) - 2  # Subtract 2 for start and end points
    
    # Initialize new data array
    new_data = np.zeros_like(self.data)
    new_data[:,0] = self.data[:,0]  # Copy dates
    new_data[:,4:8] = self.data[:,4:8]  # Copy uncertainties
    
    # Process each component
    for comp in components:
        if comp not in ['N', 'E', 'U']:
            continue
            
        comp_idx = {'N': 1, 'E': 2, 'U': 3}[comp]
        
        # Get breakpoints for this component
        bp_dates = bp[comp][0]
        bp_values = bp[comp][1]
        
        if len(bp_dates) <= 2:  # Only start and end points
            # No breakpoints to simplify, just interpolate
            new_data[:, comp_idx] = np.interp(new_data[:,0], bp_dates, bp_values)
            continue
        
        # Start with all breakpoints
        current_dates = bp_dates.copy()
        current_values = bp_values.copy()
        
        initial_bp_count = len(current_dates) - 2  # Subtract 2 for start and end points
        DEBUG("Processing %s component: starting with %d breakpoints" % (comp, initial_bp_count))
        
        # Iteratively remove breakpoints (except first and last)
        removed_any = True
        iteration = 0
        while removed_any and len(current_dates) > 2:
            removed_any = False
            iteration += 1
            
            # Try removing each breakpoint (except first and last)
            for i in range(1, len(current_dates) - 1):
                # Create test breakpoints without the current breakpoint
                test_dates = np.delete(current_dates, i)
                test_values = np.delete(current_values, i)
                
                # Interpolate test model on original dates
                test_interp = np.interp(self.data[:,0], test_dates, test_values)
                
                # Compare with original L1-trend values
                original_values = self.data[:, comp_idx]
                max_diff = np.max(np.abs(test_interp - original_values))
                
                # If within tolerance, remove this breakpoint
                if max_diff <= tolerance_m:
                    current_dates = test_dates
                    current_values = test_values
                    removed_any = True
                    DEBUG("Iteration %d: removed breakpoint %d from %s component (max_diff: %.3f mm)" % 
                           (iteration, i, comp, max_diff * 1000.0))
                    break  # Start over with the new set of breakpoints
        
        # Interpolate final simplified breakpoints on original dates
        new_data[:, comp_idx] = np.interp(new_data[:,0], current_dates, current_values)
        
        final_bp_count = len(current_dates) - 2  # Subtract 2 for start and end points
        VERBOSE("Initial/Final breakpoints for %s component: %d -> %d" % (comp, initial_bp_counts[comp], final_bp_count))
    
    # Create new Gts object
    simplified_gts = self.copy()
    simplified_gts.data = new_data
    
    return simplified_gts
  

def simplify_l1trend_with_fisher_test(self, rawts, components='ENU', alpha=0.05):
    """
    Simplify an L1-trend filtered time series by comparing to the unfiltered time series using Fisher-Snedecor tests.
    
    This function iteratively removes breakpoints and tests their usefulness by comparing
    the fit quality with and without each breakpoint using Fisher-Snedecor tests.
    
    Parameters
    ----------
    rawts : pyacs.gts.Gts.Gts
        The unfiltered (raw) time series to compare against.
    components : str
        Components to process. Default is 'ENU'.
    alpha : float
        Significance level for the Fisher-Snedecor test. Default is 0.05.
        
    Returns
    -------
    pyacs.gts.Gts.Gts
        Simplified Gts object with unnecessary breakpoints removed.
    """
    VERBOSE("Simplifying l1trend with Fisher test for %s using raw time series %s" % (self.code, rawts.code))
    VERBOSE("Components: %s, Significance level: %.3f" % (components, alpha))
    
    # Check that self and rawts have the same data dimension
    if self.data.shape != rawts.data.shape:
        ERROR("l1trend and rawts have different data dimensions for %s" % (self.code), exit=True)
    
    # Check that dates match
    if not np.allclose(self.data[:, 0], rawts.data[:, 0], rtol=1e-10):
        ERROR("l1trend and rawts have different dates for %s" % (self.code), exit=True)
    
    # Get breakpoints from the current L1-trend
    bp = self.l1trend_to_breakpoints(tol='auto')
    
    # Count initial breakpoints for each component
    initial_bp_counts = {}
    for comp in components:
        if comp in ['N', 'E', 'U']:
            initial_bp_counts[comp] = len(bp[comp][0]) - 2  # Subtract 2 for start and end points
    
    # Initialize new data array
    new_data = np.zeros_like(self.data)
    new_data[:, 0] = self.data[:, 0]  # Copy dates
    new_data[:, 4:8] = self.data[:, 4:8]  # Copy uncertainties
    
    # Process each component
    for comp in components:
        if comp not in ['N', 'E', 'U']:
            continue
            
        comp_idx = {'N': 1, 'E': 2, 'U': 3}[comp]
        
        # Get breakpoints for this component
        bp_dates = bp[comp][0]
        bp_values = bp[comp][1]
        
        if len(bp_dates) <= 2:  # Only start and end points
            # No breakpoints to simplify, just interpolate
            new_data[:, comp_idx] = np.interp(new_data[:, 0], bp_dates, bp_values)
            continue
        
        # Start with all breakpoints
        current_dates = bp_dates.copy()
        current_values = bp_values.copy()
        
        initial_bp_count = len(current_dates) - 2  # Subtract 2 for start and end points
        VERBOSE("Processing %s component: starting with %d breakpoints" % (comp, initial_bp_count))
        
        # Iteratively remove breakpoints (except first and last)
        removed_any = True
        iteration = 0
        while removed_any and len(current_dates) > 2:
            removed_any = False
            iteration += 1
            
            # Try removing each breakpoint (except first and last)
            for i in range(1, len(current_dates) - 1):
                # Get the segment defined by previous and next breakpoints
                sbp_idx = i - 1  # Start breakpoint index
                ebp_idx = i + 1  # End breakpoint index
                
                sbp_date = current_dates[sbp_idx]
                ebp_date = current_dates[ebp_idx]
                
                # Find indices for this segment in the time series
                segment_mask = (self.data[:, 0] >= sbp_date) & (self.data[:, 0] <= ebp_date)
                segment_indices = np.where(segment_mask)[0]
                
                if len(segment_indices) < 3:  # Need at least 3 points for meaningful test
                    continue
                
                # Extract data for this segment
                segment_dates = self.data[segment_indices, 0]
                segment_raw = rawts.data[segment_indices, comp_idx]
                segment_l1 = self.data[segment_indices, comp_idx]
                
                # Calculate straight line between sbp and ebp
                sbp_value = current_values[sbp_idx]
                ebp_value = current_values[ebp_idx]
                
                # Linear interpolation between sbp and ebp
                straight_line = np.interp(segment_dates, [sbp_date, ebp_date], [sbp_value, ebp_value])
                
                # Calculate chi2 for both models
                # Model 1: straight line between sbp and ebp
                chi2_straight = np.sum((segment_raw - straight_line)**2)
                
                # Model 2: current L1-trend (with the breakpoint)
                chi2_l1trend = np.sum((segment_raw - segment_l1)**2)
                
                # Calculate degrees of freedom
                n_points = len(segment_indices)
                df_straight = n_points - 2  # 2 parameters for straight line (slope, intercept)
                df_l1trend = n_points - 3   # 3 parameters for L1-trend with breakpoint (slope1, slope2, breakpoint)
                
                # Ensure degrees of freedom are positive
                if df_straight <= 0 or df_l1trend <= 0:
                    continue
                
                # Calculate F-statistic
                if chi2_l1trend > 0:
                    f_stat = (chi2_straight / df_straight) / (chi2_l1trend / df_l1trend)
                else:
                    f_stat = float('inf')
                
                # Calculate p-value using F-distribution
                from scipy.stats import f
                try:
                    p_value = 1 - f.cdf(f_stat, df_straight, df_l1trend)
                except:
                    p_value = 1.0  # Conservative approach
                
                # Test if breakpoint is useful
                breakpoint_useful = p_value < alpha
                
                VERBOSE("Breakpoint %d (%.3f): F=%.3f, p=%.3f, useful=%s" % 
                       (i, current_dates[i], f_stat, p_value, breakpoint_useful))
                
                # If breakpoint is not useful, remove it
                if not breakpoint_useful:
                    current_dates = np.delete(current_dates, i)
                    current_values = np.delete(current_values, i)
                    removed_any = True
                    VERBOSE("Iteration %d: removed breakpoint %d from %s component (p=%.3f)" % 
                           (iteration, i, comp, p_value))
                    break  # Start over with the new set of breakpoints
        
        # Interpolate final simplified breakpoints on original dates
        new_data[:, comp_idx] = np.interp(new_data[:, 0], current_dates, current_values)
        
        final_bp_count = len(current_dates) - 2  # Subtract 2 for start and end points
        VERBOSE("Initial/Final breakpoints for %s component: %d -> %d" % (comp, initial_bp_counts[comp], final_bp_count))
    
    # Create new Gts object
    simplified_gts = Gts(code=self.code, data=new_data)
    
    return simplified_gts
