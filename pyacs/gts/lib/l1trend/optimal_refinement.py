"""
Optimal piecewise linear function refinement algorithms.
"""

import numpy as np
from itertools import combinations
import pyacs.lib.astrotime as at
#from pyacs.gts.Gts import Gts
from scipy.signal import medfilt

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def preprocess_timeseries_for_optimization(t, y, target_length=50):
    """
    Preprocess time series for optimization by applying median filter and decimation.
    
    This function reduces the number of points while preserving the essential 
    structure of the time series, making optimization faster.
    
    Parameters
    ----------
    t : numpy.ndarray
        Time array
    y : numpy.ndarray
        Values array
    target_length : int
        Target length for the decimated series (default: 50)
        
    Returns
    -------
    tuple
        (t_decimated, y_decimated, decimation_factor)
        - t_decimated: decimated time array
        - y_decimated: decimated values array  
        - decimation_factor: factor used for decimation
    """
    n = len(t)
    
    if n <= target_length:
        # No need to decimate
        return t, y, 1
    
    # Apply median filter with window size 3 to reduce noise
    y_filtered = medfilt(y, kernel_size=3)
    #y_filtered = y

    # Calculate decimation factor
    decimation_factor = max(1, n // target_length)
    
    # Decimate by taking every decimation_factor-th point
    # Ensure we keep the first and last points
    indices = [0]  # Always keep first point
    for i in range(decimation_factor, n-1, decimation_factor):
        indices.append(i)
    if n-1 not in indices:  # Always keep last point
        indices.append(n-1)
    
    t_decimated = t[indices]
    y_decimated = y_filtered[indices]
    
    VERBOSE("Preprocessing: %d points -> %d points (decimation factor: %d)" % 
           (n, len(t_decimated), decimation_factor))
    
    return t_decimated, y_decimated, decimation_factor


def optimal_pwlf_refinement(t, y, y0, yn, nbp, weights=None):
    """
    Find the optimal piecewise linear function with nbp breakpoints fitting y.
    
    The function finds the optimal piecewise linear function (opwlf) fitting y with 
    nbp breakpoints according to the weighted L2 norm. The obtained optimal piecewise linear 
    function verifies opwlf[0]=y0 and opwlf[-1]=yn, that is end and starting point 
    of y are kept fixed.
    
    The problem is applied to small dimension (n<=6), so that a systematic search 
    can be performed and does not require complex optimization approaches.
    
    For large time series (>100 points), preprocessing is applied to decimate the 
    data while preserving essential structure.
    
    Parameters
    ----------
    t : numpy.ndarray
        Time array expressed in integer seconds
    y : numpy.ndarray
        Values of the function to fit
    y0 : float
        Value of the function at t[0] (fixed)
    yn : float
        Value of the function at t[-1] (fixed)
    nbp : int
        Number of breakpoints (excluding start and end points)
    weights : numpy.ndarray, optional
        Weights for each data point. If None, all weights are set to 1.0.
        
    Returns
    -------
    tuple
        (optimal_breakpoints, optimal_values, min_error)
        - optimal_breakpoints: array of breakpoint times
        - optimal_values: array of function values at breakpoints
        - min_error: minimum weighted L2 error achieved
    """
    
    n = len(t)
    
    # Initialize weights if not provided
    if weights is None:
        weights = np.ones(n)
    
    # Apply preprocessing for large time series to speed up optimization
    if n > 100:
        t_opt, y_opt, decimation_factor = preprocess_timeseries_for_optimization(t, y, target_length=50)
        # Adjust y0 and yn for the decimated series
        y0_opt = y_opt[0]
        yn_opt = y_opt[-1]
        # Decimate weights accordingly - need to use the same indices as the decimated data
        # Get the indices used in preprocessing
        n_original = len(t)
        if n_original <= 50:
            # No decimation needed
            weights_opt = weights
        else:
            # Apply median filter and decimation to weights as well
            weights_filtered = medfilt(weights, kernel_size=3)
            # Use the same decimation pattern as in preprocessing
            indices = [0]  # Always keep first point
            for i in range(decimation_factor, n_original-1, decimation_factor):
                indices.append(i)
            if n_original-1 not in indices:  # Always keep last point
                indices.append(n_original-1)
            weights_opt = weights_filtered[indices]
        use_preprocessing = True
    else:
        t_opt, y_opt = t, y
        y0_opt, yn_opt = y0, yn
        weights_opt = weights
        decimation_factor = 1
        use_preprocessing = False
    
    n_opt = len(t_opt)
    
    if n_opt <= 2:
        # No breakpoints possible, return start and end points
        if use_preprocessing:
            # Interpolate back to original time grid
            optimal_interp = np.interp(t, t_opt, np.linspace(y0_opt, yn_opt, n_opt))
            return np.array([t[0], t[-1]]), np.array([y0, yn]), np.sum(weights * (y - optimal_interp)**2)
        else:
            return np.array([t[0], t[-1]]), np.array([y0, yn]), np.sum(weights * (y - np.linspace(y0, yn, n))**2)
    
    if nbp >= n_opt - 1:
        # Too many breakpoints, use all points
        if use_preprocessing:
            # Interpolate back to original time grid
            optimal_interp = np.interp(t, t_opt, y_opt)
            return t, y, np.sum(weights * (y - optimal_interp)**2)
        else:
            return t, y, 0.0
    
    # Generate all possible combinations of breakpoint positions
    # We need to choose nbp positions from the interior points (excluding first and last)
    interior_indices = list(range(1, n_opt-1))
    
    if len(interior_indices) < nbp:
        # Not enough interior points, use all available
        nbp = len(interior_indices)
    
    min_error = float('inf')
    optimal_breakpoints = None
    optimal_values = None
    
    # Try all combinations of breakpoint positions
    for bp_indices in combinations(interior_indices, nbp):
        # Sort breakpoint indices
        bp_indices = sorted(bp_indices)
        
        # Create breakpoint array including start and end points
        bp_times = np.array([t_opt[0]] + [t_opt[i] for i in bp_indices] + [t_opt[-1]])
        
        # For each combination, we need to find the optimal values at breakpoints
        # that minimize the L2 error while keeping y0 and yn fixed
        
        # The problem can be solved using linear algebra
        # We'll solve for the interior breakpoint values
        
        if nbp == 0:
            # No interior breakpoints, just linear interpolation
            pwlf_values = np.linspace(y0, yn, n)
            error = np.sum((y - pwlf_values)**2)
            
            if error < min_error:
                min_error = error
                optimal_breakpoints = bp_times
                optimal_values = np.array([y0, yn])
                
        else:
            # We need to solve for the interior breakpoint values
            # The system is overdetermined, so we'll use weighted least squares
            
            # Create the system matrix for piecewise linear function
            # For each point, we need to express it in terms of the interior breakpoint values
            A = np.zeros((n_opt, nbp))
            b = np.zeros(n_opt)
            
            # Fill the matrix and vector
            for i in range(n_opt):
                # Find which segment t_opt[i] belongs to
                segment_idx = 0
                for j in range(len(bp_times) - 1):
                    if bp_times[j] <= t_opt[i] <= bp_times[j + 1]:
                        segment_idx = j
                        break
                
                if segment_idx == 0:
                    # First segment: linear interpolation between y0 and first interior breakpoint
                    if nbp > 0:
                        # y = y0 + (bp1 - y0) * (t - t0) / (t1 - t0)
                        # y = y0 * (1 - alpha) + bp1 * alpha, where alpha = (t - t0) / (t1 - t0)
                        alpha = (t_opt[i] - bp_times[0]) / (bp_times[1] - bp_times[0])
                        A[i, 0] = alpha
                        b[i] = y0_opt * (1 - alpha)
                elif segment_idx == len(bp_times) - 2:
                    # Last segment: linear interpolation between last interior breakpoint and yn
                    if nbp > 0:
                        # y = bp_n + (yn - bp_n) * (t - t_n) / (t_{n+1} - t_n)
                        # y = bp_n * (1 - alpha) + yn * alpha, where alpha = (t - t_n) / (t_{n+1} - t_n)
                        alpha = (t_opt[i] - bp_times[-2]) / (bp_times[-1] - bp_times[-2])
                        A[i, -1] = 1 - alpha
                        b[i] = yn_opt * alpha
                else:
                    # Middle segment: linear interpolation between two interior breakpoints
                    bp_idx = segment_idx - 1  # Index in interior breakpoints
                    if bp_idx < nbp - 1:
                        # y = bp_j + (bp_{j+1} - bp_j) * (t - t_j) / (t_{j+1} - t_j)
                        # y = bp_j * (1 - alpha) + bp_{j+1} * alpha, where alpha = (t - t_j) / (t_{j+1} - t_j)
                        alpha = (t_opt[i] - bp_times[segment_idx]) / (bp_times[segment_idx + 1] - bp_times[segment_idx])
                        A[i, bp_idx] = 1 - alpha
                        A[i, bp_idx + 1] = alpha
                    else:
                        # This shouldn't happen with proper indexing, but just in case
                        A[i, bp_idx] = 1
            
            # Apply weights to the system
            # Ensure weights_opt has the same length as y_opt
            if len(weights_opt) != len(y_opt):
                VERBOSE("Warning: weights_opt length (%d) != y_opt length (%d), adjusting weights" % (len(weights_opt), len(y_opt)))
                # If lengths don't match, use ones for the missing weights
                if len(weights_opt) < len(y_opt):
                    weights_opt = np.concatenate([weights_opt, np.ones(len(y_opt) - len(weights_opt))])
                else:
                    weights_opt = weights_opt[:len(y_opt)]
            
            W = np.diag(np.sqrt(weights_opt))
            A_weighted = W @ A
            b_weighted = W @ (y_opt - b)
            
            # Solve the weighted least squares problem
            try:
                interior_values = np.linalg.lstsq(A_weighted, b_weighted, rcond=None)[0]
                
                # Reconstruct the full piecewise linear function
                pwlf_values = np.zeros(n_opt)
                for i in range(n_opt):
                    # Find which segment t_opt[i] belongs to
                    segment_idx = 0
                    for j in range(len(bp_times) - 1):
                        if bp_times[j] <= t_opt[i] <= bp_times[j + 1]:
                            segment_idx = j
                            break
                    
                    if segment_idx == 0:
                        # First segment
                        if nbp > 0:
                            pwlf_values[i] = y0_opt + (interior_values[0] - y0_opt) * (t_opt[i] - t_opt[0]) / (bp_times[1] - t_opt[0])
                        else:
                            pwlf_values[i] = y0_opt + (yn_opt - y0_opt) * (t_opt[i] - t_opt[0]) / (bp_times[1] - t_opt[0])
                    elif segment_idx == len(bp_times) - 2:
                        # Last segment
                        if nbp > 0:
                            pwlf_values[i] = interior_values[-1] + (yn_opt - interior_values[-1]) * (t_opt[i] - bp_times[-2]) / (bp_times[-1] - bp_times[-2])
                        else:
                            pwlf_values[i] = y0_opt + (yn_opt - y0_opt) * (t_opt[i] - t_opt[0]) / (bp_times[-1] - t_opt[0])
                    else:
                        # Middle segment
                        bp_idx = segment_idx - 1
                        if bp_idx < nbp - 1:
                            pwlf_values[i] = interior_values[bp_idx] + (interior_values[bp_idx + 1] - interior_values[bp_idx]) * (t_opt[i] - bp_times[segment_idx]) / (bp_times[segment_idx + 1] - bp_times[segment_idx])
                        else:
                            pwlf_values[i] = interior_values[bp_idx] + (yn_opt - interior_values[bp_idx]) * (t_opt[i] - bp_times[segment_idx]) / (bp_times[-1] - bp_times[segment_idx])
                
                # Calculate error
                if use_preprocessing:
                    # Interpolate back to original time grid for error calculation
                    pwlf_interp = np.interp(t, t_opt, pwlf_values)
                    error = np.sum(weights * (y - pwlf_interp)**2)
                else:
                    error = np.sum(weights_opt * (y_opt - pwlf_values)**2)
                
                if error < min_error:
                    min_error = error
                    if use_preprocessing:
                        # Map breakpoint times from decimated to original time grid
                        # Find the closest original time points to the decimated breakpoints
                        optimal_breakpoints = np.zeros_like(bp_times)
                        for j, bp_time in enumerate(bp_times):
                            # Find the closest original time point
                            closest_idx = np.argmin(np.abs(t - bp_time))
                            optimal_breakpoints[j] = t[closest_idx]
                        
                        # For values, we need to interpolate from the decimated solution
                        # to the original time grid
                        optimal_values = np.array([y0] + list(interior_values) + [yn])
                    else:
                        optimal_breakpoints = bp_times
                        optimal_values = np.array([y0] + list(interior_values) + [yn])
                    
                    # Debug: Print some information about the solution
                    #if nbp == 2:  # Special attention to nbp=2 case
                    #    VERBOSE("Found better solution for nbp=%d: error=%.6f, breakpoints=%s" % 
                    #           (nbp, error, [f"{bp:.3f}" for bp in bp_times]))
                    
            except np.linalg.LinAlgError:
                # Skip this combination if the system is singular
                continue
    
    return optimal_breakpoints, optimal_values, min_error


def optimal_pwlf_refinement_fast(t, y, y0, yn, weights=None):
    """
    Fast version of optimal_pwlf_refinement using l1trendi with criterion AICc.
    
    Parameters
    ----------
    t : numpy.ndarray
        Time array expressed in integer seconds
    y : numpy.ndarray
        Values of the function to fit
    y0 : float
        Value of the function at t[0] (fixed)
    yn : float
        Value of the function at t[-1] (fixed)
    weights : numpy.ndarray, optional
        Weights for each data point. If None, all weights are set to 1.0.
        
    Returns
    -------
    tuple
        (optimal_breakpoints, optimal_values, min_error)
    """

    # create a fake Gts object  
    fake_data = np.ones((len(t), 10)) * 1.E-3 # Gts object uses m
    fake_data[:,0] = at.datetime2decyear(at.seconds2datetime(t))
    fake_data[:,1] = y * 1E-3 # mm
    fake_code = ("fast refinement [%s-%s] [%.2lf - %.2lf]" % (at.seconds2datetime(t[0]).isoformat()[:10], at.seconds2datetime(t[-1]).isoformat()[:10], y0, yn))
    fake_gts = Gts(code=fake_code, data=fake_data)

    # run l1trendi with criterion AICc
    l1ts = fake_gts.l1trendi(criterion='AICc', refine=True, lcomponent='N')

    # simplify l1ts
    l1ts = l1ts.simplify_l1trend(tolerance=.5, components='N')

    # run refine_l1trend
    l1ts = l1ts.refine_l1trend(fake_gts, lcomponent='N')

    # get the breakpoints
    breakpoints = l1ts.l1trend_to_breakpoints(tol=0.5)['N']

    # get the breakpoints times in seconds
    breakpoints_times = at.decyear2seconds(breakpoints[0])

    # get the breakpoints values
    breakpoints_values = np.array(breakpoints[1]) * 1E-3 # mm

    # replace the first and last breakpoints with y0 and yn
    breakpoints_values[0] = y0
    breakpoints_values[-1] = yn

    # Initialize weights if not provided
    if weights is None:
        weights = np.ones(len(t))
    
    # interpolate breakpoints_times and breakpoints_values at t to compute min_error
    interpolated_values = np.interp(t, breakpoints_times, breakpoints_values)
    min_error = np.sum(weights * (y - interpolated_values)**2)

    # return
    return breakpoints_times, breakpoints_values, min_error
