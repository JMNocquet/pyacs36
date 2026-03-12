"""
Main refinement function for L1-trend analysis.
"""

import numpy as np
import pyacs.lib.astrotime as at
#from pyacs.gts.Gts import Gts

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG

from .check_trend import check_l1_trend
from .optimal_refinement import optimal_pwlf_refinement, optimal_pwlf_refinement_fast

########################################################################################################################
def refine_l1trend(self, rawts, lcomponent='ENU', min_samples_per_segment=10, 
                   threshold_bias_res_detect=10, threshold_vel=8, 
                   min_sample_segment_refinement=1, norm='L2', output='ts'):
########################################################################################################################
    """
    Refine results from l1trend by optimizing the date of breakpoints for suspected periods.
    The algorithm first identifies periods where the l1trend model might be improved.
    For every candidate period, the algorithm checks different metrics to decide if the model should be refined.
    Finally, for the selected periods, it optimizes the date of breakpoints.

    Parameters
    ----------
    rawts : pyacs.gts.Gts.Gts
        Raw time series, originally processed by l1trend
    lcomponent : str
        Component to refine (default: 'ENU')
    min_samples_per_segment : int
        Minimum number of samples for a segment to be considered as a candidate
        for improvement (default: 10)
    threshold_bias_res_detect : float
        Threshold to detect bias residuals (default: 10)
    threshold_vel : float
        Threshold to detect high velocities (default: 8)
    min_sample_segment_refinement : int
        Minimum number of samples for a segment in the refined model (default: 1)
    norm : str
        Norm to be minimized for improvement (default: 'L2')
    output : str
        Output type 'ts' for the refined Gts, 'info' for the periods suspected, 
        'both' for both (default: 'ts')

    Returns
    -------
    pyacs.gts.Gts.Gts or tuple
        Refined Gts object or information about suspected periods
    """

    # check that self and rawts have the same data dimension
    if self.data.shape != rawts.data.shape:
        ERROR("l1trend and rawts have different data dimensions for %s" % (self.code), exit=True)

    # too short time series
    if self.data.shape[0] < 3:
        WARNING("Time series is too short for refinement for %s. Keeping untouched" % (self.code))
        H_period, H_cp, H_cp_pb = {}, {}, {}
        if output == 'ts':
            return self
        if output == 'both':
            return self, H_period, H_cp, H_cp_pb
        if output == 'info':
            return H_period, H_cp, H_cp_pb


    VERBOSE("----------------------------------------------------------")
    VERBOSE("Refining l1trend results using the new strategy for %s" % (self.code))
    VERBOSE("----------------------------------------------------------")
    VERBOSE("Parameters used for refinement")
    VERBOSE("min_samples_per_segment: %d" % min_samples_per_segment)
    VERBOSE("threshold_bias_res_detect: %d" % threshold_bias_res_detect)
    VERBOSE("threshold_vel: %d" % threshold_vel)
    VERBOSE("min_sample_segment_refinement: %d" % min_sample_segment_refinement)
    VERBOSE("norm: %s" % norm)
    VERBOSE("----------------------------------------------------------")

    # Simplify l1trend first
    # JMN 15/08/2025 : removed, this function should refine the actual time series provided by the user
    #VERBOSE("Cleaning l1trend breakpoints for %s" % (self.code))
    #l1ts = self.simplify_l1trend(tolerance=.5, components='NEU')

    l1ts = self.copy()

    # Check for suspicious periods
    H_period, H_cp, H_cp_pb = check_l1_trend(rawts, l1ts,
                                            component=lcomponent,
                                            min_samples_per_segment=min_samples_per_segment,
                                            threshold_bias_res_detect=threshold_bias_res_detect,
                                            threshold_vel=threshold_vel)

    # Add periods with large offsets to suspicious periods
    VERBOSE("Checking for large offsets to add to suspicious periods")
    large_offset_threshold = 5.0  # mm - threshold for considering an offset as "large"
    
    for component in lcomponent:
        if component not in H_period:
            H_period[component] = []
        
        # Check if there are offsets in the time series
        if hasattr(rawts, 'offsets_values') and rawts.offsets_values is not None:
            comp_idx = {'N': 1, 'E': 2, 'U': 3}[component]
            
            # Get offset values for this component
            offset_values = rawts.offsets_values[:, comp_idx] * 1000  # Convert to mm
            
            # Find large offsets
            large_offset_indices = np.where(np.abs(offset_values) > large_offset_threshold)[0]
            
            if len(large_offset_indices) > 0:
                VERBOSE("Found %d large offsets (>%.1f mm) for component %s" % (len(large_offset_indices), large_offset_threshold, component))
                
                # Get offset dates
                offset_dates = rawts.offsets_dates
                
                for idx in large_offset_indices:
                    if idx < len(offset_dates):
                        offset_date = offset_dates[idx]
                        offset_magnitude = offset_values[idx]
                        
                        # Create a period around the offset (e.g., ±30 days)
                        offset_period_days = 30
                        offset_period_years = offset_period_days / 365.25
                        
                        period_start = offset_date - offset_period_years
                        period_end = offset_date + offset_period_years
                        
                        # Ensure period is within time series bounds
                        period_start = max(period_start, rawts.data[0, 0])
                        period_end = min(period_end, rawts.data[-1, 0])
                        
                        # Add to suspicious periods if not already present
                        new_period = [period_start, period_end]
                        
                        # Check if this period overlaps significantly with existing periods
                        period_exists = False
                        for existing_period in H_period[component]:
                            # Check for overlap
                            overlap_start = max(new_period[0], existing_period[0])
                            overlap_end = min(new_period[1], existing_period[1])
                            
                            if overlap_start < overlap_end:
                                # There is overlap, check if it's significant
                                overlap_duration = overlap_end - overlap_start
                                new_period_duration = new_period[1] - new_period[0]
                                
                                if overlap_duration / new_period_duration > 0.5:  # More than 50% overlap
                                    period_exists = True
                                    break
                        
                        if not period_exists:
                            H_period[component].append(new_period)
                            VERBOSE("Added offset period for %s component: [%.3f, %.3f] (offset: %.1f mm)" % 
                                   (component, period_start, period_end, offset_magnitude))
                        else:
                            VERBOSE("Offset period overlaps with existing suspicious period, skipping")
            else:
                VERBOSE("No large offsets found for component %s" % component)
        else:
            VERBOSE("No offset information available for component %s" % component)

    # Convert H_period decimal years to seconds
    H_period_seconds = {}
    for component in H_period.keys():
        H_period_seconds[component] = []
        for period in H_period[component]:
            start_seconds = at.decyear2seconds(period[0])
            end_seconds = at.decyear2seconds(period[1])
            H_period_seconds[component].append([start_seconds, end_seconds])

    # Convert H_cp to seconds
    H_cp_seconds = {}
    for component in H_cp.keys():
        H_cp_start_end = l1ts.data[H_cp[component],0]
        
        # Insert l1ts.data[0,0] at the first and l1ts.data[-1,0] at the last element
        H_cp_start_end = np.insert(H_cp_start_end, 0, l1ts.data[0,0])
        H_cp_start_end = np.insert(H_cp_start_end, -1, l1ts.data[-1,0])
        
        H_cp_seconds[component] = at.decyear2seconds(H_cp_start_end)

    # Convert l1ts.data[:,0] from decimal year to seconds for time step calculation
    t_seconds = at.decyear2seconds(l1ts.data[:,0])
    
    # Calculate time step for limiting search windows
    if len(t_seconds) > 1:
        time_step = np.median(np.diff(t_seconds))
        max_extension = 50 * time_step  # 50 samples worth of time
        #VERBOSE("Time step: %.1f seconds, max extension: %.1f seconds (50 samples)" % (time_step, max_extension))
    else:
        max_extension = None
    
    # Build H_period_test to add adjacent segments to the suspected periods
    H_period_test = {}
    for component in H_period_seconds.keys():
        H_period_test[component] = []
        for period in H_period_seconds[component]:
            start_seconds = period[0]
            end_seconds = period[1]
            
            # Find the previous and next elements from H_cp_seconds
            cp_seconds = H_cp_seconds[component]
            
            # Find the index of the closest breakpoint before start_seconds
            prev_idx = 0
            for i, cp_time in enumerate(cp_seconds):
                if cp_time < start_seconds:
                    prev_idx = i
                else:
                    break
            
            # Find the index of the closest breakpoint after end_seconds
            next_idx = len(cp_seconds) - 1
            for i, cp_time in enumerate(cp_seconds):
                if cp_time > end_seconds:
                    next_idx = i
                    break
            
            # Replace start and end with previous and next breakpoints
            new_start = cp_seconds[prev_idx]
            new_end = cp_seconds[next_idx]
            
            # Apply 30-sample restriction to prevent excessive search windows
            if max_extension is not None:
                # Limit new_start to max_extension before start_seconds
                min_start = start_seconds - max_extension
                if new_start < min_start:
                    new_start = min_start
                    VERBOSE("Limited new_start to 30 samples before original period start: %.1f -> %.1f seconds" % 
                           (cp_seconds[prev_idx], new_start))
                
                # Limit new_end to max_extension after end_seconds
                max_end = end_seconds + max_extension
                if new_end > max_end:
                    new_end = max_end
                    VERBOSE("Limited new_end to 30 samples after original period end: %.1f -> %.1f seconds" % 
                           (cp_seconds[next_idx], new_end))
            
            H_period_test[component].append([new_start, new_end])

    # Merge periods in H_period_test only if corresponding periods in H_period_seconds are contiguous
    H_period_for_test = {}
    for component in H_period_test.keys():
        H_period_for_test[component] = []
        
        if len(H_period_test[component]) == 0:
            continue
            
        # Get corresponding periods from H_period_seconds
        seconds_periods = H_period_seconds[component]
        
        # Sort periods by start time (both H_period_test and seconds periods)
        sorted_test_periods = sorted(H_period_test[component], key=lambda x: x[0])
        sorted_seconds_periods = sorted(seconds_periods, key=lambda x: x[0])
        
        # Initialize with the first period
        merged_periods = [sorted_test_periods[0]]
        merged_seconds_indices = [0]  # Track which seconds periods are merged
        
        # Check each subsequent period for contiguity in H_period_seconds
        for i, current_test_period in enumerate(sorted_test_periods[1:], 1):
            last_merged = merged_periods[-1]
            last_seconds_idx = merged_seconds_indices[-1]
            
            # Check if corresponding periods in H_period_seconds are contiguous
            tolerance = 1e-6  # 1 microsecond tolerance in seconds
            
            current_seconds_period = sorted_seconds_periods[i]
            last_seconds_period = sorted_seconds_periods[last_seconds_idx]
            
            if abs(current_seconds_period[0] - last_seconds_period[1]) <= tolerance:
                # Seconds periods are contiguous, merge the test periods
                merged_periods[-1] = [last_merged[0], current_test_period[1]]
                merged_seconds_indices[-1] = i  # Update to the latest merged index
            else:
                # Seconds periods are not contiguous, add as new period
                merged_periods.append(current_test_period)
                merged_seconds_indices.append(i)
        
        H_period_for_test[component] = merged_periods

    # Copy l1ts to refine_l1ts for refinement
    refine_l1ts = l1ts.copy()
    
    # Process each component
    for component in H_period_for_test.keys():
        if len(H_period_for_test[component]) == 0:
            VERBOSE("No test periods for component %s" % component)
            continue
        VERBOSE("----------------------------------------------------------")
        VERBOSE("Processing component %s with %d test periods" % (component, len(H_period_for_test[component])))
        VERBOSE("----------------------------------------------------------")
        
        # Get component index
        comp_idx = {'N': 1, 'E': 2, 'U': 3}[component]
        
        for i, period in enumerate(H_period_for_test[component]):
            start_seconds = float(period[0])  # Ensure scalar value
            end_seconds = float(period[1])    # Ensure scalar value
            
            # Convert seconds to datetime for human-readable output
            start_dt = at.seconds2datetime(start_seconds)
            end_dt = at.seconds2datetime(end_seconds)
            VERBOSE("----------------------------------------------------------")
            VERBOSE("Processing period %d/%d for component %s: [%s, %s]" % 
                   (i+1, len(H_period_for_test[component]), component, 
                    start_dt.isoformat()[:10], end_dt.isoformat()[:10]))
            VERBOSE("----------------------------------------------------------")
            
            # Extract time window corresponding to the tested period
            mask = (t_seconds >= start_seconds) & (t_seconds <= end_seconds)
            period_indices = np.where(mask)[0]
            
            if len(period_indices) < 3:
                VERBOSE("Period too short, skipping")
                continue
                
            # Extract data for the time window
            t_period = t_seconds[period_indices]
            y_period = rawts.data[period_indices, comp_idx] * 1E3  # Convert to mm
            l1_period = l1ts.data[period_indices, comp_idx] * 1E3  # Convert to mm
            
            # Get start and end values from l1ts
            y0 = l1_period[0]
            yn = l1_period[-1]
            
            #VERBOSE("Time window: %d points, start: %.3f mm, end: %.3f mm" % 
            #       (len(period_indices), y0, yn))
            
            # Find the indexes in t_seconds corresponding to bp_dates
            bp_dates = H_cp_seconds[component]
            bp_indices = []
            for bp_time in bp_dates:
                # Find the closest index in t_seconds
                closest_idx = np.argmin(np.abs(t_seconds - bp_time))
                bp_indices.append(closest_idx)

            bp_values = l1ts.data[bp_indices, comp_idx]

            # Find breakpoints within the period
            period_bp_indices = []
            for i, bp_idx in enumerate(bp_indices):
                if bp_idx >= period_indices[0] and bp_idx <= period_indices[-1]:
                    period_bp_indices.append(i)
            
            period_bp_dates = bp_dates[period_bp_indices]
            period_bp_values = bp_values[period_bp_indices]
            
            # Count interior breakpoints (excluding start and end of the period)
            # The period boundaries are defined by the previous and next breakpoints from H_cp_seconds
            period_start = t_period[0]
            period_end = t_period[-1]
            
            # Count breakpoints that are strictly within the period (not at boundaries)
            interior_bp_count = 0
            for bp_date in period_bp_dates:
                if period_start < bp_date < period_end:
                    interior_bp_count += 1
            
            nbp = interior_bp_count
            VERBOSE("Number of interior breakpoints in period: %d" % nbp)
            
            if nbp <= 0:
                VERBOSE("No breakpoints to refine, skipping")
                continue
                
            # Store original nbp for potential retry
            original_nbp = nbp
            
            # Limit nbp for performance - exhaustive search becomes too slow for n>=4
            max_nbp = 3  # Limit to 3 breakpoints for reasonable performance
            if nbp > max_nbp:
                VERBOSE("Too many breakpoints (%d), limiting to %d for performance" % (nbp, max_nbp))
                nbp = max_nbp
            
            # Create weights for error calculation based on bias detection
            # Points in segments where bias was detected get higher weight
            weights = np.ones(len(period_indices))
            
            # Check if this period overlaps with any bias-detected segments
            period_start = t_period[0]
            period_end = t_period[-1]
            
            # Get bias-detected segments for this component
            if component in H_period_seconds:
                bias_segments = H_period_seconds[component]
                
                for bias_start, bias_end in bias_segments:
                    # Check if there's overlap between the current period and bias segment
                    overlap_start = max(period_start, bias_start)
                    overlap_end = min(period_end, bias_end)
                    
                    if overlap_start < overlap_end:  # There is overlap
                        # Find indices within the overlap
                        overlap_mask = (t_period >= overlap_start) & (t_period <= overlap_end)
                        overlap_indices = np.where(overlap_mask)[0]
                        
                        # Give higher weight to points in bias-detected segments
                        weights[overlap_indices] = 5.0  # 5x higher weight
                        
                        # Handle both array and single datetime object cases
                        bias_start_dt = at.seconds2datetime(bias_start)
                        bias_end_dt = at.seconds2datetime(bias_end)
                        bias_start_str = bias_start_dt[0].isoformat()[:10] if hasattr(bias_start_dt, '__len__') else bias_start_dt.isoformat()[:10]
                        bias_end_str = bias_end_dt[0].isoformat()[:10] if hasattr(bias_end_dt, '__len__') else bias_end_dt.isoformat()[:10]
                        
                        VERBOSE("Period overlaps with bias segment (%s to %s), applying 5x weight to %d points" % 
                               (bias_start_str, bias_end_str, len(overlap_indices)))
            
            # Function to run optimization with given number of breakpoints
            def run_optimization(n_breakpoints):
                if n_breakpoints <= 3 and len(period_indices) < 200:
                    VERBOSE("Running optimal_pwlf_refinement (exhaustive search) with %d breakpoints" % n_breakpoints)
                    return optimal_pwlf_refinement(t_period, y_period, y0, yn, n_breakpoints, weights)
                else:
                    VERBOSE("Running optimal_pwlf_refinement_fast because exhaustive search would be too slow with %d breakpoints and %d points" % (n_breakpoints, len(period_indices)))
                    return optimal_pwlf_refinement_fast(t_period, y_period, y0, yn, weights)
            
            # First optimization attempt with original nbp
            try:
                optimal_bp, optimal_vals, min_error = run_optimization(nbp)
            except Exception as e:
                VERBOSE("Error in optimal_pwlf_refinement: %s" % str(e))
                import sys; sys.exit()
                continue
                
            original_error = np.sum(weights * (y_period - l1_period)**2)
            
            VERBOSE("Original error: %.6f, Optimal error with %d breakpoints: %.6f" % (original_error, nbp, min_error))
            
            # Test if optimal_pwlf_refinement significantly improves fit
            improvement_ratio = original_error / min_error if min_error > 0 else float('inf')
            significant_improvement = improvement_ratio > 1.10  # 10% improvement threshold
            # Show breakpoint dates in human readable format
            bp_dates_readable = []
            for bp in optimal_bp:
                bp_dt = at.seconds2datetime(bp)
                bp_str = bp_dt[0].isoformat()[:10] if hasattr(bp_dt, '__len__') else bp_dt.isoformat()[:10]
                bp_dates_readable.append(bp_str)
            VERBOSE("Optimal breakpoints for %s component: %s" % (component, bp_dates_readable))
            
            if significant_improvement:
                VERBOSE("Significant improvement detected (ratio: %.2f), updating refine_l1ts" % improvement_ratio)
                
                
                # Interpolate optimal solution back to the period indices
                optimal_interp = np.interp(t_period, optimal_bp, optimal_vals)
                
                # Insert optimal solution into refine_l1ts
                refine_l1ts.data[period_indices, comp_idx] = optimal_interp / 1E3  # Convert back to m
                
                VERBOSE("Successfully updated refine_l1ts for period %d" % (i+1))
            else:
                VERBOSE("No significant improvement with %d breakpoints (ratio: %.2f), trying with %d breakpoints" % (nbp, improvement_ratio, nbp + 1))
                
                # Try with one more breakpoint if we haven't reached the limit
                if nbp + 1 <= max_nbp and nbp + 1 < len(period_indices) - 1:
                    try:
                        optimal_bp_plus1, optimal_vals_plus1, min_error_plus1 = run_optimization(nbp + 1)
                        
                        improvement_ratio_plus1 = original_error / min_error_plus1 if min_error_plus1 > 0 else float('inf')
                        significant_improvement_plus1 = improvement_ratio_plus1 > 1.20
                        
                        VERBOSE("Error with %d breakpoints: %.6f (ratio: %.2f)" % (nbp + 1, min_error_plus1, improvement_ratio_plus1))
                        # Show breakpoint dates in human readable format
                        bp_dates_readable = []
                        for bp in optimal_bp_plus1:
                            bp_dt = at.seconds2datetime(bp)
                            bp_str = bp_dt[0].isoformat()[:10] if hasattr(bp_dt, '__len__') else bp_dt.isoformat()[:10]
                            bp_dates_readable.append(bp_str)
                        VERBOSE("Optimal breakpoints for %s component: %s" % (component, bp_dates_readable))
                        
                        if significant_improvement_plus1 and improvement_ratio_plus1 > improvement_ratio:
                            VERBOSE("Better improvement with %d breakpoints (ratio: %.2f), updating refine_l1ts" % (nbp + 1, improvement_ratio_plus1))
                            
                            
                            # Use the solution with nbp + 1 breakpoints
                            optimal_interp = np.interp(t_period, optimal_bp_plus1, optimal_vals_plus1)
                            refine_l1ts.data[period_indices, comp_idx] = optimal_interp / 1E3  # Convert back to m
                            
                            VERBOSE("Successfully updated refine_l1ts for period %d with %d breakpoints" % (i+1, nbp + 1))
                        else:
                            VERBOSE("No better improvement with %d breakpoints, keeping original" % (nbp + 1))
                    except Exception as e:
                        VERBOSE("Error in optimal_pwlf_refinement with %d breakpoints: %s" % (nbp + 1, str(e)))
                        VERBOSE("Keeping original solution")
                else:
                    VERBOSE("Cannot try %d breakpoints (limit reached or too many points), keeping original" % (nbp + 1))


    # Track which segments have been improved
    improved_segments = set()
    
    # Simplify using Fisher test only on improved segments
    VERBOSE("Simplifying using Fisher test on improved segments")
    
    # Create a copy of refine_l1ts for simplification
    simplified_ts = refine_l1ts.copy()
    
    # For each component, identify segments that were improved and apply Fisher test
    for comp in ['N', 'E', 'U']:
        comp_idx = {'N': 1, 'E': 2, 'U': 3}[comp]
        
        # Get breakpoints for this component
        bp = simplified_ts.l1trend_to_breakpoints(tol='auto')
        bp_dates = bp[comp][0]
        bp_values = bp[comp][1]
        
        if len(bp_dates) <= 2:  # Only start and end points
            continue
            
        # For each breakpoint (except first and last), check if it's in an improved segment
        for i in range(1, len(bp_dates) - 1):
            # Get the segment defined by previous and next breakpoints
            sbp_date = bp_dates[i - 1]
            ebp_date = bp_dates[i + 1]
            
            # Find indices for this segment in the time series
            segment_mask = (simplified_ts.data[:, 0] >= sbp_date) & (simplified_ts.data[:, 0] <= ebp_date)
            segment_indices = np.where(segment_mask)[0]
            
            if len(segment_indices) < 3:  # Need at least 3 points for meaningful test
                continue
            
            # Check if this segment was improved during refinement
            # We can do this by comparing the current values with the original l1ts values
            segment_current = simplified_ts.data[segment_indices, comp_idx]
            segment_original = self.data[segment_indices, comp_idx]
            
            # Calculate if there was significant change in this segment
            segment_diff = np.sum((segment_current - segment_original)**2)
            segment_improvement = segment_diff > 1e-12  # Small threshold to detect changes
            
            if segment_improvement:
                improved_segments.add((comp, sbp_date, ebp_date))
                VERBOSE("Segment %s (%.3f to %.3f) was improved, will apply Fisher test" % (comp, sbp_date, ebp_date))
    
    # Apply Fisher test only on improved segments
    improved_segments = False
    if improved_segments:
        VERBOSE("Applying Fisher test on %d improved segments" % len(improved_segments))
        
        # For each improved segment, create a temporary Gts and apply Fisher test
        for comp, sbp_date, ebp_date in improved_segments:
            comp_idx = {'N': 1, 'E': 2, 'U': 3}[comp]
            
            # Find indices for this segment
            segment_mask = (simplified_ts.data[:, 0] >= sbp_date) & (simplified_ts.data[:, 0] <= ebp_date)
            segment_indices = np.where(segment_mask)[0]
            
            if len(segment_indices) < 3:
                continue
            
            # Create temporary Gts for this segment
            segment_l1ts = simplified_ts.copy()
            segment_l1ts.data = simplified_ts.data[segment_indices, :]
            
            segment_rawts = rawts.copy()
            segment_rawts.data = rawts.data[segment_indices, :]
            
            # Apply Fisher test on this segment
            try:
                segment_simplified = segment_l1ts.simplify_l1trend_with_fisher_test(
                    rawts=segment_rawts,
                    components=comp,
                    alpha=0.05
                )
                
                # Update the main simplified_ts with the results from this segment
                simplified_ts.data[segment_indices, comp_idx] = segment_simplified.data[:, comp_idx]
                
                VERBOSE("Fisher test applied to segment %s (%.3f to %.3f)" % (comp, sbp_date, ebp_date))
                
            except Exception as e:
                VERBOSE("Error applying Fisher test to segment %s (%.3f to %.3f): %s" % (comp, sbp_date, ebp_date, str(e)))
                VERBOSE("Keeping original segment values")
    else:
        VERBOSE("No improved segments found, skipping Fisher test")

    VERBOSE("Refinement completed for %s" % (self.code))

    if output == 'ts':
        return simplified_ts
    if output == 'both':
        return simplified_ts, H_period, H_cp, H_cp_pb
    if output == 'info':
        return H_period, H_cp, H_cp_pb
