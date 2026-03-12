"""
Optimal L1-trend workflow for GPS time series analysis.

This module provides a complete workflow that combines L1-trend filtering,
simplification, and refinement for optimal results.
"""

import time
import numpy as np
import pyacs.message.verbose_message as VERBOSE
import logging
import os


def l1trend_optimal_workflow(self, criterion='AICc', component='NEU', log_dir=None):
    """
    Optimal workflow for l1 trend analysis.
    
    This workflow combines multiple L1-trend processing steps for optimal results:
    1. Initial L1-trend filtering with AICc criterion
    2. Breakpoint simplification to remove redundant breakpoints
    3. Local refinement using original data for better fit
    4. Final simplification to clean up the refined result
    
    Parameters
    ----------
    criterion : str, optional
        Criterion for L1-trend optimization ('AICc', 'BIC', 'Cp')
    component : str, optional
        Components to process ('N', 'E', 'U', 'NE', 'NU', 'EU', 'NEU'). Default is 'NEU'
    log_dir : str, optional
        Directory for Gts-specific logging. If provided, creates a log file named {self.code}.log
    
    Returns
    -------
    Gts
        Final optimized L1-trend time series
    """

    # import
    import pyacs.lib.astrotime as at

    start_time = time.time()

    # Setup Gts-specific logging if log_dir is provided
    gts_logger = None
    if log_dir is not None:
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Create Gts-specific log file
        log_file = os.path.join(log_dir, f"{self.code}.log")
        
        # Configure Gts-specific logger
        gts_logger = logging.getLogger(f'pyacs.gts.{self.code}')
        if not gts_logger.handlers:
            gts_logger.setLevel(logging.INFO)
            # Prevent propagation to parent loggers to avoid duplicate logging
            gts_logger.propagate = False
            file_handler = logging.FileHandler(log_file, mode='a')  # append mode
            file_handler.setLevel(logging.INFO)
            # Custom formatter with [Gts.code] prefix
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - [%(name)s] %(message)s')
            file_handler.setFormatter(formatter)
            gts_logger.addHandler(file_handler)
        
        # Log workflow start with [Gts.code] prefix
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Starting workflow")
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Parameters: criterion={criterion}, component={component}")

    # Validate component parameter and handle equivalent combinations
    valid_components = ['N', 'E', 'U', 'NE', 'EN', 'NU', 'UN', 'EU', 'UE', 'NEU', 'NUE', 'ENU', 'EUN', 'UNE', 'UEN']
    if component not in valid_components:
        raise ValueError(f"Invalid component '{component}'. Must be one of: {valid_components}")
    
    # Convert component string to list of individual components and sort for consistency
    component_list = sorted(list(set(component)))  # Remove duplicates and sort
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Processing components: {component_list}")
    
    # Calculate time series statistics
    n_points = len(self.data)
    time_span_years = self.data[-1, 0] - self.data[0, 0]
    
    # Convert start and end dates to human readable format
    start_date_decyear = self.data[0, 0]
    end_date_decyear = self.data[-1, 0]
    start_date_readable = at.decyear2datetime(start_date_decyear).strftime('%Y-%m-%d')
    end_date_readable = at.decyear2datetime(end_date_decyear).strftime('%Y-%m-%d')
    
    VERBOSE("Starting optimal L1trend_optimal_workflow for %s %d points, %.2f years span (%s to %s)" % 
            (self.code, n_points, time_span_years, start_date_readable, end_date_readable))
    
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] {n_points} epochs, {time_span_years:.2f} years span ({start_date_readable} to {end_date_readable})")
    
    # Step 1: remove large uncertainty points that are likely outliers
    VERBOSE("Step 1: Removing large uncertainty points that are likely outliers")
    #if gts_logger:
    #    gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 1: Removing large uncertainty points")
    ts_no_large_uncertainty = self.find_large_uncertainty().remove_outliers()
    outliers_removed = self.data.shape[0] - ts_no_large_uncertainty.data.shape[0]
    VERBOSE("Removed %d outliers" % outliers_removed)
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 1: Removed {outliers_removed} large uncertainty outliers")

    # Step 2: detrend the time series
    VERBOSE("Step 2: Detrending the time series")
    #if gts_logger:
    #   gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 2: Detrending time series")
    if ts_no_large_uncertainty.data[-1,0] - ts_no_large_uncertainty.data[0,0] > 1:
        try:
            detrended_ts = ts_no_large_uncertainty.detrend_median()
            if gts_logger:
                gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 2: Used detrend_median method")
        except:
            VERBOSE("detrend_median failed, using detrend method")
            detrended_ts = ts_no_large_uncertainty.detrend()
    else:            
        detrended_ts = ts_no_large_uncertainty.detrend()
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 2: Used detrend method (fallback)")


    # Step 3: l1trendi with criterion 
    step3_start_time = time.time()
    VERBOSE("Step 3: Initial L1-trend filtering with criterion %s" % criterion)
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Initial L1-trend filtering with criterion {criterion}")
    
    # Create a copy for processing only specified components
    detrended_ts_processed = detrended_ts.copy()
    
    # Zero out components that are not being processed
    component_indices = {'N': 1, 'E': 2, 'U': 3}
    for comp in ['N', 'E', 'U']:
        if comp not in component_list:
            detrended_ts_processed.data[:, component_indices[comp]] = 0.0
    
    l1_aicc = detrended_ts_processed.l1trendi(criterion=criterion, refine=False, component=component)
    step3_time = time.time() - step3_start_time
    VERBOSE("Step 3 completed in %.2f seconds" % step3_time)
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Completed in {step3_time:.2f} seconds")
    
    # Step 3 statistics
    VERBOSE("Step 3 - Initial L1-trend Statistics:")
    #if gts_logger:
    #   gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Computing statistics")
    
    # Number of breakpoints per component in a single line

    bp_step3 = l1_aicc.l1trend_to_breakpoints(tol='auto')

    bp_counts = []
    for comp in component_list:
        if comp in bp_step3:
            n_bp = len(bp_step3[comp][0]) - 2  # Subtract start and end points
            bp_counts.append(f"{comp}:{n_bp}")
        else:
            bp_counts.append(f"{comp}:0")
    VERBOSE("  Breakpoints: %s" % ", ".join(bp_counts))
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Breakpoints found: {', '.join(bp_counts)}")
    
    # Offset dates
    if hasattr(l1_aicc, 'offsets_dates') and l1_aicc.offsets_dates:
        VERBOSE("  Offset dates: %s" % [f"{date:.3f}" for date in l1_aicc.offsets_dates])
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Offset dates: {[f'{date:.3f}' for date in l1_aicc.offsets_dates]}")
        # Convert to human readable format
        offset_dates_readable = []
        for date in l1_aicc.offsets_dates:
            try:
                date_dt = at.decyear2datetime(date)
                if hasattr(date_dt, '__len__'):
                    date_str = date_dt[0].isoformat()[:10]
                else:
                    date_str = date_dt.isoformat()[:10]
                offset_dates_readable.append(date_str)
            except:
                offset_dates_readable.append(f"{date:.3f}")
        VERBOSE("  Offset dates (readable): %s" % offset_dates_readable)
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: Offset dates: {offset_dates_readable}")
    else:
        VERBOSE("  No offset dates found")
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 3: No offset dates found")

    

    # Step 4: remove obvious outliers
    VERBOSE("Step 4: Removing obvious outliers using flag_outliers_using_l1trend")
    #if gts_logger:
    #    gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 4: Flagging outliers using L1-trend")
    
    detrended_ts_outliers_flagged = detrended_ts_processed.flag_outliers_using_l1trend(l1_aicc, threshold=20)
    
    # Count outliers by component
    if detrended_ts_outliers_flagged.outliers:
        outlier_indices = np.array(detrended_ts_outliers_flagged.outliers)
        outlier_counts = []
        for comp in component_list:
            comp_idx = {'N': 1, 'E': 2, 'U': 3}[comp]
            # Count outliers where this component exceeds threshold
            component_outliers = np.sum(np.abs(detrended_ts_outliers_flagged.data[outlier_indices, comp_idx]) > 0.020)  # 20mm threshold
            outlier_counts.append(f"{comp}:{component_outliers}")
        VERBOSE("  Outliers removed by component: %s" % ", ".join(outlier_counts))
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 4: Outliers removed by component: {', '.join(outlier_counts)}")
    else:
        VERBOSE("  No outliers found")
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 4: No outliers found")
    
    dts_no_outliers = detrended_ts_outliers_flagged.remove_outliers()
    l1_aicc_no_outliers = l1_aicc.copy()
    l1_aicc_no_outliers.outliers = detrended_ts_outliers_flagged.outliers
    l1_aicc_no_outliers = l1_aicc_no_outliers.remove_outliers()

    
    # Step 5: simplify l1 trend
    VERBOSE("Step 5: Simplifying L1-trend breakpoints")
    #if gts_logger:
    #   gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 5: Simplifying L1-trend breakpoints")
    l1_simplified = l1_aicc_no_outliers.simplify_l1trend()
    
    # Get breakpoints after simplification
    bp_step5 = l1_simplified.l1trend_to_breakpoints(tol='auto')
    
    # Count simplified breakpoints by component
    bp_counts_step5 = []
    for comp in component_list:
        if comp in bp_step5:
            n_bp = len(bp_step5[comp][0]) - 2  # Subtract start and end points
            bp_counts_step5.append(f"{comp}:{n_bp}")
        else:
            bp_counts_step5.append(f"{comp}:0")
    
    VERBOSE("  Simplified breakpoints: %s" % ", ".join(bp_counts_step5))
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 5: New breakpoints: {', '.join(bp_counts_step5)}")



    # Step 6: refine l1 trend with the original data
    step6_start_time = time.time()
    VERBOSE("Step 6: Refining L1-trend using original data")
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 6: Refining L1-trend using original data")
    
    # Instead of splitting into segments, process the entire time series
    # but ensure the refine_l1trend function handles offsets properly
    VERBOSE("Step 6.1: Refining entire time series with offset awareness")
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 6.1: Refining entire time series with offset awareness")
    
    # Refine the entire simplified L1-trend using the full detrended time series
    l1_refined, H_period, H_cp, H_cp_pb = l1_simplified.refine_l1trend(dts_no_outliers, min_samples_per_segment=4,output='both')
    
    # Get refinement information to display refined periods
    #H_period, H_cp, H_cp_pb = l1_simplified.refine_l1trend(dts_no_outliers, min_samples_per_segment=4, output='info')
    
    # Display refined periods in human-readable format
    if H_period:
        VERBOSE("Refined periods:")
        for comp in component_list:
            if comp in H_period and H_period[comp]:
                VERBOSE("  %s component:" % comp)
                for i, period in enumerate(H_period[comp]):
                    start_date = period[0]
                    end_date = period[1]
                    try:
                        start_dt = at.decyear2datetime(start_date)
                        end_dt = at.decyear2datetime(end_date)
                        if hasattr(start_dt, '__len__'):
                            start_str = start_dt[0].isoformat()[:10]
                        else:
                            start_str = start_dt.isoformat()[:10]
                        if hasattr(end_dt, '__len__'):
                            end_str = end_dt[0].isoformat()[:10]
                        else:
                            end_str = end_dt.isoformat()[:10]
                        VERBOSE("    Period %d: %s to %s" % (i+1, start_str, end_str))
                    except:
                        VERBOSE("    Period %d: %.3f to %.3f" % (i+1, start_date, end_date))
            else:
                VERBOSE("  %s component: No periods refined" % comp)
    else:
        VERBOSE("No periods were refined")
    
    VERBOSE("Refinement completed on entire time series")
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 6: Refinement completed on entire time series")

    step6_time = time.time() - step6_start_time
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 6: Completed in {step6_time:.2f} seconds")
    
    # Step 7: simplify again
    VERBOSE("Step 7: Final simplification of refined L1-trend")
    #if gts_logger:
    #    gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 7: Final simplification of refined L1-trend")
    
    # Get breakpoints before simplification
    bp_before = l1_refined.l1trend_to_breakpoints(tol='auto')
    
    l1_final_detrended = l1_refined.simplify_l1trend()
    
    # Get breakpoints after simplification
    bp_after = l1_final_detrended.l1trend_to_breakpoints(tol='auto')
    
    # Count breakpoints removed by component
    bp_removed = []
    for comp in component_list:
        if comp in bp_before and comp in bp_after:
            n_before = len(bp_before[comp][0]) - 2  # Subtract start and end points
            n_after = len(bp_after[comp][0]) - 2
            n_removed = n_before - n_after
            bp_removed.append(f"{comp}:{n_removed}")
        elif comp in bp_before:
            n_before = len(bp_before[comp][0]) - 2
            bp_removed.append(f"{comp}:{n_before}")
        else:
            bp_removed.append(f"{comp}:0")
    
    VERBOSE("  Breakpoints removed by component: %s" % ", ".join(bp_removed))
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 7: Breakpoints removed by component: {', '.join(bp_removed)}")
    

    # Step 8: apply back the velocity removed at step #0
    VERBOSE("Step 8: Applying back the velocity removed at step #0")
    #if gts_logger:
    #    gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 8: Applying back the velocity removed at step #0")
    velocity = detrended_ts.velocity

    l1_final = l1_final_detrended.remove_velocity(velocity)

    # compute a median bias with respect to the original data and remove it
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 8.1: Computing median bias")
    median_bias = np.median(l1_final.substract_ts_daily(self).data[:, 1:4],axis=0)
    l1_final.data[:, 1:4] -= median_bias
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 8: Median bias removed: [{median_bias[0]*1000:.3f}, {median_bias[1]*1000:.3f}, {median_bias[2]*1000:.3f}] mm")

    # Step 9 interpolate the final l1 trend to the original time series
    VERBOSE("Step 9: Interpolating the final L1-trend to the original time series")
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Step 9: Interpolating the final L1-trend to the original time series")
    #from pyacs.gts.Gts import Gts
    #l1_final_interpolated = Gts(code=self.code)
    l1_final_interpolated = self.copy()
    l1_final_interpolated.data = np.zeros((self.data.shape[0], 10))
    l1_final_interpolated.data[:,4:7] = 0.001
    l1_final_interpolated.data[:, 0] = self.data[:, 0]
    
    # Interpolate only the processed components, preserve original data for unprocessed components
    component_indices = {'N': 1, 'E': 2, 'U': 3}
    for comp in ['N', 'E', 'U']:
        comp_idx = component_indices[comp]
        if comp in component_list:
            # Interpolate the processed component
            l1_final_interpolated.data[:, comp_idx] = np.interp(self.data[:, 0], l1_final.data[:, 0], l1_final.data[:, comp_idx])
        else:
            # Preserve original data for unprocessed components
            l1_final_interpolated.data[:, comp_idx] = self.data[:, comp_idx]
    
    # Calculate and display final statistics
    VERBOSE("----------------------------------------------------------")
    VERBOSE("FINAL STATISTICS")
    VERBOSE("----------------------------------------------------------")
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final statistics")
    
    # Number of breakpoints per component (from the L1-trend, not the velocity-adjusted version)
    bp_final = l1_final_detrended.l1trend_to_breakpoints(tol='auto')
    for comp in component_list:
        if comp in bp_final:
            n_bp = len(bp_final[comp][0]) - 2  # Subtract start and end points
            VERBOSE("Number of breakpoints for %s component: %d" % (comp, n_bp))
            if gts_logger:
                gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: Number of breakpoints for {comp} component: {n_bp}")
        else:
            VERBOSE("Number of breakpoints for %s component: 0" % comp)
            if gts_logger:
                gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: Number of breakpoints for {comp} component: 0")
    
    # Offset dates (from the L1-trend, not the velocity-adjusted version)
    if hasattr(l1_final_detrended, 'offsets_dates') and l1_final_detrended.offsets_dates:
        VERBOSE("Offset dates: %s" % [f"{date:.3f}" for date in l1_final_detrended.offsets_dates])
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: Offset dates: {[f'{date:.3f}' for date in l1_final_detrended.offsets_dates]}")
        # Convert to human readable format
        offset_dates_readable = []
        for date in l1_final_detrended.offsets_dates:
            try:
                date_dt = at.decyear2datetime(date)
                if hasattr(date_dt, '__len__'):
                    date_str = date_dt[0].isoformat()[:10]
                else:
                    date_str = date_dt.isoformat()[:10]
                offset_dates_readable.append(date_str)
            except:
                offset_dates_readable.append(f"{date:.3f}")
        VERBOSE("Offset dates (readable): %s" % offset_dates_readable)
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: Offset dates (readable): {offset_dates_readable}")
    else:
        VERBOSE("No offset dates found")
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: No offset dates found")
    
    # Median of fit (residuals between final L1-trend and original data)
    residuals = l1_final_interpolated.substract_ts_daily(self).data[:, 1:4] * 1000  # Convert to mm
    
    # Calculate statistics only for processed components
    component_indices = {'N': 0, 'E': 1, 'U': 2}  # 0-based indices for residuals array
    processed_residuals = []
    processed_comp_names = []
    
    for comp in ['N', 'E', 'U']:
        if comp in component_list:
            comp_idx = component_indices[comp]
            processed_residuals.append(residuals[:, comp_idx])
            processed_comp_names.append(comp)
    
    if processed_residuals:
        processed_residuals = np.column_stack(processed_residuals)
        median_residuals = np.median(processed_residuals, axis=0)
        mad_residuals = np.median(np.abs(processed_residuals - median_residuals), axis=0)
        rms_residuals = np.sqrt(np.mean(processed_residuals**2, axis=0))
        
        # Create formatted strings for display
        median_str = ", ".join([f"{comp}:{val:.2f}" for comp, val in zip(processed_comp_names, median_residuals)])
        mad_str = ", ".join([f"{comp}:{val:.2f}" for comp, val in zip(processed_comp_names, mad_residuals)])
        rms_str = ", ".join([f"{comp}:{val:.2f}" for comp, val in zip(processed_comp_names, rms_residuals)])
        
        VERBOSE("Median residuals (%s): [%s] mm" % (", ".join(processed_comp_names), median_str))
        VERBOSE("MAD residuals (%s): [%s] mm" % (", ".join(processed_comp_names), mad_str))
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: MAD residuals ({', '.join(processed_comp_names)}): [{mad_str}] mm")
        
        VERBOSE("RMS residuals (%s): [%s] mm" % (", ".join(processed_comp_names), rms_str))
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: RMS residuals ({', '.join(processed_comp_names)}): [{rms_str}] mm")
    else:
        VERBOSE("No components processed, no residual statistics calculated")
        if gts_logger:
            gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Final: No components processed, no residual statistics calculated")
    
    total_time = time.time() - start_time
    VERBOSE("Optimal L1-trend workflow completed")
    VERBOSE("Total workflow time: %.2f seconds" % total_time)
    if gts_logger:
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Workflow completed successfully")
        gts_logger.info(f"[{self.code}] [optimal_l1trend_workflow] Total workflow time: {total_time:.2f} seconds")
    #VERBOSE("Step 6 (refinement) completed in %.2f seconds" % step6_time)
    
    return l1_final_interpolated
