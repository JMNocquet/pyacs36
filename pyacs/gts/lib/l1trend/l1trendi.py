"""
Main L1-trend filtering function.
"""

import numpy as np
try:
    from trendfilter import trend_filter
except:
    from pathlib import Path
    theme_path = 'bokeh_theme.yaml'
    if not theme_path.exists():
        theme_path.write_text(
            "attrs:\n"
            "    Axis:\n"
            "        major_label_text_font_size: '12pt'\n"
            "        axis_label_text_font_size: '14pt'\n"
            "    Line:\n"
            "        line_width: 3\n"
            "    Title:\n"
            "        text_font_size: '12pt'\n"
            "    Legend:\n"
            "        background_fill_alpha: 0.8\n",
            encoding='utf-8'
        )
    from trendfilter import trend_filter

import pyacs.lib.astrotime as at

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG
from time import time
import os
from datetime import datetime



def l1trendi(self,
            alpha='auto',
            ndays_mf=1,
            criterion='BIC',
            lcomponent='ENU',
            algo='golden',
            pre_process=4500,
            bounds=[-2, 1],
            tol=0.01,
            refine=True,
            min_samples_per_segment=5,
            threshold_bias_res_detect=5,
            threshold_vel=3,
            min_sample_segment_refinement=1,
            norm='L2',
            verbose=False,
            log_dir=None,
            component=None):
    """
    Performs l1 trend filter on GPS time series. For each component, optimal filter parameter is searched using specified criterion.
    Uses https://pypi.org/project/trendfilter/

    Parameters
    ----------
    alpha : float or str
        Either hyperparameter as float or 'auto' to find optimal filtering using criterion
    ndays_mf : int
        Parameter for a median filter to be applied prior to l1-trend filtering
    criterion : str
        Criterion to select automatic optimal l1-trend filtering. Options: 'BIC' (default), 'AICc', 'Cp'
    lcomponent : str
        List of components to be filtered. Other components will remain unchanged
    algo : str
        Algorithm to find optimal alpha. Options: 'golden' (default) or 'custom'
    pre_process : float
        Pre-process time series to avoid convergence problems in l1trend. It corrects for unrealistically large
        velocity using the threshold provided (mm/yr). Correction is then added again after l1_trend.
        Default is 3600 (3.6 m/yr)
    bounds : list
        Bounds for golden section search algorithm. Default is [-2, 1]
    tol : float
        Tolerance for golden section search algorithm. Default is 0.01
    refine : bool
        If True, will run refinement after l1 trend filtering
    min_samples_per_segment : int
        Minimum number of samples per segment for refinement
    threshold_bias_res_detect : float
        Threshold for bias detection in refinement
    threshold_vel : float
        Threshold for velocity detection in refinement
    min_sample_segment_refinement : int
        Minimum number of samples per segment in refinement
    norm : str
        Norm to use for refinement. Options: 'L2' (default)
    verbose : bool
        Enable verbose output for debugging
    log_dir : str or None
        If provided, directory where the log file will be created. The log file will be named 'SITE_l1trendi.log'
        where SITE is the site code. If None, no file logging is performed.
    component : str or None
        Components to process ('N', 'E', 'U', 'NE', 'EN', 'NU', 'UN', 'EU', 'UE', 'NEU', etc.).
        If None, uses lcomponent parameter. Used for component-specific statistics calculation.

    Returns
    -------
    Gts
        l1-trend filtered time series as new Gts instance

    Raises
    ------
    ValueError
        If invalid parameters are provided
    RuntimeError
        If processing fails
    """
    from .preprocessing import pre_process_test
    from .optimization import best_l1trend_golden, best_l1trend_custom

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
        
        # Log method start with [Gts.code] prefix
        gts_logger.info(f"[{self.code}] [l1trendi] Starting method")
        gts_logger.info(f"[{self.code}] [l1trendi] Parameters: criterion={criterion}, alpha={alpha}, algo={algo}")

    # Log start of processing
    str_sdate = at.decyear2datetime(self.data[0,0]).isoformat()[:10]
    str_edate = at.decyear2datetime(self.data[-1,0]).isoformat()[:10]
    # Get time period string for logging
    str_period = f"[{str_sdate} - {str_edate}]"
    VERBOSE(f"Starting l1trend processing for {self.code} {str_period} {lcomponent} {criterion} at {datetime.now().isoformat()}")

    # Validate input parameters
    if criterion not in ['BIC', 'AICc', 'Cp']:
        error_msg = f"Invalid criterion: {criterion}. Must be one of ['BIC', 'AICc', 'Cp']"
        ERROR(error_msg)
        raise ValueError(error_msg)
    if algo not in ['golden', 'custom']:
        error_msg = f"Invalid algorithm: {algo}. Must be one of ['golden', 'custom']"
        ERROR(error_msg)
        raise ValueError(error_msg)
    if not isinstance(ndays_mf, int) or ndays_mf < 1:
        error_msg = f"Invalid ndays_mf: {ndays_mf}. Must be a positive integer"
        ERROR(error_msg)
        raise ValueError(error_msg)
    if not isinstance(pre_process, (int, float)) or pre_process <= 0:
        error_msg = f"Invalid pre_process: {pre_process}. Must be positive"
        ERROR(error_msg)
        raise ValueError(error_msg)

    # Component mapping
    COMPONENT_MAP = {
        1: {'full': 'North', 'short': 'N'},
        2: {'full': 'East', 'short': 'E'},
        3: {'full': 'Up', 'short': 'U'}
    }
    
    # Determine which components to process for statistics
    if component is not None:
        # Validate component parameter
        valid_components = ['N', 'E', 'U', 'NE', 'EN', 'NU', 'UN', 'EU', 'UE', 'NEU', 'NUE', 'ENU', 'EUN', 'UNE', 'UEN']
        if component not in valid_components:
            raise ValueError(f"Invalid component '{component}'. Must be one of: {valid_components}")
        stats_components = sorted(list(set(component)))  # Remove duplicates and sort
    else:
        stats_components = sorted(list(set(lcomponent)))  # Use lcomponent if component is None

    # Criterion index mapping
    CRITERION_INDEX = {
        'BIC': -1,
        'AICc': -2,
        'Cp': -3
    }

    # Start timing
    start_time = time()

    try:
        # Create working copy
        gts = self.copy()
        
        # Apply median filter if requested
        if ndays_mf > 1:
            VERBOSE(f"Applying median filter with window size {ndays_mf}")
            gts = self.median_filter(ndays_mf)

        # Create output Gts
        ngts = gts.copy()
        x = ngts.data[:, 0]

        # Validate data
        if ngts.data.shape[0] <= 2:
            warning_msg = f"Too few data points ({ngts.data.shape[0]}) in Gts {ngts.code}. Returning original Gts {str_period}"
            WARNING(warning_msg)
            return ngts

        # Process each component
        for i in [1, 2, 3]:
            if COMPONENT_MAP[i]['short'] not in lcomponent:
                continue

            comp_name = COMPONENT_MAP[i]['full']
            VERBOSE(f"Processing {comp_name} component for {gts.code} {str_period}")

            # Handle user-defined alpha
            if isinstance(alpha, float):
                VERBOSE(f"Using fixed alpha={alpha} for {comp_name} component")
                if gts_logger:
                    gts_logger.info(f"[{self.code}] [l1trendi] Using fixed alpha={alpha} for {comp_name} component")
                ngts.data[:, i] = trend_filter(x, gts.data[:, i] * 1.E3, l_norm=1, alpha_2=alpha)['y_fit'] * 1.E-3
                continue

            # Automatic alpha search
            if alpha == 'auto':
                criterion_idx = CRITERION_INDEX[criterion]
                VERBOSE(f"Searching optimal alpha using {criterion} criterion for {comp_name} component")
                if gts_logger:
                    gts_logger.info(f"[{self.code}] [l1trendi] Searching optimal alpha using {criterion} criterion for {comp_name} component")

                # Pre-process if requested
                if pre_process > 0:
                    VERBOSE(f"Optimal alpha: Running pre-process detection for large offsets in {comp_name} component")
                    if gts_logger:
                        gts_logger.info(f"[{self.code}] [l1trendi] Optimal alpha: Running pre-process detection for large offsets in {comp_name} component")
                    lidx_offset = pre_process_test(gts, COMPONENT_MAP[i]['short'], threshold=pre_process)
                    
                    if len(lidx_offset) > 0:
                        # Handle time series split
                        idx_offset = lidx_offset[0]
                        date_offset = at.decyear2datetime(gts.data[idx_offset,0]).isoformat()[:19]
                        warning_msg = (f"Optimal alpha: Splitting Gts {self.code} {comp_name} component at {date_offset}")
                        if gts_logger:
                            gts_logger.info(f"[{self.code}] [l1trendi] {warning_msg}")
                        else:
                            VERBOSE(warning_msg)

                        
                        # Process each segment
                        segments = [
                            (gts.data[:idx_offset, :], "before"),
                            (gts.data[idx_offset:, :], "after")
                        ]
                        
                        results = []
                        for segment_data, position in segments:
                            segment_gts = gts.copy()
                            segment_gts.data = segment_data
                            
                            result = segment_gts.l1trendi(
                                alpha=alpha, ndays_mf=ndays_mf, criterion=criterion,
                                lcomponent=COMPONENT_MAP[i]['short'], algo=algo,
                                pre_process=pre_process, bounds=bounds, tol=tol,
                                refine=False,
                                min_samples_per_segment=min_samples_per_segment,
                                threshold_bias_res_detect=threshold_bias_res_detect,
                                threshold_vel=threshold_vel,
                                min_sample_segment_refinement=min_sample_segment_refinement,
                                norm=norm,
                                log_dir=log_dir
                            )
                            results.append(result)
                        
                        # Merge results
                        ngts.data[:, i] = np.hstack((results[0].data[:, i], results[1].data[:, i]))
                        date_offset_dec_year = (ngts.data[idx_offset,0]+ngts.data[idx_offset-1,0])/2
                        date_offset_str = at.decyear2datetime(date_offset_dec_year).isoformat()[:19]
                        ngts.offsets_dates.append(date_offset_dec_year)
                        continue

                # Process without splitting
                idx_offset = []
                try:
                    # Create component mask for this specific component
                    current_comp = COMPONENT_MAP[i]['short']
                    component_mask = np.ones(len(x), dtype=bool)  # Default: all points active
                    
                    # If this component is not in the stats_components, create a mask that excludes zero values
                    if current_comp not in stats_components:
                        # Only consider points where this component has non-zero values
                        component_mask = np.abs(gts.data[:, i]) > 1e-10
                    
                    if algo == 'golden':
                        ngts.data[:, i] = best_l1trend_golden(x, gts.data[:, i] * 1.E3, criterion_idx, bounds=bounds, tol=tol, component_mask=component_mask) * 1.E-3
                        VERBOSE(f"Optimal alpha: Successfully processed {comp_name} component with algo golden")
                        if gts_logger:
                            gts_logger.info(f"[{self.code}] [l1trendi] Optimal alpha: Successfully processed {comp_name} component with algo golden")
                    else:
                        ngts.data[:, i] = best_l1trend_custom(x, gts.data[:, i] * 1.E3, criterion_idx, component_mask=component_mask)[0] * 1.E-3
                        VERBOSE(f"Optimal alpha: Successfully processed {comp_name} component with algo custom")
                        if gts_logger:
                            gts_logger.info(f"[{self.code}] [l1trendi] Optimal alpha: Successfully processed {comp_name} component with algo custom")

                except Exception as e:
                    warning_msg = f"Failed in algo {algo} search to process {comp_name} component: {str(e)}"
                    if gts_logger:
                        gts_logger.error(f"[{self.code}] [l1trendi] {warning_msg}")
                    ERROR(warning_msg, exit=True)
                    continue

        # print offsets_dates
        ngts.offsets_dates = sorted(set(ngts.offsets_dates))
        mesg = f"List of offsets detected for {ngts.code}"
        if gts_logger:
            gts_logger.info(f"[{self.code}] [l1trendi] {mesg}")
        else:
            VERBOSE(mesg)

        for date_offset in ngts.offsets_dates:
            mesg = f"{date_offset:.5f} {at.decyear2datetime(date_offset).isoformat()[:19]}"
            if gts_logger:
                gts_logger.info(f"[{self.code}] [l1trendi] {mesg}")
            VERBOSE(mesg)

        # Apply refinement if requested
        if refine:
            VERBOSE(f"Applying refinement for {ngts.code} {lcomponent} {str_period}")
            if gts_logger:
                gts_logger.info(f"[{self.code}] [l1trendi] Applying refinement for {ngts.code} {lcomponent} {str_period}")
            ngts = ngts.refine_l1trend(
                self,
                min_samples_per_segment=min_samples_per_segment,
                threshold_bias_res_detect=threshold_bias_res_detect,
                threshold_vel=threshold_vel,
                min_sample_segment_refinement=min_sample_segment_refinement,
                norm=norm,
                lcomponent=lcomponent
            )

        # Log completion
        elapsed = time() - start_time
        VERBOSE(f"Completed l1trendi processing for {self.code} in {elapsed:.1f} seconds")
        if gts_logger:
            gts_logger.info(f"[{self.code}] [l1trendi] Completed processing in {elapsed:.1f} seconds")
            gts_logger.info(f"[{self.code}] [l1trendi] Method completed successfully")
        
        return ngts

    except Exception as e:
        error_msg = f"Failed to process Gts {self.code}: {str(e)}"
        if gts_logger:
            gts_logger.error(f"[{self.code}] [l1trendi] {error_msg}")
        raise RuntimeError(f"l1trend processing failed: {str(e)}")
        ERROR(error_msg, exit=True)
