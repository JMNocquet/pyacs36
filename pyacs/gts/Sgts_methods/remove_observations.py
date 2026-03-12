###############################################################################
def remove_observations(self, observations_to_remove, date_tolerance=0.25, verbose=False):
###############################################################################
    """Remove observations given as (code, date) pairs.

    Parameters
    ----------
    observations_to_remove : list of tuple
        (code, date) with date in decimal year.
    date_tolerance : float, optional
        Date matching tolerance in days. Default is 0.25.
    verbose : bool, optional
        Verbose mode. Default is False.

    Returns
    -------
    Sgts
        New Sgts with specified observations removed.

    Notes
    -----
    Date matching uses the same tolerance as other pyacs methods.
    """
    
    # import
    import numpy as np
    import copy
    
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    
    import inspect
    
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)
    
    # Create a deep copy of the current Sgts object
    new_sgts = copy.deepcopy(self)
    
    # Group observations by site code for efficiency
    observations_by_site = {}
    for code, date in observations_to_remove:
        if code not in observations_by_site:
            observations_by_site[code] = []
        observations_by_site[code].append(date)
    
    # Process each site
    total_removed = 0
    for code, dates_to_remove in observations_by_site.items():
        if code not in new_sgts.__dict__:
            WARNING("Site code %s not found in Sgts object" % code)
            continue
            
        gts = new_sgts.__dict__[code]
        if gts.data is None or gts.data.shape[0] == 0:
            WARNING("Site %s has no data" % code)
            continue
        
        # Get indices of observations to remove
        from pyacs.gts.Gts import get_index_from_dates
        indices_to_remove = get_index_from_dates(dates_to_remove, gts.data, tol=date_tolerance)
        
        if len(indices_to_remove) == 0:
            if verbose:
                VERBOSE("No observations found to remove for site %s" % code)
            continue
        
        # Create mask for observations to keep
        mask = np.ones(gts.data.shape[0], dtype=bool)
        mask[indices_to_remove] = False
        
        # Remove observations
        gts.data = gts.data[mask]
        
        # Update other data arrays if they exist
        if gts.data_xyz is not None:
            gts.data_xyz = gts.data_xyz[mask]
        if gts.data_corr_neu is not None:
            gts.data_corr_neu = gts.data_corr_neu[mask]
        if gts.data_corr_xyz is not None:
            gts.data_corr_xyz = gts.data_corr_xyz[mask]
        
        # Update outliers list if it exists
        if hasattr(gts, 'outliers') and gts.outliers:
            # Adjust outlier indices after removal
            new_outliers = []
            for outlier_idx in gts.outliers:
                # Count how many removed observations come before this outlier
                removed_before = sum(1 for idx in indices_to_remove if idx < outlier_idx)
                new_idx = outlier_idx - removed_before
                if new_idx >= 0 and new_idx < gts.data.shape[0]:
                    new_outliers.append(new_idx)
            gts.outliers = new_outliers
        
        total_removed += len(indices_to_remove)
        
        if verbose:
            VERBOSE("Removed %d observations from site %s" % (len(indices_to_remove), code))
    
    MESSAGE("Total observations removed: %d" % total_removed)
    
    return new_sgts
