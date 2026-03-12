###############################################################################
def get_dates(self, fmt='isoformat'):
###############################################################################
    """Return start and end dates of all time series in this Sgts.

    Parameters
    ----------
    fmt : str, optional
        Output format: 'isoformat', 'datetime', 'decyear', 'mjd', 'dayno', 'cal'. Default is 'isoformat'.

    Returns
    -------
    tuple
        (start_date, end_date) in the requested format.

    Notes
    -----
    Formats: isoformat (ISO 8601), datetime, decyear, mjd, dayno (dayno, year), cal (day, month, year).
    """
    
    # import
    import numpy as np
    import pyacs.lib.astrotime as at
    
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    
    import inspect
    
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)
    
    # Check if Sgts has any time series
    if len(self.lcode()) == 0:
        ERROR("Sgts object contains no time series", exit=True)
    
    # Collect all start and end dates
    start_dates = []
    end_dates = []
    
    for code in self.lcode():
        gts = self.__dict__[code]
        
        # Check if Gts has data
        if gts.data is None or gts.data.shape[0] == 0:
            WARNING("Site %s has no data, skipping" % code)
            continue
        
        # Get start and end dates (first and last decimal years)
        start_date = gts.data[0, 0]   # First observation date
        end_date = gts.data[-1, 0]    # Last observation date
        
        start_dates.append(start_date)
        end_dates.append(end_date)
        
        DEBUG("Site %s: %f to %f" % (code, start_date, end_date))
    
    if len(start_dates) == 0:
        ERROR("No valid time series found in Sgts object", exit=True)
    
    # Convert to numpy arrays and find global min/max
    start_dates = np.array(start_dates)
    end_dates = np.array(end_dates)
    
    global_start = np.min(start_dates)
    global_end = np.max(end_dates)
    
    VERBOSE("Global date range: %f to %f (decimal years)" % (global_start, global_end))
    
    # Convert to requested format
    if fmt == 'decyear':
        return (global_start, global_end)
    
    elif fmt == 'datetime':
        start_datetime = at.decyear2datetime(global_start)
        end_datetime = at.decyear2datetime(global_end)
        return (start_datetime, end_datetime)
    
    elif fmt == 'isoformat':
        start_datetime = at.decyear2datetime(global_start)
        end_datetime = at.decyear2datetime(global_end)
        return (start_datetime.isoformat(), end_datetime.isoformat())
    
    elif fmt == 'mjd':
        start_mjd = at.decyear2mjd(global_start)
        end_mjd = at.decyear2mjd(global_end)
        return (start_mjd, end_mjd)
    
    elif fmt == 'dayno':
        start_dayno = at.decyear2dayno(global_start)
        end_dayno = at.decyear2dayno(global_end)
        return (start_dayno, end_dayno)
    
    elif fmt == 'cal':
        start_cal = at.decyear2cal(global_start)
        end_cal = at.decyear2cal(global_end)
        return (start_cal, end_cal)
    
    else:
        ERROR("Unknown format '%s'. Available formats: 'isoformat', 'datetime', 'decyear', 'mjd', 'dayno', 'cal'" % fmt, exit=True)
